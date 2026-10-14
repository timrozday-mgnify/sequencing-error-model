import os
import shutil
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from sequencing_error_model import bias, recovery, sites
from sequencing_error_model import generate as gen
from sequencing_error_model.fit import kmer
from sequencing_error_model.fit.quality import Array
from sequencing_error_model.observations import CountTable
from sequencing_error_model.sources import bam, skiver_analyze
from sequencing_error_model.spec import Component


def _aln(template: str, read: str, start: int = 0, reverse: bool = False) -> bam.Alignment:
    """One end-to-end aligned record, already in read orientation (as `bam._aligned` yields them)."""
    r = gen.Read(read, "I" * len(read), f"{len(read)}M", np.zeros(0, np.int64))
    return ("g", start, start + len(template), reverse, template, r, 1)


def test_alleles_are_reported_in_reference_orientation() -> None:
    ref = "ACGTACGTAC"
    # Forward: a substitution at reference site 3 reads C.
    assert sites._alleles(0, len(ref), False, ref, _aln(ref, "ACGCACGTAC")[5]) == [(3, "C")]
    # Reverse: the record's template and read are reverse-complemented, so template index t is site n-1-t and
    # the read base is complemented back.
    rc = gen._revcomp(ref)
    assert sites._alleles(0, len(ref), True, rc, _aln(rc, rc[:6] + "C" + rc[7:])[5]) == [(3, "G")]


def test_linked_needs_the_same_allele_pair_on_several_reads() -> None:
    ref = "ACGTACGTACGTACGTACGT"
    both = ref[:4] + "C" + ref[5:12] + "G" + ref[13:]  # the same two alleles, on one haplotype
    apart = [ref[:4] + "C" + ref[5:], ref[:12] + "G" + ref[13:]]  # each allele alone: no pair to recur
    assert sites.linked([_aln(ref, both)] * 3) == {"g": {4, 12}}
    assert sites.linked([_aln(ref, both)] * 2) == {}  # below min_pairs
    assert sites.linked([_aln(ref, s) for s in apart] * 5) == {}
    # A different partner allele each time is what independent errors look like.
    scattered = [_aln(ref, ref[:4] + "C" + ref[5:i] + "G" + ref[i + 1 :]) for i in range(8, 18)]
    assert sites.linked(scattered) == {}


def test_conservative_drops_low_depth_multi_copy_and_linked_sites() -> None:
    ref = "ACGTACGTACGTACGTACGTAAGGCCTTAACG"
    linked_read = ref[:4] + "C" + ref[5:12] + "G" + ref[13:]
    # 13 reads (3 with one linked allele pair) over all but the last four bases, and 20 more over the first
    # eight only - so those sit at 33x against a contig median of 13x, and the tail has no depth at all.
    alignments = (
        [_aln(ref[:28], ref[:28])] * 10 + [_aln(ref[:28], linked_read[:28])] * 3 + [_aln(ref[:8], ref[:8])] * 20
    )
    counts = bam.count_alleles(alignments, {"g": len(ref)})
    masked, report = sites.conservative(counts, sites.linked(alignments), lambda c: ref, min_depth=5.0, copy_ratio=2.0)
    mask = masked["g"]
    assert mask[:8].all()  # 33x against a median of 13: two collapsed copies
    assert mask[28:].all()  # no depth
    assert mask[12] and not mask[9]  # the linked pair's second site, and an ordinary site beside it
    row = report[0]
    assert row["dropped_low_depth"] == 5  # the four uncovered bases and the trailing bookkeeping site
    assert row["dropped_multi_copy"] == 8
    assert row["dropped_linked"] == 1  # site 4 is inside the doubled region, so linkage is charged site 12
    assert row["dropped"] == int(mask.sum())
    assert abs(sum(row["kept_composition"]["context"].values()) - 1) < 1e-3  # rounded shares
    # Linkage alone, on a contig with flat depth, drops exactly the two variant sites.
    flat = [_aln(ref, ref)] * 10 + [_aln(ref, linked_read)] * 6
    masked, report = sites.conservative(
        bam.count_alleles(flat, {"g": len(ref)}), sites.linked(flat), lambda c: ref, min_depth=5.0
    )
    assert np.flatnonzero(masked["g"][: len(ref)]).tolist() == [4, 12]
    assert report[0]["dropped_linked"] == 2
    # What it dropped is what the variant reads carry: a C and a G, both at a reference A.
    assert report[0]["dropped_composition"]["ops"] | {} == {
        "alt_A": 0.0,
        "alt_C": 0.5,
        "alt_G": 0.5,
        "alt_T": 0.0,
        "deletion": 0.0,
        "insertion": 0.0,
    }
    assert report[0]["dropped_composition"]["context"] == {"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0}


def _reference_ratio(tmp_path: Path, ani: float, *, site_mask: str | None) -> tuple[float, float | None]:
    """The `reference` mode's rate ratio on one population, with and without a site model."""
    truth = recovery.scale_error_rate(recovery.example_spec(), 0.2)
    doc = bias.scenario(
        truth,
        tmp_path / f"ani{ani}-mask{site_mask}",
        ani=ani,
        n_haplotypes=1 if ani >= 1.0 else 2,
        minor_fraction=0.5,
        coverage=18.0,
        genome_length=6000,
        read_length=120,
        insert_mean=180.0,
        insert_sd=15.0,
        modes=("reference",),
        q_reads=400,
        seed=1,
        site_mask=site_mask,
    )
    mode = doc["modes"]["reference"]
    assert "skipped" not in mode, mode
    return mode["rate_ratio"], bias._masked_fraction(mode["site_mask"])


def test_the_conservative_mask_removes_most_of_the_reference_modes_variant_inflation(tmp_path: Path) -> None:
    if shutil.which("minibwa") is None:  # as tests/test_bias.py: CI installs the aligners
        if os.environ.get("REQUIRE_ALIGNERS") == "1":
            pytest.fail("minibwa is not on PATH")
        pytest.skip("minibwa is not on PATH")
    plain, none = _reference_ratio(tmp_path, 0.95, site_mask=None)
    masked, share = _reference_ratio(tmp_path, 0.95, site_mask="conservative")
    assert none is None
    # ~2.5% root-to-tip divergence on a ~4% error rate, so the unmasked mode reads ~1.5x the truth (phase 7).
    assert plain > 1.35, plain
    assert masked < 1.15, (plain, masked)
    assert share is not None and share < 0.25, share  # and it pays for it with a modest share of the sites
    # The clonal point is where the mask must cost nothing: it drops low-depth and multi-copy sites there too.
    clonal, clonal_share = _reference_ratio(tmp_path, 1.0, site_mask="conservative")
    assert abs(clonal - 1) < 0.1, clonal
    assert clonal_share is not None and clonal_share < 0.2, clonal_share


_JOINT_TOKENS = ["QualityWindow(0)", "Context(1,1)"]


def _pileup(
    rng: np.random.Generator, ref: str, depth: int, rate: float, variants: dict[int, int]
) -> list[bam.Alignment]:
    """`depth` end-to-end records over `ref` with independent errors at `rate`, and, at each site of
    `variants`, that many records carrying one alternative allele - a variant of known frequency."""
    alt = {s: "ACGT"[("ACGT".index(ref[s]) + 1) % 4] for s in variants}
    out = []
    for i in range(depth):
        read = np.array(list(ref))
        hit = rng.random(len(ref)) < rate
        read[hit] = rng.choice(list("ACGT"), size=int(hit.sum()))
        for s, carriers in variants.items():
            if i < carriers:
                read[s] = alt[s]
        out.append(_aln(ref, "".join(read)))
    return out


def test_joint_finds_a_variant_its_own_head_e_cannot_explain_and_costs_the_clonal_case_nothing() -> None:
    rng = np.random.default_rng(0)
    ref = "".join(rng.choice(list("ACGT"), size=800))
    alignments = _pileup(rng, ref, 60, 0.02, {400: 24})  # one site at 40%, on a 2% error rate
    counts = bam.count_alleles(alignments, {"g": len(ref)})
    masked, report = sites.joint(alignments, counts, lambda c: ref, _JOINT_TOKENS, flank=(1, 1), m=0)
    assert masked["g"][400]
    row = report[0]
    assert row["dropped_variant"] <= 2  # the variant, and at most one site where errors piled up by chance
    assert row["dropped"] == int(masked["g"].sum())
    # The clonal pileup is where head E's expected counts *are* the observation, so nothing is called a variant.
    clonal = _pileup(rng, ref, 60, 0.02, {})
    _, clonal_report = sites.joint(
        clonal, bam.count_alleles(clonal, {"g": len(ref)}), lambda c: ref, _JOINT_TOKENS, flank=(1, 1), m=0
    )
    assert clonal_report[0]["dropped_variant"] == 0
    assert clonal_report[0]["dropped_multi_copy"] == 0


def test_joint_scores_the_count_against_head_e_not_against_a_frequency_threshold() -> None:
    """A minor allele below the error rate is not separable (§6.6 signal 2), and the model says so rather than
    dropping every site with a non-reference base."""
    rng = np.random.default_rng(1)
    ref = "".join(rng.choice(list("ACGT"), size=800))
    # 1 read in 60 carries the allele: the same count two errors at that site would make.
    alignments = _pileup(rng, ref, 60, 0.02, {400: 1})
    counts = bam.count_alleles(alignments, {"g": len(ref)})
    masked, report = sites.joint(alignments, counts, lambda c: ref, _JOINT_TOKENS, flank=(1, 1), m=0)
    assert not masked["g"][400]
    assert report[0]["dropped_variant"] <= 1


def test_the_joint_model_removes_the_reference_modes_variant_inflation(tmp_path: Path) -> None:
    if shutil.which("minibwa") is None:  # as tests/test_bias.py: CI installs the aligners
        if os.environ.get("REQUIRE_ALIGNERS") == "1":
            pytest.fail("minibwa is not on PATH")
        pytest.skip("minibwa is not on PATH")
    plain, _ = _reference_ratio(tmp_path, 0.95, site_mask=None)
    joint, share = _reference_ratio(tmp_path, 0.95, site_mask="joint")
    assert plain > 1.35, plain
    assert joint < 1.15, (plain, joint)
    assert share is not None and share < 0.25, share
    clonal, clonal_share = _reference_ratio(tmp_path, 1.0, site_mask="joint")
    assert abs(clonal - 1) < 0.1, clonal
    assert clonal_share is not None and clonal_share < 0.2, clonal_share


def test_keys_drops_a_planted_variant_and_leaves_the_clonal_fixture_alone() -> None:
    """The `kmer` key test on the committed skiver fixture (713 keys at ~44x, clonal, 1.1% error): it calls no
    key a variant as it stands, and with a 40% allele planted in 30 keys - the count a minor strain would leave -
    it calls those 30 and carries their excess into the level."""
    a = skiver_analyze.read_analyze(Path(__file__).parent / "fixtures" / "skiver-v0.3.2" / "analyze")
    tokens = ["QualityWindow(0)", "Context(1,1)"]
    _, clean = sites.keys(a, tokens)
    assert clean["dropped_variant"] == 0
    assert abs(clean["level"] - 1) < 0.01

    planted = {}
    for i, r in enumerate(a.kvmer[:30]):
        moved = int(0.4 * r.consensus_count)
        op = f"{r.consensus_value[6]}>{'C' if r.consensus_value[6] == 'A' else 'A'}"
        planted[i] = replace(
            r, consensus_count=r.consensus_count - moved, op_counts={**r.op_counts, op: r.op_counts[op] + moved}
        )
    variant = replace(a, kvmer=[planted.get(i, r) for i, r in enumerate(a.kvmer)])
    tested, report = sites.keys(variant, tokens)
    called = {i for i, r in enumerate(tested.kvmer) if not r.passes_filter}
    assert set(planted) <= called, sorted(set(planted) - called)
    assert report["dropped_variant"] <= 30 + 2
    assert report["level"] < 0.9  # the planted keys' excess, removed from summary_phred.csv's level


def test_the_multiplicity_flag_is_free_on_clonal_keys_and_fires_on_far_values() -> None:
    """Phase 6b item 1's reference-free multiplicity flag: a key whose observations often sit two or more edits
    from its consensus is another locus, since two independent errors cost the rate squared.

    The committed fixture is a random genome, so every key is single-locus at k = 11 and the flag must cost
    nothing: its far mass tops out below the default threshold. Planting one key's mass at a far value fires it,
    and `far_fraction=1.0` switches it off.
    """
    a = skiver_analyze.read_analyze(Path(__file__).parent / "fixtures" / "skiver-v0.3.2" / "analyze")
    tokens = ["QualityWindow(0)", "Context(1,1)"]
    _, clean = sites.keys(a, tokens)
    assert clean["dropped_multi_locus"] == 0
    assert clean["far_fraction_median"] == 0.0
    assert abs(clean["level"] - 1) < 0.01  # no level cost, because no key is dropped

    # A far value on one key: a third of its observations two or more edits off, which no error rate explains.
    r = a.kvmer[0]
    far = replace(r, total_count=r.total_count + r.consensus_count // 2)
    planted = replace(a, kvmer=[far, *a.kvmer[1:]])
    _, report = sites.keys(planted, tokens)
    assert report["dropped_multi_locus"] == 1
    # One key of 713 is a small share of a clonal run's far observations, which is the point: clonal far mass
    # is diffuse and low per key (median 0, max 0.077 here), while a multi-locus key concentrates it.
    assert 0 < report["far_observations_dropped"] < 0.1
    assert sites.keys(planted, tokens, far_fraction=1.0)[1]["dropped_multi_locus"] == 0
    # It selects on the other values, not on the single-edit outcome head E is fitted on, so the level stands.
    assert abs(report["level"] - 1) < 0.02


def _plant_allele(
    a: skiver_analyze.SkiverAnalyze, keys: set[int], pos: int, frac: float = 0.4
) -> skiver_analyze.SkiverAnalyze:
    """Move `frac` of each chosen key's consensus observations onto one alt base at value position `pos`.

    A minor allele as `kvmer.csv` records one: the op count rises, and because the running counts are a
    *survival* series, every `consensus_count_up_to_v{t}` from `pos` on falls by the same number.
    """
    rows = []
    for i, r in enumerate(a.kvmer):
        if i not in keys:
            rows.append(r)
            continue
        moved = int(frac * r.consensus_count)
        base = r.consensus_value[pos - 1]
        op = f"{base}>{'C' if base == 'A' else 'A'}"
        rows.append(
            replace(
                r,
                consensus_count=r.consensus_count - moved,
                neighbor_count=r.neighbor_count + moved,
                op_counts={**r.op_counts, op: r.op_counts[op] + moved},
                consensus_count_up_to_v=tuple(
                    c - moved if t >= pos else c for t, c in enumerate(r.consensus_count_up_to_v, 1)
                ),
            )
        )
    return replace(a, kvmer=rows)


def _kvmer_tables(a: skiver_analyze.SkiverAnalyze) -> tuple[CountTable, CountTable]:
    tables = {t.source.split(":")[1]: t for t in skiver_analyze.tables(a, use_all=True)}
    return tables["kvmer"], tables["kvmer_survival"]


def _context_log_odds(head: Sequence[Component]) -> Array:
    w = next(c for c in head if c.name == "Context").params["weights"]
    return np.asarray(w[..., 1:].sum(axis=-1) - w[..., 0])


def test_censoring_a_variant_position_costs_a_tenth_of_discarding_the_key() -> None:
    """Phase 8: censor the variant position, do not discard the key.

    A 40% allele planted at one value position of 60 keys. Censoring must find those positions, and it must
    cost an order of magnitude less exposure than dropping the keys whole - a key is k + v bases, so dropping
    it for one site removes the other positions' observations too. The distortion it leaves on head E's
    context is measured against the same fit on the clean fixture.
    """
    pos = 6
    a = skiver_analyze.read_analyze(Path(__file__).parent / "fixtures" / "skiver-v0.3.2" / "analyze")
    tokens = ["QualityWindow(0)", "Context(1,1)"]
    planted = _plant_allele(a, set(range(60)), pos)

    def head(x: skiver_analyze.SkiverAnalyze, censored: dict[str, list[int]] | None = None) -> tuple[Component, ...]:
        kv, sv = _kvmer_tables(x)
        return kmer.fit(kv, tokens, flank=(1, 1), survival=sv, censored=censored)

    censored, report = sites.censor(planted, tokens)
    planted_loci = {r.key + r.consensus_value for r in planted.kvmer[:60]}
    assert {k for k, v in censored.items() if pos in v} >= planted_loci
    # Every call is at the planted position: nothing smeared over the key's other positions.
    assert {p for k, v in censored.items() if k in planted_loci for p in v} == {pos}
    assert report["censored_exposure_fraction"] * 8 < report["whole_key_exposure_fraction"]

    truth = _context_log_odds(head(a))
    naive = np.abs(_context_log_odds(head(planted)) - truth).mean()
    with_censoring = np.abs(_context_log_odds(head(planted, censored)) - truth).mean()
    # Censoring removes most of the allele's pull on the context shape, at a tenth of the exposure.
    assert with_censoring < naive / 3, (with_censoring, naive)


def test_censoring_is_nearly_free_on_clonal_keys() -> None:
    a = skiver_analyze.read_analyze(Path(__file__).parent / "fixtures" / "skiver-v0.3.2" / "analyze")
    censored, report = sites.censor(a, ["QualityWindow(0)", "Context(1,1)"])
    assert report["positions"] <= 3  # 713 keys x 13 positions, clonal at ~1.1% error
    assert report["censored_exposure_fraction"] < 0.001
    # Whole-key discard would cost an order of magnitude more for the same calls.
    assert report["whole_key_exposure_fraction"] > 8 * report["censored_exposure_fraction"]
    assert sum(len(v) for v in censored.values()) == report["positions"]


def test_q_diagnostic_and_fit_follow_miscalibrated_q_through_variation(tmp_path: Path) -> None:
    """The "Q is never a label" guard with variants (§10, phase 8): the truth's Q dependence reversed, so the top
    reported Q carries the worst error, then ANI 95% and the joint site model. The fitted rate rises with Q as the
    injected truth does - a fit reading 10^(-Q/10) would fall - and so does the diagnostic's kept-site curve,
    while the dropped sites read flat and high, as real alleles do."""
    if shutil.which("minibwa") is None:
        if os.environ.get("REQUIRE_ALIGNERS") == "1":
            pytest.fail("minibwa is not on PATH")
        pytest.skip("minibwa is not on PATH")
    # scale 0.05 puts the reversed truth's marginal (~4.6%) at the calibrated grid's ~4.2%
    truth = recovery.miscalibrate(recovery.scale_error_rate(recovery.example_spec(), 0.05))
    doc = bias.scenario(
        truth,
        tmp_path,
        ani=0.95,
        n_haplotypes=2,
        minor_fraction=0.5,
        coverage=18.0,
        genome_length=6000,
        read_length=120,
        insert_mean=180.0,
        insert_sd=15.0,
        modes=("reference",),
        q_reads=400,
        seed=1,
        site_mask="joint",
    )
    mode = doc["modes"]["reference"]
    assert abs(mode["rate_ratio"] - 1) < 0.1, mode["rate_ratio"]
    assert (np.diff(mode["rate_by_q_fit"]) > 0).all(), mode["rate_by_q_fit"]
    # Q >= 12 within 20%; the lowest-error bin (now Q2) carries the alignment's mismatch floor, as the clonal run
    # does at 1.7x, and is bounded rather than held to the tolerance.
    ratio = mode["rate_by_q_ratio"]
    assert max(abs(r - 1) for r in ratio[1:]) < 0.2 and ratio[0] < 2.5, ratio

    kept, dropped = (mode["q_diagnostic"][side] for side in ("kept", "dropped"))
    assert (np.diff([kept[q]["mismatch_rate"] for q in sorted(kept)]) > 0).all(), kept
    low, high = min(dropped), max(dropped)
    assert dropped[low]["mismatch_rate"] > 10 * kept[low]["mismatch_rate"], (kept, dropped)
    assert dropped[high]["mismatch_rate"] > 2 * kept[high]["mismatch_rate"], (kept, dropped)
