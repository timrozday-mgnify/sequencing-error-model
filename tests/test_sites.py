import os
import shutil
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from sequencing_error_model import bias, recovery, sites
from sequencing_error_model import generate as gen
from sequencing_error_model.sources import bam, skiver_analyze


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
