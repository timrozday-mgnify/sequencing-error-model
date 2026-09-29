import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from sequencing_error_model import bias, recovery, sites
from sequencing_error_model import generate as gen
from sequencing_error_model.sources import bam


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


def _reference_ratio(tmp_path: Path, ani: float, *, site_mask: bool) -> tuple[float, float | None]:
    """The `reference` mode's rate ratio on one population, with and without the mask."""
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
    plain, none = _reference_ratio(tmp_path, 0.95, site_mask=False)
    masked, share = _reference_ratio(tmp_path, 0.95, site_mask=True)
    assert none is None
    # ~2.5% root-to-tip divergence on a ~4% error rate, so the unmasked mode reads ~1.5x the truth (phase 7).
    assert plain > 1.35, plain
    assert masked < 1.15, (plain, masked)
    assert share is not None and share < 0.25, share  # and it pays for it with a modest share of the sites
    # The clonal point is where the mask must cost nothing: it drops low-depth and multi-copy sites there too.
    clonal, clonal_share = _reference_ratio(tmp_path, 1.0, site_mask=True)
    assert abs(clonal - 1) < 0.1, clonal
    assert clonal_share is not None and clonal_share < 0.2, clonal_share
