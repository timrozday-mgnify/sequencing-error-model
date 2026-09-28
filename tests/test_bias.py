import json
import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from sequencing_error_model import bias, recovery
from sequencing_error_model import generate as gen
from sequencing_error_model import variation as var
from sequencing_error_model.spec import ErrorModelSpec


def _require(tool: str) -> None:
    """CI sets REQUIRE_ALIGNERS=1 and installs the aligners; elsewhere a missing one skips."""
    if shutil.which(tool) is None:
        if os.environ.get("REQUIRE_ALIGNERS") == "1":
            pytest.fail(f"{tool} is not on PATH")
        pytest.skip(f"{tool} is not on PATH")


def _quiet() -> ErrorModelSpec:
    """The example truth at a fifth of its error rate: overlaps that noisy fail placement (phase 4)."""
    return recovery.scale_error_rate(recovery.example_spec(), 0.2)


def test_scale_error_rate_shifts_only_the_non_match_logits() -> None:
    truth = recovery.example_spec()
    scaled = recovery.scale_error_rate(truth, 0.25)
    a, b = (s.error_head[0].params["bias"] for s in (truth, scaled))
    assert b[0] == a[0]
    assert np.allclose(b[1:] - a[1:], np.log(0.25))


def test_labels_size_the_problem_and_a_clonal_population_has_no_variant() -> None:
    truth = _quiet()
    rng = np.random.default_rng(0)
    genome = "".join(rng.choice(list("ACGT"), size=4000))
    for ani, n, density in ((1.0, 1, 0.0), (0.98, 2, 0.0), (1.0, 1, 0.01)):
        pop = bias.population([("g", genome)], rng, ani=ani, n_haplotypes=n, minor_density=density)
        frags = var.fragments(pop, pop.samples[0], gen.insert_sizes(180, 20), 300, 120, rng)
        templates = [t for f in frags for t in f.templates]
        mates = [1, 2] * len(frags)
        got = bias.labels(frags, gen.generate(truth, templates, mates, rng))
        # One label per *read* base, so indels move the total off the template length by the indel rate.
        assert abs(got["bases"] / (2 * len(frags) * 120) - 1) < 0.02
        assert abs(sum(got[k] for k in ("match", "variant", "error", "adapter")) - 1) < 1e-9
        assert got["error"] > 0
        # Variation is the only thing that makes a read base differ from the consensus without an error.
        assert (got["variant"] > 0) == (ani < 1.0 or density > 0)


def test_pe_overlap_is_variant_immune_where_the_reference_mode_is_not(tmp_path: Path) -> None:
    _require("minibwa")
    truth = _quiet()
    rates = {}
    for ani in (1.0, 0.95):
        doc = bias.scenario(
            truth,
            tmp_path / f"ani{ani}",
            ani=ani,
            n_haplotypes=1 if ani == 1.0 else 2,
            minor_fraction=0.5,
            coverage=12.0,
            genome_length=6000,
            read_length=120,
            insert_mean=180.0,
            insert_sd=15.0,
            modes=("pe-overlap", "reference"),
            q_reads=400,
            seed=1,
        )
        assert doc["labels"]["variant"] > 0.02 if ani < 1.0 else doc["labels"]["variant"] == 0
        rates[ani] = {m: doc["modes"][m]["rate_ratio"] for m in ("pe-overlap", "reference")}
    # Both modes recover the clonal point; only `reference` counts the strains' real differences as errors.
    assert all(abs(r - 1) < 0.1 for r in rates[1.0].values()), rates
    assert abs(rates[0.95]["pe-overlap"] - 1) < 0.12, rates
    # ~2.5% root-to-tip divergence on top of a ~3% error rate, so the reference mode reads ~1.5x the truth.
    assert rates[0.95]["reference"] > 1.4, rates
    assert rates[0.95]["reference"] - rates[0.95]["pe-overlap"] > 0.3, rates


def test_grid_records_a_missing_skiver_and_still_tabulates(tmp_path: Path) -> None:
    truth = _quiet()
    points = bias.grid(
        truth,
        tmp_path,
        ani=(1.0,),
        coverage=(8.0,),
        modes=("pe-overlap", "kmer"),
        genome_length=3000,
        read_length=100,
        insert_mean=150.0,
        insert_sd=12.0,
        q_reads=200,
        skiver=None,
    )
    assert len(points) == 1 and points[0]["clonal"]
    rows = {row["mode"]: row for row in bias.table(points)}
    assert set(rows) == {"pe-overlap", "kmer"}
    assert "skiver" in rows["kmer"]["skipped"]
    assert rows["pe-overlap"]["rate_ratio"] is not None
    assert json.dumps(points)  # the report is JSON, so the grid can be written and diffed


def test_pair_flags_sets_the_mate_bits_from_the_name_suffix(tmp_path: Path) -> None:
    sam = tmp_path / "a.sam"
    sam.write_text(
        "@HD\tVN:1.6\nr1:m1\t0\tg\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\nr1:m2\t16\tg\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
    )
    bias._pair_flags(sam)
    rows = [line.split("\t") for line in sam.read_text().splitlines() if not line.startswith("@")]
    assert [r[0] for r in rows] == ["r1", "r1"]
    assert [int(r[1]) for r in rows] == [0x1 | 0x40, 16 | 0x1 | 0x80]


def test_reference_choices_are_not_the_reads_own_sequence() -> None:
    rng = np.random.default_rng(2)
    genome = "".join(rng.choice(list("ACGT"), size=3000))
    pop = bias.population([("g", genome)], rng, ani=0.98, n_haplotypes=2)
    seqs = {which: bias._reference(pop, which, 0.95, rng)[0][1] for which in bias.REFERENCES}
    assert seqs["consensus"] == genome
    assert seqs["majority"] == pop.haplotypes["hap1"][0].sequence
    assert seqs["relative"] != genome
    with pytest.raises(ValueError, match="unknown reference"):
        bias._reference(pop, "assembly", 0.95, rng)
