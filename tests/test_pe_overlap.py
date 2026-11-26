from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from sequencing_error_model import generate as gen
from sequencing_error_model import recovery
from sequencing_error_model import spec as spec_io
from sequencing_error_model.fit import error
from sequencing_error_model.sources import pe_overlap
from sequencing_error_model.spec import Component


def test_place_finds_offset_read_through_and_indels() -> None:
    rng = np.random.default_rng(0)
    for size in (30, 45, 60, 75):
        frag = "".join(rng.choice(list("ACGT"), size=size))
        r1, r2 = (
            (s + a)[:50].ljust(50, "N") for s, a in ((frag, gen.ADAPTERS[0]), (gen._revcomp(frag), gen.ADAPTERS[1]))
        )
        assert pe_overlap.place(r1, r2) == (size - 50, False)
        mid = 49 - (max(0, size - 50) + min(50, size)) // 2 + (size - 50)
        mutated = r2[:mid] + r2[mid + 1 :] + "A"  # delete a base mid-overlap
        s, indel = pe_overlap.place(r1, mutated)
        assert s is None or indel, size
    assert pe_overlap.place("ACGT" * 10, "") == (None, False)
    # A tandem repeat overlaps equally well at several offsets: ambiguous, so not placed.
    frag = "".join(rng.choice(list("ACGT"), size=15)) * 5
    assert pe_overlap.place(frag[:50], gen._revcomp(frag)[:50]) == (None, False)
    # Unrelated mates (long inserts) are neither placed nor mistaken for indels.
    for _ in range(50):
        r1, r2 = ("".join(rng.choice(list("ACGT"), size=125)) for _ in range(2))
        assert pe_overlap.place(r1, r2) == (None, False)


def test_cli_fits_spec_from_fastq(tmp_path: Path) -> None:
    model = recovery.example_spec()
    rng = np.random.default_rng(1)
    templates, mates = recovery.paired_draw(40, 45, 8, genome_length=5000)(rng, 400)
    reads = gen.generate(model, templates, mates, rng, error_rate_scale=0.2)
    for mate in (1, 2):
        records = (f"@p{i}/{mate}\n{r.sequence}\n+\n{r.quality}\n" for i, r in enumerate(reads[mate - 1 :: 2]))
        (tmp_path / f"r{mate}.fastq").write_text("".join(records))
    out = tmp_path / "spec"
    args = [str(tmp_path / "r1.fastq"), str(tmp_path / "r2.fastq"), "--output", str(out), "--iterations", "2"]
    assert pe_overlap.main(args) == 0
    spec = spec_io.load(out)
    assert spec.provenance["stats"]["pairs"] == 200 and spec.provenance["identified_ops"] == ["substitution"]
    assert np.isneginf(spec.error_head[0].params["bias"][5:]).all()


def test_pe_overlap_recovers_substitution_head() -> None:
    example = recovery.example_spec()
    head0 = example.error_head[0]
    # ~3.5% errors, not the example's 14%: overlaps that noisy fail placement. (At ~1% the per-Q rates need
    # ~10k pairs to settle within 10%, too slow for CI.) Q37, the lowest-rate bin, needs 5,000 pairs here.
    bias = head0.params["bias"] - np.r_[0, np.full(9, 1.5)]
    quiet = (replace(head0, params={**head0.params, "bias": bias}), *example.error_head[1:])
    truth = replace(example, error_head=quiet)
    fitted, report = recovery.recover(
        truth,
        n_reads=10000,
        seed=4,
        mode=pe_overlap.mode,
        draw=recovery.paired_draw(40, 45, 8),
        substitutions_only=True,
    )
    stats = fitted.provenance["stats"]
    assert stats["disagree"] and stats["indel"] and stats["no_overlap"] < stats["pairs"] / 10, stats
    assert fitted.provenance["mode"] == "pe-overlap"
    assert report.failures() == [], (report.scalars, report.curves["rate_by_q_true"], report.curves["rate_by_q_fit"])


def test_long_inserts_give_no_overlap_evidence_and_say_so() -> None:
    """Inserts longer than twice the read length place nowhere, so head E is not identified (phase 7 grid)."""
    truth = recovery.scale_error_rate(recovery.example_spec(), 0.2)
    rng = np.random.default_rng(0)
    templates, mates = recovery.paired_draw(40, 200, 10, genome_length=5000)(rng, 80)
    reads = gen.generate(truth, templates, mates, rng)
    with pytest.raises(ValueError, match="no pair overlapped"):
        pe_overlap.mode(templates, reads, mates, truth)


def test_attribution_weighs_each_explanation_by_head_q() -> None:
    """A low Q on G is likelier if the template base there is G, so a disagreement where mate 1 reads G at Q2 and
    mate 2 reads A at Q37 leans further to mate 1 with head Q than with head E alone (the reversed-Q truth's
    failure, plan phase 8)."""
    truth = recovery.example_spec()
    rng = np.random.default_rng(0)
    template = "".join(rng.choice(list("ACGT"), size=40))
    template = template[:20] + "A" + template[21:]
    r1, q1 = template[:20] + "G" + template[21:], [37] * 20 + [2] + [37] * 19
    ev = pe_overlap.collect([(r1, q1, gen._revcomp(template), [37] * 40)], 1, (2, 2))
    ev.alphabet = set(truth.quality_alphabet)
    groups = list(ev.disputed.items())
    assert len(groups) == 1
    alone = pe_overlap._posterior(ev, groups, truth.error_head)
    with_q = pe_overlap._posterior(ev, groups, truth.error_head, truth.quality_head)
    assert with_q[0] > alone[0]


@pytest.mark.parametrize(
    ("name", "want"),
    [
        ("@HWI-C00124:212:C9T75ANXX:4:1101:3159:2229 1:N:0:", "4_1101"),  # CASAVA 1.8+, seven fields
        ("HS24_135:4:1101:10001:74476/1", "4_1101"),  # the older five-field layout
        ("@M05310:18:000000000-CB339:1:1101:9946:1807", "1_1101"),
        ("@SRR24523812.1 1/1", None),  # SRA-renamed: the tile is gone
        ("@ERR14365452.1/1", None),
        ("", None),
    ],
)
def test_parse_tile_reads_both_illumina_layouts(name: str, want: str | None) -> None:
    assert pe_overlap.parse_tile(name) == want


def test_read_pairs_refuses_names_without_a_tile(tmp_path: Path) -> None:
    for mate in (1, 2):
        (tmp_path / f"r{mate}.fastq").write_text(f"@SRR1.1 {mate}/1\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError, match="carries no flowcell tile"):
        list(pe_overlap.read_pairs(tmp_path / "r1.fastq", tmp_path / "r2.fastq", tiles=True))


def test_tile_recovers_a_known_per_tile_rate_ratio() -> None:
    """Phase 4's `Tile`: the one read-level covariate that is observed, so an overlap can fit it.

    Pairs drawn with a 3x error-odds difference between two tiles. The fitted head must reproduce each tile's
    rate on the rows it saw, and a head without the term must still reproduce the pooled rate - otherwise the
    term is absorbing level rather than being identified.
    """
    truth_head = recovery.scale_error_rate(recovery.example_spec(), 0.1)
    weights = np.zeros((2, 10))
    weights[1, 1:] = np.log(3.0)
    truth = replace(
        truth_head,
        error_head=(*truth_head.error_head, Component("Tile(2)", {"weights": weights, "prior": np.ones(2)})),
    )
    rng = np.random.default_rng(3)
    templates, mates = recovery.paired_draw(40, 45, 8)(rng, 800)
    tiles = np.array([i // 2 % 2 for i in range(len(templates))])  # both mates of a pair share a tile
    reads = gen.generate(truth, templates, mates, rng, tiles=tiles)
    pairs = [
        (a.sequence, [ord(c) - 33 for c in a.quality], b.sequence, [ord(c) - 33 for c in b.quality], f"1_{t}")
        for a, b, t in zip(reads[::2], reads[1::2], tiles[::2], strict=True)
    ]
    ev = pe_overlap.collect(pairs, 1, (2, 2))
    assert ev.tiles == {"1_0": 0, "1_1": 1}
    quality_tokens = ["QualityMarkov(1)", "Position(4)", "Mate"]
    fitted = pe_overlap.fit(ev, ["QualityWindow(1)", "Context(1,1)", "Tile(2)"], quality_tokens, iterations=3)
    flat = pe_overlap.fit(ev, ["QualityWindow(1)", "Context(1,1)"], quality_tokens, iterations=3)

    table = pe_overlap.table(ev, fitted.error_head, fitted.quality_head)
    n = np.array(list(table.counts.values()), float)
    tile_at = table.fields.index("tile")

    def rate(spec: spec_io.ErrorModelSpec, sel: "np.ndarray[Any, Any] | None" = None) -> float:
        p = error.probabilities(list(spec.error_head), spec.quality_alphabet, table)
        rows = np.ones(len(n), bool) if sel is None else sel
        return float(n[rows] @ (1 - p[rows, 0]) / n[rows].sum())

    for tile in (0, 1):
        sel = np.array([k[tile_at] == tile for k in table.counts])
        # Against the truth on the same rows, as `compare.components` does: the overlap's exposure is its own.
        # Tolerance over seeds 3/7/11: 0.89-0.97. Fitting no tile effect at all would read ~0.46 here.
        assert rate(fitted, sel) / rate(truth, sel) == pytest.approx(1.0, abs=0.2), tile
    hot = np.array([k[tile_at] == 1 for k in table.counts])
    ratio = rate(fitted, hot) / rate(fitted, ~hot)
    assert ratio / (rate(truth, hot) / rate(truth, ~hot)) == pytest.approx(1.0, abs=0.15), ratio
    # Dropping the term leaves the pooled rate where it was, so it is identified, not soaking up the level.
    assert rate(flat) == pytest.approx(rate(fitted), rel=0.02)


def test_cli_tiles_fits_one_level_per_tile(tmp_path: Path) -> None:
    model = recovery.example_spec()
    rng = np.random.default_rng(1)
    templates, mates = recovery.paired_draw(40, 45, 8, genome_length=5000)(rng, 400)
    reads = gen.generate(model, templates, mates, rng, error_rate_scale=0.2)
    for mate in (1, 2):
        records = (
            # Real MiSeq names, two tiles: the CLI must read the tile out of them.
            f"@M05310:18:000000000-CB339:1:{1101 + i % 2}:{i}:{i} {mate}:N:0:\n{r.sequence}\n+\n{r.quality}\n"
            for i, r in enumerate(reads[mate - 1 :: 2])
        )
        (tmp_path / f"r{mate}.fastq").write_text("".join(records))
    out = tmp_path / "spec"
    args = [str(tmp_path / "r1.fastq"), str(tmp_path / "r2.fastq"), "--output", str(out), "--iterations", "2",
            "--tiles"]  # fmt: skip
    assert pe_overlap.main(args) == 0
    spec = spec_io.load(out)
    assert spec.provenance["tiles"] == ["1_1101", "1_1102"]
    assert next(c for c in spec.error_head if c.name == "Tile").token == "Tile(2)"

    # A Tile token passed by hand as well is a mistake, and so is a run whose names carry no tile.
    assert pytest.raises(SystemExit, pe_overlap.main, [*args, "--error-tokens", "QualityWindow(1)", "Tile(2)"])
    plain = tmp_path / "plain.fastq"
    plain.write_text("@SRR1.1 1/1\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError, match="carries no flowcell tile"):
        pe_overlap.main([str(plain), str(plain), "--output", str(tmp_path / "s2"), "--tiles"])
