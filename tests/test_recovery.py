import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from sequencing_error_model import generate as gen
from sequencing_error_model import recovery


def _distance(a: str, b: str) -> int:
    row = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        prev, row[0] = row[0], i
        for j, cb in enumerate(b, 1):
            prev, row[j] = row[j], min(row[j] + 1, row[j - 1] + 1, prev + (ca != cb))
    return row[-1]


def test_realign_is_optimal_and_left_aligned() -> None:
    def read(seq: str, cigar: str, quality: str | None = None) -> gen.Read:
        return gen.Read(seq, quality or "I" * len(seq), cigar, np.zeros(0, np.int64))

    scores = recovery.ALIGNER_SCORES["minibwa"]
    assert gen.realign("ACGGGGT", read("ACGGGT", "6M1D"), scores)[0].cigar == "2M1D4M"
    assert gen.realign("ACGGGT", read("ACGGGGT", "6M1I"), scores)[0].cigar == "2M1I4M"
    assert gen.realign("AACCGGTTAC", read("AACTTAC", "7M3D"), scores)[0].cigar == "3M3D4M"  # one gap, not three
    clipped = gen.realign("ACGT", read("ACGTA", "4M1I", "ABCDE"), (0, 1, 0, 1))[0]
    assert (clipped.cigar, clipped.sequence, clipped.quality, clipped.clipped) == ("4M", "ACGT", "ABCD", ("", "E"))
    # Free template ends: the read placed inside a longer window.
    placed, start, end = gen.realign("TTTTACGTACGGGG", read("ACGTAC", "6M8D"), scores, free_template_ends=True)
    assert (placed.cigar, start, end) == ("6M", 4, 10)

    # Unit scores make the optimum the minimum edit count, and the band never cuts it off.
    rng = np.random.default_rng(0)
    templates = ["".join(rng.choice(list("ACGT"), size=int(rng.integers(30, 41)))) for _ in range(200)]
    reads = gen.generate(recovery.example_spec(), templates, [1] * len(templates), rng)
    hidden = 0
    for template, true in zip(templates, reads, strict=True):
        best = gen.realign(template, true, (0, 1, 0, 1))[0]
        edits = gen._edits(template, best) + len(best.clipped[1])
        assert edits == _distance(template, true.sequence)
        hidden += gen._edits(template, true) - edits
    assert hidden > 0  # the generator's CIGARs carry edits no alignment recovers


def test_observable_left_aligns_on_the_reference() -> None:
    # A reverse-strand read with one G of the GGG run deleted: the gap sits at the run's first base on the
    # forward strand (6M1D7M there), so at its last base in read orientation.
    window, forward = "GATCACGGGTCATG", "GATCACGGTCATG"
    read = gen.Read(gen._revcomp(forward), "ABCDEFGHIJKLM", "", np.zeros(0, np.int64))
    best, start, end = recovery.realign_observable(window, read, True, recovery.ALIGNER_SCORES["minibwa"])
    assert (best.cigar, best.sequence, best.quality, start, end) == ("7M1D6M", read.sequence, read.quality, 0, 14)
    # An insertion before the window's first base (1I14M forward) ends the read: clipped, as `sources.bam` does.
    read = gen.Read(gen._revcomp("T" + window), "ABCDEFGHIJKLMNO", "", np.zeros(0, np.int64))
    best = recovery.realign_observable(window, read, True, recovery.ALIGNER_SCORES["minibwa"])[0]
    assert (best.cigar, best.sequence, best.clipped) == ("14M", read.sequence[:-1], ("", "O"))


def test_cigar_mode_recovers_example_spec() -> None:
    truth = recovery.example_spec()
    fitted, report = recovery.recover(truth, n_reads=3000, seed=3)
    assert [c.token for c in fitted.error_head] == [c.token for c in truth.error_head]
    assert report.failures() == [], report.scalars
    assert abs(report.scalars["q_lag1_fit"] - report.scalars["q_lag1_true"]) < 0.05, report.scalars


def test_indel_components_isolate_homopolymer_and_context() -> None:
    truth = recovery.example_spec()
    rng = np.random.default_rng(0)
    templates = ["".join(rng.choice(list("ACGT"), size=40)) for _ in range(2000)]
    mates = [1] * len(templates)
    table = recovery._tuples(truth, templates, gen.generate(truth, templates, mates, rng), mates)
    assert set(recovery.indel_components(truth, truth, table).values()) == {1.0}

    head = list(truth.error_head)
    runs, ctx = (head[i].params["weights"].copy() for i in (2, 1))
    runs[2:, 5:] = 0.0  # no homopolymer effect on indels
    ctx[0, 0, 9] = 1.0  # deletions after A
    head[2] = replace(head[2], params={"weights": runs})
    head[1] = replace(head[1], params={"weights": ctx})
    moved = recovery.indel_components(truth, replace(truth, error_head=tuple(head)), table)
    assert moved["D run 3+"] < 0.5 and moved["I run 3+"] < 0.5 and moved["I run 1"] > 0.9
    assert moved["D -1A"] > 1.5 * moved["D -1C"] and 0.9 < moved["I -1A"] / moved["I -1C"] < 1.1


def test_cli_writes_report(tmp_path: Path) -> None:
    out = tmp_path / "report.json"
    code = recovery.main(["--reads", "400", "--seed", "1", "--output", str(out)])
    doc = json.loads(out.read_text())
    assert code == (1 if doc["failures"] else 0)
    assert {"op_tv", "rate_ratio", "q_position_tv_max"} <= doc["scalars"].keys()
    assert len(doc["curves"]["rate_by_q_true"]) == len(doc["curves"]["reported_q"])
