import re
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from sequencing_error_model import generate as gen
from sequencing_error_model.fit import error, quality
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.spec import Component, ErrorModelSpec

ALPHABET = (2, 12, 23, 37)
K = len(error.CATEGORIES)


def spec() -> ErrorModelSpec:
    errors = np.r_[0, np.ones(K - 1)]
    lags = np.zeros((1, 5, 4))
    lags[0, :4] = 2.5 * np.eye(4)
    q_context = np.zeros((3, 5, 4))
    q_context[1, 2] = [1.5, 0.5, -0.5, -1.5]  # low Q on G
    window = np.zeros((3, 5, K))
    window[1, :4] = np.outer([2.0, 1.0, 0.0, -1.0], errors)
    window[2, :4] = np.outer([0.7, 0.3, 0.0, -0.3], errors)
    e_context = np.zeros((3, 5, K))
    e_context[0, 2, 1:5] = 1.2  # substitutions after G
    runs = np.zeros((5, K))
    runs[2:, 5:] = 1.5  # indels in runs of 3+
    return ErrorModelSpec(
        ALPHABET,
        {"mode": "reference"},
        (
            Component("QualityMarkov(1)", {"bias": np.array([-1.0, 0.0, 0.5, 1.5]), "lags": lags}),
            Component(
                "Position(3)",
                {
                    "knots": np.linspace(0, np.log(40), 3),
                    "start": np.zeros((3, 4)),
                    "end": np.array([[1.5, 0.5, -0.5, -1.5], [0] * 4, [0] * 4]),
                },
            ),
            Component("Mate", {"weights": np.array([[0, 0, 0, 0.5], [0.5, 0, 0, -0.5]])}),
            Component("Context(1,1)", {"weights": q_context}),
        ),
        (
            Component(
                "QualityWindow(1)",
                {"bias": np.array([0, -4.0, -4.0, -4.0, -4.0, -6.5, -6.5, -6.5, -6.5, -6.0]), "window": window},
            ),
            Component("Context(1,1)", {"weights": e_context}),
            Component("Homopolymer", {"weights": runs}, meta={"flank": [2, 2]}),
            Component("Mate", {"weights": np.outer([0, 0.3], errors)}),
        ),
    )


def templates(n: int, seed: int) -> tuple[list[str], list[int]]:
    rng = np.random.default_rng(seed)
    return ["".join(rng.choice(list("ACGT"), size=rng.integers(30, 41))) for _ in range(n)], [
        int(m) for m in rng.integers(1, 3, size=n)
    ]


def test_reads_are_alignment_consistent() -> None:
    model = spec()
    seqs, mates = templates(300, seed=1)
    seqs[0] = "ACGTNNNNACGT"
    reads = gen.generate(model, seqs, mates, np.random.default_rng(2), max_ins_run=3, error_rate_scale=8)
    ops: Counter[str] = Counter()
    for template, read in zip(seqs, reads, strict=True):
        n = Counter[str]()
        for count, op in re.findall(r"(\d+)([MID])", read.cigar):
            n[op] += int(count)
        ops += n
        assert n["M"] + n["D"] == len(template) == len(read.q_track)
        assert n["M"] + n["I"] == len(read.sequence) == len(read.quality)
        assert {ord(c) - 33 for c in read.quality} <= set(ALPHABET)
        _, tq = gen.align(template, read)
        deleted = [t for t, q in enumerate(tq) if q is None]
        assert all(q == read.q_track[t] for t, q in enumerate(tq) if q is not None)
        assert len(deleted) == n["D"]
    assert "NNNN" in reads[0].sequence  # non-ACGT bases pass through as matches, with no insertions between
    assert ops["I"] and ops["D"], ops


def centred(w: np.ndarray) -> np.ndarray:
    w = w - w.mean(axis=-1, keepdims=True)
    return (w - w.mean(axis=-2, keepdims=True)).ravel()


def test_refit_from_generated_reads_recovers_spec() -> None:
    model = spec()
    seqs, mates = templates(3000, seed=3)
    reads = gen.generate(model, seqs, mates, np.random.default_rng(4))
    table = gen.observations(zip(seqs, reads, mates, strict=True), flank=(2, 2), m=1)

    fitted_e = error.fit(table, [c.token for c in model.error_head], ALPHABET)
    n = np.array(list(table.counts.values()))
    q = np.array([k[0] for k in table.counts])
    err = {
        name: n * (1 - error.probabilities(c, ALPHABET, table)[:, 0])
        for name, c in (("true", model.error_head), ("fit", fitted_e))
    }
    assert abs(err["fit"].sum() / err["true"].sum() - 1) < 0.05
    by_q = [err["fit"][q == a].sum() / err["true"][q == a].sum() for a in ALPHABET]
    assert np.all(np.abs(np.array(by_q) - 1) < 0.1), by_q

    oi = table.fields.index("op")
    final: Counter[Key] = Counter()
    for key, c in table.counts.items():
        if key[oi] == "=" or re.fullmatch(r"[ACGT]>[ACGT]", str(key[oi])):
            final[key] += c
    base = CountTable("generator", table.fields, "base", True, final, table.meta)
    q_table = base.marginal("q-1", "pos_start", "pos_end", "mate", "context", "q")
    fitted_q = quality.fit(q_table, [c.token for c in model.quality_head], ALPHABET)
    for want, got, param, sl in (
        (model.quality_head[0], fitted_q[0], "lags", np.s_[:, :4]),
        (model.quality_head[3], fitted_q[3], "weights", np.s_[:]),
    ):
        r = np.corrcoef(centred(want.params[param][sl]), centred(got.params[param][sl]))[0, 1]
        assert r > 0.9, (want.token, r)


def test_cli_matches_contract(tmp_path: Path) -> None:
    model = spec()
    model.save(tmp_path / "spec")
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">0/1\nACGTACGTACGTACGTACGT\n>1/2\nGGGCCCAAATTTGGGCCC\n>2/1\n\n")
    outputs: list[Any] = []
    for i in range(2):
        out = tmp_path / f"out{i}.fastq"
        args = ["--model", str(tmp_path / "spec"), "--input", str(fasta), "--output", str(out), "--paired"]
        assert gen.main([*args, "--seed", "7"]) == 0
        outputs.append(out.read_text())
    assert outputs[0] == outputs[1]
    lines = outputs[0].splitlines()
    assert len(lines) == 12 and [lines[i].split()[0] for i in (0, 4, 8)] == ["@0/1", "@1/2", "@2/1"]
    assert all(re.fullmatch(r"@\S+ cigar:(\d+[MID])*", lines[i]) for i in (0, 4, 8))
