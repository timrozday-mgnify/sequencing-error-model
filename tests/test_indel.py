import re
from dataclasses import replace

import numpy as np

from sequencing_error_model import generate as gen
from sequencing_error_model import recovery
from sequencing_error_model.fit import indel
from sequencing_error_model.observations import CountTable
from sequencing_error_model.spec import Component, ErrorModelSpec


def _truth() -> ErrorModelSpec:
    """The example spec with 4x the indel rate and longer indels in homopolymer runs of 3+."""
    base = recovery.example_spec()
    run = np.zeros((2, indel.MAX_RUN, 4))
    run[:, 2:, 1:] = [1.0, 1.5, 2.0]
    lengths = Component(
        "IndelLength(4)",
        {"bias": np.tile([0.0, -1.5, -3.0, -4.0], (2, 1)), "run": run, "q": np.zeros((2, 4, 4))},
        meta={"max_run": indel.MAX_RUN},
    )
    head0 = base.error_head[0]
    bias = head0.params["bias"].copy()
    bias[5:] += np.log(4)
    return replace(
        base, error_head=(replace(head0, params={**head0.params, "bias": bias}), *base.error_head[1:], lengths)
    )


def _columns(table: CountTable, alphabet: tuple[int, ...]) -> tuple[np.ndarray, ...]:
    keys = list(table.counts)
    kind = np.array([indel.KINDS.index(str(k[0])) for k in keys])
    length, run = (np.array([int(str(k[i])) for k in keys]) for i in (1, 2))
    q = np.searchsorted(alphabet, [int(str(k[3])) for k in keys])
    return kind, np.minimum(length, 4), run, q, np.array(list(table.counts.values()), float)


def test_generated_lengths_follow_truth() -> None:
    truth, rng = _truth(), np.random.default_rng(0)
    seqs = ["".join(rng.choice(list("ACGT"), size=int(rng.integers(60, 80)))) for _ in range(3000)]
    mates = [1] * len(seqs)
    reads = gen.generate(truth, seqs, mates, rng)
    for template, read in zip(seqs, reads, strict=True):
        ops = {op: sum(int(c) for c, o in re.findall(r"(\d+)([MID])", read.cigar) if o == op) for op in "MID"}
        assert ops["M"] + ops["D"] == len(template)
        assert ops["M"] + ops["I"] == len(read.sequence) == len(read.quality)
    kind, length, run, q, n = _columns(gen.indel_events(zip(seqs, reads, mates, strict=True)), truth.quality_alphabet)
    p = indel.probabilities_at(truth.error_head[-1], kind, run, q)
    for k in (0, 1):
        for long_run in (False, True):
            sel = (kind == k) & ((run >= 3) == long_run)
            observed = np.bincount(length[sel] - 1, n[sel], 4) / n[sel].sum()
            expected = n[sel] @ p[sel] / n[sel].sum()
            assert n[sel].sum() > 300 and 0.5 * np.abs(observed - expected).sum() < 0.05, (k, long_run)
    # Per-event head E rows: one row per indel event, so as many indel rows as events.
    table = gen.observations(zip(seqs, reads, mates, strict=True), (2, 2), 1, per_event=True)
    rows = sum(c for key, c in table.counts.items() if "-" in str(key[-1]))
    assert rows == n.sum()


def test_recovers_indel_lengths() -> None:
    truth = _truth()
    fitted, report = recovery.recover(truth, 3000, seed=1)
    assert not report.failures(), report.failures()
    assert [c.token for c in fitted.error_head] == [c.token for c in truth.error_head]
    rng = np.random.default_rng(5)
    seqs = ["".join(rng.choice(list("ACGT"), size=40)) for _ in range(2000)]
    reads = gen.generate(truth, seqs, [1] * len(seqs), rng)
    kind, _, run, q, n = _columns(
        gen.indel_events(zip(seqs, reads, [1] * len(seqs), strict=True)), truth.quality_alphabet
    )
    p_true, p_fit = (indel.probabilities_at(s.error_head[-1], kind, run, q) for s in (truth, fitted))
    assert n @ (0.5 * np.abs(p_true - p_fit).sum(axis=1)) / n.sum() < 0.05
