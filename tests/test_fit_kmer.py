"""Synthetic recovery for the latent-position head E fit from `kvmer.csv` (plan phase 6).

Loci are simulated straight from a known head E: every value position draws a category from the model, and
the locus counts keep only the observations that matched throughout or carried exactly one edit, as skiver
does. An insertion is counted as one edit and leaves its position matching, which is the exposure
`fit.kmer` assumes; the generator's redraw after an insertion is a second-order correction here.
"""

from collections import Counter
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pytest

from sequencing_error_model.fit import error, kmer
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.sources import skiver_analyze
from sequencing_error_model.spec import Component

K = len(error.CATEGORIES)
ALPHABET = (kmer.DUMMY_Q,)
FLANK = (2, 2)
FIXTURES = Path(__file__).parent / "fixtures" / "skiver-v0.3.2"


def truth() -> tuple[Component, ...]:
    rng = np.random.default_rng(1)
    bias = np.array([0.0, *[-5.5] * 4, *[-7.0] * 4, -7.0])
    context = rng.normal(scale=0.25, size=(5, 5, K))
    context[1, 2, 1:5] += 1.2  # substitutions after G
    context[3, 0, 5:] += 1.0  # indels before A
    return (
        Component("QualityWindow(0)", {"bias": bias, "window": np.zeros((1, 2, K))}),
        Component("Context(2,2)", {"weights": context}),
    )


def windows(locus: str, k: int, v: int, flank: tuple[int, int]) -> list[str]:
    """The context window of every value position, "." past the locus (independent of `fit.kmer`'s own)."""
    left, right = flank
    padded = locus + "." * right
    return [padded[k + t - left : k + t + right + 1] for t in range(v)]


def context_table(contexts: Sequence[str]) -> CountTable:
    counts: Counter[Key] = Counter({(c, kmer.DUMMY_Q): 1 for c in dict.fromkeys(contexts)})
    return CountTable("sim", ("context", "q"), "base", False, counts, {"flank": list(FLANK)})


def loci_probabilities(
    components: Sequence[Component], n_loci: int, seed: int, k: int = 11, v: int = 13
) -> tuple[list[str], np.ndarray]:
    """Random loci and, per locus, P(category) [v, K] at every value position under `components`."""
    rng = np.random.default_rng(seed)
    loci = ["".join(rng.choice(list("ACGT"), size=k + v)) for _ in range(n_loci)]
    per_locus = [windows(locus, k, v, FLANK) for locus in loci]
    table = context_table([c for w in per_locus for c in w])
    index = {key[0]: i for i, key in enumerate(table.counts)}
    p = error.probabilities(components, ALPHABET, table)
    return loci, np.stack([p[[index[c] for c in ctx]] for ctx in per_locus])


def truncation_ratio(p: np.ndarray) -> float:
    """Dropping 2+-edit values biases any kvmer fit down; this is the error rate it leaves, over the true one.

    Only values matching the consensus or one edit away are counted, so a position's exposure is
    P(0 or 1 edit in the value) and its errors are P(exactly 1 edit here) (plan §6.2).
    """
    pe = 1 - p[:, :, 0]
    survive = np.prod(1 - pe, axis=1)
    one = survive * (pe / (1 - pe)).sum(axis=1)
    return float(one.sum() / (p.shape[1] * (survive + one).sum()) / pe.mean())


def simulate(
    components: Sequence[Component], n_loci: int, observations: int, seed: int, k: int = 11, v: int = 13
) -> CountTable:
    rng = np.random.default_rng(seed + 1)
    loci, p = loci_probabilities(components, n_loci, seed, k, v)

    counts: Counter[Key] = Counter()
    for locus, rows in zip(loci, p, strict=True):
        cum = rows.cumsum(axis=1)
        draws = rng.random((v, observations))
        cats = np.argmax(draws[:, :, None] < cum[:, None, :], axis=-1)
        single = np.flatnonzero(np.count_nonzero(cats, axis=0) == 1)  # skiver drops 2+-edit values
        counts[(locus, "=")] += observations - int(np.count_nonzero(cats.any(axis=0)))
        for t, cat in zip(*np.nonzero(cats[:, single]), strict=True):
            name, base = error.CATEGORIES[cats[t, single[cat]]], locus[k + t]
            counts[(locus, name if name.startswith("->") else base + ">-" if name == "-" else base + name)] += 1
    return CountTable("skiver_analyze:kvmer", ("locus", "op"), "value", True, +counts, {"k": k, "v": v})


def log_odds(w: np.ndarray) -> np.ndarray:
    """Error-vs-match log-odds per offset and base, centred over bases (their mean sits in the bias)."""
    lo: np.ndarray = w[:, :4, 1:] - w[:, :4, :1]
    return lo - lo.mean(axis=1, keepdims=True)


def error_rate(components: Sequence[Component], table: CountTable) -> float:
    return float(1 - error.probabilities(components, ALPHABET, table)[:, 0].mean())


def test_recovers_context_through_the_latent_position() -> None:
    true = truth()
    table = simulate(true, n_loci=3000, observations=40, seed=2)
    fitted = kmer.fit(table, ["QualityWindow(0)", "Context(2,2)"])
    assert [c.token for c in fitted] == [c.token for c in true]

    # Context effects on substitutions, away from the masked centre offset.
    want, got = (log_odds(c[1].params["weights"]) for c in (true, fitted))
    sel = np.s_[[0, 1, 3, 4], :, :4]
    r = np.corrcoef((want - want.mean(axis=0))[sel].ravel(), (got - got.mean(axis=0))[sel].ravel())[0, 1]
    assert r > 0.9, r

    # The error rate on a held-out context exposure, biased down by exactly the mass skiver's 2+-edit
    # values carry away (the intercept is refitted against summary_phred.csv under the FASTQ exposure).
    loci, p = loci_probabilities(true, 300, seed=5)
    held_out = context_table([c for locus in loci for c in windows(locus, 11, 13, FLANK)])
    rate = [error_rate(c, held_out) for c in (true, fitted)]
    expected = truncation_ratio(p)
    assert abs(rate[1] / rate[0] / expected - 1) < 0.05, (rate, expected)
    indels = [error.probabilities(c, ALPHABET, held_out)[:, 5:].sum(axis=1).mean() for c in (true, fitted)]
    share = [i / r for i, r in zip(indels, rate, strict=True)]
    assert abs(share[1] / share[0] - 1) < 0.15, share


def test_fits_a_real_analyze_run() -> None:
    kvmer = next(
        t
        for t in skiver_analyze.tables(skiver_analyze.read_analyze(FIXTURES / "analyze"))
        if t.source.split(":")[1] == "kvmer"
    )
    fitted = kmer.fit(kvmer, ["QualityWindow(0)", "Context(1,1)", "Homopolymer"], flank=(2, 2))
    assert [c.token for c in fitted] == ["QualityWindow(0)", "Context(1,1)", "Homopolymer"]
    assert all(np.all(np.isfinite(a)) for c in fitted for a in c.params.values())


def test_rejects_tokens_and_tables_it_cannot_fit() -> None:
    table = CountTable(
        "skiver_analyze:kvmer", ("locus", "op"), "value", True, Counter({("A" * 24, "="): 1}), {"k": 11, "v": 13}
    )
    with pytest.raises(ValueError, match="QualityWindow\\(0\\)"):
        kmer.fit(table, ["QualityWindow(1)", "Context(1,1)"])
    with pytest.raises(ValueError, match="QualityWindow\\(0\\)"):
        kmer.fit(table, ["QualityWindow(0)", "Strand"])
    with pytest.raises(ValueError, match="does not fit in a 11-base key"):
        kmer.fit(table, ["QualityWindow(0)", "Context(12,0)"])
