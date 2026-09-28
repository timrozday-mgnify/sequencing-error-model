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


REAL_ALPHABET = (2, 12, 23, 37)


def truth_with_q() -> tuple[Component, ...]:
    """The same context effects, plus a centre-Q effect on every error category."""
    intercept, context = truth()
    window = np.zeros((1, len(REAL_ALPHABET) + 1, K))
    window[0, :4] = np.outer([2.0, 1.0, 0.0, -1.0], np.r_[0, np.ones(K - 1)])
    return Component("QualityWindow(0)", {"bias": intercept.params["bias"], "window": window}), context


def exposure_table(n_reads: int, seed: int, length: int = 60) -> CountTable:
    """A FASTQ-like (observed context, centre Q) exposure where low Q and the hard context co-occur.

    Qualities are drawn low after a G, which is where `truth` puts its extra substitutions, so a model that
    multiplies the two marginals into every cell over-counts those cells.
    """
    rng = np.random.default_rng(seed)
    counts: Counter[Key] = Counter()
    for _ in range(n_reads):
        padded = ".." + "".join(rng.choice(list("ACGT"), size=length)) + ".."
        for t in range(length):
            hard = padded[t + 1] == "G"
            q = rng.choice(REAL_ALPHABET, p=[0.30, 0.30, 0.25, 0.15] if hard else [0.05, 0.15, 0.30, 0.50])
            counts[(padded[t : t + 5], int(q))] += 1
    return CountTable("fastq_quality:context", ("context", "q"), "base", False, counts, {"flank": list(FLANK)})


def cells(exposure: CountTable) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """(count, context index, Q index, contexts) in `exposure.counts` order."""
    contexts: dict[str, int] = {}
    keys = [(str(c), int(str(q))) for c, q in exposure.counts]
    ctx = np.array([contexts.setdefault(c, len(contexts)) for c, _ in keys])
    q = np.array([REAL_ALPHABET.index(v) for _, v in keys])
    return np.array(list(exposure.counts.values()), float), ctx, q, list(contexts)


def pooled_context_head(
    true: Sequence[Component], exposure: CountTable, n: np.ndarray, ctx: np.ndarray, contexts: list[str], level: float
) -> tuple[Component, ...]:
    """What `kvmer.csv` gives: P(op | context) with Q pooled over the exposure, at `level` times the true rate.

    Pooling is what makes marginal matching necessary: the context margin has no Q in it, and the exposure's
    Q tilt per context is exactly the overlap the two margins would otherwise double-count. `level` stands in
    for the truncation bias of a kvmer-only fit (§6.2).
    """
    per_context = np.zeros((len(contexts), K))
    np.add.at(per_context, ctx, n[:, None] * error.probabilities(true, REAL_ALPHABET, exposure))
    counts: Counter[Key] = Counter()
    for c, row in enumerate(per_context):
        counts[(contexts[c], kmer.DUMMY_Q, "=")] = row[0] + (1 - level) * row[1:].sum()
        for cat, errors in enumerate(row[1:], start=1):
            if errors > 0:
                counts[(contexts[c], kmer.DUMMY_Q, kmer._op(cat, contexts[c][FLANK[0]]))] = level * errors
    table = CountTable("kvmer-like", ("context", "q", "op"), "base", True, counts, {"flank": list(FLANK)})
    return error.fit(table, ["QualityWindow(0)", "Context(2,2)"], (kmer.DUMMY_Q,), l2=1e-3)


def rate_per_cell(components: Sequence[Component], exposure: CountTable) -> np.ndarray:
    rate: np.ndarray = 1 - error.probabilities(components, REAL_ALPHABET, exposure)[:, 0]
    return rate


def group_rate(n: np.ndarray, rate: np.ndarray, key: np.ndarray) -> np.ndarray:
    return np.bincount(key, n * rate) / np.bincount(key, n)


def phred_bins(rates: np.ndarray) -> list[skiver_analyze.RateBin]:
    """`summary_phred.csv` rows for `REAL_ALPHABET`; only the rate and a non-zero exposure are read."""
    return [
        skiver_analyze.RateBin(q, q, float(r), (0.0, 0.0), 1000, 1) for q, r in zip(REAL_ALPHABET, rates, strict=True)
    ]


def test_marginal_matching_reproduces_both_skiver_marginals() -> None:
    true = truth_with_q()
    exposure = exposure_table(400, seed=7)
    n, ctx, q, contexts = cells(exposure)
    truth_rate = rate_per_cell(true, exposure)

    # skiver's two marginals under this exposure: P(error | Q) from summary_phred.csv, and a context margin
    # with Q pooled out, 18 % low as a kvmer-only fit leaves it (§6.2).
    phred = phred_bins(group_rate(n, truth_rate, q))
    head = pooled_context_head(true, exposure, n, ctx, contexts, level=0.82)
    fitted = kmer.fit_centre_q(head, phred, exposure, ["QualityWindow(0)", "Context(2,2)"], REAL_ALPHABET)

    got = rate_per_cell(fitted, exposure)
    # The Q margin is skiver's, and the level with it, despite the context head coming in 18 % low.
    np.testing.assert_allclose(group_rate(n, got, q), group_rate(n, truth_rate, q), rtol=0.01)
    assert abs((n @ got) / (n @ truth_rate) - 1) < 0.01
    # The context margin too, and with both margins right, the joint rate of every cell.
    by_context = [group_rate(n, r, ctx) for r in (truth_rate, got)]
    assert np.corrcoef(by_context[0], by_context[1])[0, 1] > 0.99, by_context

    # Raking, not multiplying the two rates into every cell: that double-counts low Q after a G. The truth here
    # is log-additive, which is the assumption marginal matching rests on and cannot itself test (§6.1).
    naive = by_context[1][ctx] * group_rate(n, got, q)[q] / ((n @ got) / n.sum())
    hard = np.array([c[1] == "G" for c in contexts])[ctx] & (q == 0)
    for label, where in (("all cells", np.ones(len(n), bool)), ("hard context at the lowest Q", hard)):
        off = [float(np.average(np.abs(r[where] / truth_rate[where] - 1), weights=n[where])) for r in (got, naive)]
        assert off[0] < 0.02 and off[1] > 0.2, (label, off)


def test_warns_about_exposure_it_cannot_match() -> None:
    """An N centre base or a Q skiver reports no rate for is dropped loudly, not silently."""
    counts: Counter[Key] = Counter({("ACGTA", 12): 100, ("ACNTA", 12): 5, ("ACGTA", 37): 3})
    exposure = CountTable("fastq_quality:context", ("context", "q"), "base", False, counts, {"flank": list(FLANK)})
    phred = [skiver_analyze.RateBin(12, 12, 0.01, (0.0, 0.0), 1000, 1)]  # nothing for Q 37
    with pytest.warns(RuntimeWarning, match="no centre-Q evidence"):
        raked = kmer.rake(truth(), phred, exposure, REAL_ALPHABET)
    assert {k[1] for k in raked.counts} == {12}
    assert sum(n for k, n in raked.counts.items() if k[2] != "=") == pytest.approx(1.0)  # 1 % of the 100 kept bases


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
