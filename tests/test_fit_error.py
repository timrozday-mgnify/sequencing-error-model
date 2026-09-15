from collections import Counter
from collections.abc import Sequence
from typing import Any

import numpy as np
import pytest

from sequencing_error_model import select
from sequencing_error_model.fit import error
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.spec import Component

ALPHABET = (2, 12, 23, 37)
FIELDS = ("q-1", "q", "q+1", "context", "pos_start", "pos_end", "mate", "strand", "gc")
K = len(error.CATEGORIES)


def truth() -> tuple[Component, ...]:
    rng = np.random.default_rng(1)
    errors = np.r_[0, np.ones(K - 1)]
    window = np.zeros((3, 5, K))
    window[1, :4] = np.outer([2.0, 1.0, 0.0, -1.0], errors)  # centre Q: more errors at low Q
    window[2, :4] = np.outer([0.7, 0.3, 0.0, -0.3], errors)  # next base's Q
    context = rng.normal(scale=0.2, size=(5, 5, K))
    context[1, 2, 1:5] += 1.2  # substitutions after G
    homopolymer = np.zeros((5, K))
    homopolymer[2:, 5:] = np.outer([0.8, 1.4, 2.0], np.ones(5))  # indels in runs
    interaction = np.zeros((3, 5, 1, K))
    interaction[1, 1, 0, 5:] = 1.5  # indels at low-Q C
    interaction[:, :4] -= interaction[:, :4].mean(axis=1, keepdims=True)
    return (
        Component(
            "QualityWindow(1)",
            {"bias": np.array([0, -4.0, -4.0, -4.0, -4.0, -6.0, -6.0, -6.0, -6.0, -5.5]), "window": window},
        ),
        Component("Context(2,2)", {"weights": context}),
        Component("Homopolymer", {"weights": homopolymer}),
        Component(
            "Position(3)",
            {"knots": np.linspace(0, np.log(40), 3), "start": np.zeros((3, K)), "end": np.outer([1.0, 0, 0], errors)},
        ),
        Component("Mate", {"weights": np.outer([0, 0.4], errors)}),
        Component("Strand", {"weights": np.outer([0, -0.2], errors)}),
        Component("GC(2)", {"knots": np.linspace(0, 100, 2), "weights": np.outer([-0.5, 0.5], errors)}),
        Component("QualityxContext(1)", {"quality": np.array([[1.0], [0.5], [0.0], [0.0]]), "context": interaction}),
    )


def exposure(n: int, seed: int) -> CountTable:
    """Base tuples from random reads, qualities and covariates, without op labels."""
    rng = np.random.default_rng(seed)
    counts: Counter[Key] = Counter()
    for _ in range(n):
        seq = "".join(rng.choice(list("ACGT"), size=rng.integers(30, 41)))
        q = [None, *rng.choice(ALPHABET, size=len(seq), p=[0.1, 0.2, 0.3, 0.4]).tolist(), None]
        mate, strand, padded = int(rng.integers(1, 3)), str(rng.choice(["+", "-"])), f"..{seq}.."
        gc = 10 * min(9, int(10 * sum(b in "GC" for b in seq) / len(seq)))
        for t in range(len(seq)):
            key = (q[t], q[t + 1], q[t + 2], padded[t : t + 5], t + 1, len(seq) - t, mate, strand, (gc, gc + 10))
            counts[key] += 1
    return CountTable("test", FIELDS, "base", False, counts, {"flank": (2, 2)})


def labelled(components: Sequence[Component], table: CountTable, seed: int) -> CountTable:
    rng = np.random.default_rng(seed)
    p = error.probabilities(components, ALPHABET, table)
    ci, counts = table.fields.index("context"), Counter[Key]()
    for (key, n), pk in zip(table.counts.items(), p, strict=True):
        base = str(key[ci])[2]
        for cat, m in zip(error.CATEGORIES, rng.multinomial(n, pk), strict=True):
            op = "=" if cat == "=" else f"{base}>-" if cat == "-" else cat if cat[0] == "-" else base + cat
            if m:
                counts[(*key, op)] += int(m)
    return CountTable("test", (*FIELDS, "op"), "base", True, counts, table.meta)


def log_odds(w: np.ndarray, symbols: int) -> np.ndarray:
    """Error-vs-match log-odds per slot, centred over symbols (their mean is absorbed by the bias)."""
    lo: np.ndarray = w[:, :symbols, 1:] - w[:, :symbols, :1]
    return lo - lo.mean(axis=1, keepdims=True)


def error_rate_by_q(components: Sequence[Component], table: CountTable) -> np.ndarray:
    n = np.array(list(table.counts.values()))
    err = n * (1 - error.probabilities(components, ALPHABET, table)[:, 0])
    q = np.array([k[1] for k in table.counts])
    return np.array([err[q == a].sum() / n[q == a].sum() for a in ALPHABET])


def test_recovers_known_head() -> None:
    true = truth()
    fitted = error.fit(labelled(true, exposure(6000, seed=2), seed=3), [c.token for c in true], ALPHABET)
    assert [c.token for c in fitted] == [c.token for c in true]
    by = {c.name: c for c in fitted}
    want = {c.name: c for c in true}

    # Dominant effects: Q window on all errors, context on substitutions (off the masked centre), runs on indels.
    effects: list[tuple[np.ndarray, np.ndarray, Any]] = [
        (log_odds(want["QualityWindow"].params["window"], 4), log_odds(by["QualityWindow"].params["window"], 4), ...),
        (
            log_odds(want["Context"].params["weights"], 4),
            log_odds(by["Context"].params["weights"], 4),
            np.s_[[0, 1, 3, 4], :, :4],
        ),
        (
            want["Homopolymer"].params["weights"][:, 5:] - want["Homopolymer"].params["weights"][:, :1],
            by["Homopolymer"].params["weights"][:, 5:] - by["Homopolymer"].params["weights"][:, :1],
            np.s_[:4],
        ),
    ]
    for w, got, sel in effects:
        r = np.corrcoef((w - w.mean(axis=0))[sel].ravel(), (got - got.mean(axis=0))[sel].ravel())[0, 1]
        assert r > 0.9, r

    # Q is only a feature: the truth is far from 10^(-Q/10), and the fit must follow the truth.
    held_out = exposure(3000, seed=4)
    rate_true, rate_fit = error_rate_by_q(true, held_out), error_rate_by_q(fitted, held_out)
    assert np.all(np.abs(rate_fit / rate_true - 1) < 0.1), (rate_true, rate_fit)
    n = np.array(list(held_out.counts.values()))
    marginal = [n @ (1 - error.probabilities(c, ALPHABET, held_out)[:, 0]) for c in (true, fitted)]
    assert abs(marginal[1] / marginal[0] - 1) < 0.05, marginal


def test_select_keeps_true_effects_only() -> None:
    by = {c.name: c for c in truth()}
    context = by["Context"].params["weights"][1:4]
    true = (by["QualityWindow"], Component("Context(1,1)", {"weights": context}))
    train, test = (labelled(true, exposure(n, seed=s), seed=s) for n, s in ((3000, 6), (1500, 7)))
    first = ("QualityWindow(0)", "QualityWindow(1)")
    result = select.select(error, train, test, ALPHABET, first, [("Context(1,1)",), ("Mate",), ("Strand",)])
    assert [c.token for c in result.components] == ["QualityWindow(1)", "Context(1,1)"], result.trace
    assert all(s.test_log_likelihood < 0 for s in result.trace)


def test_mask_and_categories() -> None:
    table = exposure(20, seed=5)
    p = error.probabilities(truth(), ALPHABET, table)
    centre = np.array(["ACGT".index(str(k[3])[2]) for k in table.counts])
    assert np.allclose(p.sum(axis=1), 1) and np.all(p[np.arange(len(p)), 1 + centre] == 0)


@pytest.mark.parametrize(
    ("tokens", "fields", "key", "meta", "match"),
    [
        (["Mate"], ("context", "mate", "op"), ("A", 1, "="), {}, "QualityWindow"),
        (["QualityWindow(0)", "QualityMarkov(1)"], ("context", "q", "op"), ("A", 2, "="), {}, "QualityWindow"),
        (["QualityWindow(0)", "GC(1)"], ("context", "q", "gc", "op"), ("A", 2, (0, 10), "="), {}, "bad arguments"),
        (["QualityWindow(1)"], ("context", "q", "op"), ("A", 2, "="), {}, "q-1"),
        (["QualityWindow(0)"], ("context", "q", "op"), ("A", 3, "="), {}, "alphabet"),
        (["QualityWindow(0)"], ("context", "q", "op"), ("A", 2, "C>G"), {}, "true base"),
        (["QualityWindow(0)"], ("context", "q", "op"), ("A", 2, "!"), {}, "exact-position"),
        (["QualityWindow(0)"], ("context", "q", "op"), ("N", 2, "="), {}, "centre base"),
        (["QualityWindow(0)", "Context(2,1)"], ("context", "q", "op"), ("ACG", 2, "="), {"flank": (1, 1)}, "flank"),
        (["QualityWindow(0)", "QualityxContext(1)"], ("context", "q", "op"), ("A", 2, "="), {}, "flank"),
    ],
)
def test_invalid(tokens: list[str], fields: tuple[str, ...], key: Key, meta: dict[str, Any], match: str) -> None:
    table = CountTable("test", fields, "base", True, Counter({key: 1}), meta)
    with pytest.raises(ValueError, match=match):
        error.fit(table, tokens, ALPHABET)


def test_needs_truth_labels() -> None:
    table = CountTable("test", ("context", "q"), "base", False, Counter({("A", 2): 1}))
    with pytest.raises(ValueError, match="truth-bearing"):
        error.fit(table, ["QualityWindow(0)"], ALPHABET)
