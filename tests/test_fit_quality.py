from collections import Counter
from collections.abc import Sequence
from typing import Any

import numpy as np
import pytest

from sequencing_error_model import select
from sequencing_error_model.fit import quality
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.spec import Component

ALPHABET = (2, 12, 23, 37)
FIELDS = ("q-1", "pos_start", "pos_end", "mate", "context", "q")


def truth() -> tuple[Component, ...]:
    rng = np.random.default_rng(1)
    lags = rng.normal(size=(1, 5, 4))
    lags[0, :4] += 2.5 * np.eye(4)  # sticky qualities
    context = rng.normal(scale=0.3, size=(3, 5, 4))
    context[1, 2] = [1.5, 0.5, -0.5, -1.5]  # low Q on centre G
    return (
        Component("QualityMarkov(1)", {"bias": np.array([-1.0, 0.0, 0.5, 1.5]), "lags": lags}),
        Component(
            "Position(3)",
            {
                "knots": np.linspace(0, np.log(40), 3),
                "start": np.zeros((3, 4)),
                "end": np.array([[1.5, 0.5, -0.5, -1.5], [0, 0, 0, 0], [0, 0, 0, 0]]),  # end-of-read decay
            },
        ),
        Component("Mate", {"weights": np.array([[0, 0, 0, 0.5], [0.5, 0, 0, -0.5]])}),
        Component("Context(1,1)", {"weights": context}),
    )


def simulate(components: Sequence[Component], n: int, seed: int) -> tuple[list[str], list[int], list[Any]]:
    rng = np.random.default_rng(seed)
    reads = ["".join(rng.choice(list("ACGT"), size=rng.integers(30, 41))) for _ in range(n)]
    mates = list(rng.integers(1, 3, size=n))
    return reads, mates, quality.sample(components, ALPHABET, reads, mates, rng)


def tuples(reads: list[str], mates: list[int], quals: list[Any]) -> CountTable:
    counts: Counter[Key] = Counter()
    for seq, mate, q in zip(reads, mates, quals, strict=True):
        padded = f".{seq}."
        for t, qt in enumerate(q):
            counts[(int(q[t - 1]) if t else None, t + 1, len(q) - t, int(mate), padded[t : t + 3], int(qt))] += 1
    return CountTable("test", FIELDS, "base", False, counts, {"flank": (1, 1)})


def centred(w: np.ndarray) -> np.ndarray:
    w = w - w.mean(axis=-1, keepdims=True)
    return (w - w.mean(axis=-2, keepdims=True)).ravel()


def position_hist(quals: list[Any], upto: int) -> np.ndarray:
    h = np.zeros((upto, len(ALPHABET)))
    for q in quals:
        np.add.at(h, (np.arange(upto), np.searchsorted(ALPHABET, q[:upto])), 1)
    return h / h.sum(axis=1, keepdims=True)


def test_recovers_known_head() -> None:
    true = truth()
    table = tuples(*simulate(true, 4000, seed=2))
    fitted = quality.fit(table, [c.token for c in true], ALPHABET)
    assert [c.token for c in fitted] == [c.token for c in true]

    by = {c.name: c for c in fitted}
    effects: list[tuple[str, str, Any]] = [
        ("QualityMarkov", "lags", np.s_[:, :4]),
        ("Context", "weights", np.s_[:]),
        ("Mate", "weights", np.s_[:]),
    ]
    for name, param, sl in effects:
        want = next(c for c in true if c.name == name).params[param][sl]
        r = np.corrcoef(centred(want), centred(by[name].params[param][sl]))[0, 1]
        assert r > 0.9, (name, r)

    reads, mates, want = simulate(true, 20000, seed=3)
    got = quality.sample(fitted, ALPHABET, reads, mates, np.random.default_rng(4))
    tv = 0.5 * np.abs(position_hist(want, 30) - position_hist(got, 30)).sum(axis=1)
    assert tv.max() < 0.05, tv


def test_select_finds_true_components() -> None:
    true = truth()
    train, test = tuples(*simulate(true, 2000, seed=6)), tuples(*simulate(true, 1000, seed=7))
    first = ("QualityMarkov(0)", "QualityMarkov(1)")
    result = select.select(quality, train, test, ALPHABET, first, [("Position(3)",), ("Mate",), ("Context(1,1)",)])
    tokens = [c.token for c in result.components]
    assert tokens[0] == "QualityMarkov(1)" and set(tokens) == {c.token for c in true}, result.trace
    assert sum(s.accepted for s in result.trace) == len(true)


def test_sample_emits_only_alphabet() -> None:
    _, _, quals = simulate(truth(), 50, seed=5)
    assert set(np.concatenate(quals)) <= set(ALPHABET)


@pytest.mark.parametrize(
    ("tokens", "fields", "key", "meta", "match"),
    [
        (["Mate"], ("mate", "q"), (1, 2), {}, "QualityMarkov"),
        (["QualityMarkov(1)", "Strand"], ("q-1", "q"), (None, 2), {}, "QualityMarkov"),
        (["QualityMarkov(1)", "Position(1)"], ("q-1", "q"), (None, 2), {}, "bad arguments"),
        (["QualityMarkov(2)"], ("q-1", "q"), (None, 2), {}, "q-2"),
        (["QualityMarkov(0)"], ("q",), (3,), {}, "alphabet"),
        (["QualityMarkov(0)", "Position(3)"], ("pos_start", "q"), (1, 2), {}, "pos_end"),
        (["QualityMarkov(0)", "Context(2,1)"], ("context", "q"), ("ACG", 2), {"flank": (1, 1)}, "flank"),
    ],
)
def test_invalid(tokens: list[str], fields: tuple[str, ...], key: Key, meta: dict[str, Any], match: str) -> None:
    table = CountTable("test", fields, "base", False, Counter({key: 1}), meta)
    with pytest.raises(ValueError, match=match):
        quality.fit(table, tokens, ALPHABET)
