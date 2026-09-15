"""Head Q: P(q_t | q_{t-1..t-m}, x_{t-L..t+R}, position, mate) as a log-linear softmax over the
quality alphabet (plan §5.2, §5.4).

Every component adds a K-vector of logits (K = alphabet size) for each feature slot a row
activates, so the model is logits = X @ theta for a sparse design X. Fitting and sampling
build X the same way, so the generator samples exactly the fitted model.

Components and their params (trailing axis K):

- `QualityMarkov(m)` (required): `bias` [K], `lags` [m, K+1, K]; row K is "before the read start".
  Fields `q-1`..`q-m`, with None for lags before the read start.
- `Position(n)`: `knots` [n] over log position, `start` and `end` [n, K]: linear-spline effects of
  the distance from the read start and end. Fields `pos_start`, `pos_end`.
- `Mate`: `weights` [2, K]. Field `mate`.
- `Context(L,R)`: `weights` [L+R+1, 5, K] over bases A, C, G, T and anything else ("." beyond a
  read end, N). Field `context`, whose table flank must cover (L, R).

Input tables count (covariates, q). Per-observation tuples condition on true bases; FASTQ counts
on observed bases (§6.2). Qualities are the modelled output here, never error evidence.
"""

import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from math import prod
from typing import Any

import numpy as np
import numpy.typing as npt
from scipy import optimize, sparse, special

from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.spec import Component

Array = npt.NDArray[Any]
COMPONENTS = ("QualityMarkov", "Position", "Mate", "Context")
_BASE = np.full(256, 4, dtype=np.int64)
_BASE[np.frombuffer(b"ACGT", np.uint8)] = np.arange(4)


@dataclass
class _Rows:
    lags: Array  # [R, m] alphabet indices, K before the read start
    pos_start: Array
    pos_end: Array
    mate: Array
    context: Array  # [R, L0+R0+1] base indices
    flank: tuple[int, int]


def _shapes(c: Component, k: int) -> dict[str, tuple[int, ...]]:
    """Fitted param shapes, in slot order, without the trailing K axis."""
    if c.name == "QualityMarkov":
        return {"bias": (), "lags": (c.args[0], k + 1)}
    if c.name == "Position":
        return {"start": (c.args[0],), "end": (c.args[0],)}
    if c.name == "Mate":
        return {"weights": (2,)}
    left, right = c.args
    return {"weights": (left + right + 1, 5)}


def _hat(pos: Array, knots: Array) -> Array:
    return np.stack([np.interp(np.log(pos), knots, e) for e in np.eye(len(knots))], axis=1)


def _design(c: Component, rows: _Rows, k: int) -> tuple[Array, Array]:
    """(slot index, value) per row for one component; both [R, J]."""
    n_rows = len(rows.mate)
    if c.name == "QualityMarkov":
        m = c.args[0]
        idx = np.hstack([np.zeros((n_rows, 1), np.int64), 1 + np.arange(m) * (k + 1) + rows.lags[:, :m]])
        return idx, np.ones(idx.shape)
    if c.name == "Position":
        knots = c.params["knots"]
        vals = np.hstack([_hat(rows.pos_start, knots), _hat(rows.pos_end, knots)])
        return np.broadcast_to(np.arange(vals.shape[1]), vals.shape), vals
    if c.name == "Mate":
        return (rows.mate - 1)[:, None], np.ones((n_rows, 1))
    left, right = c.args
    window = rows.context[:, rows.flank[0] - left : rows.flank[0] + right + 1]
    return np.arange(window.shape[1]) * 5 + window, np.ones(window.shape)


def _matrix(components: Sequence[Component], rows: _Rows, k: int) -> sparse.csr_matrix:
    blocks, offset = [], 0
    for c in components:
        idx, vals = _design(c, rows, k)
        r = np.broadcast_to(np.arange(len(idx))[:, None], idx.shape)
        blocks.append((vals.ravel(), r.ravel(), (offset + idx).ravel()))
        offset += sum(prod(s) for s in _shapes(c, k).values())
    vals, r, col = (np.concatenate(x) for x in zip(*blocks, strict=True))
    return sparse.csr_matrix((vals, (r, col)), shape=(len(rows.mate), offset))


def _theta(components: Sequence[Component], k: int) -> Array:
    return np.vstack([c.params[p].reshape(-1, k) for c in components for p in _shapes(c, k)])


def _check(components: Sequence[Component]) -> None:
    names = [c.name for c in components]
    if not names or names[0] != "QualityMarkov" or not set(names) <= set(COMPONENTS) or len(set(names)) != len(names):
        raise ValueError(f"head Q needs QualityMarkov(m) first, then distinct {COMPONENTS[1:]}, got {names}")
    arity = {"QualityMarkov": 1, "Position": 1, "Mate": 0, "Context": 2}
    for c in components:
        if len(c.args) != arity[c.name] or (c.name == "Position" and c.args[0] < 2):
            raise ValueError(f"bad arguments in {c.token!r}")


def fit(table: CountTable, tokens: Sequence[str], alphabet: Sequence[int], *, l2: float = 1.0) -> tuple[Component, ...]:
    """Fit head Q components (in `tokens` order) by penalised maximum likelihood.

    `l2` is a Gaussian prior precision on every weight; it also pins the softmax gauge.
    """
    components = [Component(t) for t in tokens]
    _check(components)
    k, qi = len(alphabet), table.fields.index("q")
    q_index: dict[Any, int] = {q: i for i, q in enumerate(alphabet)}
    keys: dict[Key, int] = {}
    r_idx, q_idx, n = [], [], []
    for key, cnt in table.counts.items():
        if key[qi] not in q_index:
            raise ValueError(f"{table.source}: Q {key[qi]!r} is not in the alphabet {tuple(alphabet)}")
        r_idx.append(keys.setdefault(key[:qi] + key[qi + 1 :], len(keys)))
        q_idx.append(q_index[key[qi]])
        n.append(cnt)
    counts = np.zeros((len(keys), k))
    np.add.at(counts, (r_idx, q_idx), n)

    m = components[0].args[0]
    need = {f"q-{i}" for i in range(1, m + 1)}
    for c in components:
        need |= {"Position": {"pos_start", "pos_end"}, "Mate": {"mate"}, "Context": {"context"}}.get(c.name, set())
    if missing := sorted(need - set(table.fields)):
        raise ValueError(f"{table.source}: head Q components {tuple(tokens)} need fields {missing}")
    flank: tuple[int, int] = tuple(table.meta.get("flank", (0, 0)))
    for c in components:
        if c.name == "Context" and (c.args[0] > flank[0] or c.args[1] > flank[1]):
            raise ValueError(f"{table.source}: context flank {flank} does not cover {c.token}")

    n_rows = len(keys)
    cols = dict(zip([f for f in table.fields if f != "q"], zip(*keys, strict=True), strict=True))
    lag_index: dict[Any, int] = q_index | {None: k}
    ones = np.ones(n_rows, np.int64)
    context = "".join(cols.get("context", "." * n_rows))
    rows = _Rows(
        np.array([[lag_index[q] for q in cols[f"q-{i}"]] for i in range(1, m + 1)], np.int64).reshape(m, n_rows).T,
        np.asarray(cols.get("pos_start", ones)),
        np.asarray(cols.get("pos_end", ones)),
        np.asarray(cols.get("mate", ones)),
        _BASE[np.frombuffer(context.encode(), np.uint8)].reshape(n_rows, -1),
        flank,
    )
    if any(c.name == "Position" for c in components):
        top = np.log(max(rows.pos_start.max(), rows.pos_end.max(), 2))
        components = [
            Component(c.token, {"knots": np.linspace(0, top, c.args[0])}) if c.name == "Position" else c
            for c in components
        ]

    x = _matrix(components, rows, k)
    totals = counts.sum(axis=1)

    def objective(flat: Array) -> tuple[float, Array]:
        theta = flat.reshape(-1, k)
        logits = x @ theta
        logz = special.logsumexp(logits, axis=1)
        nll = totals @ logz - np.sum(counts * logits) + 0.5 * l2 * np.sum(theta**2)
        grad = x.T @ (totals[:, None] * np.exp(logits - logz[:, None]) - counts) + l2 * theta
        return float(nll), grad.ravel()

    res = optimize.minimize(objective, np.zeros(x.shape[1] * k), jac=True, method="L-BFGS-B")
    if not res.success:
        warnings.warn(f"head Q fit did not converge: {res.message}", RuntimeWarning, stacklevel=2)

    out, theta, i = [], res.x.reshape(-1, k), 0
    for c in components:
        params = dict(c.params)
        for name, shape in _shapes(c, k).items():
            params[name] = theta[i : i + prod(shape)].reshape(*shape, k)
            i += prod(shape)
        out.append(Component(c.token, params))
    return tuple(out)


def sample(
    components: Sequence[Component],
    alphabet: Sequence[int],
    reads: Sequence[str],
    mates: Sequence[int],
    rng: np.random.Generator,
) -> list[Array]:
    """Sample a Q track per read, position by position, vectorised across reads."""
    _check(components)
    k, alpha = len(alphabet), np.asarray(alphabet)
    theta = _theta(components, k)
    m = components[0].args[0]
    flank = next(((c.args[0], c.args[1]) for c in components if c.name == "Context"), (0, 0))
    lengths, mate = np.array([len(s) for s in reads]), np.asarray(mates)
    top = int(lengths.max(initial=0))
    text = "".join("." * flank[0] + s.ljust(top, ".") + "." * flank[1] for s in reads)
    bases = _BASE[np.frombuffer(text.encode(), np.uint8)].reshape(len(reads), -1)
    q = np.full((len(reads), top), k)
    for t in range(top):
        act = np.flatnonzero(lengths > t)
        lags = np.stack([q[act, t - i] if t >= i else np.full(len(act), k) for i in range(1, m + 1)], axis=1)
        rows = _Rows(
            lags.reshape(len(act), m),
            np.full(len(act), t + 1),
            lengths[act] - t,
            mate[act],
            bases[act, t : t + sum(flank) + 1],
            flank,
        )
        p = special.softmax(_matrix(components, rows, k) @ theta, axis=1)
        q[act, t] = np.minimum((rng.random((len(act), 1)) > p.cumsum(axis=1)).sum(axis=1), k - 1)
    return [alpha[q[i, :n]] for i, n in enumerate(lengths)]
