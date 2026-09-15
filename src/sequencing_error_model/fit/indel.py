"""`IndelLength(n)` (plan §5.4): P(length | kind, homopolymer run, Q) per indel event, a log-linear softmax over
lengths 1..n fitted from `generate.indel_events` tables (unit `error`, truth-bearing).

Head E decides where an indel starts, with one draw per event; this component decides its length. Lengths
above n count as n, so a spec generates at most n bases per event.

Params (trailing axis n): `bias` [2, n] per kind (I, D); `run` [2, R, n] by the template homopolymer run length,
capped at R = `meta["max_run"]`; `q` [2, A, n] by the Q at the event's template base (for a deletion, the next
read base's Q, as in `observations`). Q is a categorical feature, never an error probability.
"""

import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np
from scipy import optimize, special

from sequencing_error_model.fit.quality import Array
from sequencing_error_model.observations import CountTable
from sequencing_error_model.spec import Component

KINDS = ("I", "D")
MAX_RUN = 8  # ponytail: fixed cap; make it a token argument if long-read homopolymers need more


def split(components: Sequence[Component]) -> tuple[tuple[Component, ...], Component | None]:
    """Head E components without `IndelLength`, and `IndelLength` if present."""
    lengths = next((c for c in components if c.name == "IndelLength"), None)
    return tuple(c for c in components if c.name != "IndelLength"), lengths


def _check(token: str) -> int:
    c = Component(token)
    if c.name != "IndelLength" or len(c.args) != 1 or c.args[0] < 1:
        raise ValueError(f"expected IndelLength(n) with n >= 1, got {token!r}")
    return c.args[0]


def _logits(params: dict[str, Array], kind: Array, run: Array, q: Array) -> Array:
    return np.asarray(params["bias"][kind] + params["run"][kind, np.minimum(run, MAX_RUN) - 1] + params["q"][kind, q])


def probabilities_at(c: Component, kind: Array, run: Array, q: Array) -> Array:
    """P(length = 1..n) [R, n] from kind indices (0 I, 1 D), run lengths and alphabet indices of Q."""
    return np.asarray(special.softmax(_logits(c.params, kind, run, q), axis=1))


def _rows(table: CountTable, n: int, alphabet: Sequence[int]) -> tuple[Array, Array, Array, Array]:
    """Kind, run and Q index per distinct covariate row, and counts [R, n] by capped length."""
    fields = ("indel", "indel_length", "run", "q")
    if table.unit != "error" or not set(fields) <= set(table.fields):
        raise ValueError(f"{table.source}: IndelLength needs an error-unit table with fields {fields}")
    idx = [table.fields.index(f) for f in fields]
    q_index: dict[Any, int] = {q: i for i, q in enumerate(alphabet)}
    keys: dict[tuple[int, int, int], int] = {}
    r, col, w = [], [], []
    for key, cnt in table.counts.items():
        kind, length, run, q = (key[i] for i in idx)
        if q not in q_index:
            raise ValueError(f"{table.source}: Q {q!r} is not in the alphabet {tuple(alphabet)}")
        r.append(keys.setdefault((KINDS.index(str(kind)), int(str(run)), q_index[q]), len(keys)))
        col.append(min(int(str(length)), n) - 1)
        w.append(cnt)
    counts = np.zeros((len(keys), n))
    np.add.at(counts, (r, col), w)
    cols = np.array(list(keys), np.int64).reshape(-1, 3)
    return cols[:, 0], cols[:, 1], cols[:, 2], counts


def fit(table: CountTable, token: str, alphabet: Sequence[int], *, l2: float = 1.0) -> Component:
    """Fit `IndelLength(n)` by penalised maximum likelihood; `l2` is a Gaussian prior precision on every weight."""
    n = _check(token)
    kind, run, q, counts = _rows(table, n, alphabet)
    shapes = {"bias": (2, n), "run": (2, MAX_RUN, n), "q": (2, len(alphabet), n)}
    cuts = np.cumsum([int(np.prod(s)) for s in shapes.values()])
    runi, totals = np.minimum(run, MAX_RUN) - 1, counts.sum(axis=1)

    def unpack(flat: Array) -> dict[str, Array]:
        return {k: p.reshape(s) for (k, s), p in zip(shapes.items(), np.split(flat, cuts[:-1]), strict=True)}

    def objective(flat: Array) -> tuple[float, Array]:
        logits = _logits(unpack(flat), kind, run, q)
        logp = logits - special.logsumexp(logits, axis=1, keepdims=True)
        g = totals[:, None] * np.exp(logp) - counts
        grads = {k: np.zeros(s) for k, s in shapes.items()}
        np.add.at(grads["bias"], kind, g)
        np.add.at(grads["run"], (kind, runi), g)
        np.add.at(grads["q"], (kind, q), g)
        nll = -np.sum(counts * logp) + 0.5 * l2 * flat @ flat
        return float(nll), np.concatenate([v.ravel() for v in grads.values()]) + l2 * flat

    res = optimize.minimize(objective, np.zeros(cuts[-1]), jac=True, method="L-BFGS-B")
    if not res.success:
        warnings.warn(f"IndelLength fit did not converge: {res.message}", RuntimeWarning, stacklevel=2)
    return Component(token, unpack(res.x), meta={"max_run": MAX_RUN})


def log_likelihood(c: Component, alphabet: Sequence[int], table: CountTable) -> float:
    """Log-likelihood of an indel event table's lengths under a fitted `IndelLength`."""
    kind, run, q, counts = _rows(table, _check(c.token), alphabet)
    logits = _logits(c.params, kind, run, q)
    return float(np.sum(counts * (logits - special.logsumexp(logits, axis=1, keepdims=True))))
