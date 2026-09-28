"""Head E from unmodified skiver's `kvmer.csv`, where the edit position inside the value is latent (plan §6.1).

A `kvmer.csv` row is one locus: a key (k bases) plus its consensus value (v bases), the number of observed
values matching the consensus throughout, and, per error operation, how many observed values differ from the
consensus by exactly that one edit. *Where* in the value the edit sat is not reported, so each single-edit
count is a mixture over the value positions whose consensus base the operation can act on (any position for
an insertion). Positions differ in context across loci, which is what identifies `Context(L,R)`.

EM over that latent position: the E-step weighs the positions of one (locus, operation) count by the current
model's odds of that operation there, and the M-step refits head E on those positions' contexts with those
weights, plus the match exposure every observation leaves at every other position. The likelihood of a locus
with N observations that matched or carried one edit is

    prod_t P(= | ctx_t)^N  *  prod_ops ( sum_t r_op(t) )^n_op,

with r_op(t) = P(op | ctx_t) / P(= | ctx_t) for a substitution or deletion (the edit replaces that position's
match) and r_op(t) = P(op | ctx_t) for an insertion (it sits *before* a position that still matches, as head E
draws it). Contexts at the value's right end run past the locus, and those bases are "." (unknown), the same
level head E uses beyond a read end.

`kvmer.csv` carries no quality, position, strand or GC, so this fits the context part only: the returned
`QualityWindow(0)` is an intercept over a one-value dummy alphabet, and the centre-Q term comes from
`summary_phred.csv` by marginal matching under the FASTQ exposure (plan §5.6).

ponytail: the per-locus `consensus_count_up_to_v*` columns also localise each observation's *first* mismatch,
which would pin the latent position directly; they mix in the multi-edit values skiver drops from the op
counts, so they are left out. Fold them in as an E-step prior if a recovery run shows position unidentifiable.
"""

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from sequencing_error_model.fit import error
from sequencing_error_model.fit.quality import _BASE, Array
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.spec import Component

DUMMY_Q = 0  # the one-value alphabet standing in for the centre Q `kvmer.csv` does not report
_INSERTIONS = slice(5, 9)  # error.CATEGORIES


@dataclass(frozen=True)
class _Loci:
    """Per-locus columns: context index and consensus base per value position, plus the observation counts."""

    contexts: list[str]  # unique context strings, indexed by `ctx`
    ctx: Array  # [n_loci, v] index into `contexts`
    centre: Array  # [n_loci, v] consensus base index
    total: Array  # [n_loci] observations that matched or carried one edit
    rows: dict[str, tuple[Array, Array]]  # op -> (locus indices, counts)


def _flank(tokens: Sequence[str]) -> tuple[int, int]:
    """The context window head E's `tokens` need: the widest `Context(L,R)`."""
    args = [Component(t).args for t in tokens if Component(t).name == "Context"]
    return max((a[0] for a in args), default=0), max((a[1] for a in args), default=0)


def _decode(code: int, width: int) -> str:
    return "".join("ACGT."[(code // 5**i) % 5] for i in range(width - 1, -1, -1))


def _prepare(table: CountTable, flank: tuple[int, int]) -> _Loci:
    if tuple(table.fields) != ("locus", "op"):
        raise ValueError(f"{table.source}: kvmer fits need fields ('locus', 'op'), got {table.fields}")
    k, v = int(table.meta["k"]), int(table.meta["v"])
    left, right = flank
    if left > k:
        raise ValueError(f"{table.source}: a Context left flank of {left} does not fit in a {k}-base key")
    loci: dict[str, int] = {}
    totals: Counter[int] = Counter()
    rows: dict[str, list[tuple[int, float]]] = {}
    for (locus, op), n in table.counts.items():
        if len(str(locus)) != k + v:
            raise ValueError(f"{table.source}: locus {locus!r} is not k + v = {k + v} bases")
        i = loci.setdefault(str(locus), len(loci))
        totals[i] += n
        if op != "=":
            rows.setdefault(str(op), []).append((i, n))

    bases = _BASE[np.frombuffer("".join(loci).encode(), np.uint8)].reshape(len(loci), k + v)
    centre = bases[:, k:]
    if centre.max(initial=0) > 3:
        raise ValueError(f"{table.source}: consensus values must be A, C, G or T")
    width = left + right + 1
    padded = np.concatenate([bases, np.full((len(loci), right), 4, np.int64)], axis=1)
    window = padded[:, k - left + np.arange(v)[:, None] + np.arange(width)]  # [n_loci, v, width]
    codes, ctx = np.unique(window @ 5 ** np.arange(width - 1, -1, -1), return_inverse=True)
    return _Loci(
        [_decode(int(c), width) for c in codes],
        ctx.reshape(len(loci), v),
        centre,
        np.array([totals[i] for i in range(len(loci))], float),
        {op: (np.array([i for i, _ in r]), np.array([n for _, n in r], float)) for op, r in rows.items()},
    )


def _target(op: str) -> tuple[int, int | None]:
    """(head E category, the consensus base the operation needs, None for an insertion)."""
    if op.startswith("->"):
        return error.CATEGORIES.index(op), None
    base = int(_BASE[ord(op[0])])
    if base > 3:
        raise ValueError(f"operation {op!r} does not act on a base")
    return (9 if op.endswith("-") else error.CATEGORIES.index(op[1:])), base


def _ratios(p: Array) -> Array:
    """r_op(t) per context [contexts, K]: the match odds of each operation, but an insertion's own probability."""
    r = p / np.maximum(p[:, :1], 1e-300)
    r[:, _INSERTIONS] = p[:, _INSERTIONS]
    return r


def _positions(loci: _Loci, op: str, ratios: Array) -> Array:
    """r_op(t) [rows, v] for one operation's rows, zero where its consensus base is not there."""
    cat, base = _target(op)
    r: Array = ratios[loci.ctx[loci.rows[op][0]], cat]
    return r if base is None else np.where(loci.centre[loci.rows[op][0]] == base, r, 0.0)


def _weights(loci: _Loci, ratios: Array) -> dict[str, Array]:
    """Per-(locus, op) row, the posterior [rows, v] over the latent position."""
    out = {}
    for op in loci.rows:
        r = _positions(loci, op, ratios)
        total = r.sum(axis=1, keepdims=True)
        if not total.all():
            raise ValueError(f"operation {op!r} counted at a locus with no consensus base it can act on")
        out[op] = r / total
    return out


def _table(table: CountTable, loci: _Loci, flank: tuple[int, int], weights: dict[str, Array]) -> CountTable:
    """Head E rows over contexts: errors split by the position posterior, matches as the leftover exposure."""
    counts: Counter[Key] = Counter()
    n_contexts = len(loci.contexts)
    match = np.repeat(loci.total[:, None], loci.ctx.shape[1], axis=1)
    for op, (idx, n) in loci.rows.items():
        w = n[:, None] * weights[op]
        if _target(op)[1] is not None:  # a substitution or deletion replaces that position's match
            np.add.at(match, idx, -w)
        by_context = np.bincount(loci.ctx[idx].ravel(), w.ravel(), minlength=n_contexts)
        counts.update({(loci.contexts[c], DUMMY_Q, op): x for c, x in enumerate(by_context) if x > 0})
    by_context = np.bincount(loci.ctx.ravel(), match.ravel(), minlength=n_contexts)
    counts.update({(loci.contexts[c], DUMMY_Q, "="): x for c, x in enumerate(by_context) if x > 0})
    meta = {**table.meta, "flank": list(flank)}
    return CountTable(f"{table.source}:latent-position", ("context", "q", "op"), "base", True, counts, meta)


def _log_likelihood(loci: _Loci, p: Array, ratios: Array) -> float:
    """The mixture log-likelihood above, in nats."""
    ll = float(loci.total @ np.log(np.maximum(p[loci.ctx, 0], 1e-300)).sum(axis=1))
    for op, (_, n) in loci.rows.items():
        ll += float(n @ np.log(np.maximum(_positions(loci, op, ratios).sum(axis=1), 1e-300)))
    return ll


def fit(
    table: CountTable,
    tokens: Sequence[str],
    *,
    flank: tuple[int, int] | None = None,
    iterations: int = 20,
    tol: float = 1e-4,
    l2: float = 1.0,
    seed: int = 0,
) -> tuple[Component, ...]:
    """Fit head E `tokens` from a `skiver_analyze:kvmer` table by EM over the latent edit position.

    `tokens` must start with `QualityWindow(0)` (the intercept; `kvmer.csv` has no quality) and may use only
    the context components `Context(L,R)` and `Homopolymer`. `flank` is the context window the rows carry
    (default: the widest `Context(L,R)`), so `Homopolymer` can be given a wider one. EM stops when the
    log-likelihood gains less than `tol` nats per observation.
    """
    if list(tokens[:1]) != ["QualityWindow(0)"] or not {Component(t).name for t in tokens[1:]} <= {
        "Context",
        "Homopolymer",
    }:
        raise ValueError(f"kvmer fits need QualityWindow(0) then Context(L,R)/Homopolymer, got {list(tokens)}")
    window = _flank(tokens) if flank is None else flank
    loci = _prepare(table, window)
    contexts = CountTable(
        table.source,
        ("context", "q"),
        "base",
        False,
        Counter({(c, DUMMY_Q): 1 for c in loci.contexts}),
        {"flank": list(window)},
    )
    weights = _weights(loci, _ratios(np.ones((len(loci.contexts), len(error.CATEGORIES))) / len(error.CATEGORIES)))
    components: tuple[Component, ...] | None = None
    last = -np.inf
    for _ in range(iterations):
        rows = _table(table, loci, window, weights)
        components = error.fit(rows, tokens, (DUMMY_Q,), l2=l2, seed=seed, init=components)
        p = error.probabilities(components, (DUMMY_Q,), contexts)
        ratios = _ratios(p)
        total = _log_likelihood(loci, p, ratios)
        weights = _weights(loci, ratios)
        if total - last < tol * loci.total.sum():
            break
        last = total
    assert components is not None
    return components
