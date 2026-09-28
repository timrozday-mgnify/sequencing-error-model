"""Head E for `kmer` default mode: the context from `kvmer.csv`, the centre-Q term by marginal matching (§6.1).

**`fit`: the context, through a latent position.** A `kvmer.csv` row is one locus: a key (k bases) plus its
consensus value (v bases), the number of observed values matching the consensus throughout, and, per error
operation, how many observed values differ from the consensus by exactly that one edit. *Where* in the value
the edit sat is not reported, so each single-edit count is a mixture over the value positions whose consensus
base the operation can act on (any position for an insertion). Positions differ in context across loci, which
is what identifies `Context(L,R)`.

EM over that latent position: the E-step weighs the positions of one (locus, operation) count by the current
model's odds of that operation there, and the M-step refits head E on those positions' contexts with those
weights, plus the match exposure every observation leaves at every other position. The likelihood of a locus
with N observations that matched or carried one edit is

    prod_t P(= | ctx_t)^N  *  prod_ops ( sum_t r_op(t) )^n_op,

with r_op(t) = P(op | ctx_t) / P(= | ctx_t) for a substitution or deletion (the edit replaces that position's
match) and r_op(t) = P(op | ctx_t) for an insertion (it sits *before* a position that still matches, as head E
draws it). Contexts at the value's right end run past the locus, and those bases are "." (unknown), the same
level head E uses beyond a read end.

`kvmer.csv` carries no quality, so that fit runs on a one-value dummy alphabet and its `QualityWindow(0)` is
only an intercept.

**`rake` and `fit_centre_q`: the centre-Q term.** `summary_phred.csv` gives P(error | centre Q) and the fitted
context head gives P(error | context); the FASTQ gives the exposure, the joint count of the two. `rake` spreads
error counts over those cells by IPF until both margins are reproduced, so a difficult context at low Q is not
counted twice, and `fit_centre_q` fits the log-additive head E to the result. Neighbouring-Q and Q x context
terms are not identifiable here and stay out of the model (plan §5.6, §6.2).

ponytail: the per-locus `consensus_count_up_to_v*` columns also localise each observation's *first* mismatch,
which would pin the latent position directly; they mix in the multi-edit values skiver drops from the op
counts, so they are left out. Fold them in as an E-step prior if a recovery run shows position unidentifiable.
"""

import warnings
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from sequencing_error_model.fit import error
from sequencing_error_model.fit.quality import _BASE, Array
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.sources.skiver_analyze import RateBin
from sequencing_error_model.spec import Component

DUMMY_Q = 0  # the one-value alphabet standing in for the centre Q `kvmer.csv` does not report
_INSERTIONS = slice(5, 9)  # error.CATEGORIES
_TINY = 1e-300


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


def _op(category: int, centre: str) -> str:
    """The `op` label for a head E category at a base, inverting `error._category`."""
    name = error.CATEGORIES[category]
    return name if name.startswith("->") else f"{centre}>-" if name == "-" else centre + name


def _ratios(p: Array) -> Array:
    """r_op(t) per context [contexts, K]: the match odds of each operation, but an insertion's own probability."""
    r = p / np.maximum(p[:, :1], _TINY)
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
    ll = float(loci.total @ np.log(np.maximum(p[loci.ctx, 0], _TINY)).sum(axis=1))
    for op, (_, n) in loci.rows.items():
        ll += float(n @ np.log(np.maximum(_positions(loci, op, ratios).sum(axis=1), _TINY)))
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


def _rates(phred: Sequence[RateBin]) -> dict[int, float]:
    """P(error | centre Q) per reported Q from `summary_phred.csv`, skipping bins with no exposure.

    The rate is skiver's per-Q Weibull rate, not `num_error / (num_correct + num_error)`: the counts come from
    a scan that stops at the first mismatch, so their ratio is a hazard, not a per-base rate (plan §5.6).
    """
    return {b.lo: b.per_base_error_rate for b in phred if b.num_correct + b.num_error > 0}


def _cells(
    exposure: CountTable, alphabet: Sequence[int], rates: dict[int, float]
) -> tuple[list[str], Array, Array, Array]:
    """(contexts, context index, Q index into `alphabet`, count) for the exposure cells head E can use."""
    if not {"context", "q"} <= set(exposure.fields):
        raise ValueError(f"{exposure.source}: the exposure needs `context` and `q` fields, got {exposure.fields}")
    if exposure.truth:
        raise ValueError(f"{exposure.source}: the exposure is unlabelled reads; a truth-bearing table is a mistake")
    centre = int(tuple(exposure.meta.get("flank", (0, 0)))[0])
    index = {q: i for i, q in enumerate(alphabet) if q in rates}
    contexts: dict[str, int] = {}
    rows: list[tuple[int, int, float]] = []
    kept = dropped = 0.0
    for (ctx, q), n in exposure.marginal("context", "q").counts.items():
        if str(ctx)[centre] in "ACGT" and q in index:
            rows.append((contexts.setdefault(str(ctx), len(contexts)), index[int(str(q))], n))
            kept += n
        else:
            dropped += n
    if dropped:
        warnings.warn(
            f"{exposure.source}: {dropped / (kept + dropped):.2%} of the exposure has an N centre base or a Q "
            f"outside both the alphabet and skiver's reported bins, and gets no centre-Q evidence",
            RuntimeWarning,
            stacklevel=3,
        )
    if not rows:
        raise ValueError(f"{exposure.source}: no exposure left to match skiver's Q marginal against")
    ctx_idx, q_idx, n_cell = (np.array(col) for col in zip(*rows, strict=True))
    return list(contexts), ctx_idx, q_idx, n_cell.astype(float)


def _ipf(n: Array, keys: Sequence[Array], margins: Sequence[Array], iterations: int, tol: float) -> Array:
    """Error counts per cell whose `keys` margins are `margins`, raked from the margins' product.

    Starting at the product leaves the cell odds with no interaction beyond what the margins force, which is the
    point: the alternative (multiplying both rates into every cell) double-counts a hard context at low Q.
    Both margins must total the same; the caller rescales.
    """
    e = n.copy()
    for key, margin in zip(keys, margins, strict=True):
        e *= margin[key] / np.maximum(np.bincount(key, n, minlength=len(margin)), _TINY)[key]
    e *= margins[0].sum() / max(e.sum(), _TINY)
    for _ in range(iterations):
        gap = 0.0
        for key, margin in zip(keys, margins, strict=True):
            got = np.bincount(key, e, minlength=len(margin))
            gap = max(gap, float(np.abs(margin - got).sum() / max(margin.sum(), _TINY)))
            e *= np.where(got[key] > 0, margin[key] / np.maximum(got[key], _TINY), 0.0)
        if gap < tol:
            break
    else:
        warnings.warn(f"marginal matching did not converge in {iterations} iterations", RuntimeWarning, stacklevel=3)
    return e


def rake(
    context_head: Sequence[Component],
    phred: Sequence[RateBin],
    exposure: CountTable,
    alphabet: Sequence[int],
    *,
    iterations: int = 200,
    tol: float = 1e-8,
) -> CountTable:
    """Expected head E counts over the FASTQ exposure, matching skiver's Q and context error marginals (§5.6).

    `context_head` is a head E fitted from `kvmer.csv` (`fit`), `phred` is `summary_phred.csv` and `exposure`
    is `fastq_quality:context`, the joint count of centre Q and *observed* context. Cell error counts are raked
    (IPF) until they sum to `phred`'s rate per Q and to `context_head`'s rate per context. The level comes from
    `phred` alone: the context margin is rescaled to its total, which is also what corrects the truncation bias
    of a `kvmer.csv`-only fit (§6.2). Each cell's error mass is split over the operations by `context_head`'s
    composition for that context, so the op mix stays a function of context only.

    Observed bases stand in for the true ones in the exposure, and `phred`'s margin counts skiver's value-side,
    first-error-stopping scan while the exposure is every base of every read: both are approximations this mode
    cannot avoid (§6.2).
    """
    rates = _rates(phred)
    contexts, ctx, q, n = _cells(exposure, alphabet, rates)
    table = CountTable(
        exposure.source,
        ("context", "q"),
        "base",
        False,
        Counter({(c, DUMMY_Q): 1 for c in contexts}),
        dict(exposure.meta),
    )
    p = error.probabilities(context_head, (DUMMY_Q,), table)
    by_q = np.bincount(q, n, minlength=len(alphabet)) * np.array([rates.get(a, 0.0) for a in alphabet])
    by_ctx = np.bincount(ctx, n, minlength=len(contexts)) * (1 - p[:, 0])
    by_ctx *= by_q.sum() / max(by_ctx.sum(), _TINY)
    e = _ipf(n, (q, ctx), (by_q, by_ctx), iterations, tol)
    if (over := e > n).any():
        warnings.warn(
            f"marginal matching wants more errors than bases in {int(over.sum())} cells "
            f"({n[over].sum() / n.sum():.2%} of the exposure); clipped, so the margins are approximate",
            RuntimeWarning,
            stacklevel=2,
        )
        e = np.minimum(e, n)
    composition = p[:, 1:] / np.maximum(1 - p[:, :1], _TINY)
    centre = int(tuple(exposure.meta.get("flank", (0, 0)))[0])
    counts: Counter[Key] = Counter()
    for c, qi, cell, errors in zip(ctx, q, n, e, strict=True):
        context = contexts[c]
        counts[(context, alphabet[qi], "=")] += cell - errors
        for cat, share in enumerate(composition[c], start=1):
            if share > 0:
                counts[(context, alphabet[qi], _op(cat, context[centre]))] += errors * share
    return CountTable(f"{exposure.source}:raked", ("context", "q", "op"), "base", True, +counts, dict(exposure.meta))


def fit_centre_q(
    context_head: Sequence[Component],
    phred: Sequence[RateBin],
    exposure: CountTable,
    tokens: Sequence[str],
    alphabet: Sequence[int],
    *,
    l2: float = 1e-3,
    **kwargs: float | None,
) -> tuple[Component, ...]:
    """Head E over the real quality alphabet, fitted to the raked counts (`rake`); `kwargs` go to `error.fit`.

    `tokens` must start with `QualityWindow(0)`: only the centre Q is identified in default mode, so no wider
    window may be asked for. `error.fit`'s `smooth` is worth a value on an unbinned alphabet, and it is what
    carries a Q bin skiver reports no rate for.

    `l2` is weak by default: the raked counts are expected values, not draws, so the usual shrinkage would pull
    the fit off the very margins it exists to reproduce (it showed as a few percent at the extreme Q, where the
    error mass is thinnest). It stays non-zero to pin the softmax gauge.
    """
    if list(tokens[:1]) != ["QualityWindow(0)"]:
        raise ValueError(f"default mode identifies the centre Q only, so tokens start with QualityWindow(0): {tokens}")
    table = rake(context_head, phred, exposure, alphabet)
    return error.fit(table, tokens, alphabet, l2=l2, **kwargs)  # type: ignore[arg-type]
