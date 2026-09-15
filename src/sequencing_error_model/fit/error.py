"""Head E: P(op_t | x_{t-L..t+R}, q_{t-m..t+m}, position, mate, strand, GC) as a log-linear softmax over
the error categories (plan §5.2, §5.4), fitted from exact-position tuples (`pe-overlap`, `reference`).

Categories (`CATEGORIES`, K = 10): "=" match, ">A".."">T" substitution to that base, "->A".."->T" insertion
before t, "-" deletion of t. Substituting a base by itself is impossible, so that category is masked in
the likelihood and in `probabilities` (plan §3, lesson 2): no fitted mass leaks into categories the
generator drops.

Each table row is one template base: the `context` centre is its true base (A, C, G or T) and `op` is
"=", "X>Y" with X the centre, "->Y" or "X>-".
ponytail: one outcome per template base, so an insertion and a substitution at the same base can't both
be counted; add a separate insertion slot when `reference` tuples need it.

Components and params (A = alphabet size, trailing axis K):

- `QualityWindow(m)` (required, first): `bias` [K], `window` [2m+1, A+1, K] over offsets -m..m; row A is
  "beyond the read end". Fields `q` and `q-1`..`q-m`, `q+1`..`q+m` (None beyond a read end).
- `Context(L,R)`: `weights` [L+R+1, 5, K] as in head Q. Field `context`.
- `Homopolymer`: `weights` [W, K] by the run length of the centre base inside the fitted table's context
  window (W bases wide, so runs are capped at W). Field `context`.
- `Position(n)`, `Mate`: as in head Q. `Strand`: `weights` [2, K] for "+", "-".
- `GC(n)`: `knots` [n] over GC %, `weights` [n, K]: linear spline of the `gc` bin midpoint.
- `QualityxContext(r)`: `quality` [A, r] and `context` [3, 5, r, K], a rank-r interaction of the centre Q
  with bases x_{t-1..t+1}.

Q enters only as a categorical feature, never as 10^(-Q/10): the calibration is learned.
"""

import warnings
from collections.abc import Sequence
from math import prod
from typing import Any

import numpy as np
from scipy import optimize, sparse, special

from sequencing_error_model.fit.quality import _BASE, Array, _hat
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.spec import Component

CATEGORIES = ("=", ">A", ">C", ">G", ">T", "->A", "->C", "->G", "->T", "-")
COMPONENTS = ("QualityWindow", "Context", "Homopolymer", "Position", "Mate", "Strand", "GC", "QualityxContext")
_ARITY = {"QualityWindow": 1, "Context": 2, "Homopolymer": 0, "Position": 1, "Mate": 0, "Strand": 0, "GC": 1}
_K = len(CATEGORIES)


def _check(components: Sequence[Component]) -> None:
    names = [c.name for c in components]
    if not names or names[0] != "QualityWindow" or not set(names) <= set(COMPONENTS) or len(set(names)) != len(names):
        raise ValueError(f"head E needs QualityWindow(m) first, then distinct {COMPONENTS[1:]}, got {names}")
    for c in components:
        a = c.args
        low = {"Position": 2, "GC": 2, "QualityxContext": 1}.get(c.name, 0)
        if len(a) != _ARITY.get(c.name, 1) or (low and a[0] < low):
            raise ValueError(f"bad arguments in {c.token!r}")


def _fields(c: Component) -> set[str]:
    if c.name == "QualityWindow":
        return {"q"} | {f"q{s}{i}" for i in range(1, c.args[0] + 1) for s in "+-"}
    extra = {"Position": {"pos_start", "pos_end"}, "Mate": {"mate"}, "Strand": {"strand"}, "GC": {"gc"}}
    return {"context"} | extra.get(c.name, set()) | ({"q"} if c.name == "QualityxContext" else set())


def _validate(components: Sequence[Component], fields: Sequence[str], meta: dict[str, Any], source: str) -> None:
    need = {"context"}.union(*(_fields(c) for c in components))
    if missing := sorted(need - set(fields)):
        raise ValueError(f"{source}: head E components {[c.token for c in components]} need fields {missing}")
    left, right = tuple(meta.get("flank", (0, 0)))
    for c in components:
        cover = c.args if c.name == "Context" else (1, 1) if c.name == "QualityxContext" else (0, 0)
        if cover[0] > left or cover[1] > right:
            raise ValueError(f"{source}: context flank {(left, right)} does not cover {c.token}")


def _data(fields: Sequence[str], keys: Sequence[Key], alphabet: Sequence[int], meta: dict[str, Any]) -> dict[str, Any]:
    """Per-row feature columns: alphabet indices, base indices, covariates."""
    n = len(keys)
    cols = dict(zip(fields, zip(*keys, strict=True), strict=True)) if n else {f: () for f in fields}
    q_index: dict[Any, int] = {q: i for i, q in enumerate(alphabet)}
    d: dict[str, Any] = {"k_q": len(alphabet), "flank": tuple(meta.get("flank", (0, 0)))}
    for f, col in cols.items():
        if f.startswith("q"):
            index = q_index if f == "q" else q_index | {None: len(alphabet)}
            if bad := set(col) - set(index):
                raise ValueError(f"Q values {sorted(map(str, bad))} in {f!r} are not in the alphabet {tuple(alphabet)}")
            d[f] = np.array([index[v] for v in col], np.int64)
        elif f == "context":
            d[f] = _BASE[np.frombuffer("".join(col).encode(), np.uint8)].reshape(n, -1)
        elif f == "strand":
            d[f] = np.array([s == "-" for s in col], np.int64)
        elif f == "gc":
            d[f] = np.array([(lo + hi) / 2 for lo, hi in col], float)
        else:
            d[f] = np.asarray(col)
    d["centre"] = d["context"][:, d["flank"][0]]
    if n and d["centre"].max() > 3:
        raise ValueError("head E rows need a true centre base A, C, G or T")
    return d


def _shapes(c: Component, d: dict[str, Any]) -> dict[str, tuple[int, ...]]:
    """Linear param shapes, in slot order, without the trailing K axis."""
    if c.name == "QualityWindow":
        return {"bias": (), "window": (2 * c.args[0] + 1, d["k_q"] + 1)}
    if c.name == "Context":
        return {"weights": (sum(c.args) + 1, 5)}
    if c.name == "Homopolymer":
        return {"weights": (c.params["weights"].shape[0] if "weights" in c.params else d["context"].shape[1],)}
    if c.name in ("Position",):
        return {"start": (c.args[0],), "end": (c.args[0],)}
    if c.name in ("Mate", "Strand"):
        return {"weights": (2,)}
    if c.name == "GC":
        return {"weights": (c.args[0],)}
    return {}  # QualityxContext is bilinear, handled by _interaction


def _design(c: Component, d: dict[str, Any]) -> tuple[Array, Array]:
    """(slot index, value) per row for the linear part of one component; both [R, J]."""
    n, f0 = len(d["centre"]), d["flank"][0]
    if c.name == "QualityWindow":
        m = c.args[0]
        q = np.stack([d[f"q{o:+d}" if o else "q"] for o in range(-m, m + 1)], axis=1)
        idx = np.hstack([np.zeros((n, 1), np.int64), 1 + np.arange(2 * m + 1) * (d["k_q"] + 1) + q])
        return idx, np.ones(idx.shape)
    if c.name == "Context":
        left, right = c.args
        window = d["context"][:, f0 - left : f0 + right + 1]
        return np.arange(window.shape[1]) * 5 + window, np.ones(window.shape)
    if c.name == "Homopolymer":
        same = d["context"] == d["centre"][:, None]
        run = np.cumprod(same[:, f0::-1], axis=1).sum(axis=1) + np.cumprod(same[:, f0 + 1 :], axis=1).sum(axis=1)
        return np.minimum(run, _shapes(c, d)["weights"][0])[:, None] - 1, np.ones((n, 1))
    if c.name == "Position":
        knots = c.params["knots"]
        vals = np.hstack([_hat(d["pos_start"], knots), _hat(d["pos_end"], knots)])
        return np.broadcast_to(np.arange(vals.shape[1]), vals.shape), vals
    if c.name == "Mate":
        return (d["mate"] - 1)[:, None], np.ones((n, 1))
    if c.name == "Strand":
        return d["strand"][:, None], np.ones((n, 1))
    if c.name == "GC":
        knots = c.params["knots"]
        vals = np.stack([np.interp(d["gc"], knots, e) for e in np.eye(len(knots))], axis=1)
        return np.broadcast_to(np.arange(vals.shape[1]), vals.shape), vals
    return np.zeros((n, 0), np.int64), np.zeros((n, 0))


def _matrix(components: Sequence[Component], d: dict[str, Any]) -> sparse.csr_matrix:
    blocks, offset = [], 0
    for c in components:
        idx, vals = _design(c, d)
        r = np.broadcast_to(np.arange(len(idx))[:, None], idx.shape)
        blocks.append((vals.ravel(), r.ravel(), (offset + idx).ravel()))
        offset += sum(prod(s) for s in _shapes(c, d).values())
    vals, r, col = (np.concatenate(x) for x in zip(*blocks, strict=True))
    return sparse.csr_matrix((vals, (r, col)), shape=(len(d["centre"]), offset))


def _centre(v: Array) -> Array:
    """Sum-to-zero over A, C, G, T per offset, so the interaction can't absorb a centre-Q main effect.

    A symmetric projection, so it also maps gradients.
    """
    out = v.copy()
    out[:, :4] -= v[:, :4].mean(axis=1, keepdims=True)
    return out


def _onehots(d: dict[str, Any]) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """Centre Q [R, A] and bases x_{t-1..t+1} [R, 15] as sparse indicator matrices."""
    n, f0 = len(d["centre"]), d["flank"][0]
    xq = sparse.csr_matrix((np.ones(n), (np.arange(n), d["q"])), shape=(n, d["k_q"]))
    ctx = np.arange(3) * 5 + d["context"][:, f0 - 1 : f0 + 2]
    xc = sparse.csr_matrix((np.ones(3 * n), (np.repeat(np.arange(n), 3), ctx.ravel())), shape=(n, 15))
    return xq, xc


def _interaction(u: Array, v: Array, xq: sparse.csr_matrix, xc: sparse.csr_matrix) -> tuple[Array, Array, Array]:
    """QualityxContext logits [R, K], plus the per-row factors [R, r] and [R, r, K] for its gradient."""
    rank = u.shape[1]
    uq = np.asarray(xq @ u)
    per_row = np.asarray(xc @ _centre(v).reshape(15, -1)).reshape(len(uq), rank, _K)
    return np.einsum("nj,njk->nk", uq, per_row), uq, per_row


def _category(op: Any, centre: str) -> int:
    if op == "=":
        return 0
    if op in CATEGORIES[5:9]:
        return CATEGORIES.index(op)
    if op == f"{centre}>-":
        return 9
    if isinstance(op, str) and len(op) == 3 and op[0] == centre and op[2] != "-":  # CountTable rejects "A>A"
        return CATEGORIES.index(op[1:])
    raise ValueError(f"op {op!r} is not an exact-position outcome for true base {centre!r}")


def _logits(components: Sequence[Component], d: dict[str, Any]) -> Array:
    theta = np.vstack([c.params[p].reshape(-1, _K) for c in components for p in _shapes(c, d)])
    logits = _matrix(components, d) @ theta
    for c in components:
        if c.name == "QualityxContext":
            logits = logits + _interaction(c.params["quality"], c.params["context"], *_onehots(d))[0]
    return np.asarray(logits)


def _mask(d: dict[str, Any]) -> Array:
    mask = np.zeros((len(d["centre"]), _K), bool)
    mask[np.arange(len(mask)), 1 + d["centre"]] = True
    return mask


def _labelled(
    table: CountTable, components: Sequence[Component], alphabet: Sequence[int]
) -> tuple[Array, dict[str, Any]]:
    """Validate a labelled table for `components`; return counts [R, K] and the feature columns."""
    if "op" not in table.fields:
        raise ValueError(f"{table.source}: head E needs op labels from a truth-bearing source")
    fields = [f for f in table.fields if f != "op"]
    _validate(components, fields, table.meta, table.source)
    oi, ci, f0 = table.fields.index("op"), fields.index("context"), tuple(table.meta.get("flank", (0, 0)))[0]
    keys: dict[Key, int] = {}
    r_idx, cat, n = [], [], []
    for key, cnt in table.counts.items():
        rest = key[:oi] + key[oi + 1 :]
        r_idx.append(keys.setdefault(rest, len(keys)))
        cat.append(_category(key[oi], str(rest[ci])[f0]))
        n.append(cnt)
    d = _data(fields, list(keys), alphabet, table.meta)
    counts = np.zeros((len(keys), _K))
    np.add.at(counts, (r_idx, cat), n)
    return counts, d


def log_likelihood(components: Sequence[Component], alphabet: Sequence[int], table: CountTable) -> float:
    """Log-likelihood of a labelled table's counts under fitted head E components."""
    _check(components)
    counts, d = _labelled(table, components, alphabet)
    logits = np.where(_mask(d), -np.inf, _logits(components, d))
    logp = logits - special.logsumexp(logits, axis=1, keepdims=True)
    return float(np.sum(counts[counts > 0] * logp[counts > 0]))


def fit(
    table: CountTable, tokens: Sequence[str], alphabet: Sequence[int], *, l2: float = 1.0, seed: int = 0
) -> tuple[Component, ...]:
    """Fit head E components (in `tokens` order) by penalised maximum likelihood.

    `l2` is a Gaussian prior precision on every weight; it also pins the softmax gauge. `seed` initialises
    the `QualityxContext` factors (zero is a saddle point).
    """
    components = [Component(t) for t in tokens]
    _check(components)
    counts, d = _labelled(table, components, alphabet)

    for i, c in enumerate(components):
        if c.name == "Position":
            top = np.log(max(d["pos_start"].max(), d["pos_end"].max(), 2))
            components[i] = Component(c.token, {"knots": np.linspace(0, top, c.args[0])})
        elif c.name == "GC":
            components[i] = Component(c.token, {"knots": np.linspace(0, 100, c.args[0])})
    x = _matrix(components, d)
    n_lin, k_q = x.shape[1] * _K, len(alphabet)
    rank = next((c.args[0] for c in components if c.name == "QualityxContext"), 0)
    split = n_lin + k_q * rank
    mask, totals = _mask(d), counts.sum(axis=1)
    start = np.zeros(split + 15 * rank * _K)
    start[n_lin:split] = np.random.default_rng(seed).normal(scale=0.1, size=k_q * rank)
    xq, xc = _onehots(d) if rank else (x[:, :0], x[:, :0])

    def objective(flat: Array) -> tuple[float, Array]:
        logits = np.asarray(x @ flat[:n_lin].reshape(-1, _K))
        if rank:
            u, v = flat[n_lin:split].reshape(k_q, rank), flat[split:].reshape(3, 5, rank, _K)
            inter, uq, per_row = _interaction(u, v, xq, xc)
            logits += inter
        logits[mask] = -np.inf
        top = logits.max(axis=1, keepdims=True)
        logp = logits - top - np.log(np.exp(logits - top).sum(axis=1, keepdims=True))
        logp[mask] = 0.0
        nll = -np.sum(counts * logp) + 0.5 * l2 * flat @ flat
        g = totals[:, None] * np.where(mask, 0.0, np.exp(logp)) - counts
        grads = [np.asarray(x.T @ g).ravel()]
        if rank:
            grads.append(np.asarray(xq.T @ np.einsum("njk,nk->nj", per_row, g)).ravel())
            dv = np.asarray(xc.T @ (uq[:, :, None] * g[:, None, :]).reshape(len(g), -1))
            grads.append(_centre(dv.reshape(3, 5, rank, _K)).ravel())
        return float(nll), np.concatenate(grads) + l2 * flat

    res = optimize.minimize(objective, start, jac=True, method="L-BFGS-B")
    if not res.success:
        warnings.warn(f"head E fit did not converge: {res.message}", RuntimeWarning, stacklevel=2)

    out, theta, i = [], res.x[:n_lin].reshape(-1, _K), 0
    for c in components:
        params = dict(c.params)
        for name, shape in _shapes(c, d).items():
            params[name] = theta[i : i + prod(shape)].reshape(*shape, _K)
            i += prod(shape)
        if c.name == "QualityxContext":
            v = _centre(res.x[split:].reshape(3, 5, rank, _K))
            params = {"quality": res.x[n_lin:split].reshape(k_q, rank), "context": v}
        out.append(Component(c.token, params))
    return tuple(out)


def probabilities(components: Sequence[Component], alphabet: Sequence[int], table: CountTable) -> Array:
    """P(category) [R, K] for each key of `table`, in `table.counts` order; masked categories are 0."""
    _check(components)
    fields = [f for f in table.fields if f != "op"]
    _validate(components, fields, table.meta, table.source)
    idx = [table.fields.index(f) for f in fields]
    d = _data(fields, [tuple(k[i] for i in idx) for k in table.counts], alphabet, table.meta)
    return np.asarray(special.softmax(np.where(_mask(d), -np.inf, _logits(components, d)), axis=1))
