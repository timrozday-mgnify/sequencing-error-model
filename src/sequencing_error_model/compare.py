"""Cross-mode comparison on shared support (plan §7, phase 5).

Two modes observe different things: `pe-overlap` sees only substitutions, only at overlap positions, and not errors
both mates share (PCR, library, cluster); `reference` sees every op at every aligned base, but also reference
errors and strain variation. So comparisons are restricted to what both observe.

- `evidence(a, b, by)`: two labelled head E tables (`pe_overlap.table`, `generate.observations` of BAM records)
  restricted to substitution support (match and substitution rows), marginalised onto the covariates `by`, and
  kept at covariate values both observe with at least `min_exposure` rows each. Each rate is standardised to the
  pooled exposure over those values, so differences in coverage (overlap positions, Q mix) don't read as rate
  differences. `excess` = rate_b − rate_a: with a = `pe-overlap` and b = `reference` on the same reads it
  estimates PCR/library substitutions, plus any residual variation or reference error.
- `models(a, b, table)`: both specs' head E on one table, conditioned on no indel: op TV, rates, rate by Q, and
  each spec's log-likelihood per row.
"""

from collections import Counter
from collections.abc import Sequence
from dataclasses import replace
from typing import Any

import numpy as np

from sequencing_error_model.fit import error, indel
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.spec import ErrorModelSpec


def _substitution_support(table: CountTable) -> CountTable:
    if "op" not in table.fields:
        raise ValueError(f"{table.source}: comparison needs op labels")
    oi = table.fields.index("op")
    kept: Counter[Key] = Counter({k: n for k, n in table.counts.items() if k[oi] == "=" or "-" not in str(k[oi])})
    return replace(table, counts=kept)


def _exposure(table: CountTable, by: Sequence[str]) -> dict[Key, tuple[float, float]]:
    """(rows, substitutions) per value of `by`, on substitution support."""
    if missing := sorted(set(by) - set(table.fields)):
        raise ValueError(f"{table.source}: no fields {missing}")
    t = _substitution_support(table)
    oi, idx = t.fields.index("op"), [t.fields.index(f) for f in by]
    out: dict[Key, tuple[float, float]] = {}
    for k, n in t.counts.items():
        rows, subs = out.get(key := tuple(k[i] for i in idx), (0.0, 0.0))
        out[key] = (rows + n, subs + n * (k[oi] != "="))
    return out


def evidence(a: CountTable, b: CountTable, by: Sequence[str] = ("q",), min_exposure: float = 100) -> dict[str, Any]:
    """Substitution rates of two tables on the covariate values both observe, standardised to pooled exposure."""
    ea, eb = _exposure(a, by), _exposure(b, by)
    shared = sorted((k for k in ea.keys() & eb.keys() if min(ea[k][0], eb[k][0]) >= min_exposure), key=str)
    if not shared:
        raise ValueError(f"{a.source} and {b.source} share no values of {list(by)} with {min_exposure} rows each")
    w = np.array([ea[k][0] + eb[k][0] for k in shared])
    ra, rb = (np.array([e[k][1] / e[k][0] for k in shared]) for e in (ea, eb))
    rate_a, rate_b = float(w @ ra / w.sum()), float(w @ rb / w.sum())
    return {
        "sources": [a.source, b.source],
        "by": list(by),
        "values": [str(k if len(k) != 1 else k[0]) for k in shared],
        "coverage": [sum(e[k][0] for k in shared) / sum(v[0] for v in e.values()) for e in (ea, eb)],
        "rate_a": rate_a,
        "rate_b": rate_b,
        "rate_ratio": rate_b / rate_a,
        "excess": rate_b - rate_a,
        "rates_a": ra.tolist(),
        "rates_b": rb.tolist(),
    }


def models(a: ErrorModelSpec, b: ErrorModelSpec, table: CountTable, by: str = "q") -> dict[str, Any]:
    """Both specs' head E on `table`'s substitution-support rows, each conditioned on no indel."""
    t = _substitution_support(table)
    oi, ci, f0 = t.fields.index("op"), t.fields.index("context"), t.meta.get("flank", (0, 0))[0]
    n = np.array(list(t.counts.values()), float)
    cat = np.array([error._category(k[oi], str(k[ci])[f0]) for k in t.counts])
    p = []
    for spec in (a, b):
        q = error.probabilities(indel.split(spec.error_head)[0], spec.quality_alphabet, t)[:, :5]
        p.append(q / q.sum(axis=1, keepdims=True))
    rates = [float(n @ (1 - x[:, 0]) / n.sum()) for x in p]
    key = np.array([str(k[t.fields.index(by)]) for k in t.counts])
    values = sorted(set(key), key=lambda v: (len(v), v))
    return {
        "rows": float(n.sum()),
        "op_tv": float(n @ (0.5 * np.abs(p[0] - p[1]).sum(axis=1)) / n.sum()),
        "rate_a": rates[0],
        "rate_b": rates[1],
        "rate_ratio": rates[1] / rates[0],
        "excess": rates[1] - rates[0],
        "log_likelihood_per_row": [float(n @ np.log(x[np.arange(len(n)), cat]) / n.sum()) for x in p],
        "by": by,
        "values": values,
        "rates_by": [[float(n[key == v] @ (1 - x[key == v, 0]) / n[key == v].sum()) for v in values] for x in p],
    }
