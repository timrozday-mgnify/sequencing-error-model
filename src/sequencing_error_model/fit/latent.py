"""`Latent(S)`: a read-level class shared by head E and head Q (plan §5.2, §5.4), fitted by EM.

Each read draws a class s from `prior` [S]. Both heads carry a `Latent(S)` component whose weights [S, K] shift
their logits for that class, so one class can hold reads with lower qualities and more errors at the same
qualities: the per-read heterogeneity the baseline runs found on both platforms.

EM over reads: the E-step weighs each read's classes by the joint log-likelihood of its head E and head Q rows;
the M-step refits both heads (warm-started) on the rows expanded by class with those weights, and the prior as
their mean. Classes start from quantiles of each read's observed error rate (op labels, never Q), ties broken
at random.

ponytail: one class per read (a mixture), not an HMM over positions; add transitions if long reads show regimes
changing within a read. Expanding rows by class costs S x the rows in Python; fine for thousands of reads.
"""

import warnings
from collections import Counter
from collections.abc import Sequence
from dataclasses import replace

import numpy as np
from scipy import special

from sequencing_error_model.fit import error, quality
from sequencing_error_model.fit.quality import Array
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.spec import Component


def _reads(table: CountTable) -> Array:
    ri = table.fields.index("read")
    return np.array([k[ri] for k in table.counts], np.int64)


def _expand(table: CountTable, weights: Array) -> CountTable:
    """Rows without `read`, one copy per class with a `latent` field, counts scaled by the read's class weight."""
    ri, n_classes = table.fields.index("read"), weights.shape[1]
    counts: Counter[Key] = Counter()
    for (key, n), r in zip(table.counts.items(), _reads(table), strict=True):
        rest, w = key[:ri] + key[ri + 1 :], weights[r]
        for s in range(n_classes):
            if w[s] > 1e-6:
                counts[(*rest, s)] += n * w[s]
    fields = (*table.fields[:ri], *table.fields[ri + 1 :], "latent")
    return replace(table, fields=fields, counts=counts)


def _as_class(table: CountTable, s: int) -> CountTable:
    return replace(
        table, fields=(*table.fields, "latent"), counts=Counter({(*k, s): n for k, n in table.counts.items()})
    )


def fit(
    table: CountTable,
    q_table: CountTable,
    error_tokens: Sequence[str],
    quality_tokens: Sequence[str],
    alphabet: Sequence[int],
    n_classes: int,
    *,
    iterations: int = 20,
    tol: float = 1e-4,
    seed: int = 0,
) -> tuple[tuple[Component, ...], tuple[Component, ...], Component]:
    """Head E (from `table`), head Q (from `q_table`) and `Latent(S)` with its prior, by EM.

    Both tables need a `read` field indexing the same reads. `Latent(S)` is appended to both heads' tokens. EM stops
    when the log-likelihood gains less than `tol` nats per read.
    """
    if n_classes < 2:
        raise ValueError("Latent(S) needs at least 2 classes")
    if "read" not in table.fields or "read" not in q_table.fields:
        raise ValueError("Latent(S) needs per-read rows: tables with a `read` field")
    token = f"Latent({n_classes})"
    e_tokens, q_tokens = [*error_tokens, token], [*quality_tokens, token]
    reads, oi = _reads(table), table.fields.index("op")
    n_reads = 1 + int(max(reads.max(initial=0), _reads(q_table).max(initial=0)))
    n = np.array(list(table.counts.values()), float)
    is_error = np.array([k[oi] != "=" for k in table.counts])
    draws = np.bincount(reads, n, minlength=n_reads)
    errs = np.bincount(reads, n * is_error, minlength=n_reads)
    rng = np.random.default_rng(seed)
    rank = np.lexsort((rng.random(n_reads), errs / np.maximum(draws, 1)))
    weights = np.zeros((n_reads, n_classes))
    weights[rank, np.arange(n_reads) * n_classes // n_reads] = 1.0

    e_head: tuple[Component, ...] | None = None
    q_head: tuple[Component, ...] | None = None
    last = -np.inf
    for _ in range(iterations):
        e_head = error.fit(_expand(table, weights), e_tokens, alphabet, init=e_head, seed=seed)
        q_head = quality.fit(_expand(q_table, weights), q_tokens, alphabet, init=q_head)
        prior = weights.mean(axis=0)
        ll = np.log(np.maximum(prior, 1e-300)) + np.stack(
            [
                error.read_log_likelihood(e_head, alphabet, _as_class(table, s), n_reads)
                + quality.read_log_likelihood(q_head, alphabet, _as_class(q_table, s), n_reads)
                for s in range(n_classes)
            ],
            axis=1,
        )
        total = float(special.logsumexp(ll, axis=1).sum())
        weights = special.softmax(ll, axis=1)
        if total - last < tol * n_reads:
            break
        last = total
    else:
        warnings.warn(
            f"Latent({n_classes}) EM did not converge in {iterations} iterations", RuntimeWarning, stacklevel=2
        )
    assert e_head is not None and q_head is not None
    return e_head, q_head, Component(token, {"prior": prior})
