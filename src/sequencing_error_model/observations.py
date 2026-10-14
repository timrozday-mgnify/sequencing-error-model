"""Sparse count tables: the observation schema every evidence source reduces to (plan §7).

A `CountTable` counts deduplicated keys over a subset of the shared field vocabulary.
A table with fewer fields is a marginal of the joint model: fitters marginalise the
model onto `fields` before comparing. Later sources (skiver dump, PE overlap, BAM) emit
the same tables with more fields.

Quality is never error truth (plan §1): the `op` field is only allowed on tables from
truth-bearing sources (`truth=True`), so a table built from reported qualities cannot
carry error labels.

Counts are non-negative numbers. Fractional counts are expected counts, e.g. `pe-overlap`'s EM attribution
of a mate disagreement.

Field vocabulary (positions are 1-based):

- `op`: "=" (base matched, or observation survived), "!" (error of unrecorded type),
  or one of `ERROR_OPS` ("A>C" substitution, "->A" insertion, "A>-" deletion).
- `context`: bases x_{t-L..t+R} with `meta["flank"] = (L, R)`; "-" as the centre of an
  insertion, "." beyond a read end.
- `locus`: a whole skiver key + consensus value (k + v bases); the edit position is latent.
- `q`: the centre base's reported Q; `q-1`, `q+2`, ...: Q at that offset.
- `t`: position in the skiver value.
- `pos_start`, `pos_end`: position from the read start, from the read end.
- `strand` ("+", "-"), `mate` (1, 2), `gc` ((lo, hi) GC % bin), `length` (read length).
- `indel` ("I", "D"), `indel_length` (bases in one indel event), `run` (template homopolymer run length at the
  event); error labels like `op`, so they also need a truth-bearing source.
- `read` (index of the read in its source), `latent` (read-level class of `Latent(S)`, 0-based),
  `tile` (the read's flowcell tile, as an index into `Tile(n)`'s own levels).
"""

import re
from collections import Counter
from collections.abc import Hashable, Iterable, Iterator
from dataclasses import dataclass, field, replace
from itertools import islice
from typing import Any, TypeVar

# Upstream skiver `ALL_OPERATIONS` order.
ERROR_OPS = (
    *("A>C", "A>G", "A>T", "G>A", "G>C", "G>T", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
    *("->A", "->C", "->G", "->T"),
    *("A>-", "C>-", "G>-", "T>-"),
)
OP_VALUES = frozenset(("=", "!", *ERROR_OPS))
# What one count is: a base, a skiver value observation, a read, or an error event with no exposure.
UNITS = frozenset(("base", "value", "read", "error"))
_FIELD = re.compile(
    r"op|context|locus|q|q[+-][1-9]\d*|t|pos_start|pos_end|strand|mate|gc|length"
    r"|indel|indel_length|run|read|latent|tile"
)

Key = tuple[Hashable, ...]


def count(items: Iterable[tuple[Key, int]]) -> Counter[Key]:
    """Sum counts per key, dropping zeros."""
    c: Counter[Key] = Counter()
    for key, n in items:
        c[key] += n
    return +c


_T = TypeVar("_T")


def subset(items: Iterable[_T], limit: int | None = None, every_n: int = 1, offset: int = 0) -> Iterator[_T]:
    """Keep every `every_n`th record from `offset`, at most `limit` of them.

    **Never slice a run contiguously.** At identical read counts, a head slice of a real Illumina run
    reads 2.5x the per-base error rate of a tail slice of the same file (plan phase 6b item 2:
    0.0081 against 0.0032 on ERR2935851), because flowcell position is read-level heterogeneity and
    file order approximates it. So `limit` on its own is a *biased* subsample, not merely a small one,
    and any cross-mode ratio built from two different blocks of one file measures read selection.
    `every_n` spreads the sample over the whole run instead; `offset` takes a disjoint second sample.

    Synthetic reads from `generate` are exchangeable, so `limit` alone is sound there.
    """
    if every_n < 1 or offset < 0:
        raise ValueError(f"every_n must be >= 1 and offset >= 0, got every_n={every_n}, offset={offset}")
    return islice(islice(items, offset, None, every_n), limit)


@dataclass(frozen=True)
class CountTable:
    source: str
    fields: tuple[str, ...]
    unit: str
    truth: bool  # labels come from a truth-bearing source (skiver consensus, overlap, alignment)
    counts: Counter[Key]
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        bad = [f for f in self.fields if not _FIELD.fullmatch(f)]
        if bad or len(set(self.fields)) != len(self.fields):
            raise ValueError(f"{self.source}: unknown or repeated fields in {self.fields}")
        if self.unit not in UNITS:
            raise ValueError(f"{self.source}: unit must be one of {sorted(UNITS)}, got {self.unit!r}")
        if {"op", "indel"} & set(self.fields) and not self.truth:
            raise ValueError(f"{self.source}: op labels need a truth-bearing source; quality is never error truth")
        op = self.fields.index("op") if "op" in self.fields else None
        for key, n in self.counts.items():
            if len(key) != len(self.fields) or not isinstance(n, int | float) or n < 0:
                raise ValueError(f"{self.source}: bad entry {key!r}: {n!r} for fields {self.fields}")
            if op is not None and key[op] not in OP_VALUES:
                raise ValueError(f"{self.source}: unknown op {key[op]!r}")

    def marginal(self, *fields: str) -> "CountTable":
        """Sum counts over every field not in `fields`."""
        idx = [self.fields.index(f) for f in fields]
        return replace(self, fields=fields, counts=count((tuple(k[i] for i in idx), n) for k, n in self.counts.items()))
