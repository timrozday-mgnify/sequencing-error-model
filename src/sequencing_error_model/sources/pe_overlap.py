"""`pe-overlap` evidence mode (plan §1.1, §6.4): agreement of Illumina mates in their overlap as error truth.

Per pair, streaming:

1. `place` finds the gapless offset of reverse-complemented R2 on R1 from bases alone, never Q: the offset with
   the lowest mismatch fraction over at least `min_overlap` comparable bases (then the longest). When the insert
   is shorter than the read, the adapter tails fall outside the overlap, so read-through needs no adapter list.
2. The pair is dropped and counted when an indel is implied (splitting the overlap at a breakpoint and shifting
   one side by up to `max_shift` removes at least `min_gain` mismatches), or else when the best offset's
   mismatch fraction exceeds `max_mismatch`.
3. Every overlapping base where both mates read A, C, G or T gives one head E row per mate, in that mate's own
   read orientation (Q window, context, position, strand). Where the mates agree, their base is the template
   base (both erring to the same base is ignored as a nuisance). Where they disagree it is latent: one mate
   matched and the other substituted. Other context bases are the consensus inside the overlap (N where the
   mates disagree) and the mate's own bases outside it.
4. `fit` attributes each disagreement by EM with head E over both mates' rows, starting from a coin flip and
   never from Q.

Head Q is fitted from every base of every pair, placed or not, conditioned on observed bases (as
`fastq_quality`). Fitting it only on kept overlaps would bias it: the pairs dropped for indels or mismatches are
the low-Q ones. ponytail: observed bases stand in for true ones, fine at Illumina error rates; condition on the
overlap consensus if a platform's rate makes the difference visible.

Only substitutions are identified: the fitted head E gets -inf indel logits and provenance records
`identified_ops`. Errors both mates share (PCR, library, cluster) cancel and are not in the rates. The rows are
the overlap exposure, so position effects are fitted conditional on the positions the overlap covers.
"""

import argparse
import json
import sys
from collections import Counter
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import asdict, dataclass, field, replace
from itertools import islice
from pathlib import Path
from typing import cast

import numpy as np

from sequencing_error_model.fit import error, quality
from sequencing_error_model.fit.quality import _BASE, Array
from sequencing_error_model.generate import Read, _flank, _revcomp, gc_bin
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.sources.fastq_quality import read_fastq
from sequencing_error_model.spec import Component, ErrorModelSpec

_COMP = dict(zip("ACGT", "TGCA", strict=True))
Pair = tuple[str, Sequence[int], str, Sequence[int]]  # R1 bases, R1 Q, R2 bases, R2 Q


@dataclass
class Stats:
    pairs: int = 0
    no_overlap: int = 0
    indel: int = 0
    agree: int = 0  # overlapping template bases where the mates agree
    disagree: int = 0


@dataclass
class Evidence:
    m: int
    flank: tuple[int, int]
    certain: Counter[Key] = field(default_factory=Counter)
    # (mate 1 row, mate 2 row) if mate 1 read the template base, then the same if mate 2 did
    disputed: Counter[tuple[Key, Key, Key, Key]] = field(default_factory=Counter)
    # head Q rows (q-1..q-m, pos_start, pos_end, mate, observed context, q) over every base of every pair
    quality: Counter[Key] = field(default_factory=Counter)
    alphabet: set[int] = field(default_factory=set)
    stats: Stats = field(default_factory=Stats)


def fields(m: int) -> tuple[str, ...]:
    """Row fields, in `generate.observations` order."""
    window = (f"q{s}{o}" for o in range(1, m + 1) for s in "-+")
    return ("q", *window, "context", "pos_start", "pos_end", "mate", "strand", "gc", "op")


def place(
    r1: str, r2: str, min_overlap: int = 20, max_mismatch: float = 0.2, max_shift: int = 3, min_gain: int = 3
) -> tuple[int | None, bool]:
    """(offset s, indel): reverse-complemented `r2`[k] pairs with `r1`[s + k]; s is None if nothing qualifies."""
    a, b = (_BASE[np.frombuffer(s.upper().encode(), np.uint8)] for s in (r1, _revcomp(r2.upper())))
    if not len(a) or not len(b):
        return None, False
    ok = (a[:, None] < 4) & (b[None, :] < 4)
    diff = ok & (a[:, None] != b[None, :])
    diag = (np.subtract.outer(np.arange(len(a)), np.arange(len(b))) + len(b) - 1).ravel()
    n = np.bincount(diag, ok.ravel(), len(a) + len(b) - 1)
    frac = np.where(n >= max(min_overlap, 1), np.bincount(diag, diff.ravel(), len(n)) / np.maximum(n, 1), np.inf)
    if not np.isfinite(frac.min()):
        return None, False
    best = np.flatnonzero(frac == frac.min())
    s = int(best[np.argmax(n[best])]) - len(b) + 1

    i = np.arange(max(0, s), min(len(a), s + len(b)))

    def mismatches(offset: int) -> Array:
        k = i - offset
        inside = (k >= 0) & (k < len(b))
        return np.asarray(np.r_[0, np.cumsum(np.where(inside, diff[i, np.clip(k, 0, len(b) - 1)], True))])

    c0 = mismatches(s)
    best_split = int(c0[-1])
    for shift in (*range(-max_shift, 0), *range(1, max_shift + 1)):
        c1 = mismatches(s + shift)
        split = np.minimum(c0 + c1[-1] - c1, c1 + c0[-1] - c0)  # one side at s, the other shifted
        best_split = min(best_split, int(split.min()))
    # Checked before the gapless cap (an indel mid-overlap fails every gapless offset), but the split alignment
    # must pass the cap itself: unrelated mates, ~50% mismatched at their best offset, gain from any shift.
    if c0[-1] - best_split >= min_gain and best_split <= max_mismatch * len(i):
        return s, True
    return (s if frac.min() <= max_mismatch else None), False


def collect(
    pairs: Iterable[Pair],
    m: int = 1,
    flank: tuple[int, int] = (2, 2),
    min_overlap: int = 20,
    max_mismatch: float = 0.2,
    max_shift: int = 3,
    min_gain: int = 3,
) -> Evidence:
    """Overlap rows for Q windows of ±`m` and contexts of `flank`, counted over `pairs`.

    The mismatch cap is loose on purpose: a tight one drops exactly the error-rich pairs and biases rates down.
    """
    ev, (left, right) = Evidence(m, flank), flank
    # ponytail: Python loop per base; move row building to numpy if real runs are too slow.
    for r1, q1, r2, q2 in pairs:
        ev.stats.pairs += 1
        reads, qs = (r1.upper(), r2.upper()), (list(q1), list(q2))
        ev.alphabet.update(qs[0], qs[1])
        for mate, (seq, q) in enumerate(zip(reads, qs, strict=True), 1):
            padded = "." * left + seq + "." * right
            for p in range(len(seq)):
                lags = tuple(q[p - k] if p >= k else None for k in range(1, m + 1))
                ev.quality[(*lags, p + 1, len(seq) - p, mate, padded[p : p + left + right + 1], q[p])] += 1
        s, indel = place(reads[0], reads[1], min_overlap, max_mismatch, max_shift, min_gain)
        if s is None or indel:
            ev.stats.no_overlap += s is None
            ev.stats.indel += indel
            continue
        est = [list(reads[0]), list(reads[1])]
        sites = []
        for i in range(max(0, s), min(len(reads[0]), s + len(reads[1]))):
            j = len(reads[1]) - 1 - (i - s)
            b1, b2 = reads[0][i], _COMP.get(reads[1][j], "N")
            if b1 != b2 or b1 not in _COMP:
                est[0][i] = est[1][j] = "N"
            if b1 in _COMP and b2 in _COMP:
                sites.append((i, j, b1, b2))
        mates = [("".join(e), q, gc_bin("".join(e)), mate) for mate, (e, q) in enumerate(zip(est, qs, strict=True), 1)]
        for i, j, b1, b2 in sites:
            o2 = reads[1][j]
            if b1 == b2:
                ev.stats.agree += 1
                ev.certain[_row(mates[0], i, b1, b1, m, flank)] += 1
                ev.certain[_row(mates[1], j, o2, o2, m, flank)] += 1
            else:
                ev.stats.disagree += 1
                first = (_row(mates[0], i, b1, b1, m, flank), _row(mates[1], j, _COMP[b1], o2, m, flank))
                ev.disputed[(*first, _row(mates[0], i, b2, b1, m, flank), _row(mates[1], j, o2, o2, m, flank))] += 1
    return ev


def _row(
    read: tuple[str, list[int], tuple[int, int], int],
    p: int,
    centre: str,
    observed: str,
    m: int,
    flank: tuple[int, int],
) -> Key:
    """One mate's row at read position `p` if its true base there is `centre`; `read` is (bases, Q, GC bin, mate)."""
    seq, q, gc, mate = read
    n, (left, right) = len(seq), flank
    ctx = "".join(seq[p + o] if 0 <= p + o < n else "." for o in range(-left, right + 1))
    window = tuple(q[p + o] if 0 <= p + o < n else None for k in range(1, m + 1) for o in (-k, k))
    op = "=" if centre == observed else f"{centre}>{observed}"
    return (q[p], *window, ctx[:left] + centre + ctx[left + 1 :], p + 1, n - p, mate, "+-"[mate - 1], gc, op)


def fit(
    ev: Evidence,
    error_tokens: Sequence[str],
    quality_tokens: Sequence[str],
    *,
    iterations: int = 30,
    tol: float = 1e-3,
    seed: int = 0,
    sources: Sequence[str] = ("pe-overlap",),
) -> ErrorModelSpec:
    """Fit both heads from overlap evidence, attributing disagreements by EM with head E.

    Each iteration refits head E (warm-started) on expected counts, then updates P(mate 1 read the template
    base) per disagreement; it stops when no posterior moves by `tol`, or after `iterations`.
    """
    alphabet, groups = tuple(sorted(ev.alphabet)), list(ev.disputed.items())
    first, sizes = np.full(len(groups), 0.5), np.array([n for _, n in groups], float)
    fields_, meta = fields(ev.m), {"flank": ev.flank}
    error_head: tuple[Component, ...] | None = None
    for _ in range(max(iterations, 1)):
        counts: dict[Key, float] = dict(ev.certain)
        for (rows, _n), n, w in zip(groups, sizes, first, strict=True):
            for r, c in zip(rows, (n * w, n * w, n * (1 - w), n * (1 - w)), strict=True):
                counts[r] = counts.get(r, 0.0) + float(c)
        expected = cast("Counter[Key]", Counter({k: v for k, v in counts.items() if v > 0}))
        table = CountTable("pe-overlap", fields_, "base", True, expected, meta)
        error_head = error.fit(table, error_tokens, alphabet, seed=seed, init=error_head)
        if not groups:
            break
        keys = Counter(dict.fromkeys((k for rows, _ in groups for k in rows), 1))
        p = error.probabilities(error_head, alphabet, CountTable("pe-overlap", fields_, "base", True, keys, meta))
        prob = {k: p[r, error._category(k[-1], str(k[-7])[ev.flank[0]])] for r, k in enumerate(keys)}
        like = np.array([[prob[a1] * prob[a2], prob[b1] * prob[b2]] for (a1, a2, b1, b2), _ in groups])
        update = like[:, 0] / like.sum(axis=1)
        moved, first = np.abs(update - first).max(), update
        if moved < tol:
            break
    assert error_head is not None

    lags = [f"q-{i}" for i in range(1, ev.m + 1)]
    q_table = CountTable(
        "pe-overlap:quality", (*lags, "pos_start", "pos_end", "mate", "context", "q"), "base", False, ev.quality
    )
    q_lags = lags[: Component(quality_tokens[0]).args[0]]
    quality_head = quality.fit(
        replace(q_table.marginal(*q_lags, "pos_start", "pos_end", "mate", "context", "q"), meta={"flank": ev.flank}),
        quality_tokens,
        alphabet,
    )
    bias = error_head[0].params["bias"].copy()
    bias[5:] = -np.inf  # indels are not identified by gapless overlaps
    head0 = error_head[0]
    error_head = (Component(head0.token, {**head0.params, "bias": bias}, meta=head0.meta), *error_head[1:])
    provenance = {
        "mode": "pe-overlap",
        "sources": list(sources),
        "identified_ops": ["substitution"],
        "stats": asdict(ev.stats),
    }
    return ErrorModelSpec(alphabet, provenance, quality_head, error_head)


def mode(templates: list[str], reads: list[Read], mates: list[int], truth: ErrorModelSpec) -> ErrorModelSpec:
    """Recovery-harness mode: `reads` are interleaved mate 1, mate 2 pairs; fits with the truth's tokens."""
    e, q = _flank(truth.error_head), _flank(truth.quality_head)
    m = max(truth.error_head[0].args[0], truth.quality_head[0].args[0])
    pairs = (
        (a.sequence, [ord(c) - 33 for c in a.quality], b.sequence, [ord(c) - 33 for c in b.quality])
        for a, b in zip(reads[::2], reads[1::2], strict=True)
    )
    ev = collect(pairs, m, (max(e[0], q[0], 2), max(e[1], q[1], 2)))
    return fit(ev, [c.token for c in truth.error_head], [c.token for c in truth.quality_head], sources=["generator"])


def read_pairs(r1: Path, r2: Path) -> Iterator[Pair]:
    for (s1, q1), (s2, q2) in zip(read_fastq(r1), read_fastq(r2), strict=True):
        yield s1, q1, s2, q2


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m sequencing_error_model.sources.pe_overlap", description="Fit a spec from mate overlaps."
    )
    p.add_argument("r1", type=Path, help="mate 1 FASTQ[.gz]")
    p.add_argument("r2", type=Path, help="mate 2 FASTQ[.gz], in the same order")
    p.add_argument("--output", type=Path, required=True, metavar="SPEC_DIR")
    p.add_argument("--error-tokens", nargs="+", default=["QualityWindow(1)", "Context(1,1)", "Position(4)", "Mate"])
    p.add_argument("--quality-tokens", nargs="+", default=["QualityMarkov(1)", "Position(4)", "Mate", "Context(1,1)"])
    p.add_argument("--max-pairs", type=int)
    p.add_argument("--min-overlap", type=int, default=20)
    p.add_argument("--iterations", type=int, default=30, help="maximum EM iterations")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args(argv)
    m = max(Component(args.error_tokens[0]).args[0], Component(args.quality_tokens[0]).args[0])
    ev = collect(islice(read_pairs(args.r1, args.r2), args.max_pairs), m, min_overlap=args.min_overlap)
    spec = fit(
        ev,
        args.error_tokens,
        args.quality_tokens,
        iterations=args.iterations,
        seed=args.seed,
        sources=[str(args.r1), str(args.r2)],
    )
    spec.save(args.output)
    print(json.dumps(spec.provenance["stats"], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
