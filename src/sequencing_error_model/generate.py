"""Native generator (plan §5.3): template sequences → reads with bases, qualities and a CIGAR, from a spec.

Per batch of reads, vectorised across bases:

1. sample the template-indexed Q track from head Q (conditioned on true bases);
2. draw a head E outcome for every template base from its Q window;
3. materialise the read: a match or substitution emits one base with q_t; an insertion emits the inserted
   base and draws again at the same base (at most `max_ins_run` insertions); a deletion emits nothing and
   drops q_t.

ponytail: inserted bases reuse q_t until the `InsertionQuality` sub-head exists (needs `reference` tuples).
Template bases other than A, C, G, T are copied through as N matches.

`align` and `observations` invert this: head E tuples re-derived from template, read and CIGAR, the input of
the recovery harness (and the core of `reference` mode). The CLI (`sem-generate`) keeps the fork's
`skiver-generate` contract, so genome-blender can switch by pointing its generate command at it.
"""

import argparse
import gzip
import re
import sys
from collections import Counter
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any, cast

import numpy as np

from sequencing_error_model import spec as spec_io
from sequencing_error_model.fit import error, quality
from sequencing_error_model.fit.quality import _BASE, Array
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.spec import Component, ErrorModelSpec

_K = len(error.CATEGORIES)
_ASCII = np.frombuffer(b"ACGTN", np.uint8)
_M, _I, _D = (ord(c) for c in "MID")
_CIGAR = re.compile(r"(\d+)([MID])")


@dataclass(frozen=True)
class Read:
    sequence: str
    quality: str  # Phred+33
    cigar: str
    q_track: Array  # template-indexed Q, deleted positions included


def gc_bin(template: str) -> tuple[int, int]:
    """10 % GC bin of a template, the `gc` covariate."""
    lo = 10 * min(9, int(10 * sum(b in "GCgc" for b in template) / max(len(template), 1)))
    return lo, lo + 10


def _flank(components: Sequence[Component]) -> tuple[int, int]:
    need = [(0, 0)]
    for c in components:
        if c.name == "Context":
            need.append((c.args[0], c.args[1]))
        elif c.name == "QualityxContext":
            need.append((1, 1))
        elif c.name == "Homopolymer":
            need.append((c.meta["flank"][0], c.meta["flank"][1]))
    return max(n[0] for n in need), max(n[1] for n in need)


def _cigar(ops: Array) -> str:
    if not len(ops):
        return ""
    starts = np.flatnonzero(np.r_[True, ops[1:] != ops[:-1]])
    return "".join(f"{n}{chr(c)}" for n, c in zip(np.diff(np.r_[starts, len(ops)]), ops[starts], strict=True))


def generate(
    model: ErrorModelSpec,
    templates: Sequence[str],
    mates: Sequence[int],
    rng: np.random.Generator,
    *,
    max_ins_run: int = 10,
    error_rate_scale: float = 1.0,
) -> list[Read]:
    """One read per template; mate 2 also uses the "-" strand covariate."""
    if max_ins_run < 1:
        raise ValueError("max_ins_run must be at least 1")
    templates = [t.upper() for t in templates]
    alphabet = np.asarray(model.quality_alphabet)
    tracks = quality.sample(model.quality_head, model.quality_alphabet, templates, mates, rng)
    lengths = np.array([len(t) for t in templates], np.int64)
    n = int(lengths.sum())
    starts = np.r_[0, np.cumsum(lengths)[:-1]].astype(np.int64)
    within = np.arange(n) - np.repeat(starts, lengths)
    left, right = _flank(model.error_head)
    bases = _BASE[np.frombuffer("".join("." * left + t + "." * right for t in templates).encode(), np.uint8)]
    at = np.repeat(starts + np.arange(len(templates)) * (left + right), lengths) + within
    context = bases[at[:, None] + np.arange(left + right + 1)]
    qidx = np.searchsorted(alphabet, np.concatenate([*tracks, np.zeros(0, np.int64)]))
    mate, rlen = np.repeat(np.asarray(mates, np.int64), lengths), np.repeat(lengths, lengths)
    cols: dict[str, Any] = {
        "k_q": len(alphabet),
        "flank": (left, right),
        "context": context,
        "centre": context[:, left],
        "q": qidx,
        "pos_start": within + 1,
        "pos_end": rlen - within,
        "mate": mate,
        "strand": (mate == 2).astype(np.int64),
        "gc": np.repeat([sum(gc_bin(t)) / 2 for t in templates], lengths),
    }
    for o in range(1, model.error_head[0].args[0] + 1):
        for j, name in ((within - o, f"q-{o}"), (within + o, f"q+{o}")):
            ok = (j >= 0) & (j < rlen)
            cols[name] = np.where(ok, qidx[np.where(ok, np.arange(n) + j - within, 0)], len(alphabet))

    probs = error.probabilities_at(model.error_head, cols) if n else np.zeros((0, _K))
    probs[:, 1:] *= error_rate_scale
    probs[cols["centre"] > 3] = np.eye(_K)[0]
    probs /= probs.sum(axis=1, keepdims=True)
    no_ins = probs.copy()
    no_ins[:, 5:9] = 0
    no_ins /= no_ins.sum(axis=1, keepdims=True)

    def draw(p: Array, idx: Array) -> Array:
        return np.minimum((rng.random((len(idx), 1)) > p[idx].cumsum(axis=1)).sum(axis=1), _K - 1)

    final = draw(probs, np.arange(n))
    ev_pos, ev_cat = [], []
    pending = np.flatnonzero((final >= 5) & (final <= 8))
    for r in range(max_ins_run):
        ev_pos.append(pending)
        ev_cat.append(final[pending])
        final[pending] = draw(probs if r + 1 < max_ins_run else no_ins, pending)
        pending = pending[(final[pending] >= 5) & (final[pending] <= 8)]
    # Events in read order: a base's insertions (in draw order), then its final outcome.
    pos = np.concatenate([*ev_pos, np.arange(n)])
    rank = np.concatenate([np.full(len(p), r) for r, p in enumerate(ev_pos)] + [np.full(n, max_ins_run)])
    order = np.lexsort((rank, pos))
    pos, cat = pos[order], np.concatenate([*ev_cat, final])[order]
    base = np.where(cat == 0, context[pos, left], np.where(cat < 5, cat - 1, cat - 5))
    ops = np.where(cat == 9, _D, np.where(cat >= 5, _I, _M)).astype(np.uint8)
    seq, qual, emit = _ASCII[np.minimum(base, 4)], (alphabet[qidx[pos]] + 33).astype(np.uint8), ops != _D
    bounds = cast(list[int], np.searchsorted(pos, np.r_[starts, n]).tolist())
    reads = []
    for i, track in enumerate(tracks):
        a, b = bounds[i], bounds[i + 1]
        e = emit[a:b]
        reads.append(Read(seq[a:b][e].tobytes().decode(), qual[a:b][e].tobytes().decode(), _cigar(ops[a:b]), track))
    return reads


def align(template: str, read: Read) -> tuple[list[tuple[int, str]], list[int | None]]:
    """Walk the CIGAR: (template position, op) per draw, and each template base's aligned read Q (None if
    deleted). Insertions are attributed to the following template base."""
    rows: list[tuple[int, str]] = []
    tq: list[int | None] = []
    r = 0
    for count, op in _CIGAR.findall(read.cigar):
        for _ in range(int(count)):
            t = len(tq)
            if op == "I":
                rows.append((t, "->" + read.sequence[r]))
                r += 1
                continue
            x = template[t].upper()
            if op == "M":
                y = read.sequence[r]
                rows.append((t, "=" if x == y else f"{x}>{y}"))
                tq.append(ord(read.quality[r]) - 33)
                r += 1
            else:
                rows.append((t, f"{x}>-"))
                tq.append(None)
    return rows, tq


def observations(
    records: Iterable[tuple[str, Read, int]], flank: tuple[int, int], m: int, source: str = "generator"
) -> CountTable:
    """Head E tuples from (template, read, mate) triples, one row per draw, as `fit.error` expects.

    The Q window is template-indexed from the read: a deleted base takes the Q of the next read base (the
    previous one at the read end), so windows touching a deletion differ from the generator's Q track.
    Rows at template bases other than A, C, G, T are skipped.
    """
    left, right = flank
    fields = ("q", *(f"q{s}{o}" for o in range(1, m + 1) for s in "-+"))
    fields += ("context", "pos_start", "pos_end", "mate", "strand", "gc", "op")
    counts: Counter[Key] = Counter()
    for template, read, mate in records:
        rows, tq = align(template, read)
        filled, nxt = list(tq), None
        for t in reversed(range(len(tq))):
            filled[t] = nxt = tq[t] if tq[t] is not None else nxt
        for t in range(1, len(filled)):
            filled[t] = filled[t] if filled[t] is not None else filled[t - 1]
        n, padded, gc = len(template), "." * left + template.upper() + "." * right, gc_bin(template)
        for t, op in rows:
            if template[t].upper() not in "ACGT" or filled[t] is None:
                continue
            window = tuple(filled[t + o] if 0 <= t + o < n else None for k in range(1, m + 1) for o in (-k, k))
            strand = "-" if mate == 2 else "+"
            counts[(filled[t], *window, padded[t : t + left + right + 1], t + 1, n - t, mate, strand, gc, op)] += 1
    return CountTable(source, fields, "base", True, counts, {"flank": flank})


def _open(path: Path, mode: str) -> IO[str]:
    return cast(IO[str], gzip.open(path, mode + "t")) if path.suffix == ".gz" else open(path, mode)  # noqa: SIM115


def _fasta(handle: IO[str]) -> Iterator[tuple[str, str]]:
    name: str | None = None
    parts: list[str] = []
    for line in handle:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(parts)
            name, parts = (line[1:].split() or [""])[0], []
        elif line and not line.startswith(";"):
            parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="sem-generate", description="Simulate reads (bases, qualities, CIGAR) from a spec."
    )
    p.add_argument("--model", type=Path, required=True, metavar="SPEC_DIR", help="error model spec directory")
    p.add_argument("--input", type=Path, help="template FASTA[.gz] (default: stdin)")
    p.add_argument("--output", type=Path, help="FASTQ[.gz] (default: stdout); headers carry cigar:CIGAR")
    p.add_argument("--paired", action="store_true", help="interleaved R1/R2; names ending /2 are mate 2")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--max-ins-run", type=int, default=10, help="max insertions before one template base")
    p.add_argument("--error-rate-scale", type=float, default=1.0, help="multiplier on every error probability")
    p.add_argument("--no-quality", action="store_true", help="emit '*' qualities (compatibility only)")
    p.add_argument("--batch", type=int, default=4096, help=argparse.SUPPRESS)
    args = p.parse_args(argv)
    if args.max_ins_run < 1 or args.error_rate_scale < 0:
        p.error("--max-ins-run must be >= 1 and --error-rate-scale >= 0")
    model, rng = spec_io.load(args.model), np.random.default_rng(args.seed)
    src = _open(args.input, "r") if args.input else sys.stdin
    out = _open(args.output, "w") if args.output else sys.stdout
    try:
        batch: list[tuple[str, str]] = []
        for record in (*_fasta(src), None):
            if record is not None:
                batch.append(record)
            if batch and (record is None or len(batch) >= args.batch):
                mates = [2 if args.paired and name.endswith("/2") else 1 for name, _ in batch]
                reads = generate(
                    model,
                    [s for _, s in batch],
                    mates,
                    rng,
                    max_ins_run=args.max_ins_run,
                    error_rate_scale=args.error_rate_scale,
                )
                for (name, _), read in zip(batch, reads, strict=True):
                    qual = "*" * len(read.sequence) if args.no_quality else read.quality
                    out.write(f"@{name} cigar:{read.cigar}\n{read.sequence}\n+\n{qual}\n")
                batch = []
    finally:
        for handle in (src, out):
            if handle not in (sys.stdin, sys.stdout):
                handle.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
