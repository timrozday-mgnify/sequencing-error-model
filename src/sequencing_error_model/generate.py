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

`fragments` samples standalone read pairs from genome contigs, with insert sizes from the spec's
`insert_size` marginal (row 0 sizes, row 1 probabilities; `insert_sizes` builds one). Mate 2's template
is the fragment's reverse complement, and a mate reads through into its adapter when the insert is
shorter than the read.
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
ADAPTERS = ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")  # TruSeq read-through, R1/R2
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class Read:
    sequence: str
    quality: str  # Phred+33
    cigar: str
    q_track: Array  # template-indexed Q, deleted positions included
    clipped: tuple[str, str] = ("", "")  # Phred+33 Q of unaligned read bases before and after the aligned part
    strand: str | None = None  # alignment strand; None: the generator's convention, "-" for mate 2


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


def insert_sizes(mean: float, sd: float) -> Array:
    """A discretised normal insert-size marginal over mean ± 4 sd (at least 1): row 0 sizes, row 1 probabilities."""
    if mean < 1 or sd < 0:
        raise ValueError("insert size mean must be >= 1 and sd >= 0")
    sizes = np.arange(max(1, int(mean - 4 * sd)), int(np.ceil(mean + 4 * sd)) + 1)
    p = np.exp(-0.5 * ((sizes - mean) / sd) ** 2) if sd else (sizes == round(mean)).astype(float)
    return np.vstack([sizes, p / p.sum()])


def fragments(
    contigs: Sequence[tuple[str, str]],
    insert_size: Array,
    n: int,
    read_length: int,
    rng: np.random.Generator,
    first: int = 0,
) -> list[tuple[str, str, str]]:
    """`n` fragments as (name, mate 1 template, mate 2 template), uniform over the placements that fit.

    Names are `contig:start-end#i` (1-based, inclusive; i counts from `first`). Past the fragment a mate reads
    its adapter, then N.
    """
    sizes, p = np.asarray(insert_size[0], np.int64), np.asarray(insert_size[1], float)
    if read_length < 1 or not len(sizes) or sizes.min() < 1 or p.min() < 0 or not p.sum():
        raise ValueError("need read_length >= 1 and an insert_size marginal of sizes >= 1 with probabilities")
    lengths = np.array([len(s) for _, s in contigs], np.int64)
    out = []
    # ponytail: O(n × contigs) Python loop; vectorise if many-contig standalone runs get slow.
    for i, size in enumerate(rng.choice(sizes, size=n, p=p / p.sum())):
        slots = np.clip(lengths - size + 1, 0, None)
        if not slots.sum():
            raise ValueError(f"insert size {size} is longer than every contig")
        c = int(rng.choice(len(contigs), p=slots / slots.sum()))
        start = int(rng.integers(slots[c]))
        name, frag = contigs[c][0], contigs[c][1][start : start + size].upper()
        r1, r2 = (s[:read_length].ljust(read_length, "N") for s in (frag + ADAPTERS[0], _revcomp(frag) + ADAPTERS[1]))
        out.append((f"{name}:{start + 1}-{start + size}#{first + i}", r1, r2))
    return out


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


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
    A read's `clipped` bases shift `pos_start` / `pos_end` to the read's own ends and fill Q windows past the
    aligned part; `strand` overrides the mate-based strand. Rows at template bases other than A, C, G, T are
    skipped.
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
        padded, gc = "." * left + template.upper() + "." * right, gc_bin(template)
        before, after = ([ord(c) - 33 for c in s] for s in read.clipped)
        filled = [*before, *filled, *after]
        strand = read.strand or ("-" if mate == 2 else "+")
        for t, op in rows:
            i = t + len(before)
            if template[t].upper() not in "ACGT" or filled[i] is None:
                continue
            window = tuple(
                filled[i + o] if 0 <= i + o < len(filled) else None for k in range(1, m + 1) for o in (-k, k)
            )
            key = (filled[i], *window, padded[t : t + left + right + 1], i + 1, len(filled) - i, mate, strand, gc, op)
            counts[key] += 1
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
    p.add_argument("--input", type=Path, help="template FASTA[.gz], or genome FASTA with --pairs (default: stdin)")
    p.add_argument("--output", type=Path, help="FASTQ[.gz] (default: stdout); headers carry cigar:CIGAR")
    p.add_argument("--paired", action="store_true", help="interleaved R1/R2; names ending /2 are mate 2")
    p.add_argument("--pairs", type=int, help="sample this many interleaved pairs from the input genome instead")
    p.add_argument("--read-length", type=int, help="read length for --pairs")
    p.add_argument("--insert-mean", type=float, help="normal insert size for --pairs (default: the spec's)")
    p.add_argument("--insert-sd", type=float, default=0.0, help="insert size sd with --insert-mean")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--max-ins-run", type=int, default=10, help="max insertions before one template base")
    p.add_argument("--error-rate-scale", type=float, default=1.0, help="multiplier on every error probability")
    p.add_argument("--no-quality", action="store_true", help="emit '*' qualities (compatibility only)")
    p.add_argument("--batch", type=int, default=4096, help=argparse.SUPPRESS)
    args = p.parse_args(argv)
    if args.max_ins_run < 1 or args.error_rate_scale < 0:
        p.error("--max-ins-run must be >= 1 and --error-rate-scale >= 0")
    model, rng = spec_io.load(args.model), np.random.default_rng(args.seed)
    sizes = model.marginals.get("insert_size", np.zeros((2, 0)))
    if args.pairs is not None:
        if args.paired or args.pairs < 0 or not args.read_length or args.read_length < 1:
            p.error("--pairs needs --read-length >= 1, a count >= 0, and excludes --paired")
        if args.insert_mean is not None:
            try:
                sizes = insert_sizes(args.insert_mean, args.insert_sd)
            except ValueError as e:
                p.error(str(e))
        if not sizes.shape[1]:
            p.error("--pairs needs --insert-mean or a spec with an insert_size marginal")
    src = _open(args.input, "r") if args.input else sys.stdin
    out = _open(args.output, "w") if args.output else sys.stdout

    def records() -> Iterator[tuple[str, str, int]]:
        if args.pairs is None:
            for name, seq in _fasta(src):
                yield name, seq, 2 if args.paired and name.endswith("/2") else 1
            return
        contigs = list(_fasta(src))
        for first in range(0, args.pairs, args.batch):
            n = min(args.batch, args.pairs - first)
            for name, r1, r2 in fragments(contigs, sizes, n, args.read_length, rng, first):
                yield from ((f"{name}/1", r1, 1), (f"{name}/2", r2, 2))

    try:
        batch: list[tuple[str, str, int]] = []
        for record in (*records(), None):
            if record is not None:
                batch.append(record)
            if batch and (record is None or len(batch) >= args.batch):
                reads = generate(
                    model,
                    [s for _, s, _ in batch],
                    [m for _, _, m in batch],
                    rng,
                    max_ins_run=args.max_ins_run,
                    error_rate_scale=args.error_rate_scale,
                )
                for (name, _, _), read in zip(batch, reads, strict=True):
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
