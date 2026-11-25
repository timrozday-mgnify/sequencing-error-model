"""`pfam` evidence mode (docs/pfam_mode_plan.md): coding constraint and mate overlap as error truth, for reads
placed on Pfam profiles by kfpfam + HMMER.

Biology varies freely from column to column, so a read's mismatches alone do not identify the error rate. Two
anchors do:

- **In-frame stops.** Inside a domain a stop is almost never biology, so a stop codon is a substitution error in
  one of the sense codons one base away. Rows are one per base that *could* have made a stop, and carry the
  substitutions it could show in `ops` (see `fit.error`). `A>G` and `T>C` never make a stop on either strand and
  are not identified by this anchor.
- **Mate overlap.** Mates of one fragment read the same molecule, so biology cancels. Both mates must be placed on
  the same columns at the offset `pe_overlap.place` finds, and the bases *outside* the counted codons must agree.
  Tandem repeats still get through, so families whose mates disagree too often are weighted down like the stops'.
  All twelve substitutions are identified, at the overlapped cycles. Which mate erred is latent, as in
  `pe_overlap`, with the stack below as an extra prior.

**Exposure.** A codon is counted only if the rest of the read places it *without that codon*: the profile scores of
one side reach `t1` bits along the hit's alignment, and the other side reaches `t2` along the diagonal it would
continue on. A codon's own errors cannot decide whether it is counted, and a read whose stop truncated its ORF is
placed by the half that remains. For that to hold, a hit must be reported whenever a side scores `t1`: search with
a cutoff several bits below it (kfpfam's `pairs.tsv`, E <= 1000), **never a table cut at gathering thresholds**. Raw
profile bits exceed HMMER's domain score by several bits. `t1` is calibrated by its plateau (`Evidence.plateau`):
too low counts chance homology. `add_errors` plus a rerun of the search is the end-to-end check.

**The stack.** The codon a stop (or a disputed base) came from has a prior: the sample's own codons at that column
plus `lam` codons' worth of the HMM's emissions, `lam` chosen by leave-one-out likelihood of the sense codons.

**Stops that are biology** (pseudogenes, recoding, a gene's own end) recur in every read of their locus; an error does
not. `replicates` counts the other fragments that show the same codon, or another, between the same flanking bases,
and `drop_replicated` removes stops that recur more often than errors could. On error-free reads of the ZymoBIOMICS
genomes at 3x, 109 in-frame stops (2.3e-4 per base, whatever `t1` and `t2`) fell to 34, nearly all at loci no other
read covered: it needs coverage, so a library at 1x keeps most of its floor. `provenance["replication"]` says how many
stops were checked.

**Families with a real stop rate** (pseudogenes, retroelements, recoding) are found by a spike-and-slab on the
family's stops against its exposure, and weighted down whole (rows as well as stops). It looks at the outcome, so it
biases the rate down a little; `tests/test_pfam.py` measures that on clonal reads.

Quality is a covariate and head Q's output, never a label: stops are labelled by the genetic code.
Substitutions only: the fitted head E gets -inf indel logits. `Tile` and longer reads' frameshifts are not here.
"""

import argparse
import gzip
import itertools
import json
import math
import re
import sys
import warnings
from collections import Counter, defaultdict
from collections.abc import Collection, Iterable, Iterator, Mapping, Sequence
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import cast

import numpy as np
from scipy import special

from sequencing_error_model.fit import error, quality
from sequencing_error_model.fit.quality import Array
from sequencing_error_model.generate import _revcomp
from sequencing_error_model.observations import CountTable, Key, subset
from sequencing_error_model.sources import fastq_quality, pe_overlap
from sequencing_error_model.sources.fastq_quality import read_fastq
from sequencing_error_model.sources.pe_overlap import _COMP, _q_rows
from sequencing_error_model.spec import Component, ErrorModelSpec

BASES, RESIDUES = "ACGT", "ACDEFGHIKLMNPQRSTVWY"
_TABLE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {"".join(c): a for c, a in zip(itertools.product("TCAG", repeat=3), _TABLE, strict=True)}
STOPS = frozenset(c for c, a in CODE.items() if a == "*")
# HMMER's background residue frequencies, which kfpfam's log-odds profiles are against.
BACKGROUND = np.array(
    [0.0787945, 0.0151600, 0.0535222, 0.0668298, 0.0397062, 0.0695071, 0.0229198, 0.0590092, 0.0594422, 0.0963728]
    + [0.0237718, 0.0414386, 0.0482904, 0.0395639, 0.0540978, 0.0683364, 0.0540687, 0.0673417, 0.0114135, 0.0304133]
)
GAP_OPEN, GAP_EXTEND = 6.0, 1.0  # bits, as kfpfam's banded score
UNALIGNED = -100.0  # a stop, or a codon off the model: no run of codons crosses it
PLATEAU = (25, 30, 35, 40, 45)
LAMBDAS = (0.1, 0.3, 1.0, 3.0, 10.0, 30.0)
_CIGAR = re.compile(r"(\d+)([MID])")
_SLAB = np.logspace(-3.5, -1, 11)  # a real stop rate per codon, log-uniform
SPREAD = 0.3  # relative sd of a family's error stops about head E's expectation
UNIDENTIFIED_BY_STOPS = ("A>G", "T>C")
FLANK = 12  # bases each side that make a locus
CEILING = 0.01  # the highest error rate per base at which a recurring stop is still read as biology

ReadPair = tuple[str, Sequence[int], str, Sequence[int]]
Column = tuple[int, int]  # model, match column (0-based)
Exposure = tuple[Key, tuple[int, int, int]]  # a head E row, the Qs of its codon


@dataclass(frozen=True)
class Profiles:
    """kfpfam's `data/prep`: per-column log-odds bits [states, 20+] over `RESIDUES`, models concatenated."""

    names: tuple[str, ...]
    offsets: Array  # [models + 1]
    scores: Array

    @classmethod
    def load(cls, prep: Path) -> "Profiles":
        names = tuple((prep / "accs.txt").read_text().split())
        offsets = np.fromfile(prep / "offs.u32", np.uint32).astype(np.int64)
        scores = np.memmap(prep / "profile.f32", np.float32, "r").reshape(-1, 21)
        if len(offsets) != len(names) + 1 or offsets[-1] != len(scores):
            raise ValueError(f"{prep}: accs.txt, offs.u32 and profile.f32 disagree on the models")
        return cls(names, offsets, scores)

    def length(self, model: int) -> int:
        return int(self.offsets[model + 1] - self.offsets[model])


@dataclass(frozen=True)
class Hit:
    """A HMMER hit of one ORF, on the strand it was translated from: alignment start (bases) and match column."""

    model: int
    reverse: bool
    start: int  # 0-based first base of the aligned stretch, on the hit's own strand
    column: int  # 0-based match column of the aligned stretch's first M
    cigar: tuple[tuple[int, str], ...]  # M: residue and column, I: residue only, D: column only


@dataclass(frozen=True)
class Codon:
    """An exposed codon: the read positions of its bases in coding order, forward-read coordinates."""

    column: Column
    positions: tuple[int, int, int]
    reverse: bool
    bases: str
    anchor: float


@dataclass
class Stats:
    reads: int = 0
    placed: int = 0
    codons: int = 0  # exposed at t1
    stops: int = 0
    overlaps: int = 0  # pairs whose mates share counted codons
    overlap_rejected: int = 0  # proposed by the columns, not confirmed by the bases outside
    shared_codons: int = 0
    agree: int = 0
    disagree: int = 0
    stop_one_mate: int = 0  # a sequencing error
    stop_both_mates: int = 0  # biology or a library error: the Q-free floor


@dataclass(frozen=True)
class Stop:
    """One stop codon: the candidate source codons with the head E rows and head Q rows each implies."""

    model: int
    column: Column
    observed: str
    sources: tuple[str, ...]
    erred: tuple[Key, ...]  # the substituted base's row, per source
    others: tuple[tuple[Key, ...], ...]  # the source's remaining opportunities, all "="
    quality: tuple[tuple[tuple[Key, ...], tuple[Key, ...]], ...]  # head Q rows with the source, and as observed
    anchor: int  # whole bits
    flanks: tuple[str, str] | None  # `FLANK` bases before and after, coding strand; None at a read end


@dataclass
class Evidence:
    m: int
    flank: tuple[int, int]
    t1: float
    # rows of bases that did not err, by family: the stop anchor's exposure, then the overlap anchor's
    # (with the codon's three Qs, which say how likely it is to be free of errors)
    certain: dict[int, Counter[Exposure]] = field(default_factory=lambda: defaultdict(Counter))
    overlap: dict[int, Counter[Key]] = field(default_factory=lambda: defaultdict(Counter))
    shared: Counter[int] = field(default_factory=Counter)  # codons both mates expose, per family
    stops: list[Stop] = field(default_factory=list)
    disputed: list[tuple[pe_overlap.Dispute, Column, str, str]] = field(default_factory=list)
    stack: dict[Column, Counter[str]] = field(default_factory=lambda: defaultdict(Counter))
    codons: Counter[int] = field(default_factory=Counter)  # exposed codons per family
    quality: Counter[Key] = field(default_factory=Counter)
    alphabet: set[int] = field(default_factory=set)
    stats: Stats = field(default_factory=Stats)
    # (stops, opportunities) by family and whole-bit anchor score: the plateau over t1
    by_anchor: dict[tuple[int, int], list[int]] = field(default_factory=dict)
    replication: dict[str, int] | None = None  # `drop_replicated`'s report

    def crude_rate(self, t1: float, weights: Mapping[int, float] | None = None) -> tuple[float, float]:
        """(substitutions per base if all were as likely, stops), anchor bar `t1`, families by `weights` (1)."""
        kept = [(v, (weights or {}).get(f, 1.0)) for (f, bits), v in self.by_anchor.items() if bits >= t1]
        stops, chance = (sum(v[i] * x for v, x in kept) for i in (0, 1))
        return (3 * stops / chance if chance else math.nan), stops

    def plateau(self, weights: Mapping[int, float] | None = None) -> dict[str, dict[str, float]]:
        low = min((bits for _, bits in self.by_anchor), default=0)
        rates = {t: self.crude_rate(t, weights) for t in PLATEAU if t >= low}
        return {str(t): {"rate": r, "stops": n} for t, (r, n) in rates.items()}

    def check_plateau(self, weights: Mapping[int, float] | None = None) -> str | None:
        """Why `t1` is too low, if a bar 10 bits higher gives a significantly lower rate (chance homology)."""
        (low, n_low), (high, n_high) = self.crude_rate(self.t1, weights), self.crude_rate(self.t1 + 10, weights)
        if n_low > n_high > 0:
            sd = math.sqrt(low**2 / n_low + high**2 / n_high)
            if low > 1.3 * high and low - high > 3 * sd:
                return f"the stop rate falls from {low:.2e} to {high:.2e} between t1 {self.t1:g} and {self.t1 + 10:g}"
        return None


def fields(m: int) -> tuple[str, ...]:
    """Row fields, in `pe_overlap` order without strand and GC, plus the ops a row can show."""
    window = (f"q{s}{o}" for o in range(1, m + 1) for s in "-+")
    return ("q", *window, "context", "pos_start", "pos_end", "mate", "ops", "op")


def add_errors(seq: str, rate: float, seed: Sequence[int]) -> str:
    """Uniform substitutions at `rate` (standard additions); `seed` names the read, so any step can redo them."""
    rng = np.random.default_rng(list(seed))
    hit = np.flatnonzero(rng.random(len(seq)) < rate)
    shift = rng.integers(1, 4, len(hit))
    bases = list(seq)
    for i, s in zip(hit, shift, strict=True):
        if bases[i] in _COMP:
            bases[i] = BASES[(BASES.index(bases[i]) + int(s)) % 4]
    return "".join(bases)


def orfs(read: int, seq: str, min_length: int = 15) -> Iterator[tuple[str, str]]:
    """(name, protein) of the six-frame stretches between stops; the name carries what `read_hits` needs back."""
    for strand, s in (("+", seq.upper()), ("-", _revcomp(seq.upper()))):
        for frame in range(3):
            at = frame
            for p in "".join(CODE.get(s[i : i + 3], "X") for i in range(frame, len(s) - 2, 3)).split("*"):
                if len(p) >= min_length:
                    yield f"R{read}_{strand}_{at}", p
                at += 3 * (len(p) + 1)


def read_hits(
    path: Path, profiles: Profiles, panel: Collection[str] | None = None, exclude: Collection[str] = ()
) -> dict[int, list[Hit]]:
    """Hits per read from kfpfam's `recruit/pairs.py` table, run on this module's `orfs` (which name the reads)."""
    model = {a: i for i, a in enumerate(profiles.names)}
    out: dict[int, list[Hit]] = defaultdict(list)
    with path.open() as f:
        col = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
        if missing := {"orf", "acc", "hmm_from", "ali_from", "cigar"} - set(col):
            raise ValueError(f"{path}: not a pairs.tsv, no columns {sorted(missing)}")
        for line in f:
            c = line.rstrip("\n").split("\t")
            acc = c[col["acc"]]
            if acc not in model:
                raise ValueError(f"{path}: {acc} is not in the profiles, so they are from another Pfam release")
            if acc in exclude or (panel is not None and acc not in panel):
                continue
            name = re.fullmatch(r"R(\d+)_([+-])_(\d+)", c[col["orf"]])
            if name is None:
                raise ValueError(f"{path}: ORF {c[col['orf']]!r} was not named by this module's `orfs`")
            ops = tuple((int(n), o) for n, o in _CIGAR.findall(c[col["cigar"]]))
            out[int(name[1])].append(
                Hit(
                    model[acc],
                    name[2] == "-",
                    int(name[3]) + 3 * (int(c[col["ali_from"]]) - 1),
                    int(c[col["hmm_from"]]) - 1,
                    ops,
                )
            )
    return dict(out)


def _runs(score: Array) -> tuple[Array, Array]:
    """Best local run of scores strictly left, and strictly right, of each codon."""
    out = []
    for s in (score, score[::-1]):
        best = np.zeros(len(s))
        run = top = 0.0
        for i, x in enumerate(s):
            best[i] = top
            run = max(0.0, run + float(x))
            top = max(top, run)
        out.append(best)
    return out[0], out[1][::-1]


def _path(hit: Hit, n: int) -> tuple[Array, Array, Array]:
    """Column per codon of the frame, gap bits paid before each codon, and which codons are insertions.

    Codons beyond the alignment continue its first and last diagonals.
    """
    columns, gaps, inserted = np.full(n, -1, np.int64), np.zeros(n), np.zeros(n, bool)
    i, k, gap = hit.start // 3, hit.column, 0.0
    for length, op in hit.cigar:
        if op == "M":
            columns[i : i + length] = np.arange(k, k + length)[: max(n - i, 0)]
            if i < n:
                gaps[i] += gap
            i, k, gap = i + length, k + length, 0.0
        elif op == "I":
            inserted[i : i + length] = True
            gaps[i : i + length] = -GAP_EXTEND
            if i < n:
                gaps[i] = -GAP_OPEN
            i += length
        else:
            k, gap = k + length, gap - GAP_OPEN - GAP_EXTEND * (length - 1)
    placed = np.flatnonzero(columns >= 0)
    if len(placed):
        head, tail = int(placed[0]), int(placed[-1])
        columns[:head] = columns[head] + np.arange(-head, 0)
        columns[tail + 1 :] = columns[tail] + np.arange(1, n - tail)
    return columns, gaps, inserted


def exposed(seq: str, hits: Iterable[Hit], profiles: Profiles, t1: float, t2: float) -> list[Codon]:
    """The read's exposed codons, at most one frame per base, best anchor first."""
    found: list[Codon] = []
    strands = {False: seq.upper(), True: _revcomp(seq.upper())}
    for hit in hits:
        st, frame = strands[hit.reverse], hit.start % 3
        n, length, offset = (len(st) - frame) // 3, profiles.length(hit.model), profiles.offsets[hit.model]
        triplets = [st[frame + 3 * i : frame + 3 * i + 3] for i in range(n)]
        residues = [CODE.get(t, "X") for t in triplets]
        letter, stop = np.array([RESIDUES.find(a) for a in residues]), np.array([a == "*" for a in residues], bool)
        scores = profiles.scores[offset : offset + length]

        # bound now: the closure is used only inside this iteration
        def along(cols: Array, letter: Array = letter, stop: Array = stop, scores: Array = scores) -> Array:
            on = (cols >= 0) & (cols < len(scores)) & ~stop
            known = on & (letter >= 0)  # an unknown residue on the model scores 0
            out = np.where(on, 0.0, UNALIGNED)
            out[known] = scores[cols[known], letter[known]]
            return out

        columns, gaps, inserted = _path(hit, n)
        left, right = _runs(np.where(inserted & ~stop, 0.0, along(columns)) + gaps)
        sides: dict[int, tuple[Array, Array]] = {}
        on = ~inserted & (columns >= 0) & (columns < length)
        for i in map(int, np.flatnonzero(on & ((left >= t1) | (right >= t1)))):
            if residues[i] == "X":
                continue
            d = int(columns[i]) - i
            if d not in sides:
                sides[d] = _runs(along(np.arange(n) + d))
            # The anchor follows the alignment; the far side must continue on this codon's own diagonal, since a
            # stop here would have ended the ORF the alignment was made from.
            anchor = max(left[i] if sides[d][1][i] >= t2 else 0.0, right[i] if sides[d][0][i] >= t2 else 0.0)
            if anchor >= t1:
                at = [frame + 3 * i + p for p in range(3)]
                at = [len(seq) - 1 - x for x in at] if hit.reverse else at
                column = (hit.model, int(columns[i]))
                found.append(Codon(column, (at[0], at[1], at[2]), hit.reverse, triplets[i], float(anchor)))
    found.sort(key=lambda c: (-c.anchor, c.column, c.positions))
    taken: set[int] = set()
    out = []
    for c in found:
        if not taken.intersection(c.positions):
            taken.update(c.positions)
            out.append(c)
    return out


def changes(codon: str) -> list[tuple[int, str, str]]:
    """(position, base, result) of the nine single-base changes of `codon`."""
    return [(p, b, codon[:p] + b + codon[p + 1 :]) for p in range(3) for b in BASES if b != codon[p]]


_OPPORTUNITY = {
    c: tuple("".join(b for p, b, to in changes(c) if p == at and to in STOPS) for at in range(3)) for c in CODE
}
_SOURCES = {s: tuple((p, to) for p, _, to in changes(s) if to not in STOPS) for s in STOPS}


@dataclass(frozen=True)
class _Read:
    seq: str
    q: list[int]
    mate: int


def _row(
    read: _Read, position: int, centre: str, observed: str, ops: str, m: int, flank: tuple[int, int], seq: str = ""
) -> Key:
    """One read-oriented row; `seq` replaces the read's bases as context (an overlap's consensus)."""
    s, n, (left, right) = seq or read.seq, len(read.seq), flank
    ctx = "".join(s[position + o] if 0 <= position + o < n else "." for o in range(-left, right + 1))
    window = tuple(read.q[position + o] if 0 <= position + o < n else None for k in range(1, m + 1) for o in (-k, k))
    context = ctx[:left] + centre + ctx[left + 1 :]
    op = "=" if centre == observed else f"{centre}>{observed}"
    return (read.q[position], *window, context, position + 1, n - position, read.mate, ops, op)


def _oriented(bases: str, reverse: bool) -> str:
    return "".join(sorted(_COMP[b] if reverse else b for b in bases))


def _stop_rows(read: _Read, codon: Codon, source: str, m: int, flank: tuple[int, int]) -> dict[int, Key]:
    """Rows of the codon's bases that could have made a stop from `source`, by codon position."""
    out = {}
    for p, ops in enumerate(_OPPORTUNITY[source]):
        if ops:
            true, seen = (b[p] if not codon.reverse else _COMP[b[p]] for b in (source, codon.bases))
            out[p] = _row(read, codon.positions[p], true, seen, _oriented(ops, codon.reverse), m, flank)
    return out


def flanks(seq: str, codon: Codon) -> tuple[str, str] | None:
    """The `FLANK` bases either side of a codon, on its strand; None where the read ends first."""
    strand = _revcomp(seq) if codon.reverse else seq
    at = len(seq) - 1 - codon.positions[0] if codon.reverse else codon.positions[0]
    if at < FLANK or at + 3 + FLANK > len(seq):
        return None
    return strand[at - FLANK : at], strand[at + 3 : at + 3 + FLANK]


Locus = tuple[int, str, tuple[str, str] | None]  # model, codon, flanks


def _find(strand: str, left: str, right: str) -> list[int]:
    """Where a codon lies between the flanks: one of them exact, the other with at most two mismatches.

    The read's own errors fall in its flanks too, and must not leave its stop looking alone.
    """
    out = []
    for exact, start, other, codon in ((left, FLANK + 3, right, FLANK), (right, -FLANK - 3, left, -3)):
        i = strand.find(exact)
        while i >= 0:
            lo = i + start
            if lo >= 0 and lo + FLANK <= len(strand):
                near = sum(a != b for a, b in zip(strand[lo : lo + FLANK], other, strict=True)) <= 2
                if near and i + codon not in out:
                    out.append(i + codon)
            i = strand.find(exact, i + 1)
    return out


def replicates(
    loci: Sequence[Locus], pairs: Iterable[ReadPair], hits: Mapping[int, Sequence[Hit]]
) -> list[tuple[int, int]]:
    """Fragments with the same codon, and with another codon, between each locus' flanks: its own is counted.

    Only reads with a hit to the locus' family are searched, which is what makes this cheap. Mates count once.
    """
    by_model: dict[int, list[int]] = defaultdict(list)
    for k, (model, _, around) in enumerate(loci):
        if around:
            by_model[model].append(k)
    out = [[0, 0] for _ in loci]
    for number, (s1, _, s2, _) in enumerate(pairs):
        seen: dict[int, bool] = {}
        for mate, s in enumerate((s1, s2)):
            wanted = [k for f in {h.model for h in hits.get(2 * number + mate, ())} for k in by_model.get(f, ())]
            if not wanted:
                continue
            strands = (s.upper(), _revcomp(s.upper()))
            for k in wanted:
                _, codon, around = loci[k]
                assert around is not None
                for st in strands:
                    if at := _find(st, *around):
                        seen[k] = seen.get(k, False) or any(st[c : c + 3] == codon for c in at)
        for k, same in seen.items():
            out[k][0 if same else 1] += 1
    return [(a, b) for a, b in out]


def recurs(same: int, other: int) -> bool:
    """Do `same` fragments of `same + other` share a stop more often than errors at `CEILING` would (p < 0.01)?

    ponytail: PCR copies of one fragment count as replicates, so a library error in them is dropped too.
    """
    return same >= 2 and float(special.gammainc(same - 1, (same + other - 1) * CEILING / 3)) < 0.01


def drop_replicated(ev: Evidence, pairs: Iterable[ReadPair], hits: Mapping[int, Sequence[Hit]]) -> dict[str, int]:
    """Remove the stops that recur at their locus (biology) from `ev`, and say what was found."""
    n = replicates([(s.model, s.observed, s.flanks) for s in ev.stops], pairs, hits)
    found = {"stops": len(n), "replicated": 0, "contradicted": 0, "alone": 0}
    kept = []
    for s, (same, other) in zip(ev.stops, n, strict=True):
        if recurs(same, other):
            found["replicated"] += 1
            ev.by_anchor[s.model, s.anchor][0] -= 1
            ev.codons[s.model] -= 1
            continue
        found["contradicted" if other else "alone"] += 1
        kept.append(s)
    ev.stops = kept
    ev.replication = found
    return found


def _stop(read: _Read, codon: Codon, m: int, flank: tuple[int, int]) -> Stop:
    sources, erred, others, q = [], [], [], []
    for p, source in _SOURCES[codon.bases]:
        rows = _stop_rows(read, codon, source, m, flank)
        erred.append(rows.pop(p))
        sources.append(source)
        others.append(tuple(rows.values()))
        base = _COMP[source[p]] if codon.reverse else source[p]
        seen = read.seq[codon.positions[p]]
        q.append(tuple(_q_rows(read.seq, read.q, read.mate, codon.positions[p], b, m, flank) for b in (base, seen)))
    return Stop(
        codon.column[0],
        codon.column,
        codon.bases,
        tuple(sources),
        tuple(erred),
        tuple(others),
        cast("tuple[tuple[tuple[Key, ...], tuple[Key, ...]], ...]", tuple(q)),
        int(codon.anchor),
        flanks(read.seq, codon),
    )


def _confirmed(a: _Read, b: _Read, s: int, inside: set[int], min_outside: int, max_outside: float) -> bool:
    """Do the mates agree outside the counted bases, with reverse-complemented `b`[k] at `a`[s + k]?"""
    rc = _revcomp(b.seq)
    sites = [
        (a.seq[i], rc[i - s])
        for i in range(max(0, s), min(len(a.seq), s + len(rc)))
        if i not in inside and a.seq[i] in _COMP and rc[i - s] in _COMP
    ]
    return len(sites) >= min_outside and sum(x != y for x, y in sites) <= max_outside * len(sites)


def _overlap(
    ev: Evidence,
    mates: tuple[_Read, _Read],
    found: tuple[list[Codon], list[Codon]],
    min_outside: int,
    max_outside: float,
) -> tuple[set[Codon], set[Codon]]:
    """Add the pair's overlap rows; return the codons of each mate they use, which the stop anchor must not reuse."""
    used: tuple[set[Codon], set[Codon]] = (set(), set())
    n2 = len(mates[1].seq)
    shared: dict[int, list[tuple[Codon, Codon]]] = defaultdict(list)
    second = {(c.column, c.reverse): c for c in found[1]}
    for c in found[0]:
        partner = second.get((c.column, not c.reverse))
        if partner is not None:
            # offset of reverse-complemented mate 2 on mate 1: base i pairs with n2 - 1 - (i - s)
            shared[c.positions[1] + partner.positions[1] - (n2 - 1)].append((c, partner))
    if not shared:
        return used
    s, pairs = max(shared.items(), key=lambda x: (len(x[1]), -abs(x[0])))
    if pe_overlap.place(mates[0].seq, mates[1].seq) != (s, False) or not _confirmed(
        *mates, s, {p for c, _ in pairs for p in c.positions}, min_outside, max_outside
    ):
        ev.stats.overlap_rejected += 1
        return used
    ev.stats.overlaps += 1
    # the mates' consensus as context, N where they differ, as `pe_overlap`
    est = [list(r.seq) for r in mates]
    for i in range(max(0, s), min(len(mates[0].seq), s + n2)):
        j = n2 - 1 - (i - s)
        if mates[0].seq[i] != _COMP.get(mates[1].seq[j], "N"):
            est[0][i] = est[1][j] = "N"
    texts = ["".join(e) for e in est]
    (r1, r2), (m, flank) = mates, (ev.m, ev.flank)
    for c1, c2 in pairs:
        used[0].add(c1)
        used[1].add(c2)
        ev.stats.shared_codons += 1
        ev.shared[c1.column[0]] += 1
        stops = sum(c.bases in STOPS for c in (c1, c2))
        ev.stats.stop_one_mate += stops == 1
        ev.stats.stop_both_mates += stops == 2
        if stops == 2:
            continue
        for p in range(3):
            i, j = c1.positions[p], c2.positions[p]
            b1, o2 = r1.seq[i], r2.seq[j]
            b2 = _COMP[o2]
            if b1 == b2:
                ev.stats.agree += 1
                ev.overlap[c1.column[0]][_row(r1, i, b1, b1, BASES, m, flank, texts[0])] += 1
                ev.overlap[c1.column[0]][_row(r2, j, o2, o2, BASES, m, flank, texts[1])] += 1
                continue
            ev.stats.disagree += 1
            first = (
                _row(r1, i, b1, b1, BASES, m, flank, texts[0]),
                _row(r2, j, _COMP[b1], o2, BASES, m, flank, texts[1]),
            )
            other = (
                _row(r1, i, b2, b1, BASES, m, flank, texts[0]),
                _row(r2, j, o2, o2, BASES, m, flank, texts[1]),
            )
            q_first = _q_rows(r1.seq, r1.q, 1, i, b1, m, flank) + _q_rows(r2.seq, r2.q, 2, j, _COMP[b1], m, flank)
            q_other = _q_rows(r1.seq, r1.q, 1, i, b2, m, flank) + _q_rows(r2.seq, r2.q, 2, j, o2, m, flank)
            # the codon if mate 1 read this base, or if mate 2 did; the other two positions as mate 1 has them
            right = c2.bases[p]
            ev.disputed.append(
                ((*first, *other, q_first, q_other), c1.column, c1.bases, c1.bases[:p] + right + c1.bases[p + 1 :])
            )
    return used


def collect(
    pairs: Iterable[ReadPair],
    hits: Mapping[int, Sequence[Hit]],
    profiles: Profiles,
    m: int = 1,
    flank: tuple[int, int] = (1, 1),
    t1: float = 35.0,
    t2: float = 6.0,
    min_outside: int = 20,
    max_outside: float = 0.05,
    quality_every: int = 1,
) -> Evidence:
    """Evidence from read pairs; read ids are 2 x pair for mate 1 and 2 x pair + 1 for mate 2, as `orfs` names them.

    A single-end library passes empty mate 2 reads.
    """
    ev = Evidence(m, flank, t1)
    low = min(t1, PLATEAU[0])
    # ponytail: Python per codon, about a minute per million placed reads; move exposure to numpy if real runs need it
    for number, (s1, q1, s2, q2) in enumerate(pairs):
        mates = (_Read(s1.upper(), list(q1), 1), _Read(s2.upper(), list(q2), 2))
        found: tuple[list[Codon], list[Codon]] = ([], [])
        for k, read in enumerate(mates):
            if not read.seq:
                continue
            ev.stats.reads += 1
            ev.alphabet.update(read.q)
            if number % quality_every == 0:
                ev.quality.update(fastq_quality.quality_rows(read.seq, read.q, read.mate, m, flank))
            if placed := hits.get(2 * number + k):
                for c in exposed(read.seq, placed, profiles, low, t2):
                    tally = ev.by_anchor.setdefault((c.column[0], int(c.anchor)), [0, 0])
                    tally[0] += c.bases in STOPS
                    tally[1] += 0 if c.bases in STOPS else sum(map(len, _OPPORTUNITY[c.bases]))
                    if c.anchor >= t1:
                        found[k].append(c)
                ev.stats.placed += bool(found[k])
        used: tuple[set[Codon], set[Codon]] = (set(), set())
        if found[0] and found[1]:
            used = _overlap(ev, mates, found, min_outside, max_outside)
        for read, codons, done in zip(mates, found, used, strict=True):
            for c in codons:
                ev.stats.codons += 1
                if c.bases not in STOPS:
                    ev.stack[c.column][c.bases] += 1
                if c in done:
                    continue
                ev.codons[c.column[0]] += 1
                if c.bases in STOPS:
                    ev.stats.stops += 1
                    ev.stops.append(_stop(read, c, m, flank))
                else:
                    q3 = (read.q[c.positions[0]], read.q[c.positions[1]], read.q[c.positions[2]])
                    ev.certain[c.column[0]].update((k, q3) for k in _stop_rows(read, c, c.bases, m, flank).values())
    return ev


class Prior:
    """Codon prior per column: the sample's stack, plus `lam` codons of the HMM's residues x the sample's usage."""

    def __init__(self, ev: Evidence, profiles: Profiles, lam: float | None = None) -> None:
        self._stack, self._profiles = ev.stack, profiles
        usage: Counter[str] = Counter()
        for c in ev.stack.values():
            usage.update(c)
        by_residue: Counter[str] = Counter()
        for codon, n in usage.items():
            by_residue[CODE[codon]] += n
        # +1 so a residue the sample never shows still has its codons
        synonyms = Counter(a for a in CODE.values())
        self._within = {c: (usage[c] + 1) / (by_residue[a] + synonyms[a]) for c, a in CODE.items() if a != "*"}
        self.lam = lam if lam is not None else self._choose()

    def _base(self, column: Column, codon: str) -> float:
        row = self._profiles.offsets[column[0]] + column[1]
        k = RESIDUES.index(CODE[codon])
        return float(2.0 ** self._profiles.scores[row, k] * BACKGROUND[k]) * self._within[codon]

    def _choose(self) -> float:
        """The `lam` with the best leave-one-out likelihood of the sense codons in columns seen more than once."""
        cells = [
            (n, sum(c.values()), self._base(col, codon))
            for col, c in itertools.islice(((k, v) for k, v in self._stack.items() if v.total() > 1), 20_000)
            for codon, n in c.items()
        ]
        if not cells:
            return 1.0
        n, total, base = (np.array(x, float) for x in zip(*cells, strict=True))
        return max(LAMBDAS, key=lambda lam: float(n @ np.log((n - 1 + lam * base) / (total - 1 + lam) + 1e-300)))

    def weight(self, column: Column, codon: str, drop: int = 0) -> float:
        """Unnormalised P(codon) at the column; `drop` removes that many copies (the read's own)."""
        if codon in STOPS:
            return 1e-6  # a stop the mates share would be biology; as a source it is all but excluded
        seen = self._stack[column][codon] if column in self._stack else 0
        return max(seen - drop, 0) + self.lam * self._base(column, codon)


def family_weights(
    stops: Mapping[int, float], codons: Mapping[int, float], errors: Mapping[int, float]
) -> dict[int, float]:
    """P(a family's stops are all errors), by a spike-and-slab on its stops against its exposure.

    `errors` is each family's expected error stops, from head E over its own rows: families differ in codon usage,
    Q and position, so an average rate would flag them for that. Counts are negative binomial with a relative
    spread of `SPREAD` in the error stops, so a fit that is off by that much for one family is not a real stop rate.
    The chance that a family has real stops is estimated from the families together; one without stops has weight 1.
    """
    families = [f for f, n in stops.items() if n > 0]
    if not families:
        return {}
    n = np.array([stops[f] for f in families], float)
    expected = np.maximum([errors[f] for f in families], 1e-12)
    size = np.array([codons[f] for f in families], float)

    def log_nb(mean: Array) -> Array:
        # error stops vary by SPREAD of their mean; the variance of real stops is Poisson on top, close enough here
        shape = (mean / (SPREAD * np.maximum(expected, 1e-12).reshape(-1, *[1] * (mean.ndim - 1)))) ** 2
        k = n.reshape(-1, *[1] * (mean.ndim - 1))
        return np.asarray(
            special.gammaln(k + shape)
            - special.gammaln(shape)
            + shape * np.log(shape / (shape + mean))
            + k * np.log(mean / (shape + mean))
        )

    log_null = log_nb(expected)
    log_real = special.logsumexp(log_nb(expected[:, None] + size[:, None] * _SLAB), axis=1) - math.log(len(_SLAB))
    share, real = 0.01, np.zeros(len(families))
    for _ in range(200):
        # families without stops are, to a good approximation, all in the spike
        real = special.expit(math.log(share / (1 - share)) + log_real - log_null)
        share = min(max(float(real.sum()) / max(len(codons), 1), 1e-9), 0.5)
    return {f: float(1 - r) for f, r in zip(families, real, strict=True)}


def _plain(key: Key) -> Key:
    """The row as if every substitution could be seen: compares the candidates on their own probabilities."""
    return (*key[:-2], BASES, key[-1])


def _probabilities(ev: Evidence, keys: Iterable[Key], head: Sequence[Component]) -> dict[Key, float]:
    unique = list(dict.fromkeys(keys))
    if not unique:
        return {}
    t = CountTable("pfam", fields(ev.m), "base", True, Counter(dict.fromkeys(unique, 1)), {"flank": ev.flank})
    p = error.probabilities(head, tuple(sorted(ev.alphabet)), t)
    ci, left = fields(ev.m).index("context"), ev.flank[0]
    return {k: float(p[r, error._category(k[-1], str(k[ci])[left])]) for r, k in enumerate(unique)}


def _quality_scores(ev: Evidence, keys: Iterable[Key], head: Sequence[Component]) -> dict[Key, float]:
    unique = list(dict.fromkeys(keys))
    if not unique:
        return {}
    t = CountTable(
        "pfam:quality",
        fastq_quality.quality_fields(ev.m),
        "base",
        False,
        Counter(dict.fromkeys(unique, 1)),
        {"flank": ev.flank},
    )
    return dict(zip(unique, map(float, quality.log_prob(head, tuple(sorted(ev.alphabet)), t)), strict=True))


def _expected_stops(ev: Evidence, head: Sequence[Component], families: Collection[int]) -> dict[int, float]:
    """Error stops head E expects in each family's exposed sense codons."""
    unique: dict[Key, int] = {}
    for f in families:
        for k, _ in ev.certain.get(f, ()):
            unique.setdefault(k[:-1], len(unique))
    if not unique:
        return dict.fromkeys(families, 0.0)
    t = CountTable("pfam", fields(ev.m)[:-1], "base", False, Counter(dict.fromkeys(unique, 1)), {"flank": ev.flank})
    p = 1 - error.probabilities(head, tuple(sorted(ev.alphabet)), t)[:, 0]
    return {f: sum(n * float(p[unique[k[:-1]]]) for (k, _), n in ev.certain.get(f, {}).items()) for f in families}


def _clean(ev: Evidence, head: Sequence[Component], limit: int = 100_000) -> dict[int, float]:
    """P(a base is not a substitution) by its reported Q, under head E over the exposed bases.

    A codon is one trial of "unchanged, or one change that makes a stop"; conditional on those outcomes the odds
    are the substitution's own. A sense codon read with an error is neither, but looks like a trial, and errors
    cluster in low-Q runs, so counting it biases the rate down (by 9 % at a 1 % error rate; plan P1).
    Q is head E's covariate here, with its calibration learned.
    """
    rows = itertools.islice((k for c in ev.certain.values() for k, _ in c), limit)
    t = CountTable("pfam", fields(ev.m), "base", True, Counter(map(_plain, rows)), {"flank": ev.flank})
    if not t.counts:
        return {}
    p = error.probabilities(head, tuple(sorted(ev.alphabet)), t)[:, 0]
    q, n = np.array([k[0] for k in t.counts]), np.array(list(t.counts.values()), float)
    return {int(x): float(np.average(p[q == x], weights=n[q == x])) for x in np.unique(q)}


def _placed(ev: Evidence, head: Sequence[Component]) -> dict[int, float]:
    """P(a family's overlaps are placed right): its disagreements against the errors head E expects in both mates.

    Tandem repeats put two loci on the same columns, so their mates disagree at many bases of a few families.
    """
    counts = Counter(column[0] for _, column, *_ in ev.disputed)
    keys = list(dict.fromkeys(k for f in counts for k in ev.overlap[f]))
    p = _probabilities(ev, keys, head)
    errors = {f: sum(n * (1 - p[k]) for k, n in ev.overlap[f].items()) for f in counts}
    return family_weights(counts, ev.shared, errors)


@dataclass
class _Weights:
    """The E step's output: weights per family, stop source, disputed base and exposed codon."""

    family: dict[int, float]
    source: list[Array]
    first: Array
    clean: dict[int, float]  # P(a base is not a substitution) by its Q
    placed: dict[int, float]  # P(a family's overlaps are placed right), from its disagreements


class _Expectation:
    """E step: attribute stops to source codons and disputes to a mate, under head E, head Q and the stack."""

    def __init__(self, ev: Evidence, prior: Prior, quality_head: Sequence[Component]) -> None:
        self.ev, self.codon_weight = ev, [np.array([prior.weight(s.column, c) for c in s.sources]) for s in ev.stops]
        keys = [k for s in ev.stops for pair in s.quality for rows in pair for k in rows]
        keys += [k for d, *_ in ev.disputed for rows in d[4:] for k in rows]
        score = _quality_scores(ev, keys, quality_head)
        self.stop_quality = [
            np.array([sum(score[k] for k in a) - sum(score[k] for k in b) for a, b in s.quality]) for s in ev.stops
        ]
        self.dispute_prior = np.array(
            [
                math.log(prior.weight(column, one, 1) / prior.weight(column, two, 1))
                + sum(score[k] for k in d[4])
                - sum(score[k] for k in d[5])
                for d, column, one, two in ev.disputed
            ]
        )

    def start(self) -> _Weights:
        source = [w / w.sum() for w in self.codon_weight]
        return _Weights({}, source, np.full(len(self.ev.disputed), 0.5), {}, {})

    def update(self, head: Sequence[Component]) -> _Weights:
        ev = self.ev
        keys = [_plain(k) for s in ev.stops for k in s.erred] + [_plain(k) for d, *_ in ev.disputed for k in d[:4]]
        p = _probabilities(ev, keys, head)
        source = []
        for s, prior, q in zip(ev.stops, self.codon_weight, self.stop_quality, strict=True):
            log = np.log(prior) + np.log([p[_plain(k)] + 1e-300 for k in s.erred]) + q
            source.append(np.asarray(special.softmax(log)))
        counts = Counter(s.model for s in ev.stops)
        family = family_weights(counts, ev.codons, _expected_stops(ev, head, counts)) if counts else {}
        like = [
            math.log(p[_plain(d[0])] * p[_plain(d[1])] + 1e-300) - math.log(p[_plain(d[2])] * p[_plain(d[3])] + 1e-300)
            for d, *_ in ev.disputed
        ]
        first = np.asarray(special.expit(np.array(like) + self.dispute_prior))
        return _Weights(family, source, first, _clean(ev, head), _placed(ev, head))


def expected(ev: Evidence, w: _Weights) -> CountTable:
    """Head E rows as expected counts under `w`: the table `fit.error` and `compare` take."""
    counts: dict[Key, float] = defaultdict(float)
    for family, rows in ev.certain.items():
        scale = w.family.get(family, 1.0)
        for (k, (a, b, c)), n in rows.items():
            counts[k] += n * scale * w.clean.get(a, 1.0) * w.clean.get(b, 1.0) * w.clean.get(c, 1.0)
    for family, agreed in ev.overlap.items():
        for k, n in agreed.items():
            counts[k] += n * w.placed.get(family, 1.0)
    for s, share in zip(ev.stops, w.source, strict=True):
        scale = w.family.get(s.model, 1.0)
        for x, key, rest in zip(share, s.erred, s.others, strict=True):
            for k in (key, *rest):
                counts[k] += float(x) * scale
    for (dispute, column, *_), f in zip(ev.disputed, w.first, strict=True):
        scale = w.placed.get(column[0], 1.0)
        for k, x in zip(dispute[:4], (f, f, 1 - f, 1 - f), strict=True):
            counts[k] += float(x) * scale
    kept = cast("Counter[Key]", Counter({k: n for k, n in counts.items() if n > 0}))
    return CountTable("pfam", fields(ev.m), "base", True, kept, {"flank": ev.flank})


def floor(stats: Stats) -> dict[str, float]:
    """The Q-free floor: stops both mates read, per shared codon, with a 95 % Poisson upper bound."""
    n, codons = stats.stop_both_mates, max(stats.shared_codons, 1)
    upper = float(special.gammaincinv(n + 1, 0.95))
    return {"stops_both_mates": n, "shared_codons": stats.shared_codons, "rate": n / codons, "upper95": upper / codons}


def fit(
    ev: Evidence,
    profiles: Profiles,
    error_tokens: Sequence[str],
    quality_tokens: Sequence[str],
    *,
    lam: float | None = None,
    iterations: int = 30,
    tol: float = 1e-3,
    smooth: float = 30.0,
    window_l2: float | None = 0.01,
    seed: int = 0,
    sources: Sequence[str] = ("pfam",),
) -> tuple[ErrorModelSpec, CountTable]:
    """Fit both heads; return the spec and the expected head E rows it was fitted from.

    Head Q first (from every base sampled, conditioned on observed bases, as `pe_overlap`), then head E by EM over the
    latent attributions, starting from the stack alone and never from Q.
    """
    if not ev.stops and not ev.overlap and not ev.disputed:
        raise ValueError(f"no stop or overlap evidence in {ev.stats.codons} exposed codons of {ev.stats.reads} reads")
    alphabet = tuple(sorted(ev.alphabet))
    lags = [f"q-{i}" for i in range(1, Component(quality_tokens[0]).args[0] + 1)]
    q_table = CountTable(
        "pfam:quality", fastq_quality.quality_fields(ev.m), "base", False, ev.quality, {"flank": ev.flank}
    )
    quality_head = quality.fit(
        q_table.marginal(*lags, "pos_start", "pos_end", "mate", "context", "q"), quality_tokens, alphabet
    )
    prior = Prior(ev, profiles, lam)
    step = _Expectation(ev, prior, quality_head)
    w, head = step.start(), None
    for _ in range(max(iterations, 1)):
        head = error.fit(
            expected(ev, w), error_tokens, alphabet, smooth=smooth, window_l2=window_l2, seed=seed, init=head
        )
        new = step.update(head)
        moved = (
            max(
                [float(np.abs(a - b).max()) for a, b in zip(new.source, w.source, strict=True)]
                + [float(np.abs(new.first - w.first).max(initial=0.0))]
                + [abs(new.family.get(f, 1.0) - w.family.get(f, 1.0)) for f in {*new.family, *w.family}]
                + [abs(new.placed.get(f, 1.0) - w.placed.get(f, 1.0)) for f in {*new.placed, *w.placed}]
            )
            if ev.stops or ev.disputed
            else 0.0
        )
        w = new
        if moved < tol:
            break
    assert head is not None
    table = expected(ev, w)
    if why := ev.check_plateau(w.family):
        warnings.warn(f"t1 looks too low: {why}", RuntimeWarning, stacklevel=2)
    bias = head[0].params["bias"].copy()
    bias[5:] = -np.inf  # no indel label here
    head = (Component(head[0].token, {**head[0].params, "bias": bias}, meta=head[0].meta), *head[1:])
    provenance = {
        "mode": "pfam",
        "sources": list(sources),
        "identified_ops": ["substitution"],
        # fitted by extrapolation, from the other ops into the same base, unless mates overlap
        "unidentified_ops": [] if ev.overlap or ev.disputed else list(UNIDENTIFIED_BY_STOPS),
        "stats": asdict(ev.stats),
        "t1": ev.t1,
        "lam": prior.lam,
        "plateau": ev.plateau(w.family),
        "plateau_warning": ev.check_plateau(w.family),
        "floor": floor(ev.stats),
        "replication": ev.replication,
        "downweighted_families": {profiles.names[f]: round(x, 4) for f, x in w.family.items() if x < 0.5},
        "downweighted_overlaps": {profiles.names[f]: round(x, 4) for f, x in w.placed.items() if x < 0.5},
    }
    return ErrorModelSpec(alphabet, provenance, quality_head, head), table


def read_pairs(
    r1: Path, r2: Path | None, limit: int | None = None, every_n: int = 1, dose: float = 0.0
) -> Iterator[ReadPair]:
    """Pairs in file order, strided (never sliced: `observations.subset`), with standard additions at `dose`."""
    second: Iterable[tuple[str, Sequence[int]]] = read_fastq(r2) if r2 else itertools.repeat(("", b""))
    both = subset(zip(read_fastq(r1), second, strict=bool(r2)), limit, every_n)
    for number, ((s1, q1), (s2, q2)) in enumerate(both):
        if dose:
            s1, s2 = (add_errors(s, dose, (number, mate)) for mate, s in enumerate((s1, s2)))
        yield s1, q1, s2, q2


def write_pairs(pairs: Iterable[ReadPair], prefix: Path) -> int:
    """`prefix`_1/2.fastq.gz with each pair named by its number, so an aligner's output maps back to read ids."""
    n = 0
    with gzip.open(f"{prefix}_1.fastq.gz", "wt", 3) as one, gzip.open(f"{prefix}_2.fastq.gz", "wt", 3) as two:
        for n, (s1, q1, s2, q2) in enumerate(pairs, 1):
            for out, s, q in ((one, s1, q1), (two, s2, q2)):
                out.write(f"@{n - 1}\n{s}\n+\n{''.join(chr(v + 33) for v in q)}\n")
    return n


def only(ev: Evidence, anchor: str) -> Evidence:
    """The evidence of one anchor, `stops` or `overlap`: for testing each against both (plan P2).

    Codons that `collect` gave to the overlap are not handed back to the stops.
    """
    if anchor == "overlap":
        return replace(ev, certain=defaultdict(Counter), stops=[], codons=Counter())
    return replace(ev, overlap=defaultdict(Counter), disputed=[], shared=Counter())


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m sequencing_error_model.sources.pfam", description="Fit a spec from reads placed on Pfam."
    )
    sub = p.add_subparsers(dest="command", required=True)
    c = sub.add_parser("reads", help="write the strided, dosed pairs, named by number, for an aligner")
    a = sub.add_parser("orfs", help="write the six-frame ORFs to search with kfpfam and HMMER")
    b = sub.add_parser("fit", help="fit a spec from the search's pairs.tsv")
    for x in (a, b, c):
        x.add_argument("r1", type=Path, help="mate 1 FASTQ[.gz], raw: not merged and not error-corrected")
        x.add_argument("r2", type=Path, nargs="?", help="mate 2 FASTQ[.gz]; omit for a single-end library")
        x.add_argument("--max-pairs", type=int)
        x.add_argument("--every-n", type=int, default=1, metavar="N", help="keep every Nth pair, never a head slice")
        x.add_argument("--dose", type=float, default=0.0, help="add substitutions at this rate; same in every command")
    a.add_argument("--output", type=Path, required=True, metavar="FASTA")
    c.add_argument("--output", type=Path, required=True, metavar="PREFIX", help="PREFIX_1.fastq.gz and PREFIX_2")
    b.add_argument("--pairs", type=Path, required=True, help="kfpfam recruit/pairs.py's table, E <= 1000")
    b.add_argument("--prep", type=Path, required=True, help="kfpfam data/prep of the release searched")
    b.add_argument("--panel", type=Path, help="accessions to use, one per line (default: all)")
    b.add_argument("--exclude", type=Path, help="accessions to leave out, one per line")
    b.add_argument("--output", type=Path, required=True, metavar="SPEC_DIR")
    b.add_argument("--error-tokens", nargs="+", default=["QualityWindow(1)", "Context(0,0)", "Position(4)", "Mate"])
    b.add_argument("--quality-tokens", nargs="+", default=["QualityMarkov(1)", "Position(48)", "Mate", "Context(1,1)"])
    b.add_argument("--t1", type=float, default=35.0, help="anchor bar, raw profile bits (see the plateau it reports)")
    b.add_argument("--t2", type=float, default=6.0)
    b.add_argument("--anchors", choices=("both", "stops", "overlap"), default="both")
    b.add_argument("--quality-every", type=int, default=20, metavar="N", help="head Q from every Nth pair")
    b.add_argument("--iterations", type=int, default=30)
    b.add_argument("--seed", type=int, default=0)
    args = p.parse_args(argv)
    pairs = read_pairs(args.r1, args.r2, args.max_pairs, args.every_n, args.dose)
    if args.command == "reads":
        if args.r2 is None:
            p.error("reads needs both mates")
        # through a staged link too: the output would truncate the reads it is made from
        outputs = {Path(f"{args.output}_{i}.fastq.gz").resolve() for i in (1, 2)}
        if outputs & {args.r1.resolve(), args.r2.resolve()}:
            p.error(f"--output {args.output} would overwrite its own input")
        print(write_pairs(pairs, args.output))
        return 0
    if args.command == "orfs":
        with args.output.open("w") as out:
            for number, (s1, _, s2, _) in enumerate(pairs):
                for mate, seq in enumerate((s1, s2)):
                    out.writelines(f">{name}\n{protein}\n" for name, protein in orfs(2 * number + mate, seq))
        return 0
    profiles = Profiles.load(args.prep)
    lists = [set(x.read_text().split()) if x else None for x in (args.panel, args.exclude)]
    hits = read_hits(args.pairs, profiles, lists[0], lists[1] or ())
    flank = (1, 1)
    m = max(Component(args.error_tokens[0]).args[0], Component(args.quality_tokens[0]).args[0])
    ev = collect(pairs, hits, profiles, m, flank, args.t1, args.t2, quality_every=args.quality_every)
    drop_replicated(ev, read_pairs(args.r1, args.r2, args.max_pairs, args.every_n, args.dose), hits)
    ev = ev if args.anchors == "both" else only(ev, args.anchors)
    spec, _ = fit(
        ev,
        profiles,
        args.error_tokens,
        args.quality_tokens,
        iterations=args.iterations,
        seed=args.seed,
        sources=[str(x) for x in (args.r1, args.r2, args.pairs) if x],
    )
    spec.save(args.output)
    print(json.dumps({k: v for k, v in spec.provenance.items() if k != "sources"}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
