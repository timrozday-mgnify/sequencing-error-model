"""Baseline harness (plan phase 5): the same read-level metrics on real reads and on reads from each simulator.

Every read set is aligned to one reference the same way, and each simulator is trained by its own profiler on the
same training alignment, so a metric's distance from the real reads compares the simulators, not the pipelines.
A `profile` holds, per read set:

- error rows (`generate.observations` of aligned records) by reported Q, by cycle (10-base bins of `pos_start`)
  and by trinucleotide context, split into match / substitution / deletion / insertion;
- edits (substitutions and indel bases) per aligned read, whose spread shows per-read error heterogeneity;
- Q by cycle over every base of every primary read (clips included), and lag-1 Q pairs;
- canonical k-mer multiplicities, and how many k-mer occurrences are absent from the reference.

`distances(real, sim)` gives one number per metric, 0 when identical. Rate metrics are the real-exposure-weighted
mean |log rate ratio| per group; `q_cycle_tv` the weighted mean per-cycle Q TV; `q_lag1` the absolute difference
in lag-1 autocorrelation; `kmer_spectrum_tv` the TV of distinct-k-mer multiplicity histograms (capped at
`_MAX_MULTIPLICITY`); `edits_per_read_tv` the TV of edits-per-read histograms (capped the same way); `kmer_absent` the |log ratio| of reference-absent k-mer occurrence fractions (error k-mers).
`verdict` says whether one simulator beats, matches or trails another per metric.

ponytail: no sampling-noise floor (a real-vs-real split); read sets of similar size keep the comparison fair.
"""

import argparse
import json
import sys
from collections import Counter
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from itertools import islice
from pathlib import Path
from typing import Any

import numpy as np
import pysam
from numpy.lib.stride_tricks import sliding_window_view

from sequencing_error_model.fit.quality import Array
from sequencing_error_model.generate import Read, _edits, observations
from sequencing_error_model.sources import bam

CLASSES = ("match", "substitution", "deletion", "insertion")
_MAX_MULTIPLICITY = 50
_CODE = np.full(256, 4, np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _CODE[ord(chr(_b).lower())] = _i


def _op_class(op: str) -> str:
    return "match" if op == "=" else "insertion" if op[0] == "-" else "deletion" if op[-1] == "-" else "substitution"


def kmers(sequences: Iterable[str], k: int = 21, chunk: int = 100_000) -> Array:
    """Canonical k-mer codes (2 bits per base) of every window without a non-ACGT base, in input order."""
    power = (4 ** np.arange(k - 1, -1, -1)).astype(np.uint64)
    out = []
    for seq in sequences:
        codes = _CODE[np.frombuffer(seq.encode(), np.uint8)]
        for lo in range(0, max(1, len(codes) - k + 1), chunk):  # long contigs in overlapping chunks
            part = codes[lo : lo + chunk + k - 1]
            if len(part) < k:
                continue
            w = sliding_window_view(part, k)
            w = w[(w < 4).all(axis=1)].astype(np.uint64)
            out.append(np.minimum(w @ power, (3 - w[:, ::-1]) @ power))
    return np.concatenate(out) if out else np.zeros(0, np.uint64)


@dataclass
class Profile:
    rows: Counter[tuple[str, Any, str]] = field(default_factory=Counter)  # (axis, value, class)
    edits: Counter[int] = field(default_factory=Counter)  # reads by edit count
    q_cycle: Counter[tuple[int, int]] = field(default_factory=Counter)  # (cycle, Q)
    lag: Array = field(default_factory=lambda: np.zeros(6))  # n, Σx, Σy, Σx², Σy², Σxy
    multiplicity: Array = field(default_factory=lambda: np.zeros(0, np.int64))  # per distinct k-mer
    kmers_absent: int = 0
    reads: int = 0

    def summary(self) -> dict[str, Any]:
        n = sum(v for (a, _, _), v in self.rows.items() if a == "q")
        rates = {c: sum(v for (a, _, k), v in self.rows.items() if a == "q" and k == c) / max(n, 1) for c in CLASSES}
        bases = sum(self.q_cycle.values())
        return {
            "reads": self.reads,
            "rows": n,
            "error_rate": 1 - rates.pop("match"),
            **{f"{c}_rate": r for c, r in rates.items()},
            "edits_per_read_dispersion": _dispersion(self.edits),
            "mean_q": sum(q * v for (_, q), v in self.q_cycle.items()) / max(bases, 1),
            "q_lag1": _corr(self.lag),
            "kmer_absent_fraction": self.kmers_absent / max(int(self.multiplicity.sum()), 1),
        }


def _dispersion(c: Counter[int]) -> float:
    n = np.array(list(c.values()), float)
    x = np.array(list(c), float)
    mean = x @ n / n.sum()
    return float(((x - mean) ** 2 @ n / n.sum()) / mean)


def _corr(s: Array) -> float:
    n, sx, sy, sxx, syy, sxy = s
    return float((n * sxy - sx * sy) / np.sqrt((n * sxx - sx**2) * (n * syy - sy**2)))


def profile(
    records: Iterable[tuple[str, Read, int]],
    reads: Iterable[tuple[str, Sequence[int]]],
    reference_kmers: Array,
    k: int = 21,
    batch: int = 10_000,
) -> Profile:
    """Profile aligned (template, read, mate) triples, and every read's (sequence, Q in read orientation).
    `reference_kmers` is `np.unique(kmers(contigs, k))`."""
    out = Profile()

    def counted() -> Iterator[tuple[str, Read, int]]:
        for record in records:
            out.edits[_edits(record[0], record[1])] += 1
            yield record

    table = observations(counted(), (1, 1), 0, "baseline")
    iq, ic, ip, io = (table.fields.index(f) for f in ("q", "context", "pos_start", "op"))
    for key, n in table.counts.items():
        c = _op_class(str(key[io]))
        for axis, value in (("q", key[iq]), ("cycle", (int(key[ip]) - 1) // 10 * 10 + 1), ("context", key[ic])):  # type: ignore[call-overload]
            out.rows[axis, value, c] += n
    codes = []
    it = iter(reads)
    while chunk := list(islice(it, batch)):
        out.reads += len(chunk)
        for _, q in chunk:
            out.q_cycle.update(enumerate(q, 1))
            x, y = np.asarray(q[:-1], float), np.asarray(q[1:], float)
            out.lag += (len(x), x.sum(), y.sum(), x @ x, y @ y, x @ y)
        codes.append(kmers((s for s, _ in chunk), k))
    distinct, out.multiplicity = np.unique(
        np.concatenate(codes) if codes else np.zeros(0, np.uint64), return_counts=True
    )
    absent = ~np.isin(distinct, reference_kmers, assume_unique=True)
    out.kmers_absent = int(out.multiplicity[absent].sum())
    return out


def _rate_distance(real: Profile, sim: Profile, axis: str, errors: Sequence[str], pooled: bool = False) -> float:
    def groups(p: Profile) -> dict[Any, tuple[int, int]]:
        g: dict[Any, list[int]] = {}
        for (a, v, c), n in p.rows.items():
            if a == axis:
                t = g.setdefault(None if pooled else v, [0, 0])
                t[0] += n
                t[1] += n * (c in errors)
        return {v: (t[0], t[1]) for v, t in g.items()}

    r, s = groups(real), groups(sim)
    shared = [v for v in r if v in s]
    if not shared:
        return float("nan")
    w = np.array([r[v][0] for v in shared], float)
    d = [abs(np.log((s[v][1] + 0.5) / (s[v][0] + 1)) - np.log((r[v][1] + 0.5) / (r[v][0] + 1))) for v in shared]
    return float(w @ d / w.sum())


def _histogram_tv(*samples: tuple[Array, Sequence[int] | None]) -> float:
    """TV between two histograms of non-negative integers (values, weights), values capped at `_MAX_MULTIPLICITY`."""
    a, b = (np.bincount(np.minimum(v, _MAX_MULTIPLICITY), w, _MAX_MULTIPLICITY + 1) for v, w in samples)
    return float(0.5 * np.abs(a / a.sum() - b / b.sum()).sum())


def distances(real: Profile, sim: Profile) -> dict[str, float]:
    errors = CLASSES[1:]
    out = {
        "rate": _rate_distance(real, sim, "q", errors, pooled=True),
        **{f"{c}_rate": _rate_distance(real, sim, "q", (c,), pooled=True) for c in errors},
        "rate_by_q": _rate_distance(real, sim, "q", errors),
        "rate_by_cycle": _rate_distance(real, sim, "cycle", errors),
        "substitution_by_context": _rate_distance(real, sim, "context", ("substitution",)),
    }
    cycles = sorted({c for c, _ in real.q_cycle} & {c for c, _ in sim.q_cycle})
    tv, weight = [], []
    for c in cycles:
        h = [{q: v for (cc, q), v in p.q_cycle.items() if cc == c} for p in (real, sim)]
        tot = [sum(x.values()) for x in h]
        tv.append(0.5 * sum(abs(h[0].get(q, 0) / tot[0] - h[1].get(q, 0) / tot[1]) for q in h[0].keys() | h[1].keys()))
        weight.append(tot[0])
    out["q_cycle_tv"] = float(np.average(tv, weights=weight)) if tv else float("nan")
    out["q_lag1"] = abs(_corr(real.lag) - _corr(sim.lag))
    out["kmer_spectrum_tv"] = _histogram_tv(*((p.multiplicity, None) for p in (real, sim)))
    out["edits_per_read_tv"] = _histogram_tv(*((np.array(list(p.edits)), list(p.edits.values())) for p in (real, sim)))
    a, b = (p.summary()["kmer_absent_fraction"] for p in (real, sim))
    out["kmer_absent"] = float(abs(np.log(b / a)))
    return out


def verdict(ours: float, theirs: float, margin: float = 0.1, floor: float = 0.01) -> str:
    """`ours` vs `theirs` (distances from the real reads): matches within `margin` relative or `floor` absolute."""
    if abs(ours - theirs) <= max(floor, margin * max(ours, theirs)):
        return "matches"
    return "beats" if ours < theirs else "trails"


def bam_reads(path: Path, reference: Path) -> Iterator[tuple[str, list[int]]]:
    """(sequence, Q in read orientation) of every primary record, mapped or not."""
    with pysam.AlignmentFile(str(path), reference_filename=str(reference)) as aln:
        for r in aln.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.query_sequence is None or r.query_qualities is None:
                continue
            q = list(r.query_qualities)
            yield (r.query_sequence, q[::-1]) if r.is_reverse else (r.query_sequence, q)


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m sequencing_error_model.baseline",
        description="Compare simulated read sets with real reads on read-level metrics.",
    )
    p.add_argument("reference", type=Path, help="reference FASTA every read set is aligned to")
    p.add_argument("real", type=Path, help="real reads aligned (BAM/SAM/CRAM)")
    p.add_argument("simulated", nargs="+", metavar="NAME=BAM", help="simulated read sets, aligned the same way")
    p.add_argument("--native", default="native", help="the NAME of the native generator's read set")
    p.add_argument("--output", type=Path, required=True, metavar="JSON")
    p.add_argument("-k", type=int, default=21)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--max-reads", type=int)
    p.add_argument("--unclip", nargs=4, type=int, metavar=("MATCH", "MISMATCH", "OPEN", "EXTEND"))
    args = p.parse_args(argv)
    sims = dict(s.split("=", 1) for s in args.simulated)
    unclip = tuple(args.unclip) if args.unclip else None
    with pysam.FastaFile(str(args.reference)) as fasta:
        ref_kmers = np.unique(kmers((fasta.fetch(c) for c in fasta.references), args.k))

    def run(path: Path) -> Profile:
        records = islice(bam.records(path, args.reference, args.min_mapq, unclip=unclip), args.max_reads)
        return profile(records, islice(bam_reads(path, args.reference), args.max_reads), ref_kmers, args.k)

    real = run(args.real)
    profiles = {name: run(Path(path)) for name, path in sims.items()}
    dist = {name: distances(real, prof) for name, prof in profiles.items()}
    report: dict[str, Any] = {
        "reference": str(args.reference),
        "unclip": args.unclip,
        "summary": {"real": real.summary(), **{n: pr.summary() for n, pr in profiles.items()}},
        "distances": dist,
    }
    if args.native in dist:
        ours = dist[args.native]
        report["verdict"] = {
            name: {m: verdict(ours[m], d[m]) for m in ours} for name, d in dist.items() if name != args.native
        }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"distances": dist, "verdict": report.get("verdict")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
