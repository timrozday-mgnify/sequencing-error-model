"""Prototype: a stop's source codon from the reads that contradict it (docs/pfam_mode_plan.md section 12, item 3).

The fit splits each stop over its candidate sources (sense codons one substitution away) by the column's stack prior.
Reads of the same locus (same column, same 12 bases either side) that carry a sense codon there say which one it
was, with no prior. This compares the two against the genome's codon, and the op spectra they give. Equal rates
per op throughout, so only the prior differs.

    python docs/pfam_source_proto.py R1 R2 --pairs P.tsv --prep PREP [--every-n N --dose D] [--sam S --reference F]
"""

import argparse
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pysam

from sequencing_error_model import pfam_truth
from sequencing_error_model.sources import pfam

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
OPS = [f"{a}>{b}" for a in "ACGT" for b in "ACGT" if a != b]


class Genome(dict[str, str]):
    """Contigs fetched when first asked for."""

    def __init__(self, fasta: Path) -> None:
        super().__init__()
        self._fa = pysam.FastaFile(str(fasta))

    def __missing__(self, contig: str) -> str:
        self[contig] = self._fa.fetch(contig).upper()
        return self[contig]


def op(source: str, stop: str, p: int, reverse: bool) -> str:
    """The substitution in read orientation."""
    a, b = source[p], stop[p]
    return f"{COMP[a]}>{COMP[b]}" if reverse else f"{a}>{b}"


def main() -> None:
    a = argparse.ArgumentParser(description="Stop source codons from contradicting reads")
    a.add_argument("r1", type=Path)
    a.add_argument("r2", type=Path)
    a.add_argument("--pairs", type=Path, required=True)
    a.add_argument("--prep", type=Path, required=True)
    a.add_argument("--every-n", type=int, default=1)
    a.add_argument("--dose", type=float, default=0.0)
    a.add_argument("--t1", type=float, default=35.0)
    a.add_argument("--t2", type=float, default=6.0)
    a.add_argument("--sam", type=Path)
    a.add_argument("--reference", type=Path)
    args = a.parse_args()

    profiles = pfam.Profiles.load(args.prep)
    hits = pfam.read_hits(args.pairs, profiles)
    pairs = pfam.read_pairs(args.r1, args.r2, None, args.every_n, args.dose)
    reads = pfam_truth.exposed_reads(pairs, hits, profiles, args.t1, args.t2)
    ev = pfam.Evidence(0, (0, 0), args.t1)
    loci: dict[tuple, dict[int, str]] = defaultdict(dict)  # locus -> fragment -> sense codon
    stops = []
    for rid, (s, _, codons) in reads.items():
        for c in codons:
            flank = pfam.flanks(s, c)
            if c.bases in pfam.STOPS:
                stops.append((rid, c, flank))
            elif "N" not in c.bases:
                ev.stack[c.column][c.bases] += 1
                if flank:
                    loci[(c.column, flank)][rid // 2] = c.bases
    prior = pfam.Prior(ev, profiles)
    placed = pfam_truth.placements(args.sam, reads, 20) if args.sam else {}
    genome = Genome(args.reference) if args.sam else None

    rows = []  # per stop: sources, column posterior, locus source or None, true source or None
    for rid, c, flank in stops:
        cand = pfam._SOURCES[c.bases]
        w = np.array([prior.weight(c.column, src) for _, src in cand])
        post = w / w.sum()
        others = {f: x for f, x in loci.get((c.column, flank), {}).items() if f != rid // 2} if flank else {}
        seen = set(others.values())
        here = next(iter(seen)) if len(seen) == 1 else None
        true = None
        if genome is not None and (p := placed.get(rid)) and (t := pfam_truth.true_bases(c, p, genome)):
            true = t[1]
        rows.append((c, cand, post, here, len(others), true))

    n = len(rows)
    with_locus = [r for r in rows if r[3] is not None]
    print(
        f"{n:,} stops; {len(with_locus):,} ({len(with_locus) / max(n, 1):.0%}) with other fragments that agree on a sense codon"
    )
    by_depth = Counter(min(r[4], 3) for r in with_locus)
    print(
        "  of those, other fragments: "
        + ", ".join(f"{k}{'+' if k == 3 else ''}: {v}" for k, v in sorted(by_depth.items()))
    )
    off = sum(r[3] not in {src for _, src in r[1]} for r in with_locus)
    print(f"  locus codon not one substitution from the stop: {off} (more than one change, or biology)")
    if genome is None:
        agree = sum(r[3] == r[1][int(np.argmax(r[2]))][1] for r in with_locus if r[3] in {src for _, src in r[1]})
        print(f"  locus source = column prior's best: {agree} of {len(with_locus) - off}")
        return

    errors = [r for r in rows if r[5] is not None and r[5] not in pfam.STOPS and r[5] in {src for _, src in r[1]}]
    print(f"\ntruth: {len(errors):,} stops are one substitution from the genome's sense codon (errors)")

    def report(name: str, use: list) -> None:
        if not use:
            return
        best = np.mean([r[1][int(np.argmax(r[2]))][1] == r[5] for r in use])
        mass = np.mean([r[2][[src for _, src in r[1]].index(r[5])] for r in use])
        loc = [r for r in use if r[3] is not None]
        right = np.mean([r[3] == r[5] for r in loc]) if loc else float("nan")
        print(
            f"  {name}: {len(use):,} stops; column prior best right {best:.0%}, mass on truth {mass:.0%}; "
            f"locus covers {len(loc) / len(use):.0%}, right {right:.1%}"
        )

    report("all errors", errors)
    report("with a locus source", [r for r in errors if r[3] is not None])
    report("without", [r for r in errors if r[3] is None])

    # the op spectrum each way gives, against the truth's, on the same stops
    spectra: dict[str, Counter[str]] = {k: Counter() for k in ("truth", "column", "locus+column")}
    for c, cand, post, here, _, true in errors:
        spectra["truth"][op(true, c.bases, next(p for p, s in cand if s == true), c.reverse)] += 1
        for (p, src), x in zip(cand, post, strict=True):
            spectra["column"][op(src, c.bases, p, c.reverse)] += x
        if here is not None and here in {s for _, s in cand}:
            spectra["locus+column"][op(here, c.bases, next(p for p, s in cand if s == here), c.reverse)] += 1
        else:
            for (p, src), x in zip(cand, post, strict=True):
                spectra["locus+column"][op(src, c.bases, p, c.reverse)] += x
    print("\nop spectrum of those stops (A>G, T>C cannot make a stop)")
    print(f"{'op':>5} {'truth':>7} {'column':>8} {'locus+column':>13}")
    for o in OPS:
        t = spectra["truth"][o]
        if t or spectra["column"][o]:
            print(f"{o:>5} {t:>7} {spectra['column'][o]:>8.1f} {spectra['locus+column'][o]:>13.1f}")
    for k in ("column", "locus+column"):
        l1 = sum(abs(spectra[k][o] - spectra["truth"][o]) for o in OPS) / len(errors)
        print(f"  {k}: share of stops put on the wrong op {l1 / 2:.1%}")


if __name__ == "__main__":
    main()
