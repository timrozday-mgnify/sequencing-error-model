"""Prototype: within-column contrast (docs/pfam_mode_plan.md section 13.2, Tier 3). Not part of the `pfam` fit.

At one Pfam column, biology does not depend on the cycle, Q or mate a base was read at; sequencing error does. So
for every base of every exposed codon, `mismatch with its column's consensus = b(column) + e(covariates)` with b
free per column (strain variation, paralogs, a wrong consensus) and e shared. e is identified up to a constant:
the contrast gives the shape, and an anchor (stops, overlap, locus) must give the level. Q is a covariate here,
never a label: the labels are mismatches with the column's other reads.

    python docs/pfam_contrast_proto.py table R1 R2 --pairs P.tsv --prep PREP [--every-n N --dose D] [--sam S --reference F] --output T.npz
    python docs/pfam_contrast_proto.py fit T.npz [--minus T0.npz]
"""

import argparse
from array import array
from pathlib import Path

import numpy as np
import pysam

from sequencing_error_model import pfam_truth
from sequencing_error_model.sources import pfam

BASE = {b: i for i, b in enumerate("ACGT")}
COMP = np.array([3, 2, 1, 0])  # A<->T, C<->G by index
Q_EDGES = (15, 20, 25, 30, 35, 38)
CYCLE_EDGES = (25, 50, 75, 100, 125)
MIN_DEPTH = {"group": 20, "locus": 3}  # bases per group: fewer cannot separate b from e; a locus needs a majority
FLANK = 12  # bases each side that make a locus (pfam_locus)
FIELDS = ("group", "locus", "base", "q", "cycle", "mate", "rid", "true", "length", "insert")


def table(args: argparse.Namespace) -> dict[str, np.ndarray]:
    """One row per base of every exposed codon, in read order and orientation."""
    profiles = pfam.Profiles.load(args.prep)
    hits = pfam.read_hits(args.pairs, profiles)
    keys: dict[tuple[int, int, int, bool], int] = {}
    loci: dict[tuple[int, int, int, bool, str], int] = {}
    rows = {f: array("q") for f in FIELDS}  # typed: tens of millions of bases at P2 depth
    for number, (s1, q1, s2, q2) in enumerate(pfam.read_pairs(args.r1, args.r2, None, args.every_n, args.dose)):
        for mate, (seq, qual) in enumerate(((s1, q1), (s2, q2))):
            rid = 2 * number + mate
            if rid not in hits:
                continue
            s, q = seq.upper(), list(qual)
            for c in pfam.exposed(s, hits[rid], profiles, args.t1, args.t2):
                lo, hi = min(c.positions), max(c.positions)
                flank = (
                    f"{s[lo - FLANK : lo]}|{s[hi + 1 : hi + 1 + FLANK]}" if lo >= FLANK and hi + FLANK < len(s) else ""
                )
                for idx, at in enumerate(c.positions):
                    if s[at] not in BASE:
                        continue
                    # the strand is part of the group: a column's read-orientation bases differ by strand
                    g = keys.setdefault((*c.column, idx, c.reverse), len(keys))
                    # the same column and the same bases either side: one gene in one strain (-1: a flank runs off the read)
                    locus = loci.setdefault((*c.column, idx, c.reverse, flank), len(loci)) if flank else -1
                    for f, v in zip(FIELDS, (g, locus, BASE[s[at]], q[at], at, mate, rid, -1, len(s), -1), strict=True):
                        rows[f].append(v)
    out = {f: np.frombuffer(v, np.int64).copy() for f, v in rows.items()}
    if args.sam:
        add_truth(out, args.sam, args.reference)
    return out


def add_truth(t: dict[str, np.ndarray], sam: Path, reference: Path) -> None:
    """The genome's base under each row, in read orientation, from the read's primary alignment."""
    genome: dict[str, str] = {}
    with pysam.FastaFile(str(reference)) as fa, pysam.AlignmentFile(str(sam), check_sq=False) as f:
        for r in f.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_name is None or not r.reference_name:
                continue
            rid = 2 * int(r.query_name) + int(r.is_read2)
            lo, hi = np.searchsorted(t["rid"], (rid, rid + 1))
            n = r.infer_read_length()
            if lo == hi or n is None:
                continue
            sites = pfam_truth.sites(r, n)
            if r.is_proper_pair and r.template_length:
                t["insert"][lo:hi] = abs(r.template_length)
            contig = genome.setdefault(r.reference_name, fa.fetch(r.reference_name).upper())
            for k in range(lo, hi):
                site = sites[t["cycle"][k]]
                if site >= 0 and (b := BASE.get(contig[site], -1)) >= 0:
                    t["true"][k] = COMP[b] if r.is_reverse else b


def contrast(group: np.ndarray, z: np.ndarray, y: np.ndarray, w: np.ndarray | None = None) -> np.ndarray:
    """e per covariate bin in `y = b(group) + e(z)`, by alternating means; centred, since only its shape is identified."""
    w = np.ones(len(y)) if w is None else w
    ng, nz = np.bincount(group, w), np.bincount(z, w)
    e = np.zeros(len(nz))
    for _ in range(200):
        b = np.bincount(group, w * (y - e[z]), len(ng)) / np.maximum(ng, 1e-12)
        new = np.bincount(z, w * (y - b[group]), len(nz)) / np.maximum(nz, 1e-12)
        if np.abs(new - e).max() < 1e-9:
            e = new
            break
        e = new
    return e - np.average(e[z], weights=w)


def fit(t: dict[str, np.ndarray], minus: dict[str, np.ndarray] | None, by: str = "group", boot: int = 20) -> None:
    t = {k: v[t[by] >= 0] for k, v in t.items()}
    keep = np.bincount(t[by])[t[by]] >= MIN_DEPTH[by]
    t = {k: v[keep] for k, v in t.items()}
    _, group = np.unique(t[by], return_inverse=True)
    cons = np.zeros((group.max() + 1, 4))
    np.add.at(cons, (group, t["base"]), 1)
    y = (t["base"] != cons.argmax(1)[group]).astype(float)
    known = t["true"] >= 0
    err = (t["base"] != t["true"]) & known
    print(
        f"\n=== by {by}: {len(y):,} bases in {group.max() + 1:,} groups (depth >= {MIN_DEPTH[by]}); mismatch with consensus {y.mean():.2e}"
    )
    if known.any():
        print(f"truth on {known.mean():.0%}: error {err[known].mean():.2e} per base")
    rng = np.random.default_rng(0)
    covariates = {
        "Q": (np.digitize(t["q"], Q_EDGES), ["<15", "15-19", "20-24", "25-29", "30-34", "35-37", ">=38"]),
        "cycle": (np.digitize(t["cycle"], CYCLE_EDGES), ["1-25", "26-50", "51-75", "76-100", "101-125", ">125"]),
        "mate": (t["mate"], ["R1", "R2"]),
    }
    for name, (z, labels) in covariates.items():
        _, z = np.unique(z, return_inverse=True)
        e = contrast(group, z, y)
        sd = np.std([contrast(group, z, y, rng.poisson(1.0, group.max() + 1)[group]) for _ in range(boot)], 0)
        naive = np.bincount(z, y) / np.bincount(z)
        # the level is not the contrast's: take the truth's (the anchors' in real use), so the shape can be compared
        level = err[known].mean() if known.any() else 0.0
        print(
            f"\n{name}: shape + level{'' if known.any() else ' (no truth: level 0)'}; naive = consensus mismatch, no column effect"
        )
        print(f"{'bin':>8} {'bases':>10} {'truth':>9} {'contrast':>18} {'naive':>9}")
        for i, lab in enumerate(labels[: len(e)]):
            n = int((z == i).sum())
            truth = err[(z == i) & known].mean() if known.any() else float("nan")
            print(f"{lab:>8} {n:>10,} {truth:>9.2e} {level + e[i]:>9.2e} ± {sd[i]:.0e} {naive[i]:>9.2e}")
    if minus is not None:
        dose(t, minus, group, y, by)


def dose(t: dict[str, np.ndarray], t0: dict[str, np.ndarray], group: np.ndarray, y: np.ndarray, by: str) -> None:
    """Standard additions add error flat in Q: the dosed shape minus the undosed one should be flat."""
    # ponytail: the two runs are fitted apart, so a column's consensus can differ between them; fine at 1e-2
    t0 = {k: v[t0[by] >= 0] for k, v in t0.items()}
    t0 = {k: v[np.bincount(t0[by])[t0[by]] >= MIN_DEPTH[by]] for k, v in t0.items()}
    _, g0 = np.unique(t0[by], return_inverse=True)
    c0 = np.zeros((g0.max() + 1, 4))
    np.add.at(c0, (g0, t0["base"]), 1)
    y0 = (t0["base"] != c0.argmax(1)[g0]).astype(float)
    _, z = np.unique(np.digitize(t["q"], Q_EDGES), return_inverse=True)
    _, z0 = np.unique(np.digitize(t0["q"], Q_EDGES), return_inverse=True)
    diff = contrast(group, z, y) - contrast(g0, z0, y0)
    print("\ndosed minus undosed, by Q bin (flat if the contrast sees added error as error; centred):")
    print("  " + "  ".join(f"{x:+.1e}" for x in diff))


def main() -> None:
    p = argparse.ArgumentParser(description="Within-column contrast prototype")
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("table")
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
    a.add_argument("--output", type=Path, required=True)
    b = sub.add_parser("fit")
    b.add_argument("table", type=Path)
    b.add_argument("--minus", type=Path, help="the undosed table, for the dose check")
    args = p.parse_args()
    if args.command == "table":
        np.savez_compressed(args.output, allow_pickle=False, **table(args))
    else:
        load = lambda x: dict(np.load(x))  # noqa: E731
        for by in ("group", "locus"):
            fit(load(args.table), load(args.minus) if args.minus else None, by)


if __name__ == "__main__":
    main()
