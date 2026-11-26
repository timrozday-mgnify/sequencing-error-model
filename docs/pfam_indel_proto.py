"""Prototype: indels by locus (docs/pfam_mode_plan.md section 13.1). Not part of the `pfam` fit.

Reads that share a Pfam column, a strand and the 12 bases before an exposed codon are one locus. Each read is compared
with the consensus of the locus' other reads over the next `W` bases: at its first mismatch, the rest either matches as
a substitution, or matches after a shift of 1-3 bases (a deletion or an insertion in the read). Sequencing indels are
mostly 1 or 2 long; strain indels in coding sequence are mostly 3, so length mod 3 separates them, as the codon effect
does for substitutions. No Q, no reference.

    python docs/pfam_indel_proto.py inject R1 R2 --rate 5e-4 --output PREFIX     # PREFIX_{1,2}.fastq.gz and PREFIX.truth.tsv
    python docs/pfam_indel_proto.py detect R1 R2 --pairs P.tsv --prep PREP [--every-n N] [--truth PREFIX.truth.tsv]
"""

import argparse
import gzip
from collections import defaultdict
from pathlib import Path

import numpy as np

from sequencing_error_model.sources import pfam
from sequencing_error_model.sources.fastq_quality import read_fastq

FLANK, W, T, SHIFTS = 12, 24, 10, (1, 2, 3)  # left key, window compared, tail that must match exactly, shifts tried
LAST = W - T - max(SHIFTS)  # the last window offset at which an indel can still be told from a substitution
LENGTHS, SHARE = (1, 2, 3), (0.6, 0.2, 0.2)
CODE = np.full(256, 4, np.uint8)
for i, b in enumerate(b"ACGT"):
    CODE[b] = i


def inject(args: argparse.Namespace) -> None:
    """Indels at `rate` per base, at known positions in the output reads; inserted bases copy their neighbour's Q."""
    rng = np.random.default_rng(1)
    truth = []
    number = -1
    with gzip.open(f"{args.output}_1.fastq.gz", "wt", 3) as o1, gzip.open(f"{args.output}_2.fastq.gz", "wt", 3) as o2:
        outs = (o1, o2)
        for number, pair in enumerate(zip(read_fastq(args.r1), read_fastq(args.r2), strict=True)):
            for mate, (seq, qual) in enumerate(pair):
                q, s, qs, at = list(qual), [], [], 0
                for j in np.flatnonzero(rng.random(len(seq)) < args.rate):
                    if j < at:
                        continue  # inside the last deletion
                    s.append(seq[at:j])
                    qs.extend(q[at:j])
                    n = int(rng.choice(LENGTHS, p=SHARE))
                    pos = sum(map(len, s))
                    if rng.random() < 0.5:
                        s.append("".join(rng.choice(list("ACGT"), n)))
                        qs.extend([q[j]] * n)
                        truth.append((number, mate, pos, "ins", n))
                        at = j
                    else:
                        truth.append((number, mate, pos, "del", n))
                        at = j + n
                s.append(seq[at:])
                qs.extend(q[at:])
                text = "".join(s)
                outs[mate].write(f"@{number}\n{text}\n+\n{''.join(chr(x + 33) for x in qs[: len(text)])}\n")
    with open(f"{args.output}.truth.tsv", "w") as f:
        f.write("pair\tmate\tpos\tkind\tlength\n")
        f.writelines("\t".join(map(str, t)) + "\n" for t in truth)
    print(f"{len(truth):,} indels in {number + 1:,} pairs")


def windows(args: argparse.Namespace) -> tuple[dict[tuple, list[tuple[int, int]]], dict[int, np.ndarray]]:
    """Per locus key, (read, window start); and each read's bases as codes."""
    profiles = pfam.Profiles.load(args.prep)
    hits = pfam.read_hits(args.pairs, profiles)
    groups: dict[tuple, list[tuple[int, int]]] = defaultdict(list)
    texts: dict[int, str] = {}
    for number, (s1, _, s2, _) in enumerate(pfam.read_pairs(args.r1, args.r2, None, args.every_n)):
        for mate, seq in enumerate((s1, s2)):
            rid = 2 * number + mate
            if rid not in hits:
                continue
            s = texts[rid] = seq.upper()
            for c in pfam.exposed(s, hits[rid], profiles, args.t1, args.t2):
                lo = min(c.positions)
                if lo >= FLANK and lo + W <= len(s):
                    groups[(*c.column, c.reverse, s[lo - FLANK : lo])].append((rid, lo))
    if args.recruit == "flank":
        # Exposure needs the hit to run through the window, which a frameshift in it ends: a selection on the very
        # error counted. So a locus, once placed by >= 2 exposed reads, takes every read with its 12-base key.
        keys: dict[str, tuple | None] = {}
        for g, members in groups.items():
            if len(members) >= 2:
                keys[g[-1]] = None if g[-1] in keys else g  # a key two loci share places neither
        if args.scan:
            # standard additions with no second search: the loci are placed by these reads, and a dosed copy of
            # them (same numbering) is what is scanned and judged
            groups = defaultdict(list)
            dosed = pfam.read_pairs(args.scan[0], args.scan[1], None, args.every_n)
            texts = {
                2 * n + m: s.upper() for n, p in enumerate(dosed) for m, s in enumerate(p[::2]) if 2 * n + m in texts
            }
        for rid, s in texts.items():
            for i in range(FLANK, len(s) - W + 1):
                if g := keys.get(s[i - FLANK : i]):
                    groups[g].append((rid, i))
        groups = {g: list(dict.fromkeys(m)) for g, m in groups.items()}
    seqs = {rid: CODE[np.frombuffer(texts[rid].encode(), np.uint8)] for m in groups.values() for rid, _ in m}
    return groups, seqs


def call(read: np.ndarray, cons: np.ndarray, ok: np.ndarray) -> tuple[int, str, int] | None:
    """(offset, kind, length) of the read's first difference from the consensus, or None where it agrees to `LAST`."""
    differ = ok & (read != cons)
    if not differ[: LAST + 1].any():
        return None
    m = int(np.argmax(differ))

    def agree(
        r: int, c: int
    ) -> bool:  # read from r and consensus from c, for T bases        # every compared consensus base must be confident: a gap matching anything lets a shift fit by chance
        return bool(ok[c : c + T].all() and np.array_equal(read[r : r + T], cons[c : c + T]))

    if agree(m + 1, m + 1):
        return m, "sub", 0
    for d in SHIFTS:
        if agree(m, m + d):  # the read skips d consensus bases
            return m, "del", d
        if agree(m + d, m):  # the read has d extra bases
            return m, "ins", d
    return m, "other", 0


def detect(args: argparse.Namespace) -> None:
    groups, seqs = windows(args)
    covered = {rid: np.zeros(len(s), bool) for rid, s in seqs.items()}
    events: set[tuple[int, int, str, int]] = set()
    loci = reads = 0
    for members in groups.values():
        if len(members) < 3:
            continue
        loci += 1
        win = np.stack([seqs[r][lo : lo + W] for r, lo in members])
        counts = np.stack([(win == b).sum(0) for b in range(5)])  # ACGT and N, x W
        for k, (rid, lo) in enumerate(members):
            own = counts.copy()
            own[win[k], np.arange(W)] -= 1  # leave the read out of its own consensus
            cons, top = own[:4].argmax(0), own[:4].max(0)
            ok = (top >= 2) & (2 * top > len(members) - 1)
            found = call(win[k], cons, ok)
            end = LAST if found is None else min(found[0], LAST)
            covered[rid][lo : lo + end + 1] |= ok[: end + 1]
            reads += 1
            if found and found[1] in ("del", "ins"):
                events.add((rid, lo + found[0], found[1], found[2]))
    n = sum(int(c.sum()) for c in covered.values())
    by = defaultdict(int)
    for _, _, kind, d in events:
        by[(kind, d)] += 1
    print(f"{loci:,} loci (>= 3 reads), {reads:,} windows; {n:,} read positions where an indel would be seen")
    print("detected per covered base: " + ", ".join(f"{k}{d} {c} ({c / n:.1e})" for (k, d), c in sorted(by.items())))
    shift = sum(c for (k, d), c in by.items() if d % 3)
    print(f"frameshift (1, 2) {shift} = {shift / n:.2e} per base; in frame (3) {by[('del', 3)] + by[('ins', 3)]}")
    if args.truth:
        score(args.truth, events, covered)


def score(path: Path, events: set[tuple[int, int, str, int]], covered: dict[int, np.ndarray]) -> None:
    true = defaultdict(list)
    for line in path.read_text().splitlines()[1:]:
        pair, mate, pos, kind, length = line.split()
        true[2 * int(pair) + int(mate)].append((int(pos), kind, int(length)))
    seen = [(rid, p, k, d) for rid, c in covered.items() for p, k, d in true.get(rid, ()) if p < len(c) and c[p]]
    n = sum(int(c.sum()) for c in covered.values())
    hit = lambda rid, p, k, d: any(  # noqa: E731
        (rid, q, k, d) in events for q in range(p, p + 5)
    )  # a repeat moves the first mismatch to the run's end
    found = [e for e in seen if hit(*e)]
    real = sum(
        any(abs(q - p) <= 4 and k == kk and d == dd for q, kk, dd in true.get(rid, ())) for rid, p, k, d in events
    )
    print(
        f"truth: {len(seen):,} injected indels at covered positions ({len(seen) / n:.2e} per base); {len(found):,} found ({len(found) / max(len(seen), 1):.0%})"
    )
    print(f"       {len(events):,} called, {real:,} of them injected ({real / max(len(events), 1):.0%})")
    for d in SHIFTS:
        s = [e for e in seen if e[3] == d]
        f = [e for e in s if hit(*e)]
        print(f"       length {d}: {len(f)} of {len(s)} found")


def main() -> None:
    p = argparse.ArgumentParser(description="Indels by locus prototype")
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("inject")
    a.add_argument("r1", type=Path)
    a.add_argument("r2", type=Path)
    a.add_argument("--rate", type=float, required=True)
    a.add_argument("--output", type=Path, required=True)
    b = sub.add_parser("detect")
    b.add_argument("r1", type=Path)
    b.add_argument("r2", type=Path)
    b.add_argument("--pairs", type=Path, required=True)
    b.add_argument("--prep", type=Path, required=True)
    b.add_argument("--every-n", type=int, default=1)
    b.add_argument("--t1", type=float, default=35.0)
    b.add_argument("--t2", type=float, default=6.0)
    b.add_argument("--truth", type=Path)
    b.add_argument("--recruit", choices=("flank", "exposure"), default="flank", help="how a read joins a locus")
    b.add_argument("--scan", type=Path, nargs=2, metavar=("R1", "R2"), help="a dosed copy to judge, placed by R1 R2")
    args = p.parse_args()
    inject(args) if args.command == "inject" else detect(args)


if __name__ == "__main__":
    main()
