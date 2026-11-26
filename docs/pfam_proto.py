"""Prototype behind docs/pfam_mode_plan.md: substitution errors from in-frame stops inside Pfam domains.

A read is placed on a Pfam profile by a HMMER hit. Each codon of the read is *exposed* when the rest of the
read, without that codon, still places it: a run of codons on one side scores >= T1 bits and one on the other
side >= T2. Whether a codon is exposed therefore cannot depend on its own errors. Among exposed codons an
in-frame stop is (almost) never biology, so stops / opportunities estimates the substitution rate.

    python docs/pfam_proto.py orfs R1.fq.gz R2.fq.gz OUT [--n N] [--dose D]        # OUT.faa, OUT.coords.tsv
    ... kfpfam filter + recruit/pairs.py on OUT.faa -> pairs.tsv ...
    python docs/pfam_proto.py exposure R1 R2 PREP pairs.tsv coords.tsv OUT.npy [--n N] [--dose D]
    python docs/pfam_proto.py fit PREP OUT.npy [BASELINE.npy]
    python docs/pfam_proto.py check                                                # self-check, no data
    python docs/pfam_proto.py identify                                             # what one or two reads identify

PREP is kfpfam's data/prep (accs.txt, offs.u32, profile.f32). Read names R<i> are kfpfam's seeded shuffle of
R1 + R2 (measure/real_orfs.py). `--dose D` adds uniform substitutions at rate D, the same ones in both steps.
"""

import argparse
import gzip
import itertools
import random
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

BASES, AAS = "ACGT", "ACDEFGHIKLMNPQRSTVWY"
COMP = str.maketrans("ACGTN", "TGCAN")
_TABLE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {"".join(c): a for c, a in zip(itertools.product("TCAG", repeat=3), _TABLE, strict=True)}
CLASSES = [a + b for a in BASES for b in BASES if a != b]
HMMER_BG = np.array(
    [0.0787945, 0.0151600, 0.0535222, 0.0668298, 0.0397062, 0.0695071, 0.0229198, 0.0590092, 0.0594422, 0.0963728]
    + [0.0237718, 0.0414386, 0.0482904, 0.0395639, 0.0540978, 0.0683364, 0.0540687, 0.0673417, 0.0114135, 0.0304133]
)


def reads(r1: str, r2: str, n: int | None, dose: float) -> list[tuple[str, str, int]]:
    """(sequence, qualities, mate) in kfpfam's order, substitutions added at `dose`."""
    out: list[tuple[str, str, int]] = []
    for mate, path in enumerate((r1, r2), 1):
        with gzip.open(path, "rt") as f:
            lines = [x.rstrip("\n") for x in f]
        out += [(s, q, mate) for s, q in zip(lines[1::4], lines[3::4], strict=True)]
    random.Random(5).shuffle(out)
    out = out[:n]
    if dose:
        rng = np.random.default_rng(1)
        for i, (s, q, mate) in enumerate(out):
            hit = np.flatnonzero(rng.random(len(s)) < dose)
            if len(hit):
                b = list(s)
                for j, k in zip(hit, rng.integers(0, 3, len(hit)), strict=True):
                    if b[j] in BASES:
                        b[j] = [x for x in BASES if x != b[j]][k]
                out[i] = ("".join(b), q, mate)
    return out


def orfs(rs: list[tuple[str, str, int]], out: str, min_len: int = 15) -> None:
    with open(f"{out}.faa", "w") as faa, open(f"{out}.coords.tsv", "w") as co:
        co.write("orf\tstrand\tstart\tread_len\n")
        for r, (s, _, _) in enumerate(rs):
            j = 0
            for sign, st in (("+", s), ("-", s[::-1].translate(COMP))):
                for frame in range(3):
                    at = 0
                    for p in "".join(CODE.get(st[i : i + 3], "X") for i in range(frame, len(st) - 2, 3)).split("*"):
                        if len(p) >= min_len:
                            faa.write(f">R{r}_{j}\n{p}\n")
                            co.write(f"R{r}_{j}\t{sign}\t{frame + 3 * at}\t{len(s)}\n")
                            j += 1
                        at += len(p) + 1


def outside(score: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Best local run of codon scores strictly left, and strictly right, of each codon."""

    def left(s: np.ndarray) -> np.ndarray:
        out, run, best = np.zeros(len(s)), 0.0, 0.0
        for i, x in enumerate(s):
            out[i] = best
            run = max(0.0, run + x)
            best = max(best, run)
        return out

    return left(score), left(score[::-1])[::-1]


def exposure(
    rs: list[tuple[str, str, int]], prep: str, pairs: str, coords: str, t1: float, t2: float, hit: float
) -> np.ndarray:
    """One row per exposed codon: read, mate, strand (1 = reverse), model, column, codon (0-63), stop, pos, Q x 3."""
    offs = np.fromfile(f"{prep}/offs.u32", np.uint32).astype(np.int64)
    prof = np.memmap(f"{prep}/profile.f32", np.float32, "r").reshape(-1, 21)
    model = {x: i for i, x in enumerate(Path(f"{prep}/accs.txt").read_text().split())}
    hits = []
    with open(pairs) as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            # ponytail: ungapped hits only (four in five at 125 bp); walk the CIGAR for longer reads
            if (
                float(c[4]) >= hit
                and c[14][:-1].isdigit()
                and c[14][-1] == "M"
                and int(c[0][1:].split("_")[0]) < len(rs)
            ):
                hits.append((float(c[4]), c[0], model[c[1]], int(c[6]), int(c[8])))
    need, place = {h[1] for h in hits}, {}
    with open(coords) as f:
        for line in f:
            if line[: line.find("\t")] in need:
                c = line.split("\t")
                place[c[0]] = (c[1], int(c[2]))
    best: dict[int, tuple[float, int, int, str, int]] = {}
    for s, orf, m, hmm_from, ali_from in hits:
        r, (strand, start) = int(orf[1:].split("_")[0]), place[orf]
        if r not in best or best[r][0] < s:
            best[r] = (s, m, hmm_from - 1 - (start // 3 + ali_from - 1), strand, start % 3)
    rows = []
    for r, (_, m, off, strand, frame) in best.items():
        seq, qual, mate = rs[r]
        n, rev = (len(seq) - frame) // 3, strand == "-"
        st = seq[::-1].translate(COMP) if rev else seq
        cods = [st[frame + 3 * i : frame + 3 * i + 3] for i in range(n)]
        aas, length = [CODE.get(c, "X") for c in cods], offs[m + 1] - offs[m]
        col = np.arange(n) + off
        score = np.array(
            [
                -100.0
                if a == "*" or not 0 <= k < length
                else (0.0 if a == "X" else float(prof[offs[m] + k, AAS.index(a)]))
                for a, k in zip(aas, col, strict=True)
            ]
        )
        lo, hi = outside(score)
        ok = (np.maximum(lo, hi) >= t1) & (np.minimum(lo, hi) >= t2) & (col >= 0) & (col < length)
        for i in np.flatnonzero(ok):
            if aas[i] != "X":
                at = [frame + 3 * i + p for p in range(3)]
                at = [len(seq) - 1 - x for x in at] if rev else at
                code = sum(BASES.index(b) * 4 ** (2 - p) for p, b in enumerate(cods[i]))
                rows.append((r, mate, rev, m, col[i], code, aas[i] == "*", at[1], *(ord(qual[x]) - 33 for x in at)))
    return np.array(rows, dtype=np.int32)


def codon(code: int) -> str:
    return BASES[code // 16] + BASES[code // 4 % 4] + BASES[code % 4]


def changes(c: str) -> list[tuple[str, str]]:
    return [(c[:p] + x + c[p + 1 :], c[p] + x) for p in range(3) for x in BASES if x != c[p]]


def on_read(cls: str, rev: bool) -> int:
    return CLASSES.index(cls.translate(COMP) if rev else cls)


def fit(t: np.ndarray, prep: str, lam: float = 1.0) -> tuple[np.ndarray, np.ndarray, int]:
    """Substitution rates by class, in read orientation, with opportunities and the number of stops.

    EM over the codon a stop came from. Its prior is the sample's own codons at that column (the stack), plus
    `lam` codons' worth of the HMM's emissions x the sample's codon usage.
    """
    offs = np.fromfile(f"{prep}/offs.u32", np.uint32).astype(np.int64)
    prof = np.memmap(f"{prep}/profile.f32", np.float32, "r").reshape(-1, 21)
    stop = t[:, 6] == 1
    chance, usage = np.zeros(12), Counter()
    for (rev, code), k in Counter(map(tuple, t[~stop][:, [2, 5]].tolist())).items():
        usage[codon(code)] += k
        for to, cls in changes(codon(code)):
            chance[on_read(cls, bool(rev))] += k * (CODE[to] == "*")
    by_aa = Counter()
    for c, k in usage.items():
        by_aa[CODE[c]] += k
    stack: dict[int, Counter[str]] = defaultdict(Counter)
    cols = {int(offs[m] + k) for m, k in t[stop][:, 3:5].tolist()}
    for m, k, code in t[~stop][:, 3:6].tolist():
        if int(offs[m] + k) in cols:
            stack[int(offs[m] + k)][codon(code)] += 1
    events = []
    for rev, m, k, code in t[stop][:, 2:6].tolist():
        row = int(offs[m] + k)
        emit = 2.0 ** prof[row, :20] * HMMER_BG
        cand = []
        for src, _ in changes(codon(code)):
            if CODE[src] != "*":
                back = next(cls for to, cls in changes(src) if to == codon(code))
                w = stack[row][src] + lam * emit[AAS.index(CODE[src])] * usage[src] / max(by_aa[CODE[src]], 1)
                cand.append((on_read(back, bool(rev)), w))
        events.append(cand)
    e = np.where(chance > 0, 1e-4, 0.0)
    for _ in range(200):
        n = np.zeros(12)
        for cand in events:
            w = np.array([e[c] * p for c, p in cand])
            if w.sum() > 0:
                for (c, _), x in zip(cand, w / w.sum(), strict=True):
                    n[c] += x
        e = np.where(chance > 0, n / np.maximum(chance, 1), 0.0)
    return e, chance, int(stop.sum())


def report(t: np.ndarray, prep: str, base: np.ndarray | None) -> None:
    e, chance, n = fit(t, prep)
    total = 3 * n / chance.sum()
    print(f"{len(t)} exposed codons, {n} stops, substitutions per base {total:.2e} +- {total / max(n, 1) ** 0.5:.1e}")
    if base is not None:
        e0, c0, n0 = fit(base, prep)
        print(f"baseline {3 * n0 / c0.sum():.2e}, increment {total - 3 * n0 / c0.sum():.3e}")
        e = e - e0
    for c, x, k in zip(CLASSES, e, chance, strict=True):
        print(f"{c[0]}>{c[1]}", "unidentified" if k == 0 else f"{x:.2e}", int(k))
    q = t[:, 8:11].min(1)
    for lo, hi in ((0, 20), (20, 30), (30, 36), (36, 99)):
        k = (q >= lo) & (q < hi)
        print(f"lowest Q in codon {lo}-{hi - 1}: {int(t[k, 6].sum())} stops in {int(k.sum())} codons")


def check() -> None:
    sources = {cls for c, a in CODE.items() if a != "*" for to, cls in changes(c) if CODE[to] == "*"}
    both = {CLASSES[on_read(x, rev)] for x in sources for rev in (False, True)}
    assert set(CLASSES) - both == {"AG", "TC"}, both  # no transition into A>G creates a stop on either strand
    lo, hi = outside(np.array([5.0, 5.0, -100.0, 7.0]))
    assert lo.tolist() == [0, 5, 10, 10] and hi.tolist() == [7, 7, 7, 0]
    assert codon(0) == "AAA" and codon(63) == "TTT" and on_read("CT", True) == CLASSES.index("GA")


def identify(two: bool, diverse: bool, start: int, columns: int = 40, steps: int = 4000) -> tuple[float, float]:
    """(true, estimated) error rate by EM on exact counts: the matrix E shared, each column's bases free.

    One view is a read; two views are both mates over one base. `diverse=False` gives every column the same bases.
    No column has a forbidden base, so nothing but the two views can separate error from biology.
    """
    e = np.array([[0, 4, 10, 3], [3, 0, 2, 12], [11, 2, 0, 4], [3, 9, 4, 0]]) / 1000
    e += np.diag(1 - e.sum(1))
    pi = np.random.default_rng(0).dirichlet([0.5] * 4, columns)
    pi = pi if diverse else np.tile(pi[:1], (columns, 1))
    counts = np.einsum("kc,ca,cb->kab", pi, e, e) if two else np.einsum("kc,ca->ka", pi, e)
    rng = np.random.default_rng(start)
    f = rng.uniform(0.0005, 0.06, (4, 4))
    np.fill_diagonal(f, 0)
    f += np.diag(1 - f.sum(1))
    q = rng.dirichlet([1] * 4, columns)
    for _ in range(steps):
        joint = np.einsum("kc,ca,cb->kabc", q, f, f) if two else np.einsum("kc,ca->kac", q, f)
        w = joint / joint.sum(-1, keepdims=True) * counts[..., None]
        q = w.sum((1, 2)) if two else w.sum(1)
        f = w.sum((0, 2)).T + w.sum((0, 1)).T if two else w.sum(0).T
        q, f = q / q.sum(1, keepdims=True), f / f.sum(1, keepdims=True)
    return 1 - float(np.diag(e).mean()), 1 - float(np.diag(f).mean())


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("orfs")
    b = sub.add_parser("exposure")
    for x in (a, b):
        x.add_argument("r1")
        x.add_argument("r2")
        x.add_argument("--n", type=int)
        x.add_argument("--dose", type=float, default=0.0)
    a.add_argument("out")
    for name in ("prep", "pairs", "coords", "out"):
        b.add_argument(name)
    b.add_argument("--t1", type=float, default=25.0)
    b.add_argument("--t2", type=float, default=6.0)
    b.add_argument("--hit", type=float, default=0.0, help="lowest HMMER score a hit may have; far below t1")
    c = sub.add_parser("fit")
    c.add_argument("prep")
    c.add_argument("table")
    c.add_argument("baseline", nargs="?")
    sub.add_parser("check")
    sub.add_parser("identify")
    v = p.parse_args()
    if v.cmd == "orfs":
        orfs(reads(v.r1, v.r2, v.n, v.dose), v.out)
    elif v.cmd == "exposure":
        np.save(v.out, exposure(reads(v.r1, v.r2, v.n, v.dose), v.prep, v.pairs, v.coords, v.t1, v.t2, v.hit))
    elif v.cmd == "fit":
        report(np.load(v.table), v.prep, np.load(v.baseline) if v.baseline else None)
    elif v.cmd == "identify":
        for name, two, diverse in (("one read", False, True), ("two mates", True, True), ("two, same", True, False)):
            got = [identify(two, diverse, i) for i in range(4)]
            print(f"{name:10s} truth {got[0][0]:.4f} from 4 starts", *(f"{x:.4f}" for _, x in got))
    else:
        check()


if __name__ == "__main__":
    main()
