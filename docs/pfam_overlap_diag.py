"""P2 diagnosis: how much of the genomes' error can the mate-overlap anchor see? (docs/pfam_mode_plan.md, §7 P2)

The overlap anchor counts a base as an error only where the mates disagree. An error already in the template (library
PCR, cluster amplification, DNA damage) is read the same way by both mates, so the anchor cannot see it; the genomes
and the stops can. This splits the errors at reference sites both mates of a pair cover:

- `one`: one mate differs from the genome, the other matches it (the anchor sees it);
- `shared`: both mates carry the same non-genome base (it cannot);
- `split`: both differ from the genome, with different bases (seen as a disagreement).

It also gives the error rate of every aligned base by insert size (Q >= 36, 10 or more bases from either alignment end),
since only pairs whose insert is shorter than two reads overlap at all.

Sites where these reads carry another base at >= 20 % are masked as `pfam_truth` masks them (strain variation).

    python docs/pfam_overlap_diag.py aligned.sam[.gz] genomes.fna [--output diag.json]

`aligned.sam` is P2's `keep/<run>/aligned.sam` (minimap2 -ax sr of the undosed pairs; query names are pair numbers).
"""

import argparse
import gzip
import json
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pysam

from sequencing_error_model.sources import bam

BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")
Q_BINS = ((0, 20), (20, 30), (30, 36), (36, 99))
Q_EDGES = (20, 30, 36)  # the same groups, for np.digitize
INSERT_STEP = 50
INSERT_BINS = tuple(f"{k}+" for k in range(0, 750, INSERT_STEP))  # the last holds 700 and over
SHORT, LONG = 200, 350  # inserts compared in the substitution spectrum
CODE4 = {b: i for i, b in enumerate("ACGT")}
CODE = np.full(256, -1, np.int64)
for _i, _b in enumerate(BASES):
    CODE[ord(_b)] = _i


def genome(path: Path) -> dict[str, str]:
    with pysam.FastxFile(str(path)) as f:
        return {r.name: r.sequence.upper() for r in f if r.sequence}


def reads(sam: Path, min_mapq: int) -> Iterator[pysam.AlignedSegment]:
    """Primary, mapped alignments at `min_mapq` or more."""
    with pysam.AlignmentFile(str(sam), check_sq=False) as f:
        for r in f.fetch(until_eof=True):
            if not (r.is_secondary or r.is_supplementary or r.is_unmapped) and r.mapping_quality >= min_mapq:
                yield r


def variant_masks(sam: Path, ref: dict[str, str], min_mapq: int, max_alt_freq: float) -> dict[str, np.ndarray]:
    counts: dict[str, np.ndarray] = {}
    for r in reads(sam, min_mapq):
        at = np.array(r.get_aligned_pairs(matches_only=True), np.int64).reshape(-1, 2)
        base = CODE[np.frombuffer(r.query_sequence.encode(), np.uint8)[at[:, 0]]]
        table = counts.get(r.reference_name)
        if table is None:
            table = counts[r.reference_name] = np.zeros((len(ref[r.reference_name]) + 1, 6), np.int32)
        on = base >= 0
        np.add.at(table, (at[on, 1], base[on]), 1)
    return bam.site_masks(counts, lambda c: ref[c], max_alt_freq)[0]


def q_bin(q: int) -> str:
    return next(f"q{lo}-{hi - 1}" for lo, hi in Q_BINS if lo <= q < hi)


def _sites(r: pysam.AlignedSegment) -> dict[int, tuple[str, int]]:
    """Reference site -> (base on the forward strand, Q)."""
    seq, qual = r.query_sequence, r.query_qualities
    return {rp: (seq[qp], qual[qp]) for qp, rp in r.get_aligned_pairs(matches_only=True)}


def _op(g: str, b: str, reverse: bool) -> str:
    """The substitution in the read's own orientation."""
    return (f"{g}>{b}".translate(COMP)) if reverse else f"{g}>{b}"


def _rates(bases: np.ndarray, errors: np.ndarray) -> list[float | None]:
    return [float(e / n) if n else None for n, e in zip(bases.tolist(), errors.tolist(), strict=True)]


def by_insert(strata: np.ndarray) -> dict[str, object]:
    """Error by insert: per Q group, standardised over Q x cycle quarter x mate, and the long-insert spectrum."""
    match = np.eye(4, dtype=bool)
    bases = strata.sum((4, 5))  # insert x q x quarter x mate
    errors = (strata * ~match).sum((4, 5))
    # standardised: each insert's per-stratum rates, weighted by the strata's share of all bases
    weight = bases.sum(0) / max(bases.sum(), 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        rate = errors / bases
    seen = bases > 0
    standard = [
        float((rate[i][seen[i]] * weight[seen[i]]).sum() / weight[seen[i]].sum()) if seen[i].any() else None
        for i in range(len(INSERT_BINS))
    ]
    q_names = [f"q{lo}-{hi - 1}" for lo, hi in Q_BINS]
    spectrum = {}
    short, long = slice(0, SHORT // INSERT_STEP), slice(LONG // INSERT_STEP, None)
    for mate in range(2):
        for name, part in (("short", short), ("long", long)):
            t = strata[part, :, :, mate].sum((0, 1, 2))  # genome base x read base
            spectrum[f"mate{mate + 1}_{name}"] = {
                f"{'ACGT'[a]}>{'ACGT'[b]}": float(t[a, b] / t[a].sum()) if t[a].sum() else None
                for a in range(4)
                for b in range(4)
                if a != b
            }
    return {
        "bins": list(INSERT_BINS),
        "bases": bases.sum((1, 2, 3)).tolist(),
        "crude": _rates(bases.sum((1, 2, 3)), errors.sum((1, 2, 3))),
        "standardised": standard,
        "by_q": {q: _rates(bases[:, i].sum((1, 2)), errors[:, i].sum((1, 2))) for i, q in enumerate(q_names)},
        "spectrum": spectrum,
        "short_long": [f"< {SHORT}", f">= {LONG}"],
    }


def diagnose(sam: Path, ref: dict[str, str], min_mapq: int = 20, max_alt_freq: float = 0.2) -> dict[str, object]:
    masks = variant_masks(sam, ref, min_mapq, max_alt_freq)
    n: Counter[str] = Counter()
    ops: dict[str, Counter[str]] = {k: Counter() for k in ("one", "shared", "split")}
    qs: dict[str, Counter[str]] = {k: Counter() for k in ("match", "one", "shared", "split")}
    # every aligned base of a concordant pair, inside the overlap or outside it: [bases, errors] by Q bin
    region: dict[str, dict[str, list[int]]] = {"inside": {}, "outside": {}}
    # a `one` error: the Q of the mate that erred against the Q of the mate that did not
    pairs_q: Counter[str] = Counter()
    # every aligned base of a proper pair, 10 or more from either alignment end, in read orientation:
    # insert bin x Q group x cycle quarter x mate x genome base x read base
    strata = np.zeros((len(INSERT_BINS), len(Q_EDGES) + 1, 4, 2, 4, 4), np.int64)
    mates: dict[str, pysam.AlignedSegment] = {}

    def mask_of(x: pysam.AlignedSegment) -> np.ndarray | None:
        return masks.get(x.reference_name or "")

    for r in reads(sam, min_mapq):
        other = mates.pop(r.query_name, None)
        if other is None:
            mates[r.query_name] = r
            continue
        r1, r2 = (r, other) if r.is_read1 else (other, r)
        if r1.reference_name != r2.reference_name or r1.is_reverse == r2.is_reverse:
            n["pairs_discordant"] += 1
            continue
        n["pairs"] += 1
        if r1.is_proper_pair and r1.template_length:
            size = min(abs(r1.template_length) // INSERT_STEP, len(INSERT_BINS) - 1)
            for mate, x in enumerate((r1, r2)):
                seq, qual, mask = x.query_sequence or "", x.query_qualities or [], mask_of(x)
                start, stop, length = x.query_alignment_start, x.query_alignment_end, len(seq)
                contig = ref[x.reference_name or ""]
                for qp, rp in x.get_aligned_pairs(matches_only=True):
                    g, b = contig[rp], seq[qp]
                    if g not in BASES or b not in BASES or min(qp - start, stop - 1 - qp) < 10:
                        continue
                    if mask is not None and mask[rp]:
                        continue
                    cycle = length - 1 - qp if x.is_reverse else qp
                    if x.is_reverse:
                        g, b = g.translate(COMP), b.translate(COMP)
                    strata[
                        size, int(np.digitize(qual[qp], Q_EDGES)), 4 * cycle // length, mate, CODE4[g], CODE4[b]
                    ] += 1
        s1, s2 = _sites(r1), _sites(r2)
        both = s1.keys() & s2.keys()
        n["pairs_overlapping"] += bool(both)
        contig, mask = ref[r1.reference_name], masks.get(r1.reference_name)
        for sites in (s1, s2):
            for site, (b, q) in sites.items():
                g = contig[site]
                if g in BASES and b in BASES and not (mask is not None and mask[site]):
                    cell = region["inside" if site in both else "outside"].setdefault(q_bin(q), [0, 0])
                    cell[0] += 1
                    cell[1] += b != g
        for site in both:
            g = contig[site]
            (b1, q1), (b2, q2) = s1[site], s2[site]
            if g not in BASES or b1 not in BASES or b2 not in BASES:
                continue
            if mask is not None and mask[site]:
                n["sites_masked"] += 1
                continue
            n["sites"] += 1
            if b1 == g and b2 == g:
                qs["match"][q_bin(min(q1, q2))] += 1
            elif b1 == b2:
                n["shared"] += 1
                ops["shared"][_op(g, b1, r1.is_reverse)] += 1
                qs["shared"][q_bin(min(q1, q2))] += 1
            elif b1 != g and b2 != g:
                n["split"] += 1
                ops["split"][_op(g, b1, r1.is_reverse)] += 1
                qs["split"][q_bin(min(q1, q2))] += 1
            else:
                n["one"] += 1
                b, q, rev, q_right = (b1, q1, r1.is_reverse, q2) if b1 != g else (b2, q2, r2.is_reverse, q1)
                ops["one"][_op(g, b, rev)] += 1
                qs["one"][q_bin(q)] += 1
                pairs_q[f"{q_bin(q)} erred, {q_bin(q_right)} right"] += 1
    # errors per read base: a shared or split site is an error in both reads, `one` in one of them
    total = n["one"] + 2 * (n["shared"] + n["split"])
    seen = n["one"] + 2 * n["split"]
    return {
        "counts": dict(n),
        "error_per_read_base": total / max(2 * n["sites"], 1),
        "seen_by_overlap": seen / max(total, 1),
        "shared_per_site": n["shared"] / max(n["sites"], 1),
        "ops": {k: dict(v.most_common()) for k, v in ops.items()},
        "q": {k: dict(sorted(v.items())) for k, v in qs.items()},
        "rate_by_q": {
            k: {b: {"bases": x[0], "errors": x[1], "rate": x[1] / x[0]} for b, x in sorted(v.items())}
            for k, v in region.items()
        },
        "one_by_q_pair": dict(pairs_q.most_common()),
        "insert": by_insert(strata),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("sam", type=Path)
    p.add_argument("genomes", type=Path)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    sam = a.sam
    if sam.suffix == ".gz":  # pysam reads SAM text, not gzipped SAM
        sam = sam.with_suffix("")
        if not sam.exists():
            with gzip.open(a.sam, "rb") as f, sam.open("wb") as g:
                while chunk := f.read(1 << 24):
                    g.write(chunk)
    out = diagnose(sam, genome(a.genomes), a.min_mapq)
    text = json.dumps(out, indent=1)
    if a.output:
        a.output.write_text(text)
    c = out["counts"]
    print(
        f"{c.get('sites', 0):,} overlapped sites; errors per read base {out['error_per_read_base']:.2e}; "
        f"one {c.get('one', 0):,}, shared {c.get('shared', 0):,}, split {c.get('split', 0):,}; "
        f"seen by the overlap {out['seen_by_overlap']:.2f}"
    )
    for region, rows in out["rate_by_q"].items():
        print(f"{region:>8}: " + ", ".join(f"{b} {x['rate']:.2e}" for b, x in rows.items()))
    ins = out["insert"]

    def row(values: list[float | None]) -> str:
        return " ".join(f"{v:.1e}" if v is not None else "   -   " for v in values)

    print("insert            " + " ".join(f"{b:>7}" for b in ins["bins"]))
    print("bases             " + " ".join(f"{n:>7.0e}" for n in ins["bases"]))
    print("crude             " + row(ins["crude"]))
    print("std Q*cycle*mate  " + row(ins["standardised"]))
    for q, values in ins["by_q"].items():
        print(f"{q:<18}" + row(values))
    for mate in ("mate1", "mate2"):
        short, long = ins["spectrum"][f"{mate}_short"], ins["spectrum"][f"{mate}_long"]
        fold = {op: long[op] / short[op] for op in short if short[op] and long[op]}
        top = sorted(fold, key=lambda op: -fold[op])[:4]
        print(f"{mate} long/short insert, by op: " + ", ".join(f"{op} {fold[op]:.1f}x ({long[op]:.1e})" for op in top))


if __name__ == "__main__":
    main()
