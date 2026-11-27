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


def diagnose(sam: Path, ref: dict[str, str], min_mapq: int = 20, max_alt_freq: float = 0.2) -> dict[str, object]:
    masks = variant_masks(sam, ref, min_mapq, max_alt_freq)
    n: Counter[str] = Counter()
    ops: dict[str, Counter[str]] = {k: Counter() for k in ("one", "shared", "split")}
    qs: dict[str, Counter[str]] = {k: Counter() for k in ("match", "one", "shared", "split")}
    # every aligned base of a concordant pair, inside the overlap or outside it: [bases, errors] by Q bin
    region: dict[str, dict[str, list[int]]] = {"inside": {}, "outside": {}}
    # a `one` error: the Q of the mate that erred against the Q of the mate that did not
    pairs_q: Counter[str] = Counter()
    insert: dict[int, list[int]] = {}
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
            cell = insert.setdefault(min(abs(r1.template_length) // 50 * 50, 700), [0, 0])
            for x in (r1, r2):
                seq, qual, start, stop = (
                    x.query_sequence or "",
                    x.query_qualities or [],
                    x.query_alignment_start,
                    x.query_alignment_end,
                )
                for qp, rp in x.get_aligned_pairs(matches_only=True):
                    if (
                        qual[qp] >= 36
                        and min(qp - start, stop - 1 - qp) >= 10
                        and not (mask_of(x) is not None and mask_of(x)[rp])
                    ):
                        g, b = ref[x.reference_name][rp], seq[qp]
                        if g in BASES and b in BASES:
                            cell[0] += 1
                            cell[1] += b != g
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
        "q36_by_insert": {f"{k}+": {"bases": v[0], "rate": v[1] / max(v[0], 1)} for k, v in sorted(insert.items())},
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
    print("Q>=36 by insert: " + ", ".join(f"{k} {x['rate']:.1e}" for k, x in out["q36_by_insert"].items()))


if __name__ == "__main__":
    main()
