"""`reference` mode (plan §1.1, phase 5): BAM/CRAM + reference FASTA → head E tuples and a fitted spec.

Each primary, mapped, non-duplicate read with qualities and MAPQ ≥ `min_mapq` becomes a (template, read, mate)
triple in the read's own orientation (reverse-strand records are reverse-complemented back). The triples go
through `generate.observations`, the path the recovery harness uses on the generator's own CIGARs, so a BAM
written from generated reads reproduces those tuples exactly.

Soft-clipped bases (and an insertion after the last aligned base) produce no rows, but they count toward read
positions and fill Q windows at the aligned ends; `strand` is the alignment's. Hard clips are invisible, so
positions of hard-clipped reads start at the first retained base.

Optional site masks and contig filters come from a first pass over the same records (`pileup`, `masks`). A
site is masked when one non-reference allele (base, deletion, or insertion before the site) is seen at least
twice and at a frequency ≥ `max_alt_freq`; this catches reference errors (frequency ≈ 1) and minor alleles.
Masking selects on the outcome, so it can drop true error hotspots (§6.6); it is off by default until its
clonal cost is measured. A contig is dropped when its mean depth is below `min_contig_depth`. `masks` also
returns a per-contig report (depth, raw error rate, masked sites) so outlier contigs show.

ponytail: records with N or P ops are skipped. Per-read trajectories are still to come (phase 5).
"""

import argparse
import csv
import json
import sys
from collections import Counter
from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import replace
from itertools import islice
from pathlib import Path
from typing import Any

import numpy as np
import pysam

from sequencing_error_model.fit import error, indel, quality
from sequencing_error_model.fit.quality import Array
from sequencing_error_model.generate import Read, _revcomp, align, indel_events, observations
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.spec import Component, ErrorModelSpec

_OPS = "MIDNSHP=X"
_NO_TRACK = np.zeros(0, np.int64)
_BASE = {"A": 0, "C": 1, "G": 2, "T": 3}  # pileup columns; 4 deletion, 5 insertion before the site
_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def _aligned(bam: Path, reference: Path, min_mapq: int) -> Iterator[tuple[str, int, int, bool, str, Read, int]]:
    """(contig, start, end, reverse, template, read, mate) per usable record, in read orientation."""
    with (
        pysam.AlignmentFile(str(bam), reference_filename=str(reference)) as aln,
        pysam.FastaFile(str(reference)) as fasta,
    ):
        for r in aln.fetch(until_eof=True):
            ops, seq, qual, end = r.cigartuples, r.query_alignment_sequence, r.query_qualities, r.reference_end
            contig = r.reference_name
            if (
                contig is None
                or r.is_unmapped
                or r.is_secondary
                or r.is_supplementary
                or r.is_qcfail
                or r.is_duplicate
                or r.mapping_quality < min_mapq
                or ops is None
                or seq is None
                or qual is None
                or end is None
            ):
                continue
            cigar = [("M" if o > 6 else _OPS[o], n) for o, n in ops if o not in (4, 5)]  # =/X → M, drop clips
            if any(o in "NP" for o, _ in cigar):
                continue
            template = fasta.fetch(contig, r.reference_start, end).upper()
            a, b, full = r.query_alignment_start, r.query_alignment_end, list(qual)
            before, q, after = full[:a], full[a:b], full[b:]
            if r.is_reverse:
                template, seq, cigar = _revcomp(template), _revcomp(seq), cigar[::-1]
                before, q, after = after[::-1], q[::-1], before[::-1]
            if cigar and cigar[-1][0] == "I":
                n = cigar.pop()[1]
                seq, q, after = seq[:-n], q[:-n], q[-n:] + after
            read = Read(
                seq,
                _phred(q),
                "".join(f"{n}{o}" for o, n in cigar),
                _NO_TRACK,
                (_phred(before), _phred(after)),
                "-" if r.is_reverse else "+",
            )
            yield contig, r.reference_start, end, r.is_reverse, template, read, 2 if r.is_read2 else 1


def records(
    bam: Path, reference: Path, min_mapq: int = 20, masked: dict[str, Array] | None = None
) -> Iterator[tuple[str, Read, int]]:
    """(template, read, mate) per usable record, in file order. With `masked` (a boolean site mask per kept
    contig, from `masks`), reads on other contigs are skipped and masked sites produce no rows."""
    return apply_masks(_aligned(bam, reference, min_mapq), masked)


Alignment = tuple[str, int, int, bool, str, Read, int]  # contig, start, end, reverse, template, read, mate


def apply_masks(alignments: Iterable[Alignment], masked: dict[str, Array] | None) -> Iterator[tuple[str, Read, int]]:
    """`records` from alignments in read orientation (a BAM's, or generated reads placed on a genome)."""
    for contig, start, end, reverse, template, read, mate in alignments:
        if masked is not None:
            if contig not in masked:
                continue
            sites = np.flatnonzero(masked[contig][start:end])
            # ponytail: a masked site drops the draws at its template base, so an insertion variant is masked
            # on the strand whose insertions attach to that base; mask a gap by both flanking bases if it matters.
            read = replace(read, masked=frozenset((end - start - 1 - sites if reverse else sites).tolist()))
        yield template, read, mate


def _phred(q: Sequence[int]) -> str:
    return "".join(chr(v + 33) for v in q)


def pileup(bam: Path, reference: Path, min_mapq: int = 20) -> dict[str, Array]:
    """Per contig, (length + 1) x 6 allele counts in reference orientation: A, C, G, T, deletion, insertion
    before the site. Uses the same records and CIGAR walk as `records`."""
    with pysam.FastaFile(str(reference)) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths, strict=True))
    return count_alleles(_aligned(bam, reference, min_mapq), lengths)


def count_alleles(alignments: Iterable[Alignment], lengths: dict[str, int]) -> dict[str, Array]:
    """`pileup` from alignments in read orientation and contig lengths."""
    counts: dict[str, Array] = {}
    for contig, start, end, reverse, template, read, _ in alignments:
        n = end - start
        sites, cols = [], []
        for t, op in align(template, read)[0]:
            if op[0] == "-":  # insertion before template base t: the gap before (forward) or after (reverse) it
                site, col = (start + n - t if reverse else start + t), 5
            else:
                base = template[t] if op == "=" else op[-1]
                site, col = (start + n - 1 - t if reverse else start + t), 4 if base == "-" else _BASE.get(base, -1)
                if col < 0:
                    continue
                if reverse and col < 4:
                    col = _BASE[_COMPLEMENT["ACGT"[col]]]
            sites.append(site)
            cols.append(col)
        table = counts.setdefault(contig, np.zeros((lengths[contig] + 1, 6), np.int64))
        np.add.at(table, (sites, cols), 1)
    return counts


def masks(
    counts: dict[str, Array], reference: Path, max_alt_freq: float | None = None, min_contig_depth: float = 0.0
) -> tuple[dict[str, Array], list[dict[str, Any]]]:
    """Boolean site masks for kept contigs, and one report row per contig with aligned reads."""
    with pysam.FastaFile(str(reference)) as fasta:
        return site_masks(counts, fasta.fetch, max_alt_freq, min_contig_depth)


def site_masks(
    counts: dict[str, Array],
    sequence: Callable[[str], str],
    max_alt_freq: float | None = None,
    min_contig_depth: float = 0.0,
) -> tuple[dict[str, Array], list[dict[str, Any]]]:
    """`masks` with contig sequences from `sequence(contig)`."""
    kept: dict[str, Array] = {}
    report = []
    for contig, table in counts.items():
        ref = np.array([_BASE.get(b, -1) for b in sequence(contig).upper()] + [-1])
        depth = table[:, :5].sum(axis=1)
        alt = table.copy()
        has_ref = ref >= 0
        alt[np.flatnonzero(has_ref), ref[has_ref]] = 0
        events = int(alt.sum())
        top = alt.max(axis=1)
        site_mask = np.zeros(len(ref), bool) if max_alt_freq is None else (top >= 2) & (top >= max_alt_freq * depth)
        mean_depth = float(depth.sum()) / (len(ref) - 1)
        keep = mean_depth >= min_contig_depth
        if keep:
            kept[contig] = site_mask
        report.append(
            {
                "contig": contig,
                "length": len(ref) - 1,
                "mean_depth": round(mean_depth, 3),
                "error_rate": round(events / max(1, int(depth.sum() + table[:, 5].sum())), 6),
                "masked_sites": int(site_mask.sum()),
                "kept": keep,
            }
        )
    return kept, report


def window(error_tokens: Sequence[str], quality_tokens: Sequence[str]) -> tuple[tuple[int, int], int]:
    """Context flank (at least 2 each side, for `Homopolymer`) and Q window half-width both heads need."""
    contexts = [c.args for c in map(Component, (*error_tokens, *quality_tokens)) if c.name == "Context"]
    flank = (max([2, *(a[0] for a in contexts)]), max([2, *(a[1] for a in contexts)]))
    return flank, max(Component(error_tokens[0]).args[0], Component(quality_tokens[0]).args[0])


def fit(
    table: CountTable,
    alphabet: Sequence[int],
    error_tokens: Sequence[str],
    quality_tokens: Sequence[str],
    provenance: dict[str, Any],
    indels: CountTable | None = None,
) -> ErrorModelSpec:
    """Both heads from exact-position tuples: head E from every draw, head Q from bases that emit a read base.
    An `IndelLength` token is fitted from `indels` (`generate.indel_events`); `table` then needs per-event rows."""
    body = [t for t in error_tokens if Component(t).name != "IndelLength"]
    error_head = error.fit(table, body, alphabet)
    if lengths := [t for t in error_tokens if Component(t).name == "IndelLength"]:
        if indels is None:
            raise ValueError(f"{lengths[0]} needs an indel event table")
        error_head += (indel.fit(indels, lengths[0], alphabet),)
    oi = table.fields.index("op")
    # One row per template base with an emitted read base: matches and substitutions.
    final: Counter[Key] = Counter(
        {k: n for k, n in table.counts.items() if (op := str(k[oi])) == "=" or "-" not in (op[0], op[2])}
    )
    lags = [f"q-{i}" for i in range(1, Component(quality_tokens[0]).args[0] + 1)]
    q_table = replace(table, counts=final).marginal(*lags, "pos_start", "pos_end", "mate", "context", "q")
    quality_head = quality.fit(q_table, quality_tokens, alphabet)
    return ErrorModelSpec(tuple(alphabet), provenance, quality_head, error_head)


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m sequencing_error_model.sources.bam", description="Fit a spec from reads aligned to a reference."
    )
    p.add_argument("bam", type=Path, help="BAM/SAM/CRAM (any order)")
    p.add_argument("reference", type=Path, help="reference FASTA (indexed with .fai, or writable to index)")
    p.add_argument("--output", type=Path, required=True, metavar="SPEC_DIR")
    p.add_argument("--error-tokens", nargs="+", default=["QualityWindow(1)", "Context(1,1)", "Homopolymer", "Mate"])
    p.add_argument("--quality-tokens", nargs="+", default=["QualityMarkov(1)", "Position(48)", "Mate", "Context(1,1)"])
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--max-reads", type=int)
    p.add_argument("--mask-alt-freq", type=float, help="mask sites with a non-reference allele at this frequency")
    p.add_argument("--min-contig-depth", type=float, default=0.0, help="drop contigs below this mean depth")
    p.add_argument("--contig-report", type=Path, metavar="TSV", help="per-contig depth, raw error rate, masks")
    args = p.parse_args(argv)
    provenance: dict[str, Any] = {
        "mode": "reference",
        "sources": [str(args.bam), str(args.reference)],
        "min_mapq": args.min_mapq,
    }
    masked = None
    if args.mask_alt_freq is not None or args.min_contig_depth > 0 or args.contig_report:
        counts = pileup(args.bam, args.reference, args.min_mapq)
        masked, report = masks(counts, args.reference, args.mask_alt_freq, args.min_contig_depth)
        provenance |= {
            "mask_alt_freq": args.mask_alt_freq,
            "min_contig_depth": args.min_contig_depth,
            "masked_sites": sum(r["masked_sites"] for r in report if r["kept"]),
            "contigs_kept": len(masked),
            "contigs_dropped": len(report) - len(masked),
        }
        if args.contig_report:
            with args.contig_report.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(report[0]) if report else ["contig"], delimiter="\t")
                writer.writeheader()
                writer.writerows(report)
    flank, m = window(args.error_tokens, args.quality_tokens)
    reads = islice(records(args.bam, args.reference, args.min_mapq, masked), args.max_reads)
    per_event = any(Component(t).name == "IndelLength" for t in args.error_tokens)
    table = observations(reads, flank, m, "bam", per_event)
    indels = None
    if per_event:  # a second pass over the same records
        indels = indel_events(islice(records(args.bam, args.reference, args.min_mapq, masked), args.max_reads), "bam")
    if not table.counts:
        p.error("no usable aligned reads")
    # Q windows can reach clipped bases, whose Q never appears at a centre.
    qs = {v for k in table.counts for f, v in zip(table.fields, k, strict=True) if f.startswith("q") and v is not None}
    alphabet = sorted(int(v) for v in qs)  # type: ignore[call-overload]
    spec = fit(table, alphabet, args.error_tokens, args.quality_tokens, provenance, indels)
    spec.save(args.output)
    print(json.dumps({"rows": sum(table.counts.values()), "quality_alphabet": alphabet}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
