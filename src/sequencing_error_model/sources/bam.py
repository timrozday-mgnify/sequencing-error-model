"""`reference` mode (plan §1.1, phase 5): BAM/CRAM + reference FASTA → head E tuples and a fitted spec.

Each primary, mapped, non-duplicate read with qualities and MAPQ ≥ `min_mapq` becomes a (template, read, mate)
triple in the read's own orientation (reverse-strand records are reverse-complemented back). The triples go
through `generate.observations`, the path the recovery harness uses on the generator's own CIGARs, so a BAM
written from generated reads reproduces those tuples exactly.

ponytail: soft-clipped bases are dropped, so read positions count from the first aligned base and clipped
ends contribute no rows; records with N or P ops are skipped, as is an insertion after the last aligned base.
Site masking, coverage filters, per-contig reports and per-read trajectories are still to come (phase 5).
"""

import argparse
import json
import sys
from collections import Counter
from collections.abc import Iterator, Sequence
from dataclasses import replace
from itertools import islice
from pathlib import Path
from typing import Any

import numpy as np
import pysam

from sequencing_error_model.fit import error, quality
from sequencing_error_model.generate import Read, _revcomp, observations
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.spec import Component, ErrorModelSpec

_OPS = "MIDNSHP=X"
_NO_TRACK = np.zeros(0, np.int64)


def records(bam: Path, reference: Path, min_mapq: int = 20) -> Iterator[tuple[str, Read, int]]:
    """(template, read, mate) per usable record, in file order."""
    with (
        pysam.AlignmentFile(str(bam), reference_filename=str(reference)) as aln,
        pysam.FastaFile(str(reference)) as fasta,
    ):
        for r in aln.fetch(until_eof=True):
            ops, seq, qual, end = (
                r.cigartuples,
                r.query_alignment_sequence,
                r.query_alignment_qualities,
                r.reference_end,
            )
            if (
                r.is_unmapped
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
            template, q = fasta.fetch(r.reference_name, r.reference_start, end).upper(), list(qual)
            if r.is_reverse:
                template, seq, q, cigar = _revcomp(template), _revcomp(seq), q[::-1], cigar[::-1]
            if cigar and cigar[-1][0] == "I":
                n = cigar.pop()[1]
                seq, q = seq[:-n], q[:-n]
            quals = "".join(chr(v + 33) for v in q)
            yield template, Read(seq, quals, "".join(f"{n}{o}" for o, n in cigar), _NO_TRACK), 2 if r.is_read2 else 1


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
) -> ErrorModelSpec:
    """Both heads from exact-position tuples: head E from every draw, head Q from bases that emit a read base."""
    error_head = error.fit(table, error_tokens, alphabet)
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
    args = p.parse_args(argv)
    flank, m = window(args.error_tokens, args.quality_tokens)
    table = observations(islice(records(args.bam, args.reference, args.min_mapq), args.max_reads), flank, m, "bam")
    if not table.counts:
        p.error("no usable aligned reads")
    alphabet = sorted({int(k[0]) for k in table.counts})  # type: ignore[call-overload]
    provenance = {"mode": "reference", "sources": [str(args.bam), str(args.reference)], "min_mapq": args.min_mapq}
    spec = fit(table, alphabet, args.error_tokens, args.quality_tokens, provenance)
    spec.save(args.output)
    print(json.dumps({"rows": sum(table.counts.values()), "quality_alphabet": alphabet}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
