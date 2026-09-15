"""Cross-mode comparison on shared support (plan §7, phase 5).

Two modes observe different things: `pe-overlap` sees only substitutions, only at overlap positions, and not errors
both mates share (PCR, library, cluster); `reference` sees every op at every aligned base, but also reference
errors and strain variation. So comparisons are restricted to what both observe.

- `evidence(a, b, by)`: two labelled head E tables (`pe_overlap.table`, `generate.observations` of BAM records)
  restricted to substitution support (match and substitution rows), marginalised onto the covariates `by`, and
  kept at covariate values both observe with at least `min_exposure` rows each. Each rate is standardised to the
  pooled exposure over those values, so differences in coverage (overlap positions, Q mix) don't read as rate
  differences. `excess` = rate_b − rate_a: with a = `pe-overlap` and b = `reference` on the same reads it
  estimates PCR/library substitutions, plus any residual variation or reference error.
- `models(a, b, table)`: both specs' head E on one table, conditioned on no indel: op TV, rates, rate by Q, and
  each spec's log-likelihood per row.

CLI (`python -m sequencing_error_model.compare`): one run's R1/R2 and its BAM, with the specs the two source CLIs
fitted from them, give a JSON report of both levels (models on each mode's rows).
"""

import argparse
import json
import sys
from collections import Counter
from collections.abc import Sequence
from dataclasses import asdict, replace
from itertools import islice
from pathlib import Path
from typing import Any

import numpy as np

from sequencing_error_model import spec as spec_io
from sequencing_error_model.fit import error, indel
from sequencing_error_model.generate import _flank, observations
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.sources import bam, pe_overlap
from sequencing_error_model.spec import ErrorModelSpec


def _substitution_support(table: CountTable) -> CountTable:
    if "op" not in table.fields:
        raise ValueError(f"{table.source}: comparison needs op labels")
    oi = table.fields.index("op")
    kept: Counter[Key] = Counter({k: n for k, n in table.counts.items() if k[oi] == "=" or "-" not in str(k[oi])})
    return replace(table, counts=kept)


def _exposure(table: CountTable, by: Sequence[str]) -> dict[Key, tuple[float, float]]:
    """(rows, substitutions) per value of `by`, on substitution support."""
    if missing := sorted(set(by) - set(table.fields)):
        raise ValueError(f"{table.source}: no fields {missing}")
    t = _substitution_support(table)
    oi, idx = t.fields.index("op"), [t.fields.index(f) for f in by]
    out: dict[Key, tuple[float, float]] = {}
    for k, n in t.counts.items():
        rows, subs = out.get(key := tuple(k[i] for i in idx), (0.0, 0.0))
        out[key] = (rows + n, subs + n * (k[oi] != "="))
    return out


def evidence(a: CountTable, b: CountTable, by: Sequence[str] = ("q",), min_exposure: float = 100) -> dict[str, Any]:
    """Substitution rates of two tables on the covariate values both observe, standardised to pooled exposure."""
    ea, eb = _exposure(a, by), _exposure(b, by)
    shared = sorted((k for k in ea.keys() & eb.keys() if min(ea[k][0], eb[k][0]) >= min_exposure), key=str)
    if not shared:
        raise ValueError(f"{a.source} and {b.source} share no values of {list(by)} with {min_exposure} rows each")
    w = np.array([ea[k][0] + eb[k][0] for k in shared])
    ra, rb = (np.array([e[k][1] / e[k][0] for k in shared]) for e in (ea, eb))
    rate_a, rate_b = float(w @ ra / w.sum()), float(w @ rb / w.sum())
    return {
        "sources": [a.source, b.source],
        "by": list(by),
        "values": [str(k if len(k) != 1 else k[0]) for k in shared],
        "coverage": [sum(e[k][0] for k in shared) / sum(v[0] for v in e.values()) for e in (ea, eb)],
        "rate_a": rate_a,
        "rate_b": rate_b,
        "rate_ratio": rate_b / rate_a,
        "excess": rate_b - rate_a,
        "rates_a": ra.tolist(),
        "rates_b": rb.tolist(),
    }


def models(a: ErrorModelSpec, b: ErrorModelSpec, table: CountTable, by: str = "q") -> dict[str, Any]:
    """Both specs' head E on `table`'s substitution-support rows, each conditioned on no indel."""
    t = _substitution_support(table)
    oi, ci, f0 = t.fields.index("op"), t.fields.index("context"), t.meta.get("flank", (0, 0))[0]
    n = np.array(list(t.counts.values()), float)
    cat = np.array([error._category(k[oi], str(k[ci])[f0]) for k in t.counts])
    p = []
    for spec in (a, b):
        q = error.probabilities(indel.split(spec.error_head)[0], spec.quality_alphabet, t)[:, :5]
        p.append(q / q.sum(axis=1, keepdims=True))
    rates = [float(n @ (1 - x[:, 0]) / n.sum()) for x in p]
    key = np.array([str(k[t.fields.index(by)]) for k in t.counts])
    values = sorted(set(key), key=lambda v: (len(v), v))
    return {
        "rows": float(n.sum()),
        "op_tv": float(n @ (0.5 * np.abs(p[0] - p[1]).sum(axis=1)) / n.sum()),
        "rate_a": rates[0],
        "rate_b": rates[1],
        "rate_ratio": rates[1] / rates[0],
        "excess": rates[1] - rates[0],
        "log_likelihood_per_row": [float(n @ np.log(x[np.arange(len(n)), cat]) / n.sum()) for x in p],
        "by": by,
        "values": values,
        "rates_by": [[float(n[key == v] @ (1 - x[key == v, 0]) / n[key == v].sum()) for v in values] for x in p],
    }


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m sequencing_error_model.compare", description="Compare pe-overlap and reference on one run."
    )
    p.add_argument("r1", type=Path, help="mate 1 FASTQ[.gz]")
    p.add_argument("r2", type=Path, help="mate 2 FASTQ[.gz], in the same order")
    p.add_argument("bam", type=Path, help="the same reads aligned (BAM/SAM/CRAM)")
    p.add_argument("reference", type=Path, help="reference FASTA")
    p.add_argument("--overlap-spec", type=Path, required=True, help="spec from sources.pe_overlap on R1/R2")
    p.add_argument("--reference-spec", type=Path, required=True, help="spec from sources.bam on the BAM")
    p.add_argument("--output", type=Path, required=True, metavar="JSON")
    p.add_argument("--by", nargs="+", default=["q", "mate"], help="covariates for the evidence-level comparison")
    p.add_argument("--min-exposure", type=float, default=100)
    p.add_argument("--min-overlap", type=int, default=20)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--max-pairs", type=int)
    p.add_argument("--max-reads", type=int)
    p.add_argument("--mask-alt-freq", type=float, help="mask sites with a non-reference allele at this frequency")
    args = p.parse_args(argv)
    a, b = spec_io.load(args.overlap_spec), spec_io.load(args.reference_spec)
    # One window for both tables, wide enough for both specs' head E.
    heads = [a.error_head, b.error_head]
    flanks = [_flank(h) for h in heads]
    flank = (max(2, *(f[0] for f in flanks)), max(2, *(f[1] for f in flanks)))
    m = max(h[0].args[0] for h in heads)

    ev = pe_overlap.collect(
        islice(pe_overlap.read_pairs(args.r1, args.r2), args.max_pairs), m, flank, min_overlap=args.min_overlap
    )
    overlap = pe_overlap.table(ev, a.error_head)
    masked = None
    if args.mask_alt_freq is not None:
        masked = bam.masks(bam.pileup(args.bam, args.reference, args.min_mapq), args.reference, args.mask_alt_freq)[0]
    reads = islice(bam.records(args.bam, args.reference, args.min_mapq, masked), args.max_reads)
    reference = observations(reads, flank, m, "reference")
    if not overlap.counts or not reference.counts:
        p.error("no overlap rows" if not overlap.counts else "no usable aligned reads")

    report = {
        "sources": {"overlap_spec": str(args.overlap_spec), "reference_spec": str(args.reference_spec)},
        "overlap_stats": asdict(ev.stats),
        "mask_alt_freq": args.mask_alt_freq,
        "evidence": evidence(overlap, reference, args.by, args.min_exposure),
        "models": {"reference_rows": models(a, b, reference), "overlap_rows": models(a, b, overlap)},
    }
    args.output.write_text(json.dumps(report, indent=2))
    e = report["evidence"]
    print(json.dumps({k: e[k] for k in ("rate_a", "rate_b", "excess", "coverage")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
