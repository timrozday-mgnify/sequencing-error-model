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
- `components(a, b, table)`: two specs' head E compared per component on one table - predicted error rate by Q,
  position, GC, strand and mate, the slope of one's context log-odds on the other's, and op composition. The modes
  fit different component sets, so what is compared is what the heads predict on shared rows, not parameters.
- `generated(a, b, reference, ...)`: phase 3 metrics of reads generated from both specs, via `recovery.compare`.
- `skiver_evidence(analyze, table)`: the marginals `skiver analyze` reports - P(error | Q), the read-position
  curves, GC and the trinucleotide spectrum - recomputed from a truth-bearing table, so skiver's own evidence is
  checked against `pe-overlap` and `reference` truth independently of the `kmer` fitters (plan phase 6).

CLI (`python -m sequencing_error_model.compare`): one run's R1/R2 and its BAM, with the specs the two source CLIs
fitted from them, give a JSON report of both levels (models on each mode's rows), plus the skiver comparison
with `--skiver <analyze prefix>` and, with `--kmer-spec`, the `kmer` default spec against both of them per
component, on each mode's rows, and (with `--generated-reads`) on generated reads.
"""

import argparse
import json
import sys
from collections import Counter
from collections.abc import Hashable, Iterable, Sequence
from dataclasses import asdict, replace
from itertools import islice
from pathlib import Path
from typing import Any, cast

import numpy as np
import pysam

from sequencing_error_model import recovery
from sequencing_error_model import spec as spec_io
from sequencing_error_model.fit import error, indel
from sequencing_error_model.generate import _flank, observations
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.sources import bam, pe_overlap, skiver_analyze
from sequencing_error_model.sources.skiver_analyze import SkiverAnalyze
from sequencing_error_model.spec import ErrorModelSpec


def _substitution_support(table: CountTable, indels: bool = False) -> CountTable:
    """Matches and substitutions, or every op with `indels` (for a source that labels them, e.g. skiver)."""
    if "op" not in table.fields:
        raise ValueError(f"{table.source}: comparison needs op labels")
    oi = table.fields.index("op")
    kept: Counter[Key] = Counter(
        {k: n for k, n in table.counts.items() if indels or k[oi] == "=" or "-" not in str(k[oi])}
    )
    return replace(table, counts=kept)


def _exposure(table: CountTable, by: Sequence[str], indels: bool = False) -> dict[Key, tuple[float, float]]:
    """(rows, errors) per value of `by`, on the support `indels` selects."""
    if missing := sorted(set(by) - set(table.fields)):
        raise ValueError(f"{table.source}: no fields {missing}")
    t = _substitution_support(table, indels)
    oi, idx = t.fields.index("op"), [t.fields.index(f) for f in by]
    out: dict[Key, tuple[float, float]] = {}
    for k, n in t.counts.items():
        rows, subs = out.get(key := tuple(k[i] for i in idx), (0.0, 0.0))
        out[key] = (rows + n, subs + n * (k[oi] != "="))
    return out


def evidence(
    a: CountTable,
    b: CountTable,
    by: Sequence[str] = ("q",),
    min_exposure: float = 100,
    indels: bool = False,
) -> dict[str, Any]:
    """Substitution rates of two tables on the covariate values both observe, standardised to pooled exposure."""
    ea, eb = _exposure(a, by, indels), _exposure(b, by, indels)
    shared = sorted((k for k in ea.keys() & eb.keys() if min(ea[k][0], eb[k][0]) >= min_exposure), key=str)
    if not shared:
        raise ValueError(f"{a.source} and {b.source} share no values of {list(by)} with {min_exposure} rows each")
    w = np.array([ea[k][0] + eb[k][0] for k in shared])
    ra, rb = (np.array([e[k][1] / e[k][0] for k in shared]) for e in (ea, eb))
    rate_a, rate_b = float(w @ ra / w.sum()), float(w @ rb / w.sum())
    return {
        "sources": [a.source, b.source],
        "by": list(by),
        "support": "all ops" if indels else "substitutions",
        "values": [str(k if len(k) != 1 else k[0]) for k in shared],
        "coverage": [sum(e[k][0] for k in shared) / sum(v[0] for v in e.values()) for e in (ea, eb)],
        "rate_a": rate_a,
        "rate_b": rate_b,
        "rate_ratio": rate_b / rate_a,
        "excess": rate_b - rate_a,
        "rates_a": ra.tolist(),
        "rates_b": rb.tolist(),
    }


def _counts(items: Iterable[tuple[Key, float]]) -> Counter[Key]:
    """Sum possibly fractional counts per key, dropping non-positive ones."""
    out: dict[Key, float] = {}
    for key, n in items:
        out[key] = out.get(key, 0.0) + n
    return cast("Counter[Key]", Counter({k: v for k, v in out.items() if v > 0}))


def _gc10(lo: int, hi: int) -> tuple[int, int]:
    """skiver's GC bin as the 10 % bin `generate.gc_bin` labels reads with."""
    if lo // 10 != (hi - 1) // 10:
        raise ValueError(f"GC bin [{lo}, {hi}) crosses a 10 % boundary, so it cannot be pooled into a 10 % bin")
    return (lo // 10 * 10, lo // 10 * 10 + 10)


def skiver_marginals(a: SkiverAnalyze) -> dict[str, CountTable]:
    """skiver's rate marginals as head-E-shaped tables, at the rate skiver *reports* (plan §5.6, phase 6).

    `summary_phred.csv` and `summary_gc_content.csv` carry a per-Q/per-bin Weibull rate beside scan counts whose
    ratio is a hazard, so their counts are rebuilt as (matches, errors) at the reported rate over the same
    exposure; the GC bins are pooled into the 10 % bins `generate.gc_bin` uses. `summary_read_position.csv` has
    no fitted rate, so its own counts stand: the exposure at position p is the values that survived to p, and
    the ratio is the error rate there *given survival*, which is the per-base rate when errors do not cluster.

    First-error stopping itself cannot be emulated from a per-base table, which has no value grouping, so the
    comparison is of rates (this is the "otherwise compare rates" case of the phase 6 evidence check).
    """

    def table(name: str, field: str, rows: Iterable[tuple[Hashable, float, float]]) -> CountTable:
        items = (((k, op), n) for k, rate, e in rows for op, n in (("=", e * (1 - rate)), ("!", e * rate)))
        return CountTable(f"skiver_analyze:{name}", (field, "op"), "base", True, _counts(items), {"k": a.k, "v": a.v})

    out = {
        "phred": table("phred", "q", ((b.lo, b.per_base_error_rate, b.num_correct + b.num_error) for b in a.phred)),
        "gc_content": table(
            "gc_content",
            "gc",
            ((_gc10(b.lo, b.hi), b.per_base_error_rate, b.num_correct + b.num_error) for b in a.gc_content),
        ),
    }
    for name, field, start in (("read_position_start", "pos_start", True), ("read_position_end", "pos_end", False)):
        rows = [r for r in a.read_position if r.from_start is start]
        out[name] = table(
            name,
            field,
            ((r.index, r.num_error / n, float(n)) for r in rows if (n := r.num_correct + r.num_error)),
        )
    return out


def _trinucleotide(t: CountTable) -> Counter[Key]:
    """(trinucleotide context, op) error counts, the table's context trimmed to skiver's one base each side.

    An insertion row is keyed at its insertion point, as skiver keys it: the table centres the row on the
    template base *after* the insertion, so its trinucleotide is (previous base, "-", that base).
    """
    left, right = t.meta.get("flank", (0, 0))
    if left < 1 or right < 1:
        raise ValueError(f"{t.source}: needs one flanking base each side to compare with skiver's spectrum")
    ci, oi = t.fields.index("context"), t.fields.index("op")
    return _counts(
        (
            ((s[left - 1] + "-" + s[left] if str(k[oi]).startswith("-") else s[left - 1 : left + 2], k[oi]), n)
            for k, n in t.counts.items()
            if k[oi] != "=" and (s := str(k[ci]))
        )
    )


def spectrum(a: SkiverAnalyze, t: CountTable, min_count: float = 0.0) -> dict[str, Any]:
    """`summary_error_spectrum.csv` against a truth-bearing table: shares of the error mass, no exposure.

    Restricted to the ops both label (`pe-overlap` identifies substitutions only) and to trinucleotide context,
    which is all skiver reports; a wider table's context is trimmed to it.
    """
    sk = _counts(((skiver_analyze.context(r.op, r.prev_base, r.next_base), r.op), float(r.total)) for r in a.spectrum)
    ours = _trinucleotide(t)
    ops = {str(k[1]) for k in sk} & {str(k[1]) for k in ours}
    keys = sorted(
        (k for k in sk.keys() | ours.keys() if str(k[1]) in ops and sk.get(k, 0.0) + ours.get(k, 0.0) >= min_count),
        key=str,
    )
    if not keys:
        raise ValueError(f"{t.source} and summary_error_spectrum.csv share no error operations")
    na, nb = (np.array([c.get(k, 0.0) for k in keys]) for c in (sk, ours))
    pa, pb = na / na.sum(), nb / nb.sum()
    return {
        "sources": ["skiver_analyze:spectrum", t.source],
        "ops": sorted(ops),
        "keys": [f"{c}:{op}" for c, op in keys],
        "coverage": [float(n.sum() / sum(c.values())) for n, c in ((na, sk), (nb, ours))],
        "tv": float(0.5 * np.abs(pa - pb).sum()),
        "r": float(np.corrcoef(pa, pb)[0, 1]) if len(keys) > 1 else 1.0,
        "shares_a": pa.tolist(),
        "shares_b": pb.tolist(),
    }


def _hazard(a: SkiverAnalyze, t: CountTable) -> dict[str, Any]:
    """skiver's fitted rate, Weibull and op proportions beside the table's own, as `fit.kmer.hazard` carries them.

    skiver's hazard is over positions inside a k-mer value, an exposure a read table cannot reproduce, and
    beta < 1 is clustering default mode does not identify (§6.2), so the two are reported side by side rather
    than compared.
    """
    oi = t.fields.index("op")
    n = sum(t.counts.values())
    kinds = {"substitution": 0.0, "insertion": 0.0, "deletion": 0.0}
    for k, c in t.counts.items():
        op = str(k[oi])
        if op != "=":
            kinds["insertion" if op.startswith("-") else "deletion" if op.endswith("-") else "substitution"] += c
    errors = sum(kinds.values())
    r = a.error_rate
    return {
        "skiver_per_base_error_rate": r.per_base_error_rate,
        "skiver_mean_hazard_rate": r.mean_hazard_rate,
        "skiver_weibull": [r.lambda_, r.beta],
        "skiver_op_proportions": {
            "substitution": r.substitution_error_proportion,
            "insertion": r.insertion_error_proportion,
            "deletion": r.deletion_error_proportion,
        },
        "table_rate": errors / n if n else 0.0,
        "table_op_proportions": {k: (v / errors if errors else 0.0) for k, v in kinds.items()},
        "clustering": "not comparable: skiver's hazard is per value position, and beta is clustering no mode here fits",
    }


def skiver_evidence(a: SkiverAnalyze, t: CountTable, min_exposure: float = 100) -> dict[str, Any]:
    """Every marginal skiver reports, recomputed from a truth-bearing table (plan phase 6, evidence check).

    Tests skiver's evidence independently of the `kmer` fitters: the same reads seen through `pe-overlap` or
    `reference` truth should reproduce skiver's P(error | Q), position curves, GC curve and error spectrum.
    Indels are kept on both sides, so against a substitutions-only table (`pe-overlap`) skiver's rates carry
    indel errors the table cannot see; `support` and `coverage` record what each comparison stood on. A
    marginal with no shared values is reported as `skipped` rather than failing the whole report.
    """
    out: dict[str, Any] = {}
    for name, m in skiver_marginals(a).items():
        try:
            out[name] = evidence(m, t, (m.fields[0],), min_exposure, indels=True)
        except ValueError as e:
            out[name] = {"skipped": str(e)}
    try:
        out["spectrum"] = spectrum(a, t)
    except ValueError as e:
        out["spectrum"] = {"skipped": str(e)}
    out["hazard"] = _hazard(a, t)
    return out


def models(a: ErrorModelSpec, b: ErrorModelSpec, table: CountTable, by: str = "q") -> dict[str, Any]:
    """Both specs' head E on `table`'s substitution-support rows, each conditioned on no indel."""
    t = _substitution_support(table)
    oi, ci, f0 = t.fields.index("op"), t.fields.index("context"), t.meta.get("flank", (0, 0))[0]
    n = np.array(list(t.counts.values()), float)
    cat = np.array([error._category(k[oi], str(k[ci])[f0]) for k in t.counts])
    p = _heads(a, b, t, indels=False)
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


def _values(key: "np.ndarray[Any, Any]") -> list[str]:
    return sorted(set(key.tolist()), key=lambda v: (len(v), v))


def _rates(p: "list[np.ndarray[Any, Any]]", n: "np.ndarray[Any, Any]", sel: "np.ndarray[Any, Any]") -> list[float]:
    """Exposure-weighted P(error) under each spec over the selected rows."""
    return [float(n[sel] @ (1 - x[sel, 0]) / n[sel].sum()) for x in p]


def _heads(a: ErrorModelSpec, b: ErrorModelSpec, t: CountTable, indels: bool) -> "list[np.ndarray[Any, Any]]":
    """Each spec's head E on `t`'s rows, renormalised over the categories the comparison keeps."""
    out = []
    for s in (a, b):
        q = error.probabilities(indel.split(s.error_head)[0], s.quality_alphabet, t)
        q = q if indels else q[:, :5]
        out.append(q / q.sum(axis=1, keepdims=True))
    return out


def components(
    a: ErrorModelSpec,
    b: ErrorModelSpec,
    table: CountTable,
    by: Sequence[str] = ("q", "pos_start", "pos_end", "gc", "strand", "mate"),
    min_exposure: float = 100,
    indels: bool = False,
) -> dict[str, Any]:
    """Two specs' head E compared per component on one table's rows (plan phase 6, model check).

    Not a parameter comparison: the modes fit different component sets (the `kmer` default head has no mate
    term, a `pe-overlap` head no indels), and a shared covariate can be reached by different tokens. What is
    compared is what the heads *predict* on the same rows, exposure-weighted by the table and marginalised onto
    each covariate, so a ratio away from 1 is that covariate's effect moved, as in `recovery.indel_components`.

    - each field of `by` present on the table: predicted error rate per value, and b over a;
    - `context`: the least-squares slope and correlation of b's error log-odds on a's over contexts with at
      least `min_exposure` rows, centred, as `recovery._slope` does for `Context` parameters (1 is unbiased,
      below 1 a shrunk effect), plus the mean absolute log ratio;
    - `op`: expected shares of the error mass over `error.CATEGORIES`, and their total variation.
    """
    t = _substitution_support(table, indels)
    if not t.counts:
        raise ValueError(f"{table.source}: no rows on the comparison's support")
    n = np.array(list(t.counts.values()), float)
    p = _heads(a, b, t, indels)
    out: dict[str, Any] = {
        "sources": [str(s.provenance.get("mode", "?")) for s in (a, b)],
        "rows": float(n.sum()),
        "support": "all ops" if indels else "substitutions",
    }
    for f in by:
        if f not in t.fields:
            continue
        key = np.array([str(k[t.fields.index(f)]) for k in t.counts])
        values = [v for v in _values(key) if n[key == v].sum() >= min_exposure]
        if not values:
            continue
        pairs = [_rates(p, n, key == v) for v in values]
        ra, rb = [r[0] for r in pairs], [r[1] for r in pairs]
        out[f] = {
            "values": values,
            "rate_a": ra,
            "rate_b": rb,
            "ratio": [y / x if x else float("nan") for x, y in zip(ra, rb, strict=True)],
        }
    out["context"] = _context_shape(p, n, t, min_exposure)
    mass = [n @ y[:, 1:] for y in p]
    shares = [m / m.sum() for m in mass]
    out["op"] = {
        "categories": list(error.CATEGORIES[1 : 1 + p[0].shape[1] - 1]),
        "shares_a": shares[0].tolist(),
        "shares_b": shares[1].tolist(),
        "tv": float(0.5 * np.abs(shares[0] - shares[1]).sum()),
    }
    return out


def _context_shape(
    p: "list[np.ndarray[Any, Any]]", n: "np.ndarray[Any, Any]", t: CountTable, min_exposure: float
) -> dict[str, Any]:
    """Slope, correlation and mean absolute log ratio of b's context error log-odds on a's, both centred."""
    key = np.array([str(k[t.fields.index("context")]) for k in t.counts])
    values = [v for v in _values(key) if n[key == v].sum() >= min_exposure]
    if len(values) < 2:
        return {"skipped": f"{t.source}: fewer than 2 contexts with {min_exposure} rows"}
    rates = np.array([_rates(p, n, key == v) for v in values])  # [contexts, 2]
    lo = np.log(np.clip(rates, 1e-12, 1 - 1e-12) / np.clip(1 - rates, 1e-12, 1))
    x, y = (lo[:, i] - lo[:, i].mean() for i in (0, 1))
    return {
        "contexts": len(values),
        "log_odds_slope": float(x @ y / (x @ x)) if x @ x else float("nan"),
        "r": float(np.corrcoef(x, y)[0, 1]),
        "mean_abs_log_ratio": float(np.abs(np.log(rates[:, 1] / rates[:, 0])).mean()),
    }


def generated(
    a: ErrorModelSpec,
    b: ErrorModelSpec,
    reference: Path,
    n_reads: int,
    read_length: int,
    seed: int = 0,
    q_reads: int = 20_000,
) -> dict[str, Any]:
    """Phase 3 metrics of reads generated from both specs, on templates sliced from the real reference.

    `recovery.compare` treats `a` as the truth, so `rate_ratio` and the curve pairs read as b over a. Both
    heads Q come from the same FASTQ, so the Q-track metrics are a check that the generator applies each head,
    not a difference between the modes; the error metrics are the comparison.
    """
    with pysam.FastaFile(str(reference)) as fasta:
        genome = max((fasta.fetch(c) for c in fasta.references), key=len).upper()
    if len(genome) < read_length + 1:
        raise ValueError(f"{reference}: longest contig is shorter than one {read_length} bp read")
    rng = np.random.default_rng(seed)
    starts = rng.integers(0, len(genome) - read_length, n_reads)
    templates = [genome[s : s + read_length] for s in starts]
    mates = [1, 2] * (n_reads // 2) + [1] * (n_reads % 2)
    report = recovery.compare(a, b, templates, mates, rng, q_reads=q_reads, substitutions_only=True)
    return {"reads": n_reads, "read_length": read_length, **vars(report)}


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
    p.add_argument("--skiver", type=Path, metavar="PREFIX", help="`skiver analyze -o` prefix on the same reads")
    p.add_argument("--kmer-spec", type=Path, help="spec from fit.kmer on the same run's skiver outputs")
    p.add_argument("--generated-reads", type=int, default=0, help="with --kmer-spec: phase 3 metrics on this many")
    p.add_argument("--read-length", type=int, default=150, help="with --generated-reads")
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
    if args.kmer_spec is not None:
        k = spec_io.load(args.kmer_spec)
        report["kmer"] = {
            name: {
                "components": components(k, other, t, min_exposure=args.min_exposure),
                "models": models(k, other, t),
            }
            for name, other, t in (("overlap", a, overlap), ("reference", b, reference))
        }
        if args.generated_reads:
            report["kmer"]["generated"] = {
                name: generated(k, other, args.reference, args.generated_reads, args.read_length)
                for name, other in (("overlap", a), ("reference", b))
            }
    if args.skiver is not None:
        analyze = skiver_analyze.read_analyze(args.skiver)
        report["skiver"] = {
            name: skiver_evidence(analyze, t, args.min_exposure)
            for name, t in (("overlap", overlap), ("reference", reference))
        }
    args.output.write_text(json.dumps(report, indent=2))
    e = report["evidence"]
    print(json.dumps({k: e[k] for k in ("rate_a", "rate_b", "excess", "coverage")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
