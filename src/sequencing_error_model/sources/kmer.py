"""`kmer` default mode: one spec from unmodified `skiver analyze` outputs plus the raw FASTQ (plan §6.1).

The head E pieces are in `fit.kmer`; this composes them and fits head Q beside them:

1. the context head from `kvmer.csv` by EM over the latent edit position (`fit.kmer.fit`);
2. the centre-Q term by marginal matching against `summary_phred.csv` under the FASTQ's (observed context,
   centre Q) exposure (`fit.kmer.fit_centre_q`), which is also what sets the level;
3. `Position(n)`, `Strand` and `GC(n)` from their own skiver marginals, log-additively (`fit.kmer.position`,
   `strand`, `gc`), each shape-only so the level stays with step 2;
4. skiver's global Weibull and op proportions into the spec's `marginals` (`fit.kmer.hazard`), reports only;
5. head Q from every FASTQ base, conditioned on *observed* bases (`fastq_quality.quality_table`).

Which of step 3's components are fitted follows `error_tokens`, and only the centre Q is identifiable, so head E
must start with `QualityWindow(0)` (§5.6). `filter_stats` reports what skiver's outlier filter removed, and
`use_all` fits without it, as skiver's own `--use-all` does.
"""

import argparse
import json
import sys
import warnings
from collections import Counter
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any

from sequencing_error_model.fit import kmer, quality
from sequencing_error_model.observations import CountTable
from sequencing_error_model.sources import fastq_quality, skiver_analyze
from sequencing_error_model.sources.fastq_quality import QualityProfile
from sequencing_error_model.sources.skiver_analyze import RateBin, SkiverAnalyze
from sequencing_error_model.spec import Component, ErrorModelSpec

# Lesson §3.1, measured on the fork's amplicon recovery: v=13 recovers the source rate, v=6 is ~20% low, and
# v=1 observes no errors at all ("reject inputs with tiny v").
MIN_V = 4
GOOD_V = 13
# Phase 6b: keys must be (nearly) single-locus. On a 5 Mb isolate, 90% of key observations sit at multi-locus
# keys at k = 11, 5.4% at 15, 2.6% at 17 and 2.1% at 21 (`multiplicity.py`). skiver's own default is 21.
MIN_K = 15
GOOD_K = 17
DEFAULT_K = 21

_CONTEXT = frozenset(("QualityWindow", "Context", "Homopolymer"))
_MARGINAL = frozenset(("Position", "Strand", "GC"))
# Identified only by per-mate skiver runs, and only in the raked fit: `kvmer.csv` has no mate (phase 6b item 3).
_PER_MATE = frozenset(("Mate",))


def filter_stats(a: SkiverAnalyze) -> dict[str, Any]:
    """What skiver's outlier filter removes from `kvmer.csv`: keys, distinct consensus values, and error mass.

    On clonal reads there is no variation to remove, so everything here is cost (plan phase 6, exit criteria).
    """
    kept = [r for r in a.kvmer if r.passes_filter]

    def mass(rows: Sequence[Any]) -> float:
        return float(sum(sum(r.op_counts.values()) for r in rows))

    total, dropped = mass(a.kvmer), mass(a.kvmer) - mass(kept)
    return {
        "keys": len(a.kvmer),
        "keys_dropped": len(a.kvmer) - len(kept),
        "consensus_values": len({r.consensus_value for r in a.kvmer}),
        "consensus_values_dropped": len({r.consensus_value for r in a.kvmer}) - len({r.consensus_value for r in kept}),
        "error_mass": total,
        "error_mass_dropped": dropped,
        "error_mass_dropped_fraction": dropped / total if total else 0.0,
        "observations": float(sum(r.total_count for r in a.kvmer)),
        "observations_dropped": float(sum(r.total_count for r in a.kvmer if not r.passes_filter)),
    }


def _check_v(v: int) -> None:
    """Refuse a value length too short to observe errors, and warn below the length that recovers them.

    An error is only visible where the observed value differs from the consensus, so a short value sees a
    shrinking share of them (lesson §3.1). This is the guard that phase 6's synthetic loop reproduces.
    """
    if v < MIN_V:
        raise ValueError(
            f"skiver ran with v={v}: too short to observe errors (v=1 observes none, v=6 is ~20% low). "
            f"Rerun `skiver analyze` with v >= {GOOD_V}."
        )
    if v < GOOD_V:
        warnings.warn(
            f"skiver ran with v={v} < {GOOD_V}: the error rate will be biased low (~20% at v=6, lesson §3.1)",
            stacklevel=2,
        )


def _check_k(k: int, min_k: int = MIN_K) -> None:
    """Refuse a key length at which a real genome's keys repeat, and warn below the length that keeps them unique.

    A key at several loci has several true values: skiver reads the minority ones as errors or filters the key
    (phase 6b: at k = 11, 54% of the error mass passing skiver's filter on SRR24523812 came from repeated keys).
    Synthetic random genomes are single-locus at any k, so the harnesses pass `min_k` to measure small k.
    """
    if k < min_k:
        raise ValueError(
            f"skiver ran with k={k}: keys repeat across a real genome at this length (90% of observations at "
            f"k=11 on a 5 Mb isolate). Rerun `skiver analyze` with -k {DEFAULT_K} (skiver's default), or pass "
            f"min_k for a genome known to be single-locus at k={k}."
        )
    if k < GOOD_K:
        warnings.warn(
            f"skiver ran with k={k} < {GOOD_K}: multi-locus keys will inflate the error rate (5% of a 5 Mb "
            f"isolate's observations at k=15); complex metagenomes may need more than k={DEFAULT_K}",
            stacklevel=2,
        )


def _split(error_tokens: Sequence[str], per_mate: bool = False) -> tuple[list[str], list[str], list[Component]]:
    """(the kvmer fit's tokens, the raked fit's tokens, the marginal components), checked against default mode.

    The two token lists differ only by the per-mate terms: `kvmer.csv` carries no mate, so a `Mate` term is
    fitted on the raked counts alone, and only when per-mate skiver runs were given.
    """
    if list(error_tokens[:1]) != ["QualityWindow(0)"]:
        raise ValueError(f"default mode identifies the centre Q only, so head E starts with QualityWindow(0): "
                         f"{list(error_tokens)}")  # fmt: skip
    context: list[str] = []
    raked: list[str] = []
    marginal: list[Component] = []
    for t in error_tokens:
        c = Component(t)
        if c.name in _CONTEXT:
            context.append(t)
            raked.append(t)
        elif c.name in _MARGINAL:
            marginal.append(c)
        elif c.name in _PER_MATE:
            if not per_mate:
                raise ValueError(
                    f"{t} needs one `skiver analyze` run per mate (plan phase 6b item 3): pass `mates`, or "
                    f"--mate-analyze on the command line"
                )
            raked.append(t)
        else:
            raise ValueError(
                f"{t} is not identifiable in `kmer` default mode (§6.2); "
                f"usable: {sorted(_CONTEXT | _MARGINAL | _PER_MATE)}"
            )
    return context, raked, marginal


def fit(
    a: SkiverAnalyze,
    profile: QualityProfile,
    q_table: CountTable,
    error_tokens: Sequence[str],
    quality_tokens: Sequence[str],
    provenance: dict[str, Any],
    *,
    use_all: bool = False,
    flank: tuple[int, int] | None = None,
    level: float = 1.0,
    min_k: int = MIN_K,
    mates: Sequence[tuple[SkiverAnalyze, QualityProfile]] = (),
) -> ErrorModelSpec:
    """Both heads for `kmer` default mode. `flank` widens the context window the kvmer fit carries, and `level`
    scales `summary_phred.csv`'s rates, which is how a key filter other than skiver's own reaches the level
    (`sites.keys`). `min_k` lowers the key-length refusal (`_check_k`) for genomes known to be single-locus.

    `mates` is one (`skiver analyze`, FASTQ profile) pair per mate, in mate order, from skiver runs on that mate
    alone. It is what identifies a `Mate` term in `error_tokens`: each mate's centre Q is raked against its own
    reported rates and exposure (`fit.kmer.rake_mates`), while context, position, strand and GC stay with the
    pooled run `a` (plan phase 6b item 3)."""
    _check_v(a.v)
    _check_k(a.k, min_k)
    context_tokens, raked_tokens, marginal = _split(error_tokens, per_mate=bool(mates))
    exposure = next(t for t in fastq_quality.tables(profile) if t.source.endswith(":context"))
    kvmer = next(t for t in skiver_analyze.tables(a, use_all) if t.source.endswith(":kvmer"))
    alphabet = profile.alphabet
    if not alphabet:
        raise ValueError("the FASTQ profile has no qualities, so there is no alphabet to fit over")
    context_head = kmer.fit(kvmer, context_tokens, flank=flank)

    def scaled(bins: Sequence[RateBin]) -> list[RateBin]:
        return [replace(b, per_base_error_rate=b.per_base_error_rate * level) for b in bins]

    phred = scaled(a.phred)
    per_mate = []
    for mate_a, mate_profile in mates:
        _check_v(mate_a.v)
        _check_k(mate_a.k, min_k)
        mate_exposure = next(t for t in fastq_quality.tables(mate_profile) if t.source.endswith(":context"))
        per_mate.append((scaled(mate_a.phred), mate_exposure))
    head = list(kmer.fit_centre_q(context_head, phred, exposure, raked_tokens, alphabet, mates=per_mate))
    for c in marginal:
        if c.name == "Position":
            head.append(kmer.position(a.read_position, profile.lengths, c.args[0]))
        elif c.name == "Strand":
            head.append(kmer.strand(a.spectrum))
        else:
            head.append(kmer.gc(a.gc_content, c.args[0]))
    quality_head = quality.fit(q_table, quality_tokens, alphabet)
    return ErrorModelSpec(
        tuple(alphabet),
        {
            "k": a.k,
            "v": a.v,
            "use_all": use_all,
            "outlier_filter": filter_stats(a),
            "per_mate_skiver": [{"k": m.k, "v": m.v, "reads": p.n_reads} for m, p in mates],
            **provenance,
        },
        quality_head,
        tuple(head),
        marginals=kmer.hazard(a.error_rate),
    )


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m sequencing_error_model.sources.kmer",
        description="Fit a spec from unmodified `skiver analyze` outputs and the reads they came from.",
    )
    p.add_argument("analyze", type=Path, metavar="PREFIX", help="`skiver analyze -o` prefix")
    p.add_argument("fastq", type=Path, nargs="+", help="the same reads as FASTQ[.gz]")
    p.add_argument("--output", type=Path, required=True, metavar="SPEC_DIR")
    p.add_argument(
        "--error-tokens", nargs="+", default=["QualityWindow(0)", "Context(1,1)", "Position(4)", "Strand", "GC(3)"]
    )
    p.add_argument("--quality-tokens", nargs="+", default=["QualityMarkov(1)", "Position(48)", "Context(1,1)"])
    p.add_argument("--use-all", action="store_true", help="keep kvmer keys skiver's outlier filter rejected")
    p.add_argument(
        "--mate-analyze",
        type=Path,
        nargs="+",
        metavar="PREFIX",
        help="one `skiver analyze -o` prefix per mate, in mate order, from skiver run on that mate alone; "
        "what a Mate error token needs (plan phase 6b item 3). The FASTQ arguments must then be one per mate.",
    )
    p.add_argument("--flank", type=int, nargs=2, metavar=("LEFT", "RIGHT"), default=[2, 2], help="context window")
    p.add_argument("--max-reads", type=int)
    p.add_argument("--min-k", type=int, default=MIN_K, help="refuse skiver runs with a shorter key (phase 6b)")
    args = p.parse_args(argv)
    flank = (args.flank[0], args.flank[1])
    m = Component(args.quality_tokens[0]).args[0]
    a = skiver_analyze.read_analyze(args.analyze)
    profile = fastq_quality.profile_fastq(*args.fastq, order=m, flank=flank, max_reads=args.max_reads)
    mates = []
    if args.mate_analyze:
        if len(args.mate_analyze) != len(args.fastq):
            p.error(f"--mate-analyze takes one prefix per mate FASTQ: {len(args.mate_analyze)} and {len(args.fastq)}")
        mates = [
            (
                skiver_analyze.read_analyze(prefix),
                fastq_quality.profile_fastq(fastq, order=m, flank=flank, max_reads=args.max_reads),
            )
            for prefix, fastq in zip(args.mate_analyze, args.fastq, strict=True)
        ]
    # Head Q's `mate` field is only right when the files are known to be one per mate; pooled, every read
    # would be labelled mate 1, so a Mate quality token is refused there rather than fitted on a wrong label.
    if mates:
        q_tables = [
            fastq_quality.quality_table(f, m=m, flank=flank, mate=i, max_reads=args.max_reads)
            for i, f in enumerate(args.fastq, start=1)
        ]
        q_table = replace(q_tables[0], counts=+sum((t.counts for t in q_tables), Counter()))
    else:
        if any(Component(t).name == "Mate" for t in args.quality_tokens):
            p.error("a Mate quality token needs one FASTQ per mate, with --mate-analyze")
        q_table = fastq_quality.quality_table(*args.fastq, m=m, flank=flank, max_reads=args.max_reads)
    provenance: dict[str, Any] = {
        "mode": "kmer",
        "skiver_build": "default",
        "sources": [str(args.analyze), *map(str, args.fastq)],
        "mate_sources": [str(x) for x in (args.mate_analyze or ())],
        "identified_ops": ["substitution", "insertion", "deletion"],  # kvmer.csv reports all three, 1 bp only
    }
    spec = fit(
        a,
        profile,
        q_table,
        args.error_tokens,
        args.quality_tokens,
        provenance,
        use_all=args.use_all,
        flank=flank,
        min_k=args.min_k,
        mates=mates,
    )
    spec.save(args.output)
    print(json.dumps({"reads": profile.n_reads, "quality_alphabet": spec.quality_alphabet, **filter_stats(a)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
