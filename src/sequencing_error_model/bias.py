"""Phase 7 problem-size grid (plan §9 phase 7): run each evidence mode *unmodified* over reads that carry
biological variation, and report its bias against the clonal truth per head E component.

Test infrastructure only (§1 scope). `variation` builds the population and the per-base truth; this module
generates reads from it with a known error spec, runs the three modes on them, and compares each fitted head E
with the truth on held-out **clonal** templates - so every deviation is the variation's cost, not the fitter's:

- `pe-overlap`: mate overlaps of the same pairs, which see a minor allele in both mates and so stay blind to it;
- `reference`: single-end alignment to a reference the run does not own - the population consensus, the
  majority strain, or an external relative at `relative_ani` - then a head E refit from the alignments;
- `kmer`: a released `skiver analyze` with its own outlier filter, profiled from the same FASTQ.

`--site-mask` turns on phase 8's conservative mask (`sites.conservative`) inside the `reference` mode, so the
same grid measures what separation buys over the unmodified run.

The grid's clonal point (`--ani 1 --haplotypes 1 --minor-density 0`) is the phase 4-6 recovery run, which is
what makes this a bias table rather than a pile of numbers. `labels` sizes the problem independently of any fit:
the share of read bases that really differ from the consensus, beside the share the generator got wrong.

Deferred: a self-assembly reference (needs an assembler; `--reference majority` is the same coordinate mistake
without the dependency) and external haplotypes from a VCF.
"""

import argparse
import json
import subprocess
import sys
import tempfile
from collections import Counter
from collections.abc import Sequence
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np

from sequencing_error_model import compare as comparison
from sequencing_error_model import recovery, sites, variation
from sequencing_error_model import spec as spec_io
from sequencing_error_model.generate import Read, generate, insert_sizes
from sequencing_error_model.sources import bam, fastq_quality, pe_overlap, skiver_analyze
from sequencing_error_model.sources import kmer as kmer_source
from sequencing_error_model.spec import ErrorModelSpec

MODES = ("pe-overlap", "reference", "kmer")
REFERENCES = ("consensus", "majority", "relative")


def population(
    consensus: Sequence[tuple[str, str]],
    rng: np.random.Generator,
    *,
    ani: float = 0.99,
    n_haplotypes: int = 2,
    tree: str = "star",
    minor_fraction: float = 0.2,
    minor_density: float = 0.0,
    minor_frequency: tuple[float, float] = (0.01, 0.5),
    repeat: tuple[int, float] | None = None,
    ti_tv: float = 2.0,
    indel_fraction: float = 0.1,
) -> variation.Population:
    """One scenario's population. `ani=1`, `n_haplotypes=1`, `minor_density=0` is the clonal point."""
    base = [(name, seq.upper()) for name, seq in consensus]
    repeats: tuple[variation.Repeat, ...] = ()
    if repeat is not None:
        base, one = variation.repeat_copy(base, int(repeat[0]), repeat[1], rng, ti_tv=ti_tv)
        repeats = (one,)
    haps = variation.haplotypes(base, n_haplotypes, ani, rng, tree=tree, ti_tv=ti_tv, indel_fraction=indel_fraction)
    minor = variation.minor_alleles(base, minor_density, rng, frequency=minor_frequency, ti_tv=ti_tv)
    abundance = variation.abundances(list(haps), [minor_fraction if len(haps) > 1 else 0.0])
    return variation.Population(tuple(base), haps, abundance, minor, repeats)


def labels(fragments: Sequence[variation.Fragment], reads: Sequence[Read]) -> dict[str, float]:
    """The share of read bases that match, really differ from the consensus, or were mis-read, from the CIGARs
    and each mate's consensus counterpart - the size of the problem, before any mode is run."""
    counts = Counter[str]()
    for i, frag in enumerate(fragments):
        for mate in (0, 1):
            counts.update(variation.truth(frag.templates[mate], frag.consensus[mate], reads[2 * i + mate]))
    total = max(counts.total(), 1)
    return {"bases": float(counts.total()), **{k: counts[k] / total for k in ("match", "variant", "error", "adapter")}}


def _fastq(path: Path, fragments: Sequence[variation.Fragment], reads: Sequence[Read], mate_suffix: bool) -> None:
    with path.open("w") as out:
        for i, frag in enumerate(fragments):
            for mate in (1, 2):
                r = reads[2 * i + mate - 1]
                name = f"{frag.name}:m{mate}" if mate_suffix else f"{frag.name}/{mate}"
                out.write(f"@{name}\n{r.sequence}\n+\n{r.quality}\n")


def _pair_flags(sam: Path) -> None:
    """Set the paired / first / second flag bits from the `:m1`, `:m2` name suffix `_fastq` wrote.

    Both aligners map single-end, so the mate a read came from survives only in its name, and head E's `Mate`
    term is fitted from the flags. The suffix is `:mN` rather than `/N` because bwa-style aligners strip a
    trailing `/1` or `/2` from the query name.
    """
    lines = []
    for line in sam.read_text().splitlines():
        if line.startswith("@"):
            lines.append(line)
            continue
        f = line.split("\t")
        if f[0].endswith((":m1", ":m2")):
            f[1] = str(int(f[1]) | 0x1 | (0x80 if f[0].endswith(":m2") else 0x40))
            f[0] = f[0][:-3]
        lines.append("\t".join(f))
    sam.write_text("\n".join(lines) + "\n")


def _reference(
    pop: variation.Population, which: str, relative_ani: float, rng: np.random.Generator
) -> list[tuple[str, str]]:
    """The FASTA the `reference` mode aligns to, none of which is the reads' own sequence."""
    if which == "consensus":
        return list(pop.consensus)
    if which == "majority":
        major = max(pop.abundance[pop.samples[0]], key=lambda h: pop.abundance[pop.samples[0]][h])
        return [(c.reference, c.sequence) for c in pop.haplotypes[major]]
    if which == "relative":
        rel = variation.haplotypes(list(pop.consensus), 1, relative_ani, rng, prefix="rel")
        return [(c.reference, c.sequence) for c in next(iter(rel.values()))]
    raise ValueError(f"unknown reference {which!r}; use one of {REFERENCES}")


def _error_failures(report: recovery.Report) -> list[str]:
    return [f for f in report.failures() if not f.startswith("per-position Q")]


def _summary(truth: ErrorModelSpec, fitted: ErrorModelSpec, report: recovery.Report, table: Any) -> dict[str, Any]:
    """One mode's bias per head E component: the marginal rate, op composition, context shape, Q and position."""
    indels = not np.isneginf(fitted.error_head[0].params["bias"][5:]).all()
    parts = comparison.components(truth, fitted, table, indels=indels)
    with np.errstate(divide="ignore", invalid="ignore"):
        by_q = np.array(report.curves["rate_by_q_fit"]) / np.array(report.curves["rate_by_q_true"])
    # Position comes from `components`, not from the report's curve: the curve has one point per exact read
    # position, and its last few carry almost no exposure, so its worst ratio is sampling noise.
    by_pos = np.array(parts.get("pos_start", {}).get("ratio", [float("nan")]))
    return {
        "rate_ratio": report.scalars["rate_ratio"],
        "rate_true": report.scalars["rate_true"],
        "rate_fit": report.scalars["rate_fit"],
        "op_tv": report.scalars["op_tv"],
        "op_shares_tv": parts["op"]["tv"],
        "context_log_odds_slope": parts["context"].get("log_odds_slope"),
        "context_r": parts["context"].get("r"),
        "rate_by_q_ratio": [round(float(x), 3) for x in by_q],
        "position_ratio_max": float(np.nanmax(np.abs(by_pos - 1))),
        "failures": _error_failures(report),
        "components": parts,
    }


def scenario(
    truth: ErrorModelSpec,
    workdir: Path,
    *,
    genome_length: int = 50_000,
    ani: float = 0.99,
    n_haplotypes: int = 2,
    tree: str = "star",
    minor_fraction: float = 0.2,
    minor_density: float = 0.0,
    minor_frequency: tuple[float, float] = (0.01, 0.5),
    repeat: tuple[int, float] | None = None,
    coverage: float = 30.0,
    read_length: int = 150,
    insert_mean: float = 220.0,
    insert_sd: float = 30.0,
    modes: Sequence[str] = MODES,
    reference: str = "consensus",
    relative_ani: float = 0.95,
    aligner: str = "minibwa",
    preset: str = "sr",
    min_mapq: int = 20,
    unclip: bool = True,
    site_mask: bool = False,
    skiver: str | None = None,
    k: int = 11,
    v: int = kmer_source.GOOD_V,
    c: int = 8,
    use_all: bool = False,
    skiver_args: Sequence[str] = (),
    q_reads: int = 20_000,
    seed: int = 0,
) -> dict[str, Any]:
    """One grid point: build the population, generate reads from `truth` over it, run `modes`, report the bias.

    Coverage is per genome copy, so `coverage * genome_length / (2 * read_length)` pairs are drawn and shared out
    by abundance. Every mode is compared against `truth` on held-out clonal templates sliced from the consensus,
    so the numbers are each mode's bias under variation and nothing else.
    """
    rng = np.random.default_rng(seed)
    workdir.mkdir(parents=True, exist_ok=True)
    genome = "".join(rng.choice(list("ACGT"), size=genome_length))
    pop = population(
        [("g", genome)],
        rng,
        ani=ani,
        n_haplotypes=n_haplotypes,
        tree=tree,
        minor_fraction=minor_fraction,
        minor_density=minor_density,
        minor_frequency=minor_frequency,
        repeat=repeat,
    )
    sample = pop.samples[0]
    n_pairs = max(int(coverage * genome_length / (2 * read_length)), 2)
    frags = variation.fragments(pop, sample, insert_sizes(insert_mean, insert_sd), n_pairs, read_length, rng)
    templates = [t for f in frags for t in f.templates]
    mates = [1, 2] * len(frags)
    reads = generate(truth, templates, mates, rng)

    consensus = pop.sequences["g"]
    n_held = max(n_pairs, 2)
    starts = rng.integers(0, len(consensus) - read_length, n_held)
    held = [consensus[s : s + read_length] for s in starts]
    held_mates = [1, 2] * (n_held // 2) + [1] * (n_held % 2)

    doc: dict[str, Any] = {
        "ani": ani,
        "achieved_ani": variation.ani(pop, *list(pop.haplotypes)[:2]) if n_haplotypes > 1 else 1.0,
        "haplotypes": n_haplotypes,
        "tree": tree,
        "minor_fraction": minor_fraction if n_haplotypes > 1 else 0.0,
        "minor_density": minor_density,
        "repeat": list(repeat) if repeat else None,
        "coverage": coverage,
        "read_length": read_length,
        "insert_mean": insert_mean,
        "genome_length": len(consensus),
        "pairs": len(frags),
        "sites": len(pop.sites()),
        "labels": labels(frags, reads),
        "modes": {},
    }

    for name in modes:
        # A mode that cannot run here (no skiver binary, no aligner, an unusable output) is recorded and the
        # other modes still land: the table is the deliverable, and a hole in one column is legible.
        try:
            fitted, extra = _fit(
                name,
                truth,
                pop,
                frags,
                reads,
                templates,
                mates,
                workdir / name,
                rng,
                reference=reference,
                relative_ani=relative_ani,
                aligner=aligner,
                preset=preset,
                min_mapq=min_mapq,
                unclip=unclip,
                site_mask=site_mask,
                skiver=skiver,
                k=k,
                v=v,
                c=c,
                use_all=use_all,
                skiver_args=skiver_args,
                read_length=read_length,
            )
        except (ValueError, OSError, subprocess.CalledProcessError) as e:
            doc["modes"][name] = {"skipped": f"{type(e).__name__}: {e}"}
            continue
        if fitted is None:
            doc["modes"][name] = extra
            continue
        substitutions_only = name == "pe-overlap"
        report = recovery.compare(
            truth, fitted, held, held_mates, rng, q_reads=q_reads, substitutions_only=substitutions_only
        )
        table = recovery._tuples(truth, held, generate(truth, held, held_mates, rng), held_mates, fitted)
        doc["modes"][name] = {**_summary(truth, fitted, report, table), **extra}
    return doc


def _fit(
    name: str,
    truth: ErrorModelSpec,
    pop: variation.Population,
    frags: Sequence[variation.Fragment],
    reads: Sequence[Read],
    templates: Sequence[str],
    mates: Sequence[int],
    workdir: Path,
    rng: np.random.Generator,
    **kw: Any,
) -> tuple[ErrorModelSpec | None, dict[str, Any]]:
    """One mode's fitted spec, plus what only that mode reports. `(None, {...})` when the mode cannot run."""
    workdir.mkdir(parents=True, exist_ok=True)
    if name == "pe-overlap":
        fitted = pe_overlap.mode(list(templates), list(reads), list(mates), truth)
        return fitted, {"stats": fitted.provenance["stats"]}
    if name == "reference":
        fasta, fastq, sam = workdir / "ref.fa", workdir / "reads.fastq", workdir / "aligned.sam"
        contigs = _reference(pop, kw["reference"], kw["relative_ani"], rng)
        with fasta.open("w") as out:
            for contig, seq in contigs:
                out.write(f">{contig}\n{seq}\n")
        _fastq(fastq, frags, reads, mate_suffix=True)
        recovery.run_aligner(kw["aligner"], fasta, fastq, sam, kw["preset"])
        _pair_flags(sam)
        scores = recovery.ALIGNER_SCORES[kw["aligner"]] if kw["unclip"] else None
        alignments = list(bam._aligned(sam, fasta, kw["min_mapq"], scores))
        masked, site_report = None, None
        if kw.get("site_mask"):
            # One extra pass over the same alignments: the phase 8 conservative mask (§6.6), which the mode
            # otherwise runs without, counting every strain allele as an error.
            seqs = dict(contigs)
            counts = bam.count_alleles(alignments, {c: len(s) for c, s in contigs})
            masked, site_report = sites.conservative(counts, sites.linked(alignments), seqs.__getitem__)
        aligned = list(bam.apply_masks(alignments, masked))
        if not aligned:
            return None, {"skipped": "no read aligned"}
        a_templates, a_reads, a_mates = (list(x) for x in zip(*aligned, strict=True))
        table = recovery._tuples(truth, a_templates, a_reads, a_mates)
        if not table.counts:
            return None, {"skipped": "the site mask left no row"}
        fitted = recovery._refit_error(truth, table, aligned)
        return fitted, {
            "reference": kw["reference"],
            "unclip": kw["unclip"],
            "mapped_fraction": len(aligned) / len(reads),
            "site_mask": site_report,
        }
    if name == "kmer":
        if not kw.get("skiver"):
            return None, {"skipped": "no skiver binary: pass --skiver"}
        fastq, prefix = workdir / "reads.fastq", workdir / "analyze"
        _fastq(fastq, frags, reads, mate_suffix=False)
        extra = ["--use-all", *kw["skiver_args"]] if kw["use_all"] else list(kw["skiver_args"])
        recovery.run_skiver(kw["skiver"], fastq, prefix, kw["k"], kw["v"], kw["c"], extra)
        a = skiver_analyze.read_analyze(prefix)
        flank, m = bam.window(recovery._ERROR_TOKENS, recovery._QUALITY_TOKENS)
        fitted = kmer_source.fit(
            a,
            fastq_quality.profile_fastq(fastq, order=m, flank=flank),
            fastq_quality.quality_table(fastq, m=m, flank=flank),
            recovery._ERROR_TOKENS,
            recovery._QUALITY_TOKENS,
            {"mode": "kmer", "skiver_build": "default", "sources": ["variation-grid"]},
            use_all=kw["use_all"],
            flank=flank,
        )
        if tuple(fitted.quality_alphabet) != tuple(truth.quality_alphabet):
            return None, {"skipped": f"the reads carry qualities {list(fitted.quality_alphabet)}"}
        return fitted, {
            "outlier_filter": kmer_source.filter_stats(a),
            "skiver_error_rate": a.error_rate.per_base_error_rate,
        }
    raise ValueError(f"unknown mode {name!r}; use one of {MODES}")


def grid(
    truth: ErrorModelSpec,
    workdir: Path,
    *,
    ani: Sequence[float] = (1.0, 0.999, 0.99),
    minor_fraction: Sequence[float] = (0.2,),
    coverage: Sequence[float] = (30.0,),
    minor_density: Sequence[float] = (0.0,),
    repeat: Sequence[tuple[int, float] | None] = (None,),
    **kw: Any,
) -> list[dict[str, Any]]:
    """Every scenario in the product, in order; a failing point is recorded, not raised, so the table lands."""
    out = []
    for i, (a, m, cov, d, rep) in enumerate(product(ani, minor_fraction, coverage, minor_density, repeat)):
        clonal = a >= 1.0 and d == 0.0
        point: dict[str, Any] = dict(
            ani=a, minor_fraction=m, coverage=cov, minor_density=d, repeat=list(rep) if rep else None, clonal=clonal
        )
        rest = {k: v for k, v in kw.items() if k != "n_haplotypes"}
        try:
            point.update(
                scenario(
                    truth,
                    workdir / f"point{i}",
                    ani=a,
                    # The clonal point is one haplotype: at ani=1 the tree mutates nothing, but a second
                    # haplotype would still split the coverage, which is not what phases 4-6 ran.
                    n_haplotypes=1 if clonal else kw.get("n_haplotypes", 2),
                    minor_fraction=m,
                    coverage=cov,
                    minor_density=d,
                    repeat=rep,
                    **rest,
                )
            )
        except (ValueError, subprocess.CalledProcessError, FileNotFoundError) as e:
            point["error"] = f"{type(e).__name__}: {e}"
        point["clonal"] = clonal
        out.append(point)
    return out


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m sequencing_error_model.bias",
        description="Phase 7 problem size: each evidence mode's bias under biological variation.",
    )
    p.add_argument("--spec", type=Path, help="truth spec directory (default: the built-in example spec)")
    p.add_argument("--error-rate-scale", type=float, default=0.2, help="scale the truth's error rate")
    p.add_argument("--output", type=Path, help="write the full report as JSON")
    p.add_argument("--modes", nargs="+", choices=MODES, default=list(MODES))
    p.add_argument("--ani", type=float, nargs="+", default=[1.0, 0.999, 0.99])
    p.add_argument("--haplotypes", type=int, default=2)
    p.add_argument("--tree", choices=("star", "random"), default="star")
    p.add_argument("--minor-fraction", type=float, nargs="+", default=[0.2])
    p.add_argument("--coverage", type=float, nargs="+", default=[30.0])
    p.add_argument("--minor-density", type=float, nargs="+", default=[0.0])
    p.add_argument("--repeat", type=float, nargs=2, metavar=("LENGTH", "IDENTITY"), help="add a divergent repeat")
    p.add_argument("--genome-length", type=int, default=50_000)
    p.add_argument("--read-length", type=int, default=150)
    p.add_argument("--insert-mean", type=float, default=220.0)
    p.add_argument("--insert-sd", type=float, default=30.0)
    p.add_argument("--reference", choices=REFERENCES, default="consensus")
    p.add_argument("--relative-ani", type=float, default=0.95)
    p.add_argument("--aligner", choices=("minibwa", "minimap2"), default="minibwa")
    p.add_argument("--preset", default="sr", help="aligner preset (minibwa -x, minimap2 -x)")
    p.add_argument("--clip", action="store_true", help="leave soft clips as they are (default: realign end to end)")
    p.add_argument(
        "--site-mask", action="store_true", help="mask sites the phase 8 conservative mask drops (`reference` mode)"
    )
    p.add_argument("--skiver", help="path to a released skiver binary, for the `kmer` mode")
    p.add_argument("--use-all", action="store_true", help="run skiver with --use-all")
    p.add_argument("-k", type=int, default=11)
    p.add_argument("-v", type=int, default=kmer_source.GOOD_V)
    p.add_argument("-c", type=int, default=8, help="FracMinHash denominator")
    p.add_argument("--skiver-arg", action="append", default=[], metavar="ARG")
    p.add_argument("--q-reads", type=int, default=20_000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--workdir", type=Path, help="keep the intermediate files here instead of a temporary directory")
    args = p.parse_args(argv)
    truth = recovery.scale_error_rate(
        spec_io.load(args.spec) if args.spec else recovery.example_spec(), args.error_rate_scale
    )
    shared: dict[str, Any] = dict(
        modes=args.modes,
        n_haplotypes=args.haplotypes,
        tree=args.tree,
        genome_length=args.genome_length,
        read_length=args.read_length,
        insert_mean=args.insert_mean,
        insert_sd=args.insert_sd,
        reference=args.reference,
        relative_ani=args.relative_ani,
        aligner=args.aligner,
        preset=args.preset,
        unclip=not args.clip,
        site_mask=args.site_mask,
        skiver=args.skiver,
        use_all=args.use_all,
        k=args.k,
        v=args.v,
        c=args.c,
        skiver_args=args.skiver_arg,
        q_reads=args.q_reads,
        seed=args.seed,
    )
    with tempfile.TemporaryDirectory() as tmp:
        points = grid(
            truth,
            args.workdir or Path(tmp),
            ani=args.ani,
            minor_fraction=args.minor_fraction,
            coverage=args.coverage,
            minor_density=args.minor_density,
            repeat=[(int(args.repeat[0]), args.repeat[1])] if args.repeat else [None],
            **shared,
        )
    doc = {"error_rate_scale": args.error_rate_scale, "points": points}
    if args.output:
        args.output.write_text(json.dumps(doc, indent=2) + "\n")
    print(json.dumps(table(points), indent=2))
    return 0


def _masked_fraction(report: Sequence[dict[str, Any]] | None) -> float | None:
    """The share of reference sites the conservative mask dropped, or None when it did not run."""
    if not report:
        return None
    dropped = sum(r["dropped"] for r in report)
    return round(float(dropped) / max(1, sum(int(r["length"]) for r in report)), 4)


def table(points: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """The bias table: one row per grid point and mode, the scalars only."""
    rows = []
    for point in points:
        for name, m in point.get("modes", {}).items():
            rows.append(
                {
                    "ani": point["ani"],
                    "minor_fraction": point["minor_fraction"],
                    "coverage": point["coverage"],
                    "minor_density": point["minor_density"],
                    "repeat": point["repeat"],
                    "variant_fraction": point.get("labels", {}).get("variant"),
                    "mode": name,
                    "masked_site_fraction": _masked_fraction(m.get("site_mask")),
                    **{
                        k: m.get(k)
                        for k in (
                            "skipped",
                            "rate_ratio",
                            "op_tv",
                            "op_shares_tv",
                            "context_log_odds_slope",
                            "position_ratio_max",
                            "failures",
                        )
                    },
                }
            )
        if "error" in point:
            rows.append({"ani": point["ani"], "coverage": point["coverage"], "error": point["error"]})
    return rows


if __name__ == "__main__":
    sys.exit(main())
