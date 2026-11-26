"""Exporter round trips with the real simulators (plan phase 9): export a spec, run the simulator on a random
genome, align its reads back (minibwa for short reads, minimap2 for long, soft clips realigned as `reference`
mode does) and re-estimate the error rate and the Q marginal, against the spec's own values from the
exporter's fidelity report.

The exit criterion: under the chosen `--q-policy` the preserved quantity is within `TOLERANCE` of the spec -
the error rate (`preserve-errors`), or the Q histogram (TV) and mean Q (`preserve-quality`); Badread and
PBSIM3 ERRHMM must preserve the error rate (Badread also the Q marginal). The violated quantity is reported
beside the fidelity report's prediction of it. `python -m sequencing_error_model.export.roundtrip` runs every
simulator found on PATH (`iss`, `art_illumina`, `badread`, `pbsim`) and exits non-zero on a failure.
"""

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from sequencing_error_model.export import Sample, art, badread, columns, error_rate, iss, pbsim3, tv
from sequencing_error_model.fit.quality import Array
from sequencing_error_model.recovery import ALIGNER_SCORES, example_spec, run_aligner, scale_error_rate
from sequencing_error_model.sources import bam
from sequencing_error_model.spec import Component, ErrorModelSpec

TOLERANCE = 0.10
SHORT, LONG = 100, 5000  # read lengths
BINARIES = {
    "iss": "iss",
    "art": "art_illumina",
    "badread": "badread",
    "pbsim3-qshmm": "pbsim",
    "pbsim3-errhmm": "pbsim",
}


def short_read_spec() -> ErrorModelSpec:
    """The recovery example (binned Illumina Q, miscalibrated) at a tenth of its error rate (~2% substitutions)."""
    return scale_error_rate(example_spec(), 0.1)


def long_read_spec() -> ErrorModelSpec:
    """A small ONT-like truth: Q 0..50 with sticky lag-1 transitions around Q18, ~6% errors (more at low Q, after
    G and in homopolymer runs). Reported Q understates the error rate, so the two policies differ."""
    alphabet, q = tuple(range(51)), np.arange(51.0)
    k = len(alphabet)
    lags = np.zeros((1, k + 1, k))
    lags[0, :k] = 3.0 * np.exp(-np.abs(q[:, None] - q[None, :]) / 4)
    errors = np.r_[0, np.ones(9)]
    window = np.zeros((3, k + 1, 10))
    window[1, :k] = np.outer(-0.12 * (q - 20), errors)
    context = np.zeros((3, 5, 10))
    context[1, 2, 1:5] = 0.8
    runs = np.zeros((5, 10))
    runs[2:, 5:] = 1.0
    return ErrorModelSpec(
        alphabet,
        {"mode": "reference", "sources": ["example"]},
        (Component("QualityMarkov(1)", {"bias": -(((q - 18) / 7) ** 2), "lags": lags}),),
        (
            Component("QualityWindow(1)", {"bias": np.r_[0, [-5.5] * 4, [-6.0] * 4, -4.5], "window": window}),
            Component("Context(1,1)", {"weights": context}),
            Component("Homopolymer", {"weights": runs}, meta={"flank": [2, 2]}),
        ),
    )


def _run(cmd: Sequence[str], cwd: Path) -> None:
    subprocess.run(list(cmd), cwd=cwd, check=True, capture_output=True, env={**os.environ, "TERM": "dumb"})


def _fastq(paths: Sequence[Path], out: Path) -> None:
    with out.open("w") as handle:
        for p in paths:
            with gzip.open(p, "rt") if p.suffix == ".gz" else p.open() as src:
                shutil.copyfileobj(src, handle)


def _simulate(name: str, report: dict[str, Any], genome: Path, work: Path, seed: int) -> Path:
    """Run the simulator with the report's arguments; returns one FASTQ of all its reads."""
    reads = work / "reads.fq"
    if name == "iss":
        _run(
            [
                "iss",
                "generate",
                "--genomes",
                str(genome),
                "--mode",
                "kde",
                "--model",
                report["arguments"][3],
                "-n",
                "4000",
                "--cpus",
                "1",
                "--seed",
                str(seed),
                "--output",
                "r",
            ],
            work,
        )
        _fastq([work / "r_R1.fastq", work / "r_R2.fastq"], reads)
    elif name == "art":
        _run(
            [
                "art_illumina",
                "-i",
                str(genome),
                "-p",
                "-m",
                "300",
                "-s",
                "20",
                "-c",
                "2000",
                "-na",
                "-nf",
                "0",
                "-rs",
                str(seed),
                "-o",
                "o",
                *report["arguments"],
            ],
            work,
        )
        _fastq([work / "o1.fq", work / "o2.fq"], reads)
    elif name == "badread":
        with reads.open("w") as out:
            subprocess.run(
                [
                    "badread",
                    "simulate",
                    "--reference",
                    str(genome),
                    "--quantity",
                    "10x",
                    "--length",
                    f"{LONG},500",
                    "--seed",
                    str(seed),
                    *report["arguments"],
                ],
                cwd=work,
                check=True,
                stdout=out,
                stderr=subprocess.DEVNULL,
            )
    else:
        args = report["arguments"] if name == "pbsim3-qshmm" else report["errhmm_arguments"]
        _run(
            [
                "pbsim",
                "--strategy",
                "wgs",
                "--genome",
                str(genome),
                "--depth",
                "10",
                "--length-sd",
                "2000",
                "--seed",
                str(seed),
                "--prefix",
                "p",
                *args,
            ],
            work,
        )
        _fastq([work / "p_0001.fq.gz"], reads)
    return reads


def _measured(reads: Path, genome: Path, short: bool, work: Path) -> tuple[Array, Sample]:
    """(Q of every simulated read base, for the Q marginal; the aligned reads, for error rates)."""
    qual = "".join(reads.read_text().splitlines()[3::4]).encode()
    q = np.frombuffer(qual, np.uint8).astype(np.int64) - 33
    aligner = "minibwa" if short else "minimap2"
    sam = work / "aligned.sam"
    run_aligner(aligner, genome, reads, sam)
    records = list(bam.records(sam, genome, unclip=ALIGNER_SCORES[aligner]))
    return q, Sample([t for t, _, _ in records], [r for _, r, _ in records], [m for _, _, m in records])


def roundtrip(name: str, q_policy: str, work: Path, seed: int = 0) -> dict[str, Any]:
    """Export, simulate, re-estimate; the report holds the fidelity report's numbers, the measured ones and the
    failures against `TOLERANCE`."""
    work.mkdir(parents=True, exist_ok=True)
    short = name in ("iss", "art")
    rng = np.random.default_rng(seed)
    genome = work / "genome.fa"
    genome.write_text(">g\n" + "".join(rng.choice(list("ACGT"), size=100_000 if short else 200_000)) + "\n")
    model = short_read_spec() if short else long_read_spec()
    exporter: Callable[..., dict[str, Any]] = {"iss": iss.export, "art": art.export, "badread": badread.export}.get(
        name, pbsim3.export
    )
    options = (
        {"read_length": SHORT, "reads": 20000, "insert_mean": 300.0, "insert_sd": 20.0}
        if short
        else {"read_length": LONG, "reads": 400}
    )
    report = exporter(model, work, rng, q_policy=q_policy, **options)
    q, aligned = _measured(_simulate(name, report, genome, work, seed), genome, short, work)
    kind = report["error_kind"]
    spec = {"error_rate": report["spec_error_rate"], "mean_q": report["spec_mean_q"]}
    measured = {"error_rate": error_rate(columns(aligned), kind), "aligned_reads": len(aligned.reads)}
    if name != "pbsim3-errhmm":
        spec_q = np.zeros(94)
        spec_q[[int(k) for k in report["spec_q_histogram"]]] = list(report["spec_q_histogram"].values())
        measured |= {"mean_q": float(q.mean()), "q_tv": tv(np.bincount(q, minlength=94) / len(q), spec_q)}
    preserves = report["preserves"] if name != "pbsim3-errhmm" else "errors"
    failures = []
    if preserves != "qualities" and abs(measured["error_rate"] / spec["error_rate"] - 1) > TOLERANCE:
        failures.append(f"error rate {measured['error_rate']:.4g} vs spec {spec['error_rate']:.4g}")
    if preserves != "errors":
        if measured["q_tv"] > TOLERANCE:
            failures.append(f"Q histogram TV {measured['q_tv']:.3f}")
        if abs(measured["mean_q"] / spec["mean_q"] - 1) > TOLERANCE:
            failures.append(f"mean Q {measured['mean_q']:.2f} vs spec {spec['mean_q']:.2f}")
    predicted = {k: report[k] for k in ("implied_error_rate", "exported_mean_q", "q_tv") if k in report}
    return {
        "simulator": name,
        "q_policy": q_policy,
        "preserves": preserves,
        "spec": spec,
        "predicted": predicted,
        "measured": measured,
        "failures": failures,
    }


RUNS = (
    [("iss", p) for p in ("preserve-quality", "preserve-errors")]
    + [("art", p) for p in ("preserve-quality", "preserve-errors")]
    + [("badread", "preserve-quality")]
    + [("pbsim3-qshmm", p) for p in ("preserve-quality", "preserve-errors")]
    + [("pbsim3-errhmm", "preserve-quality")]
)


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="python -m sequencing_error_model.export.roundtrip", description=__doc__)
    p.add_argument("--workdir", type=Path, required=True)
    p.add_argument("--output", type=Path, help="JSON report (default: stdout)")
    p.add_argument("--require", action="store_true", help="fail when a simulator is missing instead of skipping it")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--only", nargs="+", choices=sorted(BINARIES), help="run these simulators only")
    args = p.parse_args(argv)
    results = []
    for name, policy in RUNS:
        if args.only and name not in args.only:
            continue
        if shutil.which(BINARIES[name]) is None:
            results.append(
                {
                    "simulator": name,
                    "q_policy": policy,
                    "skipped": f"{BINARIES[name]} not on PATH",
                    "failures": [f"{BINARIES[name]} missing"] if args.require else [],
                }
            )
            continue
        results.append(roundtrip(name, policy, args.workdir / f"{name}-{policy}", args.seed))
        print(json.dumps(results[-1]), file=sys.stderr)
    text = json.dumps(results, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    else:
        sys.stdout.write(text)
    return 1 if any(r["failures"] for r in results) else 0


if __name__ == "__main__":
    sys.exit(main())
