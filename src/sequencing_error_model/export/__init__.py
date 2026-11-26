"""Exporters (plan §7, phase 9): a spec projected onto other read simulators' model files.

Every exporter works from one `Sample`: reads drawn from the spec by the native generator (the reference
implementation of §5.3), with their true alignments. Each simulator's tables are the marginals of that sample
its format can hold, counted the way the simulator's own profiler would count them on real alignments, so the
export is exactly as faithful as the format allows and no closed form per component is needed.

Shared pieces:

- `sample`: random uniform templates (every k-mer equally exposed), one read each, in batches.
- `columns`: one row per alignment column (M, I or D) with read position, Q, true and read base.
- `remap` (`--q-policy preserve-errors`): for simulators that draw an error with probability 10^(-Q/10)
  (ART, InSilicoSeq, PBSIM3 QSHMM), each reported Q is replaced per (mate, read position, Q) cell by the
  empirical Q of the spec's errors in that cell, so the simulator's Q-derived errors match head E and the Q
  distribution is distorted instead. `preserve-quality` (the default, §12 decision 9) emits the spec's Q and
  the simulator's error rate deviates by the spec's calibration gap.
- `fidelity`: the machine-readable report every export writes (`fidelity.json`): which quantity the policy
  preserves, by how much the other is violated, what the format drops, and the simulator arguments to use.
- Base profile fallback: `profile_tracks` draws Q tracks from an ART-format quality profile (per position,
  independent) instead of head Q, snapped to the spec's alphabet; the fallback is recorded in the report.

Q is never evidence of error here either: the remap reads error rates off the spec's own sampled errors.
"""

import argparse
import json
import sys
from collections.abc import Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np

from sequencing_error_model import spec as spec_io
from sequencing_error_model.fit.quality import Array
from sequencing_error_model.generate import _CIGAR, Read, generate
from sequencing_error_model.spec import ErrorModelSpec

Q_POLICIES = ("preserve-quality", "preserve-errors")
_M, _I, _D = (ord(c) for c in "MID")
_GAP = ord("-")


@dataclass(frozen=True)
class Sample:
    templates: list[str]
    reads: list[Read]
    mates: list[int]


def read_art_profile(path: Path) -> tuple[Array, Array]:
    """The all-bases ('.') block of an ART quality profile: (Q values [K], per-position probabilities [L, K])."""
    rows: dict[int, list[list[int]]] = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if parts and parts[0] == ".":
            rows.setdefault(int(parts[1]), []).append([int(x) for x in parts[2:]])
    qs = sorted({q for r in rows.values() for q in r[0]})
    probs = np.zeros((len(rows), len(qs)))
    for pos, (values, cumulative) in sorted(rows.items()):
        counts = np.diff(np.r_[0, cumulative])
        probs[pos, np.searchsorted(qs, values)] = counts / counts.sum()
    return np.array(qs), probs


def profile_tracks(
    profiles: Sequence[tuple[Array, Array]],
    alphabet: Sequence[int],
    lengths: Sequence[int],
    mates: Sequence[int],
    rng: np.random.Generator,
) -> list[Array]:
    """Q tracks drawn per position from base profiles (one per mate; mate 2 reuses mate 1's if only one),
    snapped to the nearest `alphabet` value. Positions past a profile's end reuse its last position."""
    alpha = np.asarray(alphabet)
    out = []
    for n, m in zip(lengths, mates, strict=True):
        qs, probs = profiles[min(m, len(profiles)) - 1]
        p = probs[np.minimum(np.arange(n), len(probs) - 1)]
        q = qs[(rng.random((n, 1)) > p.cumsum(axis=1)).sum(axis=1).clip(max=len(qs) - 1)]
        out.append(alpha[np.abs(alpha[None, :] - q[:, None]).argmin(axis=1)])
    return out


def sample(
    model: ErrorModelSpec,
    n_reads: int,
    read_length: int,
    rng: np.random.Generator,
    paired: bool = True,
    base_profile: Sequence[Path] | None = None,
    batch: int = 2000,
) -> Sample:
    """`n_reads` reads of `read_length` template bases from uniform random templates; mates alternate 1, 2
    when `paired`. Q comes from head Q, or from `base_profile` (ART profiles for R1 [and R2]). Generated in
    batches of at most `batch` reads and 2 Mbp, to bound memory."""
    profiles = [read_art_profile(p) for p in base_profile] if base_profile else None
    batch = max(1, min(batch, 2_000_000 // read_length))
    out = Sample([], [], [])
    for first in range(0, n_reads, batch):
        n = min(batch, n_reads - first)
        codes = rng.integers(0, 4, size=(n, read_length))
        templates = ["".join(r) for r in np.array(list("ACGT"))[codes]]
        mates = [1 + (i % 2 if paired else 0) for i in range(first, first + n)]
        tracks = (
            None
            if profiles is None
            else profile_tracks(profiles, model.quality_alphabet, [read_length] * n, mates, rng)
        )
        out.templates.extend(templates)
        out.mates.extend(mates)
        out.reads.extend(generate(model, templates, mates, rng, tracks=tracks))
    return out


def columns(s: Sample) -> dict[str, Array]:
    """One row per alignment column of every read: `read`, `mate`, `op` (M, I, D as bytes), `pos` (0-based read
    index; a deletion takes the next read base's, the last one at the read end), `q` (Q at `pos`), `ref` (true
    base, "-" for an insertion) and `alt` (read base, "-" for a deletion)."""
    parts: dict[str, list[Array]] = {k: [] for k in ("read", "mate", "op", "pos", "q", "ref", "alt")}
    for i, (template, read, mate) in enumerate(zip(s.templates, s.reads, s.mates, strict=True)):
        if not read.sequence:
            continue
        found = _CIGAR.findall(read.cigar)
        ops = np.repeat(np.frombuffer("".join(k for _, k in found).encode(), np.uint8), [int(n) for n, _ in found])
        ri, ti = np.cumsum(ops != _D) - 1, np.cumsum(ops != _I) - 1
        seq, tpl = (np.frombuffer(x.encode(), np.uint8) for x in (read.sequence, template))
        qual = np.frombuffer(read.quality.encode(), np.uint8).astype(np.int64) - 33
        pos = np.where(ops == _D, np.minimum(ri + 1, len(seq) - 1), ri).clip(min=0)
        for key, value in (
            ("read", np.full(len(ops), i)),
            ("mate", np.full(len(ops), mate)),
            ("op", ops),
            ("pos", pos),
            ("q", qual[pos]),
            ("ref", np.where(ops == _I, _GAP, tpl[ti.clip(min=0)])),
            ("alt", np.where(ops == _D, _GAP, seq[ri.clip(min=0)])),
        ):
            parts[key].append(value)
    return {k: np.concatenate(v) if v else np.zeros(0, np.int64) for k, v in parts.items()}


def _errors(cols: dict[str, Array], kind: str) -> tuple[Array, Array]:
    """(exposure, error) indicators per column: `substitution` counts read bases from template bases (M) and
    their mismatches; `any` counts emitted read bases (M, I) and every edit column (a deletion included)."""
    if kind == "substitution":
        return cols["op"] == _M, (cols["op"] == _M) & (cols["ref"] != cols["alt"])
    return cols["op"] != _D, (cols["op"] != _M) | (cols["ref"] != cols["alt"])


def error_rate(cols: dict[str, Array], kind: str) -> float:
    n, e = _errors(cols, kind)
    return float(e.sum() / max(n.sum(), 1))


def implied_rate(cols: dict[str, Array], kind: str) -> float:
    """The error rate a simulator drawing errors with probability 10^(-Q/10) per exposed base produces."""
    n, _ = _errors(cols, kind)
    return float(np.mean(10.0 ** (-cols["q"][n] / 10))) if n.any() else 0.0


def remap(s: Sample, kind: str, by_position: bool, q_max: int) -> Sample:
    """`preserve-errors`: every read's Q replaced by the empirical Q of `kind` errors in its (mate, position, Q)
    cell (positions pooled unless `by_position`), clipped to [1, `q_max`]. A cell's rate is shrunk toward its
    (mate, Q) rate with a prior worth one expected error, so sparse high-Q cells don't round to noise."""
    cols = columns(s)
    n, e = _errors(cols, kind)
    top = int(cols["pos"].max(initial=0)) + 1 if by_position else 1
    shape = (3, top, 94)
    cell = (cols["mate"], cols["pos"] if by_position else np.zeros(len(n), np.int64), cols["q"])
    exposure, errors = np.zeros(shape), np.zeros(shape)
    np.add.at(exposure, cell, n)
    np.add.at(errors, cell, e)
    with np.errstate(divide="ignore", invalid="ignore"):
        pooled = errors.sum(axis=1, keepdims=True) / exposure.sum(axis=1, keepdims=True)
        prior = np.where(pooled > 0, 1 / pooled, 0)
        rate = np.where(pooled > 0, (errors + prior * pooled) / (exposure + prior), 0)
        table = np.clip(np.round(-10 * np.log10(rate)), 1, q_max)
    table = np.nan_to_num(table, nan=q_max, posinf=q_max).astype(np.int64)
    reads = []
    for read, mate in zip(s.reads, s.mates, strict=True):
        q = np.frombuffer(read.quality.encode(), np.uint8).astype(np.int64) - 33
        pos = np.arange(len(q)).clip(max=top - 1) if by_position else np.zeros(len(q), np.int64)
        reads.append(replace(read, quality=(table[mate, pos, q] + 33).astype(np.uint8).tobytes().decode()))
    return replace(s, reads=reads)


def q_histogram(cols: dict[str, Array], positions: int = 0) -> Array:
    """Q histogram of emitted read bases over 0..93, overall ([94]) or for the first `positions` read positions
    ([positions, 94])."""
    keep = cols["op"] != _D
    if not positions:
        return np.asarray(np.bincount(cols["q"][keep], minlength=94) / max(keep.sum(), 1))
    keep &= cols["pos"] < positions
    h = np.zeros((positions, 94))
    np.add.at(h, (cols["pos"][keep], cols["q"][keep]), 1)
    return np.asarray(h / h.sum(axis=1, keepdims=True).clip(min=1))


def tv(a: Array, b: Array) -> float:
    """Total variation distance (the max over rows for 2-D histograms)."""
    return float(np.max(0.5 * np.abs(a - b).sum(axis=-1)))


def fidelity(
    simulator: str,
    model: ErrorModelSpec,
    spec_sample: Sample,
    emitted: Sample,
    q_policy: str | None,
    kind: str,
    dropped: Sequence[str],
    arguments: Sequence[str],
    base_profile: Sequence[Path] | None = None,
    positions: int = 0,
    **extra: Any,
) -> dict[str, Any]:
    """The fidelity report: the spec's error rate of `kind` and Q marginals (from `spec_sample`) against what
    the simulator will produce from the exported qualities (`emitted`); `q_policy` None for simulators whose
    errors don't depend on Q."""
    before, after = columns(spec_sample), columns(emitted)
    spec_rate = error_rate(before, kind)
    out: dict[str, Any] = {
        "simulator": simulator,
        "q_policy": q_policy,
        "preserves": {None: "errors and qualities", "preserve-quality": "qualities", "preserve-errors": "errors"}[
            q_policy
        ],
        "error_kind": kind,
        "spec_error_rate": spec_rate,
        "spec_mean_q": float(before["q"][before["op"] != _D].mean()),
        "spec_q_histogram": {int(q): float(p) for q, p in enumerate(q_histogram(before)) if p},
        "exported_mean_q": float(after["q"][after["op"] != _D].mean()),
        "q_tv": tv(q_histogram(before), q_histogram(after)),
        "dropped": list(dropped),
        "unidentified": [c.token for c in (*model.quality_head, *model.error_head) if not c.identified],
        "base_profile": [str(p) for p in base_profile] if base_profile else None,
        "arguments": list(arguments),
        "sample_reads": len(spec_sample.reads),
        **extra,
    }
    if positions:
        out["q_position_tv_max"] = tv(q_histogram(before, positions), q_histogram(after, positions))
    if q_policy is not None:
        implied = implied_rate(after, kind)
        out |= {"implied_error_rate": implied, "error_rate_ratio": implied / spec_rate if spec_rate else None}
    return out


def write_report(report: dict[str, Any], out: Path) -> None:
    (out / "fidelity.json").write_text(json.dumps(report, indent=2) + "\n")


def main(argv: Sequence[str] | None = None) -> int:
    from sequencing_error_model.export import art, badread, iss, pbsim3

    exporters = {"art": art.export, "iss": iss.export, "badread": badread.export, "pbsim3": pbsim3.export}
    p = argparse.ArgumentParser(prog="sem-export", description="Export a spec to a read simulator's model files.")
    p.add_argument("simulator", choices=sorted(exporters))
    p.add_argument("--model", type=Path, required=True, metavar="SPEC_DIR")
    p.add_argument("--output", type=Path, required=True, metavar="DIR", help="model files and fidelity.json")
    p.add_argument("--q-policy", choices=Q_POLICIES, default="preserve-quality", help="ART, InSilicoSeq, PBSIM3 QSHMM")
    p.add_argument("--read-length", type=int, help="template bases per sampled read (default: per simulator)")
    p.add_argument("--reads", type=int, help="reads sampled from the spec (default: per simulator)")
    p.add_argument(
        "--base-profile",
        type=Path,
        nargs="+",
        metavar="ART_PROFILE",
        help="R1 [R2] ART quality profiles to draw Q from instead of head Q",
    )
    p.add_argument("--insert-mean", type=float, help="InSilicoSeq: fragment size (default: the spec's insert_size)")
    p.add_argument("--insert-sd", type=float, default=0.0)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args(argv)
    options = {
        k: v for k, v in vars(args).items() if k not in ("simulator", "model", "output", "seed") and v is not None
    }
    args.output.mkdir(parents=True, exist_ok=True)
    report = exporters[args.simulator](
        spec_io.load(args.model), args.output, np.random.default_rng(args.seed), **options
    )
    write_report(report, args.output)
    print(json.dumps({k: report[k] for k in ("simulator", "preserves", "arguments")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
