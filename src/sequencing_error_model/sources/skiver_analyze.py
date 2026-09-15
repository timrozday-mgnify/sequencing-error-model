"""Typed, header-validated parsers for unmodified `skiver analyze` (v0.3.1-v0.3.2) CSVs.

Default-mode input (plan §2.1). Column layouts were read from upstream skiver
v0.3.2 (`src/summary.rs`, `src/inference.rs`, `src/analyze.rs`); v0.3.1 writes the
same headers. Bootstrap CI columns vary between runs on the same input; everything
else is deterministic.

Conventions pinned from upstream source:

- **`summary_phred.csv`** (plan §5.6). For each observed value, bases are compared to
  the consensus value position by position from the first value base, and the scan
  stops at the first mismatch. Every compared base contributes its own reported Q:
  matches count as `num_correct`, the first mismatch as `num_error`, and later bases
  are not counted. Only value positions t in [1 + ignore_smallest_t, v - ignore_largest_t]
  (defaults 2 and 2) are summed, and qualities are value-side only. So the counts are
  hazard-model exposures ("all bases up to and including the first error"), and
  `per_base_error_rate` is 1 - exp(-λ_Q) from a Weibull fit with the global β, not
  num_error / (num_correct + num_error).
- **`summary_gc_content.csv`** uses the same scan and t range, keyed by read GC %.
- **`summary_read_position.csv`** uses the same scan over *all* value positions
  (no t trimming), indexed from the read start and from the read end.
- **Op names.** `kvmer.csv` columns use Rust `Debug` names (`AC`, `_A`, `A_`); the
  spectrum files use `Display` names (`A>C`, `->A`, `A>-`). Both are normalised to
  the `Display` form in `OPS`.
- **`hazard_rate.csv`** and `summary_error_spectrum_dependence_on_t.csv` index t as
  k + value position; `hazard_ratio` is `NA` when there are no candidates.
"""

import csv
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from sequencing_error_model.observations import ERROR_OPS, CountTable, Key, count

SUPPORTED_VERSIONS = ("0.3.1", "0.3.2")

OPS = ERROR_OPS  # upstream `ALL_OPERATIONS` order, as `Display` names
_KVMER_OP_COLUMNS = [op.replace(">", "").replace("-", "_") for op in OPS]

CI = tuple[float, float]


class SkiverFormatError(ValueError):
    """An analyze output is missing or does not match a supported skiver version."""


@dataclass(frozen=True)
class ErrorRate:
    per_base_error_rate: float
    per_base_error_rate_ci: CI
    mean_hazard_rate: float
    mean_hazard_rate_ci: CI
    lambda_: float
    lambda_ci: CI
    beta: float
    beta_ci: CI
    key_median_coverage: int
    key_coverage_ci: CI
    true_median_coverage: float
    true_coverage_ci: CI
    substitution_error_proportion: float
    insertion_error_proportion: float
    deletion_error_proportion: float


@dataclass(frozen=True)
class HazardRow:
    t: int
    num_candidates: int
    num_survival: int
    hazard_ratio: float | None
    ci: CI


@dataclass(frozen=True)
class SpectrumRow:
    op: str
    prev_base: str
    next_base: str
    total: int
    forward: int


@dataclass(frozen=True)
class SpectrumByTRow:
    op: str
    prev_base: str
    next_base: str
    total: int
    freq_at_t: dict[int, int]


@dataclass(frozen=True)
class RateBin:
    """A `summary_phred.csv` row (lo = hi = Q) or a `summary_gc_content.csv` row (GC % in [lo, hi))."""

    lo: int
    hi: int
    per_base_error_rate: float
    per_base_error_rate_ci: CI
    num_correct: int
    num_error: int


@dataclass(frozen=True)
class ReadPositionRow:
    index: int
    from_start: bool
    num_correct: int
    num_error: int


@dataclass(frozen=True)
class KvmerRow:
    key: str
    consensus_value: str
    passes_filter: bool
    homopolymer_length: int
    consensus_count: int
    neighbor_count: int
    total_count: int
    op_counts: dict[str, int]
    consensus_count_up_to_v: tuple[int, ...]


@dataclass(frozen=True)
class SkiverAnalyze:
    k: int
    v: int
    error_rate: ErrorRate
    hazard: list[HazardRow]
    survival: dict[int, float]
    spectrum: list[SpectrumRow]
    spectrum_by_t: list[SpectrumByTRow]
    phred: list[RateBin]
    gc_content: list[RateBin]
    read_position: list[ReadPositionRow]
    kvmer: list[KvmerRow]


def _read(path: Path, header: Sequence[str], *, prefix: bool = False) -> tuple[list[str], list[list[str]]]:
    """Return (header, rows) after checking the header equals `header`, or starts with it if `prefix`."""
    if not path.is_file():
        raise SkiverFormatError(f"{path}: missing (supported skiver versions: {', '.join(SUPPORTED_VERSIONS)})")
    with path.open(newline="") as fh:
        lines = [row for row in csv.reader(fh) if row]
    got, rows = (lines[0], lines[1:]) if lines else ([], [])
    if (got[: len(header)] if prefix else got) != list(header):
        raise SkiverFormatError(
            f"{path}: unexpected header (supported skiver versions: {', '.join(SUPPORTED_VERSIONS)})\n"
            f"  expected: {','.join(header)}{',...' if prefix else ''}\n  got:      {','.join(got)}"
        )
    for n, row in enumerate(rows, start=2):
        if len(row) != len(got):
            raise SkiverFormatError(f"{path}:{n}: expected {len(got)} fields, got {len(row)}")
    return got, rows


def _ci(s: str) -> CI:
    lo, hi = s.split("~")
    return float(lo), float(hi)


def _bool(s: str) -> bool:
    if s not in ("true", "false"):
        raise SkiverFormatError(f"expected true/false, got {s!r}")
    return s == "true"


def _op(name: str) -> str:
    if name not in OPS:
        raise SkiverFormatError(f"unknown operation {name!r}")
    return name


def read_error_rate(path: Path) -> ErrorRate:
    header = (
        "per_base_error_rate,per_base_error_rate_5-95th_percentile,mean_hazard_rate,mean_hazard_rate_5-95th_percentile,"
        "lambda,lambda_5-95th_percentile,beta,beta_5-95th_percentile,key_median_coverage,key_coverage_5-95th_percentile,"
        "true_median_coverage,true_coverage_5-95th_percentile,"
        "substitution_error_proportion,insertion_error_proportion,deletion_error_proportion"
    )
    _, rows = _read(path, header.split(","))
    if len(rows) != 1:
        raise SkiverFormatError(f"{path}: expected 1 data row, got {len(rows)}")
    r = rows[0]
    return ErrorRate(
        float(r[0]), _ci(r[1]), float(r[2]), _ci(r[3]), float(r[4]), _ci(r[5]), float(r[6]), _ci(r[7]),
        int(r[8]), _ci(r[9]), float(r[10]), _ci(r[11]), float(r[12]), float(r[13]), float(r[14]),
    )  # fmt: skip


def read_hazard_rate(path: Path) -> list[HazardRow]:
    _, rows = _read(path, ["t", "num_candidates", "num_survival", "hazard_ratio", "5th_percentile", "95th_percentile"])
    return [
        HazardRow(int(t), int(nc), int(ns), None if hr == "NA" else float(hr), (float(lo), float(hi)))
        for t, nc, ns, hr, lo, hi in rows
    ]


def read_survival_rate(path: Path) -> dict[int, float]:
    _, rows = _read(path, ["t", "survival_rate"])
    return {int(t): float(s) for t, s in rows}


def read_error_spectrum(path: Path) -> list[SpectrumRow]:
    _, rows = _read(path, ["operation", "prev_base", "next_base", "total", "forward"])
    return [SpectrumRow(_op(op), p, n, int(t), int(f)) for op, p, n, t, f in rows]


def read_error_spectrum_by_t(path: Path) -> list[SpectrumByTRow]:
    header, rows = _read(path, ["operation", "prev_base", "next_base", "total"], prefix=True)
    ts = [int(c.removeprefix("freq_at_t")) for c in header[4:]]
    return [SpectrumByTRow(_op(r[0]), r[1], r[2], int(r[3]), dict(zip(ts, map(int, r[4:]), strict=True))) for r in rows]


def read_phred(path: Path) -> list[RateBin]:
    header = "qscore,empirical_qscore,per_base_error_rate,per_base_error_rate_5-95th_percentile,num_correct,num_error"
    _, rows = _read(path, header.split(","))
    return [RateBin(int(q), int(q), float(rate), _ci(ci), int(nc), int(ne)) for q, _, rate, ci, nc, ne in rows]


def read_gc_content(path: Path) -> list[RateBin]:
    header = "gc_content_min,gc_content_max_exclusive,per_base_error_rate,per_base_error_rate_5-95th_percentile,num_correct,num_error"
    _, rows = _read(path, header.split(","))
    return [RateBin(int(lo), int(hi), float(rate), _ci(ci), int(nc), int(ne)) for lo, hi, rate, ci, nc, ne in rows]


def read_read_position(path: Path) -> list[ReadPositionRow]:
    _, rows = _read(path, ["index", "from_start", "num_correct", "num_error", "error_rate"])
    return [ReadPositionRow(int(i), _bool(fs), int(nc), int(ne)) for i, fs, nc, ne, _ in rows]


def read_kvmer(path: Path) -> list[KvmerRow]:
    fixed = [
        "key",
        "consensus_value",
        "passes_filter",
        "homopolymer_length",
        "consensus_count",
        "neighbor_count",
        "total_count",
    ]
    header, rows = _read(path, [*fixed, *_KVMER_OP_COLUMNS], prefix=True)
    ops_end = len(fixed) + len(OPS)
    tail = header[ops_end:]
    if not tail or tail != [f"consensus_count_up_to_v{i}" for i in range(1, len(tail) + 1)]:
        raise SkiverFormatError(f"{path}: unexpected consensus_count_up_to_v* columns: {','.join(tail)}")
    return [
        KvmerRow(
            r[0],
            r[1],
            _bool(r[2]),
            int(r[3]),
            int(r[4]),
            int(r[5]),
            int(r[6]),
            dict(zip(OPS, map(int, r[len(fixed) : ops_end]), strict=True)),
            tuple(map(int, r[ops_end:])),
        )  # fmt: skip
        for r in rows
    ]


def read_analyze(prefix: str | Path) -> SkiverAnalyze:
    """Read every CSV written by `skiver analyze -o <prefix>`."""
    p = str(prefix)
    for newer in ("survival_rate", "summary_gc_content"):  # added in v0.3.x
        if not Path(f"{p}.{newer}.csv").is_file():
            raise SkiverFormatError(
                f"{p}.{newer}.csv: missing; outputs look older than skiver 0.3 "
                f"(supported: {', '.join(SUPPORTED_VERSIONS)})"
            )
    kvmer = read_kvmer(Path(f"{p}.kvmer.csv"))
    if not kvmer:
        raise SkiverFormatError(f"{p}.kvmer.csv: no keys")
    return SkiverAnalyze(
        k=len(kvmer[0].key),
        v=len(kvmer[0].consensus_value),
        error_rate=read_error_rate(Path(f"{p}.summary_error_rate.csv")),
        hazard=read_hazard_rate(Path(f"{p}.hazard_rate.csv")),
        survival=read_survival_rate(Path(f"{p}.survival_rate.csv")),
        spectrum=read_error_spectrum(Path(f"{p}.summary_error_spectrum.csv")),
        spectrum_by_t=read_error_spectrum_by_t(Path(f"{p}.summary_error_spectrum_dependence_on_t.csv")),
        phred=read_phred(Path(f"{p}.summary_phred.csv")),
        gc_content=read_gc_content(Path(f"{p}.summary_gc_content.csv")),
        read_position=read_read_position(Path(f"{p}.summary_read_position.csv")),
        kvmer=kvmer,
    )


def tables(a: SkiverAnalyze) -> list[CountTable]:
    """The analyze counts as observation tables, with `t` rebased to the value position.

    `summary_error_rate.csv` and `survival_rate.csv` are fitted parameters, not counts;
    they stay on `SkiverAnalyze`. `kvmer.csv` keys failing skiver's outlier filter are dropped.
    """
    k, v = a.k, a.v
    scan = "value bases up to and including the first mismatch"
    trimmed = f"{scan}, t in [3, {v - 2}]"  # default ignore_smallest_t = ignore_largest_t = 2

    def table(
        name: str, fields: tuple[str, ...], unit: str, items: Iterable[tuple[Key, int]], **meta: Any
    ) -> CountTable:
        return CountTable(f"skiver_analyze:{name}", fields, unit, True, count(items), {"k": k, "v": v, **meta})

    def context(op: str, prev: str, nxt: str) -> str:
        return prev + ("-" if op.startswith("-") else op[0]) + nxt

    passing = [r for r in a.kvmer if r.passes_filter]
    return [
        table("phred", ("q", "op"), "base", (((b.lo, op), n) for b in a.phred for op, n in (("=", b.num_correct), ("!", b.num_error))), exposure=trimmed),
        table("gc_content", ("gc", "op"), "base", ((((b.lo, b.hi), op), n) for b in a.gc_content for op, n in (("=", b.num_correct), ("!", b.num_error))), exposure=trimmed),
        table("read_position_start", ("pos_start", "op"), "base", (((r.index, op), n) for r in a.read_position if r.from_start for op, n in (("=", r.num_correct), ("!", r.num_error))), exposure=scan),
        table("read_position_end", ("pos_end", "op"), "base", (((r.index, op), n) for r in a.read_position if not r.from_start for op, n in (("=", r.num_correct), ("!", r.num_error))), exposure=scan),
        table("hazard", ("t", "op"), "value", (((h.t - k, op), n) for h in a.hazard for op, n in (("=", h.num_survival), ("!", h.num_candidates - h.num_survival)))),
        table("spectrum", ("context", "strand", "op"), "error", (((context(r.op, r.prev_base, r.next_base), s, r.op), n) for r in a.spectrum for s, n in (("+", r.forward), ("-", r.total - r.forward))), flank=(1, 1)),
        table("spectrum_by_t", ("context", "t", "op"), "error", (((context(r.op, r.prev_base, r.next_base), t - k, r.op), n) for r in a.spectrum_by_t for t, n in r.freq_at_t.items()), flank=(1, 1)),
        # "=" counts values matching the consensus throughout; single-edit op counts have a latent t
        table("kvmer", ("locus", "op"), "value", (((r.key + r.consensus_value, op), n) for r in passing for op, n in (("=", r.consensus_count), *r.op_counts.items()))),
        table("kvmer_survival", ("locus", "t", "op"), "value", (((r.key + r.consensus_value, t, "="), n) for r in passing for t, n in enumerate(r.consensus_count_up_to_v, 1))),
    ]  # fmt: skip
