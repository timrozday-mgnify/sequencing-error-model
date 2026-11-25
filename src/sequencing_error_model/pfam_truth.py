"""`pfam` mode against reference truth (docs/pfam_mode_plan.md, P2).

`score` takes the reads a `pfam` fit used, their alignment to genomes known to be in the sample, and fitted specs:

- every counted stop is classified by the genome under it: a sequencing **error** (the genome has a sense codon
  there), a **genome** stop or a **variant** (the biological floor: pseudogenes, recoding, misplaced reads), or not
  decided (`unaligned`, `indel`);
- every base of every exposed codon is labelled by the genome, whatever its codon, and the specs' predicted errors
  are compared with the observed ones per substitution, Q, mate and read position. A fit only saw stop-adjacent
  bases, so this also tests its extrapolation to the rest of the codon.

`verdict` applies the go / no-go rule of P2 to the `score` reports of several libraries.

The genomes are used here only, never by the fit. Q is a covariate of the comparison, not a label.
"""

import argparse
import bisect
import json
import math
import sys
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, cast

import numpy as np
import pysam
from scipy import special

from sequencing_error_model import spec as spec_io
from sequencing_error_model.fit import error
from sequencing_error_model.fit.quality import Array
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.sources import bam, pfam
from sequencing_error_model.sources.pe_overlap import _COMP
from sequencing_error_model.sources.pfam import BASES, STOPS, UNIDENTIFIED_BY_STOPS, Codon, Profiles
from sequencing_error_model.spec import ErrorModelSpec

CLASSES = ("error", "genome", "variant", "indel", "unaligned")
FLOOR = ("genome", "variant")
# P2's go rule. A cell is tested when its sd is at most half its tolerance; see `verdict`.
TOTAL, CELL, MAX_FLOOR, RUNS_NEEDED = 0.15, 0.25, 2e-4, 5
_OPS = tuple(f"{a}>{b}" for a in BASES for b in BASES if a != b)
TEN = tuple(op for op in _OPS if op not in UNIDENTIFIED_BY_STOPS)
_Q_EDGES, _POSITIONS, _BATCH = (20, 30, 36), 4, 200_000


@dataclass(frozen=True)
class Placement:
    """Where a read lies on the genomes: the contig, and the site under each base in read order (-1: none)."""

    contig: str
    reverse: bool
    sites: Array


@dataclass
class Cell:
    observed: float = 0.0  # errors, by the genome
    exposure: float = 0.0  # bases
    stops: float = 0.0  # error stops in the cell: the counting error of a fit to stops
    predicted: dict[str, float] = field(default_factory=lambda: defaultdict(float))


def exposed_reads(
    pairs: Iterable[pfam.ReadPair],
    hits: Mapping[int, Sequence[pfam.Hit]],
    profiles: Profiles,
    t1: float,
    t2: float,
) -> dict[int, tuple[str, bytes, list[Codon]]]:
    """(bases, Qs, exposed codons) of the reads with any, by read id."""
    out = {}
    for number, (s1, q1, s2, q2) in enumerate(pairs):
        for mate, (s, q) in enumerate(((s1, q1), (s2, q2))):
            rid = 2 * number + mate
            if (placed := hits.get(rid)) and (codons := pfam.exposed(s, placed, profiles, t1, t2)):
                out[rid] = (s.upper(), bytes(q), codons)
    return out


def placements(sam: Path, wanted: Mapping[int, tuple[str, bytes, list[Codon]]], min_mapq: int) -> dict[int, Placement]:
    """Primary alignments of the wanted reads; the query names are pair numbers (`pfam reads`)."""
    out: dict[int, Placement] = {}
    with pysam.AlignmentFile(str(sam), check_sq=False) as f:
        for r in f.fetch(until_eof=True):
            name = r.query_name
            if name is None or not name.isdigit():
                raise ValueError(f"{sam}: read {name!r} is not a pair number: align `pfam reads` output")
            rid = 2 * int(name) + int(r.is_read2)
            contig = r.reference_name
            if r.is_secondary or r.is_supplementary or r.is_unmapped or contig is None or rid not in wanted:
                continue
            n = len(wanted[rid][0])
            if r.mapping_quality < min_mapq or r.infer_read_length() != n:
                continue
            at = np.array(r.get_aligned_pairs(matches_only=True), np.int64).reshape(-1, 2)
            sites = np.full(n, -1, np.int64)
            sites[n - 1 - at[:, 0] if r.is_reverse else at[:, 0]] = at[:, 1]
            out[rid] = Placement(contig, bool(r.is_reverse), sites)
    return out


def variant_sites(
    reads: Mapping[int, tuple[str, bytes, list[Codon]]],
    placed: Mapping[int, Placement],
    genome: Mapping[str, str],
    max_alt_freq: float,
) -> dict[str, Array]:
    """Sites where these reads carry another base at `max_alt_freq` or more (`bam.site_masks`), per contig."""
    counts: dict[str, Array] = {}
    code = np.full(256, -1, np.int64)
    for i, b in enumerate(BASES):
        code[ord(b)] = i
    for rid, p in placed.items():
        base = code[np.frombuffer(reads[rid][0].encode(), np.uint8)]
        if p.reverse:
            base = np.where(base >= 0, 3 - base, -1)  # ACGT complement
        on = (p.sites >= 0) & (base >= 0)
        table = counts.get(p.contig)
        if table is None:
            table = counts[p.contig] = np.zeros((len(genome[p.contig]) + 1, 6), np.int32)
        np.add.at(table, (p.sites[on], base[on]), 1)
    return bam.site_masks(counts, lambda c: genome[c], max_alt_freq)[0]


def true_bases(codon: Codon, p: Placement, genome: Mapping[str, str]) -> tuple[str, str] | None:
    """The genome's bases under the codon, (read orientation, coding orientation); None off it or across a gap."""
    sites = [int(p.sites[x]) for x in codon.positions]
    if min(sites) < 0 or abs(sites[1] - sites[0]) != 1 or sites[2] - sites[1] != sites[1] - sites[0]:
        return None
    read = "".join(genome[p.contig][s] for s in sites)
    if p.reverse:
        read = "".join(_COMP.get(b, "N") for b in read)
    return read, "".join(_COMP.get(b, "N") for b in read) if codon.reverse else read


def stop_class(codon: Codon, p: Placement | None, genome: Mapping[str, str], masks: Mapping[str, Array]) -> str:
    """Which of `CLASSES` a stop codon is."""
    if p is None:
        return "unaligned"
    if (true := true_bases(codon, p, genome)) is None:
        return "indel"
    if true[1] in STOPS:
        return "genome"
    differ = [int(p.sites[x]) for x, a, b in zip(codon.positions, true[1], codon.bases, strict=True) if a != b]
    return "variant" if any(masks[p.contig][s] for s in differ) else "error"


class _Cells:
    """Observed and predicted substitutions per cell, flushed in batches so rows need not fit in memory."""

    def __init__(self, specs: Mapping[str, ErrorModelSpec], m: int, flank: tuple[int, int]) -> None:
        self.specs, self.m, self.flank = specs, m, flank
        self.cells: dict[str, Cell] = defaultdict(Cell)
        self.rows: Counter[Key] = Counter()
        names = pfam.fields(m)
        self.at = {f: names.index(f) for f in ("q", "context", "pos_start", "pos_end", "mate", "op")}

    def groups(self, key: Key) -> tuple[str, str, str]:
        q, start, end = (cast("int", key[self.at[f]]) for f in ("q", "pos_start", "pos_end"))
        group = bisect.bisect(_Q_EDGES, q)
        return (
            f"q:{(0, *_Q_EDGES)[group]}+",
            f"mate:{key[self.at['mate']]}",
            f"position:{(start - 1) * _POSITIONS // (start + end)}",
        )

    def add(self, key: Key, is_stop: bool) -> None:
        self.rows[key] += 1
        op = str(key[self.at["op"]])
        if is_stop and op in TEN:
            for c in (f"op:{op}", *self.groups(key)):
                self.cells[c].stops += 1
        if len(self.rows) >= _BATCH:
            self.flush()

    def flush(self) -> None:
        if not self.rows:
            return
        keys = list(self.rows)
        table = CountTable("truth", pfam.fields(self.m), "base", True, self.rows, {"flank": self.flank})
        p = {
            name: error.probabilities(s.error_head, s.quality_alphabet, table)[:, 1:5] for name, s in self.specs.items()
        }
        ten = np.array([[f"{a}>{b}" in TEN for b in BASES] for a in BASES])
        for r, key in enumerate(keys):
            n, op = self.rows[key], str(key[self.at["op"]])
            centre = str(key[self.at["context"]])[self.flank[0]]
            for j, b in enumerate(BASES):
                if b == centre:
                    continue
                cell = self.cells[f"op:{centre}>{b}"]
                cell.exposure += n
                cell.observed += n if op == f"{centre}>{b}" else 0
                for name, x in p.items():
                    cell.predicted[name] += n * float(x[r, j])
            # the other cells sum the ten identified substitutions
            use = ten[BASES.index(centre)]
            for c in self.groups(key):
                cell = self.cells[c]
                cell.exposure += n
                cell.observed += n if op in TEN else 0
                for name, x in p.items():
                    cell.predicted[name] += n * float(x[r, use].sum())
        self.rows = Counter()


def score(
    pairs: Sequence[pfam.ReadPair] | Callable[[], Iterable[pfam.ReadPair]],
    hits: Mapping[int, Sequence[pfam.Hit]],
    profiles: Profiles,
    sam: Path,
    genome: Mapping[str, str],
    specs: Mapping[str, ErrorModelSpec],
    *,
    t1: float = 35.0,
    t2: float = 6.0,
    flank: tuple[int, int] = (1, 1),
    min_mapq: int = 20,
    max_alt_freq: float = 0.2,
) -> dict[str, Any]:
    """The truth report of one library; `specs` by name, `both` being the one under test."""
    low = min(t1, pfam.PLATEAU[0])
    again = pairs if callable(pairs) else lambda: pairs
    reads = exposed_reads(again(), hits, profiles, low, t2)
    placed = placements(sam, reads, min_mapq)
    masks = variant_sites(reads, placed, genome, max_alt_freq)
    m = max(s.error_head[0].args[0] for s in specs.values())
    cells = _Cells(specs, m, flank)
    # the stops a fit would drop because they recur (`pfam.drop_replicated`), by read and position
    at = [(rid, c) for rid, (_, _, codons) in reads.items() for c in codons if c.bases in STOPS and c.anchor >= t1]
    loci = [(c.column[0], c.bases, pfam.flanks(reads[rid][0], c)) for rid, c in at]
    again_found = pfam.replicates(loci, again(), hits)
    recurring = {(rid, c.positions) for (rid, c), n in zip(at, again_found, strict=True) if pfam.recurs(*n)}
    by_anchor: dict[int, Counter[str]] = defaultdict(Counter)
    replicated: Counter[str] = Counter()
    family: dict[int, Counter[str]] = defaultdict(Counter)
    stats: Counter[str] = Counter()
    for rid, (seq, quals, codons) in reads.items():
        p, read = placed.get(rid), pfam._Read(seq, list(quals), rid % 2 + 1)
        for c in codons:
            is_stop = c.bases in STOPS
            kind = stop_class(c, p, genome, masks) if is_stop else "sense"
            chance = 0 if is_stop else sum(map(len, pfam._OPPORTUNITY[c.bases]))
            by_anchor[int(c.anchor)][kind] += 1
            by_anchor[int(c.anchor)]["chance"] += chance
            if c.anchor < t1:
                continue
            stats["codons"] += 1
            repeated = (rid, c.positions) in recurring
            if repeated:
                replicated[kind] += 1
            else:
                family[c.column[0]][kind] += 1
                family[c.column[0]]["chance"] += chance
            # the truth keeps the errors that happened to recur: dropping them is the fit's loss
            true = None if p is None else true_bases(c, p, genome)
            if p is None or true is None:
                stats["codons_unaligned" if p is None else "codons_indel"] += 1
                continue
            for x, b in zip(c.positions, true[0], strict=True):
                if b in BASES and seq[x] in BASES and not masks[p.contig][int(p.sites[x])]:
                    cells.add(
                        pfam._row(read, x, b, seq[x], BASES, m, flank), kind == "error" and not repeated and b != seq[x]
                    )
    cells.flush()
    down = {
        name: sorted({*s.provenance.get("downweighted_families", {}), *s.provenance.get("downweighted_overlaps", {})})
        for name, s in specs.items()
    }
    index = {a: i for i, a in enumerate(profiles.names)}
    dropped = {index[a] for a in down.get("both", []) if a in index}
    kept = Counter[str]()
    for f, n in family.items():
        if f not in dropped:
            kept.update(n)
    everything = sum(family.values(), Counter[str]())
    worst = sorted(family, key=lambda f: -sum(family[f][k] for k in FLOOR))[:20]
    stats["reads_exposed"], stats["reads_aligned"] = len(reads), len(placed)
    return {
        "t1": t1,
        "stats": dict(stats),
        "stops": {k: everything[k] for k in CLASSES},
        "replicated_stops": {k: replicated[k] for k in CLASSES},  # dropped, as the fit does: not in what follows
        "floor_before": _floor(everything + replicated),
        "floor": _floor(everything),
        "floor_kept_families": _floor(kept),
        "plateau": {str(t): _rates(by_anchor, t) for t in pfam.PLATEAU if t >= low},
        "floor_families": {
            profiles.names[f]: {k: family[f][k] for k in (*CLASSES, "sense")}
            for f in worst
            if any(family[f][k] for k in FLOOR)
        },
        "cells": {
            name: {
                "observed": c.observed,
                "exposure": c.exposure,
                "stops": c.stops,
                "predicted": {k: float(v) for k, v in c.predicted.items()},
            }
            for name, c in sorted(cells.cells.items())
        },
        "fits": {
            name: {
                k: s.provenance.get(k)
                for k in (
                    "stats",
                    "floor",
                    "replication",
                    "lam",
                    "plateau_warning",
                    "unidentified_ops",
                    "downweighted_families",
                )
            }
            for name, s in specs.items()
        },
    }


def _floor(n: Mapping[str, float]) -> dict[str, float]:
    """The floor as a substitution rate per base (3 x stops per opportunity, as `Evidence.crude_rate`)."""
    stops, chance = sum(n[k] for k in FLOOR), max(n["chance"], 1)
    upper = float(special.gammaincinv(stops + 1, 0.95))
    return {"stops": stops, "rate": 3 * stops / chance, "upper95": 3 * upper / chance}


def _rates(by_anchor: Mapping[int, Mapping[str, float]], t1: float) -> dict[str, float]:
    total = Counter[str]()
    for bits, n in by_anchor.items():
        if bits >= t1:
            total.update(n)
    return {k: 3 * total[k] / max(total["chance"], 1) for k in CLASSES}


def _ratio(cells: Sequence[Mapping[str, Any]], spec: str) -> tuple[float, float]:
    """(predicted over observed, its sd) of pooled cells: counting error of the stops fitted and of the truth."""
    predicted, observed, stops = (
        sum(x) for x in zip(*((c["predicted"].get(spec, 0.0), c["observed"], c["stops"]) for c in cells), strict=True)
    )
    if not observed or not stops:
        return math.nan, math.inf
    return predicted / observed, predicted / observed * math.sqrt(1 / stops + 1 / observed)


def checks(report: Mapping[str, Any], spec: str = "both") -> list[dict[str, Any]]:
    """P2's checks on one library: `within` its tolerance, `outside` it by 2 sd or more, `unclear`, or `untested`."""
    cells = report["cells"]
    out = []

    def judge(name: str, use: Sequence[Mapping[str, Any]], tolerance: float) -> None:
        ratio, sd = _ratio(use, spec)
        off = abs(ratio - 1)
        result = (
            "untested"
            if not sd <= tolerance / 2
            else "within"
            if off <= tolerance
            else "outside"
            if off > tolerance + 2 * sd
            else "unclear"
        )
        out.append({"check": name, "ratio": ratio, "sd": sd, "tolerance": tolerance, "result": result})

    judge("total", [cells[f"op:{op}"] for op in TEN if f"op:{op}" in cells], TOTAL)
    for name, c in cells.items():
        if name.startswith("op:") and name[3:] not in TEN:
            continue
        judge(name, [c], CELL)
    f = report["floor_kept_families"]
    result = "within" if f["upper95"] <= MAX_FLOOR else "outside" if f["rate"] > MAX_FLOOR else "unclear"
    out.append({"check": "floor", "rate": f["rate"], "upper95": f["upper95"], "tolerance": MAX_FLOOR, "result": result})
    return out


def passed(found: Sequence[Mapping[str, Any]]) -> bool | None:
    """Total and floor `within`, and no cell `outside`; None when the total could not be tested."""
    by = {c["check"]: c["result"] for c in found}
    if by["total"] == "untested":
        return None
    return by["total"] == "within" and by["floor"] == "within" and "outside" not in by.values()


def recovered(report: Mapping[str, Any], dose: float, spec: str, base: str = "both") -> float:
    """Standard additions: the dosed fit's rate less the undosed one's, on the same bases, over the dose (ten ops)."""
    each = []
    for op in TEN:
        c = report["cells"].get(f"op:{op}")
        if c:
            each.append((c["predicted"][spec] - c["predicted"][base]) / c["exposure"] / (dose / 3))
    return float(np.mean(each))


def verdict(reports: Mapping[str, Mapping[str, Any]], doses: Sequence[float] = (0.002, 0.01)) -> dict[str, Any]:
    """Go when `RUNS_NEEDED` libraries pass. A library without power counts as neither pass nor fail."""
    runs: dict[str, Any] = {}
    for name, r in reports.items():
        found = checks(r)
        runs[name] = {
            "passed": passed(found),
            "checks": found,
            "stops": r["stops"],
            "standard_additions": {
                str(d): recovered(r, d, f"dose{d}")
                for d in doses
                if all(f"dose{d}" in c["predicted"] for c in r["cells"].values())
            },
            # each anchor alone, and the two substitutions only the overlap sees
            "anchors": {
                s: {op: _pair(r["cells"], op, s) for op in ("total", *UNIDENTIFIED_BY_STOPS)}
                for s in ("both", "stops", "overlap")
                if s in r["fits"]
            },
        }
    n = sum(x["passed"] is True for x in runs.values())
    tested = sum(x["passed"] is not None for x in runs.values())
    # no-go only when enough libraries failed that the rest could not make it up
    short = len(runs) >= RUNS_NEEDED and len(runs) - (tested - n) < RUNS_NEEDED
    go = "go" if n >= RUNS_NEEDED else "no-go" if short else "undecided"
    return {"verdict": go, "passed": n, "tested": tested, "libraries": len(runs), "runs": runs}


def _pair(cells: Mapping[str, Mapping[str, Any]], op: str, spec: str) -> dict[str, float]:
    use = [cells[f"op:{x}"] for x in (TEN if op == "total" else (op,)) if f"op:{x}" in cells]
    predicted, observed = sum(c["predicted"].get(spec, 0.0) for c in use), sum(c["observed"] for c in use)
    return {"ratio": predicted / observed if observed else math.nan, "observed": observed}


def markdown(v: Mapping[str, Any]) -> str:
    def lines() -> Iterator[str]:
        yield f"# P2 verdict: {v['verdict']}\n"
        yield f"{v['passed']} of {v['libraries']} libraries pass ({v['tested']} had power); {RUNS_NEEDED} are needed.\n"
        yield "| library | passed | error / genome / variant / indel / unaligned stops | total | floor (95 %) | outside | untested | additions |"
        yield "|---|---|---|---|---|---|---|---|"
        for name, r in v["runs"].items():
            by = {c["check"]: c for c in r["checks"]}
            total, floor = by["total"], by["floor"]
            named = {k: [c["check"] for c in r["checks"] if c["result"] == k] for k in ("outside", "untested")}
            stops = " / ".join(str(r["stops"][k]) for k in CLASSES)
            added = ", ".join(f"{d}: {x:.2f}" for d, x in r["standard_additions"].items())
            yield (
                f"| {name} | {r['passed']} | {stops} | {total['ratio']:.2f} ± {total['sd']:.2f} | "
                f"{floor['rate']:.1e} ({floor['upper95']:.1e}) | {', '.join(named['outside']) or '-'} | "
                f"{len(named['untested'])} | {added} |"
            )
        yield "\n## Anchors alone (predicted / observed)\n"
        yield "| library | fit | ten ops | A>G | T>C |"
        yield "|---|---|---|---|---|"
        for name, r in v["runs"].items():
            for s, a in r["anchors"].items():
                yield f"| {name} | {s} | {a['total']['ratio']:.2f} | {a['A>G']['ratio']:.2f} | {a['T>C']['ratio']:.2f} |"

    return "\n".join(lines()) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m sequencing_error_model.pfam_truth", description="Score pfam fits against known genomes."
    )
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("score", help="one library: classify its stops and compare its fits with the genomes")
    a.add_argument("r1", type=Path, help="`pfam reads` output, undosed")
    a.add_argument("r2", type=Path)
    a.add_argument("--pairs", type=Path, required=True, help="the undosed search's pairs.tsv")
    a.add_argument("--prep", type=Path, required=True)
    a.add_argument("--sam", type=Path, required=True, help="r1 and r2 aligned to the genomes (SAM or BAM)")
    a.add_argument("--reference", type=Path, required=True, help="the genomes, FASTA")
    a.add_argument("--spec", nargs="+", required=True, metavar="NAME=DIR", help="fits; `both` is the one tested")
    a.add_argument("--t1", type=float, default=35.0)
    a.add_argument("--t2", type=float, default=6.0)
    a.add_argument("--min-mapq", type=int, default=20)
    a.add_argument("--max-alt-freq", type=float, default=0.2, help="a site is a variant from this frequency on")
    a.add_argument("--output", type=Path, required=True, metavar="JSON")
    b = sub.add_parser("verdict", help="the go / no-go over the libraries' score reports")
    b.add_argument("reports", type=Path, nargs="+", help="score reports, named after their library")
    b.add_argument("--output", type=Path, required=True, metavar="PREFIX", help="PREFIX.json and PREFIX.md")
    args = p.parse_args(argv)
    if args.command == "verdict":
        v = verdict({x.stem: json.loads(x.read_text()) for x in args.reports})
        args.output.with_suffix(".json").write_text(json.dumps(v, indent=2))
        args.output.with_suffix(".md").write_text(markdown(v))
        print(markdown(v))
        return 0
    specs = {}
    for item in args.spec:
        name, _, path = item.partition("=")
        specs[name] = spec_io.load(Path(path))
    if "both" not in specs:
        p.error("--spec needs both=DIR")
    profiles = Profiles.load(args.prep)
    with pysam.FastaFile(str(args.reference)) as f:
        genome = {c: f.fetch(c).upper() for c in f.references}
    report = score(
        lambda: pfam.read_pairs(args.r1, args.r2),
        pfam.read_hits(args.pairs, profiles),
        profiles,
        args.sam,
        genome,
        specs,
        t1=args.t1,
        t2=args.t2,
        min_mapq=args.min_mapq,
        max_alt_freq=args.max_alt_freq,
    )
    args.output.write_text(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
