"""Prototype: P2's redesign, scored on P2's data as development, never as a verdict (pfam_mode_plan.md §7 P2).

P2's diagnosis kept three of the four pieces, each for one job, and dropped the mate-overlap anchor:

- **level**: the stop anchor's fitted ten-op rate per base (P2's `stops` fit, read from P2's score report);
- **shape** over Q group, read-position quarter and mate: `stops` (the default) keeps the stop fit's own curve; or the
  within-locus contrast, `mismatch with the locus'
  majority = b(locus) + e(cell)`, on the ten ops stops identify. `additive` (the contrast's own form): the level plus
  e(cell), centred over the cells' bases; `relative` (the default, better in the small high-Q cells): (mean + e) /
  mean times the level. Both carry the shape of the locus bases, which are cleaner than the rest by different amounts
  in different Q groups (P2's diagnosis), so neither is the library's shape exactly;
- **op mix**: each substitution's locus label rate (a base against the unanimous base of two or more other reads),
  scaled so the ten ops sum to the level. It is the only source of A>G / T>C here. With `--locus`, each op's labels
  are scaled by the share that P2's locus mixture calls errors (`error_share`), which takes out strains' transitions;
  not where omega >= `MAX_OMEGA`, where the mixture cannot split biology from error.

Nothing here reads the truth: the shape and the op mix come from the reads alone (`pfam_p2_cells` tables, built from
reads and hits), the level from the stop fit. The truth enters only to judge, against P2's score cells (masked
variant sites, MAPQ >= 20) under P2's rule (`pfam_truth.checks`): ten-op total within 15 %, no cell outside 25 % + 2 sd.

    python docs/pfam_redesign_proto.py SCORE.json TABLE.npz [--output OUT.json]
"""

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import pfam_contrast_proto as proto  # noqa: E402
import pfam_p2_cells as cells  # noqa: E402

from sequencing_error_model import pfam_truth  # noqa: E402

BASES = "ACGT"
TEN = np.array([[f"{a}>{b}" in pfam_truth.TEN for b in BASES] for a in BASES])  # source x read base
COMP = str.maketrans("ACGT", "TGCA")
MAX_OMEGA = 0.5  # above it, biology looks like error (missense near synonymous) and the mixture cannot split them


def error_share(locus: dict[str, Any], flank: int = 12, others: int = 2) -> dict[str, float] | None:
    """Per read-orientation substitution, the share of locus labels that are errors, from P2's mixture (`pfam_locus`).

    The mixture is on the coding strand; reads cover both strands of a gene, so a read's X>Y is the coding X>Y on
    one and its complement on the other. Each op takes its complementary pair's share: strains hit both strands
    alike, and the label's own strand asymmetry (8-oxoG G>T against C>A) is kept. None where omega is too high.
    """
    row = next(x for x in locus["loci"] if x["flank"] == flank and x["others"] == others)
    if not (row["omega"] < MAX_OMEGA) or not row["ops"]:
        return None
    ops = locus["ops"]
    errors = {op: row["ops"][op] * sum(row["chance"][i]) for i, op in enumerate(ops)}
    labels = {op: sum(row["labels"][i]) for i, op in enumerate(ops)}
    out = {}
    for op in ops:
        pair = (op, op.translate(COMP))
        n = sum(labels[x] for x in pair)
        out[op] = min(sum(errors[x] for x in pair) / n, 1.0) if n else 1.0
    return out


COLUMNS = {"locus": np.int64, "base": np.int8, "q": np.uint8, "cycle": np.int16, "length": np.int16, "mate": np.int8}


def load(path: Path, max_rows: int = 50_000_000) -> tuple[dict[str, np.ndarray], int]:
    """Rows in loci, aligned or not (nothing here needs the genome), and every k-th locus beyond `max_rows`.

    Both the contrast and the op mix use loci only, and both stay unbiased on a subset of loci: this bounds memory
    at P2's ~240 M rows a library.
    """
    with np.load(path) as z:
        locus = z["locus"]
        on = locus >= 0
        every = max(1, math.ceil(int(on.sum()) / max_rows))
        keep = on & (locus % every == 0)
        del locus, on
        return {k: z[k][keep].astype(dtype) for k, dtype in COLUMNS.items()}, every


def schemes_names() -> list[str]:
    return [
        *(f"q:{x}+" for x in (0, *pfam_truth._Q_EDGES)),
        *(f"position:{i}" for i in range(pfam_truth._POSITIONS)),
        "mate:1",
        "mate:2",
    ]


def schemes(t: dict[str, np.ndarray]) -> dict[str, tuple[np.ndarray, list[str]]]:
    """Each row's Q, position and mate cell, as `pfam_truth._Cells.groups` names them."""
    position = (t["cycle"].astype(np.int32) * pfam_truth._POSITIONS) // (t["length"].astype(np.int32) + 1)
    edges = (0, *pfam_truth._Q_EDGES)
    return {
        "q": (np.digitize(t["q"], pfam_truth._Q_EDGES).astype(np.int8), [f"q:{x}+" for x in edges]),
        "position": (position.astype(np.int8), [f"position:{i}" for i in range(pfam_truth._POSITIONS)]),
        "mate": (t["mate"], ["mate:1", "mate:2"]),
    }


def predict(
    score: dict[str, Any], t: dict[str, np.ndarray], shape: str = "stops", share: dict[str, float] | None = None
) -> dict[str, float]:
    """Predicted errors per cell of `score` (P2's cells), from the stop fit's level and the reads."""
    c = score["cells"]
    bases = c["mate:1"]["exposure"] + c["mate:2"]["exposure"]
    ten_ops = [f"op:{op}" for op in pfam_truth.TEN]
    level = sum(c[k]["predicted"]["stops"] for k in ten_ops if k in c) / bases  # ten-op errors per base
    counts = cells.locus_counts(t)
    out: dict[str, float] = {}
    if shape == "stops":  # the stop fit's own Q / position / mate curve
        out.update({name: c[name]["predicted"]["stops"] for name in schemes_names() if name in c})
    else:
        # shape: the contrast on ten-op mismatches with the locus majority, every read of a locus of three or more
        at = np.where(t["locus"] >= 0, t["locus"], 0)
        deep = (t["locus"] >= 0) & (counts.sum(1)[at] >= proto.MIN_DEPTH["locus"])
        group = t["locus"][deep]
        majority = counts.argmax(1)[group]
        y = ((t["base"][deep] != majority) & TEN[majority, t["base"][deep]]).astype(np.float64)
        del at, majority
        ybar = float(y.mean())
        for z, names in schemes(t).values():
            e = proto.contrast(group, z[deep], y)
            n = np.array([c[x]["exposure"] if x in c else 0.0 for x in names[: len(e)]])
            if shape == "additive":
                rate = level + e - n @ e / n.sum()  # centred over P2's bases in these cells
            else:
                relative = (ybar + e) / ybar
                rate = level * relative / (n @ relative / n.sum())
            for name, r, k in zip(names, rate, n, strict=False):
                out[name] = max(float(r), 0.0) * k
        del y, group, deep
    # op mix: label rates per substitution, scaled so the ten ops sum to the level
    lab, cons = cells.labels(t, counts)
    src, got = cons[lab], t["base"][lab]
    pairs = np.bincount(src.astype(np.int64) * 4 + got, minlength=16).reshape(4, 4).astype(float)
    per_source = np.bincount(src, minlength=4).astype(float)
    if share is not None:  # biology out of the labels, op by op
        pairs *= np.array([[share.get(f"{a}>{b}", 1.0) if a != b else 1.0 for b in BASES] for a in BASES])
    rate = pairs / np.maximum(per_source, 1)[:, None]  # per base of the source
    ten_label = float((pairs * TEN).sum() / per_source.sum())
    for a in range(4):
        for b in range(4):
            if a != b and f"op:{BASES[a]}>{BASES[b]}" in c:
                key = f"op:{BASES[a]}>{BASES[b]}"
                out[key] = float(rate[a, b]) * level / ten_label * c[key]["exposure"]
    out["level"] = level
    return out


def judge(score: dict[str, Any], predicted: dict[str, float]) -> dict[str, Any]:
    """P2's checks (`pfam_truth.checks`) with the prototype's predictions in place of a fit's."""
    report = {
        **score,
        "cells": {
            k: {**v, "predicted": {**v["predicted"], "redesign": predicted.get(k, v["predicted"]["stops"])}}
            for k, v in score["cells"].items()
        },
    }
    found = pfam_truth.checks(report, "redesign")
    return {"passed": pfam_truth.passed(found), "checks": found}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("score", type=Path, help="P2's score/<run>.json")
    p.add_argument("table", type=Path, help="the run's pfam_p2_cells table, .npz")
    p.add_argument("--shape", choices=("stops", "relative", "additive"), default="stops")
    p.add_argument("--max-rows", type=int, default=50_000_000, help="every k-th locus beyond this many rows")
    p.add_argument("--locus", type=Path, help="P2's locus/<run>.0.json: take biology out of the op mix by its mixture")
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    score = json.loads(a.score.read_text())
    table, every = load(a.table, a.max_rows)
    share = error_share(json.loads(a.locus.read_text())) if a.locus else None
    predicted = predict(score, table, a.shape, share)
    del table
    out = {
        "predicted": predicted,
        "redesign": judge(score, predicted),
        "shape": a.shape,
        "loci_every": every,
        "error_share": share,
    }
    for spec in ("both", "stops"):
        found = pfam_truth.checks(score, spec)
        out[spec] = {"passed": pfam_truth.passed(found), "checks": found}
    by = {s: {c["check"]: c for c in out[s]["checks"]} for s in ("both", "stops", "redesign")}
    print(f"passed: both {out['both']['passed']}, stops {out['stops']['passed']}, redesign {out['redesign']['passed']}")
    print(f"{'check':>11} {'both':>14} {'stops':>14} {'redesign':>14}")
    # P2 does not judge A>G and T>C (stops cannot identify them): reported beside the checks
    for op in pfam_truth.UNIDENTIFIED_BY_STOPS:
        cell = score["cells"].get(f"op:{op}")
        if cell and cell["observed"]:
            fits = {s: cell["predicted"][s] / cell["observed"] for s in ("both", "stops")}
            fits["redesign"] = predicted[f"op:{op}"] / cell["observed"]
            out.setdefault("unjudged", {})[op] = fits
            print(f"{op} (not judged): " + ", ".join(f"{s} {r:.2f}" for s, r in fits.items()))
    if a.output:
        a.output.write_text(json.dumps(out, indent=1, default=float))
    for name in by["redesign"]:
        row = []
        for s in ("both", "stops", "redesign"):
            x = by[s].get(name, {})
            if "ratio" in x:
                mark = {"within": " ", "outside": "*", "unclear": "?", "untested": "-"}[x["result"]]
                row.append(f"{x['ratio']:>6.2f} ± {x['sd']:.2f}{mark}")
            else:
                row.append(f"{x.get('result', ''):>14}")
        print(f"{name:>11} " + " ".join(f"{r:>14}" for r in row))


if __name__ == "__main__":
    main()
