"""P2 step 2: the locus anchor and the within-locus contrast, cell by cell, against the genomes (pfam_mode_plan.md §7 P2).

P2 judged the fits per op, Q group, mate and read-position quarter, but showed the locus anchor only as a total, and on
the bases it labels. This scores both per cell against the truth over **all** exposed bases, the truth P2's checks use:

- `truth`: errors per base of every exposed codon, by the genome under it (`pfam_contrast_proto.add_truth`);
- `truth_locus`: the same on the bases the locus labels; equal to `truth` in a cell when the locus selects only
  through the cell's covariates, unequal when its clean-flank, three-fragment rule picks cleaner reads within it;
- `locus`: the label rate, a base against the unanimous base of two or more other reads of its locus (same column,
  base, strand and 12 bases each side); errors and singleton biology, as `pfam_locus` counts before its mixture;
- `contrast`: `mismatch with the locus' majority = b(locus) + e(cell)` over loci of three or more reads
  (`pfam_contrast_proto.contrast`, its `by="locus"` rows); only e's shape is identified, so its level is set to
  `truth` over all bases, and the check is of the shape alone.

Q is a covariate throughout, never a label. Ratios use P2's rule (`pfam_truth.CELL`): `within` 25 %, `outside` by more
than 25 % + 2 sd, tested only where the sd is at most half the tolerance.

    python docs/pfam_p2_cells.py R1 R2 --pairs pairs.tsv --prep PREP --sam aligned.sam --reference genomes.fna \\
        --table run.npz --output run.json      # builds run.npz once; later runs with --table alone reuse it
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import pfam_contrast_proto as proto  # noqa: E402

from sequencing_error_model.pfam_truth import CELL  # noqa: E402

BASES = "ACGT"
Q_EDGES = (20, 30, 36)  # pfam_truth's Q groups
Q_NAMES = ("q:0+", "q:20+", "q:30+", "q:36+")
OTHERS = 2  # unanimous other reads that make a label (pfam_locus)


def labels(t: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """(labelled, the other reads' base) per row: the base of its locus' other reads, where two or more agree on it."""
    on = t["locus"] >= 0
    _, locus = np.unique(np.where(on, t["locus"], -1), return_inverse=True)
    counts = np.zeros((locus.max() + 1, 4), np.int64)
    np.add.at(counts, (locus[on], t["base"][on]), 1)
    others = counts[locus].copy()
    others[np.arange(len(locus)), t["base"]] -= 1
    depth, top = others.sum(1), others.argmax(1)
    return on & (depth >= OTHERS) & (others.max(1) == depth), top


def judge(ratio: float, sd: float, tolerance: float = CELL) -> str:
    off = abs(ratio - 1)
    if not sd <= tolerance / 2:
        return "untested"
    return "within" if off <= tolerance else "outside" if off > tolerance + 2 * sd else "unclear"


def _rate(events: float, bases: float) -> float:
    return events / bases if bases else math.nan


def score(t: dict[str, np.ndarray], boot: int = 5) -> dict[str, dict[str, object]]:
    known = t["true"] >= 0
    t = {k: v[known] for k, v in t.items()}
    err = t["base"] != t["true"]
    lab, cons = labels(t)
    mism = lab & (t["base"] != cons)
    schemes = {
        "q": (np.digitize(t["q"], Q_EDGES), Q_NAMES),
        "position": (np.minimum(4 * t["cycle"] // np.maximum(t["length"], 1), 3), [f"position:{i}" for i in range(4)]),
        "mate": (t["mate"], ["mate:1", "mate:2"]),
    }
    # the contrast, as the prototype fits it: every read of a locus of three or more against the locus' majority
    # (a leave-one-out label would leave an erring read alone in its locus, and b(locus) would absorb it)
    deep = (t["locus"] >= 0) & (
        np.bincount(np.maximum(t["locus"], 0))[np.maximum(t["locus"], 0)] >= proto.MIN_DEPTH["locus"]
    )
    _, group = np.unique(t["locus"][deep], return_inverse=True)
    majority = np.zeros((group.max() + 1 if len(group) else 0, 4))
    np.add.at(majority, (group, t["base"][deep]), 1)
    y = (t["base"][deep] != majority.argmax(1)[group]).astype(float)
    out: dict[str, dict[str, object]] = {}

    def cell(name: str, sel: np.ndarray, extra: dict[str, float] | None = None) -> None:
        n, e = int(sel.sum()), int(err[sel].sum())
        nl, el, ml = int((sel & lab).sum()), int((err & sel & lab).sum()), int((mism & sel).sum())
        truth, locus = _rate(e, n), _rate(ml, nl)
        ratio = locus / truth if truth else math.nan
        sd = ratio * math.sqrt(1 / ml + 1 / e) if ml and e else math.inf
        out[name] = {
            "bases": n,
            "truth": truth,
            "truth_locus": _rate(el, nl),
            "locus": locus,
            "locus_over_truth": ratio,
            "sd": sd,
            "result": judge(ratio, sd),
            **(extra or {}),
        }

    cell("total", np.ones(len(err), bool))
    level = err.mean()
    for z, names in schemes.values():
        shape = proto.contrast(group, z[deep], y)
        # loci resampled (Poisson weights), as the prototype: the sd of the shape, not of the level
        rng = np.random.default_rng(0)
        spread = np.std(
            [proto.contrast(group, z[deep], y, rng.poisson(1.0, group.max() + 1)[group]) for _ in range(boot)], 0
        )
        for i, name in enumerate(names):
            sel = z == i
            truth = float(err[sel].mean()) if sel.any() else math.nan
            predicted = level + float(shape[i]) if i < len(shape) else math.nan
            ratio = predicted / truth if truth else math.nan
            sd = float(spread[i]) / truth if truth and i < len(spread) and boot > 1 else math.inf
            cell(
                name,
                sel,
                {
                    "contrast": predicted,
                    "contrast_over_truth": ratio,
                    "contrast_sd": sd,
                    "contrast_result": judge(ratio, sd),
                },
            )
    # ops in read orientation: rate per base of the source base, as pfam_truth's op cells
    for a in range(4):
        for b in range(4):
            if a == b:
                continue
            src_t, src_l = t["true"] == a, lab & (cons == a)
            e, n = int((src_t & (t["base"] == b)).sum()), int(src_t.sum())
            ml, nl = int((src_l & (t["base"] == b)).sum()), int(src_l.sum())
            truth, locus = _rate(e, n), _rate(ml, nl)
            ratio = locus / truth if truth else math.nan
            sd = ratio * math.sqrt(1 / ml + 1 / e) if ml and e else math.inf
            out[f"op:{BASES[a]}>{BASES[b]}"] = {
                "truth": truth,
                "locus": locus,
                "locus_over_truth": ratio,
                "sd": sd,
                "result": judge(ratio, sd),
            }
    out["total"]["labelled_share"] = float(lab.mean())
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("r1", type=Path, nargs="?")
    p.add_argument("r2", type=Path, nargs="?")
    p.add_argument("--pairs", type=Path)
    p.add_argument("--prep", type=Path)
    p.add_argument("--sam", type=Path)
    p.add_argument("--reference", type=Path)
    p.add_argument("--t1", type=float, default=35.0)
    p.add_argument("--t2", type=float, default=6.0)
    p.add_argument("--table", type=Path, required=True, help=".npz of the per-base table: built if missing")
    p.add_argument("--output", type=Path)
    p.add_argument("--boot", type=int, default=5, help="locus resamples for the contrast's sd")
    a = p.parse_args()
    if not a.table.exists():
        if not (a.r1 and a.r2 and a.pairs and a.prep and a.sam and a.reference):
            p.error("building the table needs R1 R2 --pairs --prep --sam --reference")
        a.every_n, a.dose = 1, 0.0
        np.savez_compressed(a.table, allow_pickle=False, **proto.table(a))
    out = score(dict(np.load(a.table)), a.boot)
    if a.output:
        a.output.write_text(json.dumps(out, indent=1))
    t = out["total"]
    print(
        f"{t['bases']:,} exposed bases, {t['labelled_share']:.0%} labelled; truth {t['truth']:.2e}, "
        f"truth on labelled {t['truth_locus']:.2e}, locus {t['locus']:.2e} ({t['locus_over_truth']:.2f} ± {t['sd']:.2f})"
    )
    print(
        f"{'cell':>10} {'truth':>9} {'on locus':>9} {'locus':>9} {'ratio':>12} {'result':>9} {'contrast (shape only)':>24}"
    )
    for name, c in out.items():
        if name == "total":
            continue
        shape = (
            f"{c['contrast_over_truth']:.2f} ± {c['contrast_sd']:.2f} {c['contrast_result']}" if "contrast" in c else ""
        )
        print(
            f"{name:>10} {c['truth']:>9.2e} {c.get('truth_locus', math.nan):>9.2e} {c['locus']:>9.2e} "
            f"{c['locus_over_truth']:>5.2f} ± {c['sd']:.2f} {c['result']:>9} {shape:>24}"
        )


if __name__ == "__main__":
    main()
