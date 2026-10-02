"""Key multiplicity: how many loci of an assembly share a skiver key (plan §6.3, phase 6b).

skiver v0.3.2 keys are literal k-mers read off both strands of a read (FracMinHash on the key), and a key's
value is the v bases that follow it on that strand. So every assembly position on both strands is one key
occurrence, and a key at several loci has several *true* values: skiver reads the minority ones as errors when
they are one edit from the consensus, and its outlier filter drops the key when they are frequent. At k = 11 on a
5 Mb isolate that is 90 % of the observations; at k = 21, 2 %.

- `occurrences` / `summary`: per-k multiplicity of an assembly, each occurrence weighted by its contig's coverage
  (`cov=` in shovill/SPAdes-style headers, else 1), a proxy for the observations it contributes.
- `by_loci`: a skiver `kvmer.csv` joined key by key to those loci, by multiplicity class - what passes the
  filter, and where the single-edit ("error") mass comes from.

`python -m sequencing_error_model.multiplicity ASSEMBLY -k 11 21 [--kvmer PREFIX_K11 PREFIX_K21] --output JSON`
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from sequencing_error_model.fit.quality import Array
from sequencing_error_model.generate import _fasta, _open, _revcomp
from sequencing_error_model.sources.skiver_analyze import KvmerRow, read_kvmer

_CODE = np.full(256, 4, np.int64)
_CODE[np.frombuffer(b"ACGT", np.uint8)] = np.arange(4)
CLASSES = ("0", "1", "2", "3-4", "5+")  # loci per key; 0 = the key is not in the assembly


def encode(seq: str) -> int:
    """2-bit integer of an ACGT string, first base most significant (as `occurrences` encodes windows)."""
    return int(seq.translate(str.maketrans("ACGT", "0123")), 4)


def _windows(codes: Array, w: int) -> tuple[Array, Array]:
    win = np.lib.stride_tricks.sliding_window_view(codes, w)
    return (win & 3) @ (4 ** np.arange(w - 1, -1, -1, dtype=np.int64)), (win < 4).all(axis=1)


def occurrences(contigs: Iterable[tuple[str, str, float]], k: int, v: int) -> tuple[Array, Array, Array]:
    """(key, value, weight) per key occurrence on both strands of (name, sequence, weight) contigs, N-free only."""
    if not (1 <= k <= 31 and 1 <= v <= 31):
        raise ValueError(f"k and v must be in 1..31 to fit 2-bit int64 windows (got k={k}, v={v})")
    keys, vals, wts = [], [], []
    for _, seq, weight in contigs:
        for s in (seq.upper(), _revcomp(seq.upper())):
            codes = _CODE[np.frombuffer(s.encode(), np.uint8)]
            if len(codes) < k + v:
                continue
            key, kok = _windows(codes[: len(codes) - v], k)
            val, vok = _windows(codes[k:], v)
            ok = kok & vok
            keys.append(key[ok])
            vals.append(val[ok])
            wts.append(np.full(int(ok.sum()), weight))
    if not keys:
        empty = np.zeros(0, np.int64)
        return empty, empty, np.zeros(0)
    key, val, w = np.concatenate(keys), np.concatenate(vals), np.concatenate(wts)
    order = np.lexsort((val, key))
    return key[order], val[order], w[order]


def one_edit(a: Array, b: Array, v: int) -> Array:
    """Whether each v-base value in `a` is one substitution, deletion or insertion from its counterpart in `b`.

    A shifted value is v bases long like the consensus, so an indel at i leaves a[:i] == b[:i] and the rest
    offset by one, with one base running past the end. Equal values count too; callers pass only differing ones.
    """
    x = a ^ b
    hit = sum((((x >> (2 * j)) & 3) != 0).astype(np.int64) for j in range(v)) == 1
    for i in range(v):
        lo = v - 1 - i
        rest = (1 << (2 * lo)) - 1
        b_less = ((b >> (2 * (v - i))) << (2 * lo)) | (b & rest)  # b without its base i
        a_less = ((a >> (2 * (v - i))) << (2 * lo)) | (a & rest)
        hit |= ((a >> 2) == b_less) | (a_less == (b >> 2))  # a lacks b's base i / a has an extra base at i
    return np.asarray(hit, bool)


def summary(key: Array, val: Array, w: Array, v: int) -> dict[str, float]:
    """Weighted shares of key observations by multiplicity, over `occurrences` output (sorted by key, value)."""
    newpair = np.r_[True, (key[1:] != key[:-1]) | (val[1:] != val[:-1])]
    pid = np.cumsum(newpair) - 1
    pkey, pval = key[newpair], val[newpair]
    pw, pocc = np.bincount(pid, weights=w), np.bincount(pid)
    kid = np.cumsum(np.r_[True, pkey[1:] != pkey[:-1]]) - 1
    loci, nval = np.bincount(kid, weights=pocc), np.bincount(kid)
    heaviest = np.lexsort((-pw, kid))
    first = np.r_[True, kid[heaviest][1:] != kid[heaviest][:-1]]
    majority = np.empty(int(kid.max()) + 1, np.int64)
    majority[kid[heaviest][first]] = pval[heaviest][first]
    off = pval != majority[kid]
    edit1 = np.zeros(len(pval), bool)
    edit1[off] = one_edit(pval[off], majority[kid[off]], v)
    total = float(w.sum())
    return {
        "keys": float(len(loci)),
        "occurrences": float(len(key)),
        "multi_locus": float(pw[(loci > 1)[kid]].sum()) / total,
        "multi_value": float(pw[(nval > 1)[kid]].sum()) / total,
        "off_majority": float(pw[off].sum()) / total,
        "off_majority_one_edit": float(pw[edit1].sum()) / total,  # read by skiver as a single-edit error
        "off_majority_multi_edit": float(pw[off & ~edit1].sum()) / total,  # dropped, but flags the key
    }


def _class(loci: int) -> str:
    return CLASSES[min(loci, 2)] if loci <= 2 else CLASSES[3] if loci <= 4 else CLASSES[4]


def by_loci(rows: Sequence[KvmerRow], key: Array, val: Array) -> dict[str, Any]:
    """skiver's per-key counts by the key's loci in the assembly, from `occurrences` output at skiver's k and v.

    Per class: keys, observations (`total_count`), the share passing skiver's filter, the share whose consensus is
    one of the key's true values, single-edit and other-value counts per observation, and the class's share of
    all single-edit mass and of the single-edit mass on keys passing the filter.
    """
    rows = [r for r in rows if not set(r.key + r.consensus_value) - set("ACGT")]
    uniq, start, loci = np.unique(key, return_index=True, return_counts=True)
    q = np.array([encode(r.key) for r in rows], np.int64)
    at = np.minimum(np.searchsorted(uniq, q), max(len(uniq) - 1, 0))
    found = (uniq[at] == q) if len(uniq) else np.zeros(len(q), bool)
    acc: dict[str, Array] = {c: np.zeros(7) for c in CLASSES}
    for r, f, i in zip(rows, found, at, strict=True):
        n = int(loci[i]) if f else 0
        true = f and encode(r.consensus_value) in set(val[start[i] : start[i] + n].tolist())
        other = r.total_count - r.consensus_count - r.neighbor_count
        acc[_class(n)] += [1, r.total_count, r.passes_filter, true, r.neighbor_count, other,
                           r.neighbor_count * r.passes_filter]  # fmt: skip
    edits = sum(a[4] for a in acc.values()) or 1.0
    passing = sum(a[6] for a in acc.values()) or 1.0
    out: dict[str, Any] = {
        "keys": len(rows),
        "pass_share": sum(a[2] for a in acc.values()) / max(len(rows), 1),
        "filter_removes_single_edit_mass": 1 - passing / edits,
        "classes": {},
    }
    for c, (n, obs, ok, true, single, other, single_ok) in acc.items():
        if n:
            out["classes"][c] = {
                "keys": int(n),
                "obs": int(obs),
                "pass_share": ok / n,
                "consensus_true_share": true / n,
                "single_edit_per_obs": single / max(obs, 1),
                "other_per_obs": other / max(obs, 1),
                "single_edit_mass_share": single / edits,
                "single_edit_mass_share_passing": single_ok / passing,
            }
    return out


def _coverage(path: Path) -> dict[str, float]:
    """Contig name -> `cov=` from FASTA headers, where present."""
    with _open(path, "r") as fh:
        return {
            line[1:].split()[0]: float(m.group(1))
            for line in fh
            if line.startswith(">") and (m := re.search(r"\bcov=([0-9.eE+-]+)", line))
        }


def contigs(path: Path) -> list[tuple[str, str, float]]:
    cov = _coverage(path)
    with _open(path, "r") as fh:
        return [(name, seq, cov.get(name, 1.0)) for name, seq in _fasta(fh)]


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("assembly", type=Path)
    p.add_argument("-k", type=int, nargs="+", default=[11, 21])
    p.add_argument("-v", type=int, default=13)
    p.add_argument("--kvmer", nargs="*", default=[], help="skiver analyze prefixes (or kvmer.csv), one per -k")
    p.add_argument("--output", type=Path, help="JSON report (default: stdout)")
    args = p.parse_args(argv)
    if args.kvmer and len(args.kvmer) != len(args.k):
        p.error("--kvmer takes one prefix per -k")
    seqs = contigs(args.assembly)
    report: dict[str, Any] = {}
    for i, k in enumerate(args.k):
        key, val, w = occurrences(seqs, k, args.v)
        report[str(k)] = {"assembly": summary(key, val, w, args.v)}
        if args.kvmer:
            src = Path(args.kvmer[i])
            path = src if src.suffix == ".csv" else Path(f"{src}.kvmer.csv")
            report[str(k)]["kvmer"] = by_loci(read_kvmer(path), key, val)
    text = json.dumps(report, indent=2)
    if args.output:
        args.output.write_text(text + "\n")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
