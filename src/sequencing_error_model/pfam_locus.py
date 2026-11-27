"""Locus consensus as an error label, reported for docs/pfam_mode_plan.md section 12 (not part of the `pfam` fit).

Reads that share the bases either side of a codon (the same Pfam column, the same `w` bases before and after) are one
locus. A codon that differs from the codon all the other fragments of its locus share, `others` of them or more, is an
error or a strain's variant. The codon's effect separates the two: errors make synonymous, missense and nonsense
changes at equal rates per opportunity, biology makes synonymous >> missense >> nonsense. `mixture` fits, per
substitution on the coding strand,

    rate(effect) = E + b x (1, omega, 0)        synonymous, missense, nonsense

with omega the sample's own missense to synonymous ratio among alleles that recur, which are biology. Singleton
variants are younger and less purged than recurring ones, so omega is low and E high by an amount only truth can tell:
with the genomes of the sample (`--sam`, `--reference`) each label is an `error` (the consensus is the genome's),
`biology` (the read is) or `other`, and the omega of the labels that are biology is reported. Q is not used.
"""

import argparse
import hashlib
import itertools
import json
import math
import sys
from array import array
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pysam
from scipy.optimize import minimize

from sequencing_error_model import pfam_truth
from sequencing_error_model.fit.quality import Array
from sequencing_error_model.generate import _revcomp
from sequencing_error_model.sources import pfam
from sequencing_error_model.sources.pfam import Codon, Profiles

CODONS = tuple(pfam.CODE)
INDEX = {c: i for i, c in enumerate(CODONS)}
OPS = tuple(f"{a}>{b}" for a in "ACGT" for b in "ACGT" if a != b)
EFFECTS = ("synonymous", "missense", "nonsense")
FLANKS = (9, 12, 18)  # bases each side
LEAST = 3  # fragments of a second allele that make it biology
TRUTHS = ("error", "biology", "other")
_AMINO = np.array([pfam.CODE[c] for c in CODONS])
_LETTER = np.array([list(c) for c in CODONS])
DIST = (_LETTER[:, None] != _LETTER[None]).sum(2)
STOP = _AMINO == "*"
# index into EFFECTS of a change from a sense codon
EFFECT = np.where(STOP[None], 2, np.where(_AMINO[:, None] == _AMINO[None], 0, 1))
OP = np.full((64, 64), -1)
CHANCE = np.zeros((64, 12, 3))  # changes open to each codon, by substitution and effect
for _i, _j in zip(*np.nonzero(DIST == 1), strict=True):
    _p = int(np.flatnonzero(_LETTER[_i] != _LETTER[_j])[0])
    OP[_i, _j] = OPS.index(f"{_LETTER[_i, _p]}>{_LETTER[_j, _p]}")
    CHANCE[_i, OP[_i, _j], EFFECT[_i, _j]] += 1
# columns of `table`
CODON, PAIR, TRUE, READ, START, REVERSE, FIRST = range(7)


def _digest(text: str) -> int:
    return int.from_bytes(hashlib.blake2b(text.encode(), digest_size=8).digest(), "little", signed=True)


def table(
    pairs: Iterable[pfam.ReadPair],
    hits: Mapping[int, Sequence[pfam.Hit]],
    profiles: Profiles,
    t1: float = 35.0,
    t2: float = 6.0,
) -> Array:
    """One row per exposed codon: its index, fragment, true codon (-1), read, start on its strand, strand, and a key
    per entry of `FLANKS` (0 where the read ends inside the flank)."""
    # ponytail: 72 bytes a codon, about 5 GB for 9 M reads of a mock; shard by family if that is too much
    rows = array("q")
    for number, (s1, _, s2, _) in enumerate(pairs):
        for mate, seq in enumerate((s1, s2)):
            rid = 2 * number + mate
            if rid not in hits:
                continue
            s = seq.upper()
            for c in pfam.exposed(s, hits[rid], profiles, t1, t2):
                strand = _revcomp(s) if c.reverse else s
                at = len(s) - 1 - c.positions[0] if c.reverse else c.positions[0]
                rows.extend((INDEX.get(c.bases, -1), number, -1, rid, at, c.reverse))
                rows.extend(
                    _digest(f"{c.column}{strand[at - w : at]}|{strand[at + 3 : at + 3 + w]}")
                    if at >= w and at + 3 + w <= len(s)
                    else 0
                    for w in FLANKS
                )
    return np.array(rows, np.int64).reshape(-1, FIRST + len(FLANKS))


def add_truth(
    rows: Array, sam: Path, genome: Mapping[str, str], min_mapq: int = 0, masks: Mapping[str, Array] | None = None
) -> None:
    """Fill `TRUE` with the genome's codon, from the primary alignment of each read (named by pair number).

    MAPQ 0 is kept by default: strains that tie on a read are alike under it, and dropping them would drop a mock's
    strains. `masks` (`pfam_truth.sam_variant_masks`) leaves codons over a strain's variant site unknown, with
    `min_mapq` 20 the truth `pfam_truth` judges P2 by.
    """
    reads = rows[:, READ].copy()  # contiguous: searching the column itself would copy it every time
    with pysam.AlignmentFile(str(sam), check_sq=False) as f:
        for r in f.fetch(until_eof=True):
            name, contig = r.query_name, r.reference_name
            if r.is_secondary or r.is_supplementary or r.is_unmapped or name is None or contig is None:
                continue
            rid = 2 * int(name) + int(r.is_read2)
            lo, hi = np.searchsorted(reads, (rid, rid + 1))
            n = r.infer_read_length()
            if lo == hi or r.mapping_quality < min_mapq or n is None:
                continue
            placed = pfam_truth.Placement(contig, bool(r.is_reverse), pfam_truth.sites(r, n))
            for k in range(lo, hi):
                at, reverse = int(rows[k, START]), bool(rows[k, REVERSE])
                where = [n - 1 - at - p if reverse else at + p for p in range(3)]
                codon = Codon((0, 0), (where[0], where[1], where[2]), reverse, "", 0.0)
                mask = masks.get(contig) if masks is not None else None
                if mask is not None and any(placed.sites[x] >= 0 and mask[placed.sites[x]] for x in where):
                    continue
                if found := pfam_truth.true_bases(codon, placed, genome):
                    rows[k, TRUE] = INDEX.get(found[1], -1)


def _fragments(group: Array, codon: Array, pair: Array) -> tuple[Array, Array, Array]:
    """(row, group number, codon) with each fragment once per group."""
    o = np.lexsort((pair, group))
    o = o[np.r_[True, (np.diff(group[o]) != 0) | (np.diff(pair[o]) != 0)]]
    return o, np.unique(group[o], return_inverse=True)[1], codon[o]


def consensus(group: Array, codon: Array, pair: Array, others: int = 2) -> tuple[Array, Array]:
    """(row, the codon shared by all the other fragments of its group), where there are `others` (>= 2) of them."""
    if others < 2:
        raise ValueError("one other fragment is not a consensus: which of the two erred is open")
    if not len(group):
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    row, g, c = _fragments(group, codon, pair)
    key, n = np.unique(g * 64 + c, return_counts=True)
    s = np.lexsort((-n, key // 64))
    first = np.flatnonzero(np.r_[True, np.diff(key[s] // 64) != 0])
    depth, alleles = np.bincount(g), np.diff(np.r_[first, len(s)])
    top, top_n = key[s][first] % 64, n[s][first]
    odd = (alleles[g] == 2) & (top_n[g] == depth[g] - 1) & (c != top[g])
    use = ((alleles[g] == 1) | odd) & (depth[g] > others)
    return row[use], top[g][use]


def recurring(group: Array, codon: Array, pair: Array, least: int = LEAST) -> Array:
    """Second alleles one change from the first and seen in `least` fragments or more, by substitution and effect."""
    out = np.zeros((12, 3))
    if not len(group):
        return out
    _, g, c = _fragments(group, codon, pair)
    key, n = np.unique(g * 64 + c, return_counts=True)
    at = np.flatnonzero((np.bincount(key // 64) == 2)[key // 64])[::2]
    swap = n[at] < n[at + 1]
    major = np.where(swap, key[at + 1], key[at]) % 64
    minor = np.where(swap, key[at], key[at + 1]) % 64
    ok = (np.minimum(n[at], n[at + 1]) >= least) & (DIST[major, minor] == 1) & ~STOP[major]
    np.add.at(out, (OP[major[ok], minor[ok]], EFFECT[major[ok], minor[ok]]), 1)
    return out


def counts(want: Array, got: Array) -> Array:
    """One-change disagreements by substitution and effect."""
    out = np.zeros((12, 3))
    one = DIST[want, got] == 1
    np.add.at(out, (OP[want[one], got[one]], EFFECT[want[one], got[one]]), 1)
    return out


def omega(n: Array, chance: Array) -> float:
    """Missense over synonymous rate, pooled over the substitutions (Mantel-Haenszel); nan without synonymous ones."""
    both = chance[:, 0] + chance[:, 1]
    w = np.divide(1.0, both, out=np.zeros_like(both), where=both > 0)
    below = float((n[:, 0] * chance[:, 1] * w).sum())
    return float((n[:, 1] * chance[:, 0] * w).sum()) / below if below > 0 else float("nan")


def mixture(n: Array, chance: Array, ratio: float) -> Array:
    """E per substitution: `n` ~ Poisson(chance x (E + b x (1, ratio, 0))), by maximum likelihood."""
    w = np.array([1.0, ratio, 0.0])
    seen = chance > 0

    def loss(p: Array) -> float:
        mu = chance * (np.exp(p[:12])[:, None] + np.exp(p[12:])[:, None] * w)
        return float((mu - n * np.log(np.where(seen, mu, 1.0) + 1e-300)).sum())

    fits = (
        minimize(loss, np.log([e] * 12 + [b] * 12), method="L-BFGS-B", bounds=[(-25, 0)] * 24)
        for e, b in itertools.product((1e-4, 1e-3), repeat=2)
    )
    return np.asarray(np.exp(min(fits, key=lambda f: f.fun).x[:12]))


def _per_base(rate: Array, chance: Array) -> float:
    """Per-substitution rates as errors per base, on the substitutions the codons are open to (nine per codon)."""
    total = chance.sum()
    return float(3 * (rate * chance.sum(1)).sum() / total) if total else float("nan")


def _ratio(n: Array, chance: Array) -> Array:
    return np.asarray(np.divide(n, chance, out=np.zeros_like(n, dtype=float), where=chance > 0))


def _fitted(n: Array, chance: Array, again: Array) -> tuple[float, Array | None]:
    ratio = omega(again, chance)
    return ratio, mixture(n, chance, ratio) if n.sum() and np.isfinite(ratio) else None


def selection(rows: Array, flank: int, others: int = 2) -> dict[str, Any]:
    """How clean the codons a locus labels are, against all exposed codons, measured by the stop anchor.

    A codon is *keyed* when its flanks match a locus whose other fragments, `others` or more, share one sense codon.
    That never looks at the codon itself, and a stop is labelled by the genetic code, not by a flank, so the stop rate
    on keyed codons over that on all of them is how much the clean-flank rule selects against error: a k-mer key's
    bias (docs/pfam_mode_plan.md §12), measured without a genome. Codons whose locus agrees on a stop are biology and
    left out of both. Stops are ten ops in stop-adjacent contexts, so the factor holds for all error only if the
    selection does not depend on the op; with truth, `truth` gives the factor over every base to check that.
    """
    exposed = rows[rows[:, CODON] >= 0]
    keyed_at = np.flatnonzero(exposed[:, FIRST + flank] != 0)
    sub = exposed[keyed_at]
    row, want = consensus(sub[:, FIRST + flank], sub[:, CODON], sub[:, PAIR], others)
    # `consensus` keeps one row per fragment; both mates of one are keyed alike
    ids = np.unique(sub[:, [FIRST + flank, PAIR]], axis=0, return_inverse=True)[1].ravel()
    verdict = np.full(int(ids.max()) + 1 if len(ids) else 0, -1)
    verdict[ids[row]] = STOP[want]
    keyed, biology = np.zeros(len(exposed), bool), np.zeros(len(exposed), bool)
    keyed[keyed_at[verdict[ids] == 0]] = True
    biology[keyed_at[verdict[ids] == 1]] = True
    stop, rest = STOP[exposed[:, CODON]], ~keyed & ~biology
    sk, su = int(stop[keyed].sum()), int(stop[rest].sum())
    nk, nu = int(keyed.sum()), int(rest.sum())
    rate_k, rate_all = sk / nk if nk else math.nan, (sk + su) / (nk + nu) if nk + nu else math.nan
    ratio = rate_k / rate_all if rate_all else math.nan
    # ln ratio = ln r_k - ln(f r_k + (1 - f) r_u): both counts enter through the share of stops not keyed
    sd = ratio * (su / (sk + su)) * math.sqrt(1 / sk + 1 / su) if sk and su else math.inf
    out: dict[str, Any] = {
        "codons_keyed": nk,
        "codons": nk + nu,
        "stops_keyed": sk,
        "stops": sk + su,
        "stop_rate_keyed": rate_k,
        "stop_rate": rate_all,
        "keyed_over_all": ratio,
        "sd": sd,
    }
    true, codon = exposed[:, TRUE], exposed[:, CODON]
    known = (true >= 0) & ~biology
    if known.any():
        wrong = DIST[np.where(known, true, 0), codon]
        error_stop = stop & ~STOP[np.where(known, true, 0)]
        k = known & keyed
        n_k, n_all = max(int(k.sum()), 1), int(known.sum())
        per_k, per_all = wrong[k].sum() / (3 * n_k), wrong[known].sum() / (3 * n_all)
        stops_k, stops_all = error_stop[k].sum() / n_k, error_stop[known].sum() / n_all
        out["truth"] = {
            "per_base_keyed": float(per_k),
            "per_base": float(per_all),
            "keyed_over_all": float(per_k / per_all) if per_all else math.nan,
            "error_stops_keyed_over_all": float(stops_k / stops_all) if stops_all else math.nan,
        }
    return out


def report(rows: Array, flank: int, others: int = 2, bootstrap: int = 25, seed: int = 0) -> dict[str, Any]:
    """Counts, opportunities and the mixture for one flank (an index into `FLANKS`), with truth where it is known."""
    use = rows[(rows[:, FIRST + flank] != 0) & (rows[:, CODON] >= 0)]
    group, codon, pair = use[:, FIRST + flank], use[:, CODON], use[:, PAIR]
    row, want = consensus(group, codon, pair, others)
    row, want = row[~STOP[want]], want[~STOP[want]]
    got, true = codon[row], use[row, TRUE]
    chance = np.tensordot(np.bincount(want, minlength=64).astype(float), CHANCE, 1)
    n, again = counts(want, got), recurring(group, codon, pair)
    ratio, e = _fitted(n, chance, again)
    rng = np.random.default_rng(seed)
    # the sd of the counts alone: labels and recurring alleles, not omega's bias
    redo = [_fitted(rng.poisson(n), chance, rng.poisson(again))[1] for _ in range(bootstrap if e is not None else 0)]
    spread = [_per_base(x, chance) for x in redo if x is not None]
    nonsense = chance[:, 2].sum()
    out: dict[str, Any] = {
        "flank": FLANKS[flank],
        "others": others,
        "codons": len(row),
        "share_of_exposed": len(row) / max(len(rows), 1),
        "more_than_one_change": int((DIST[want, got] > 1).sum()),
        "labels": n.astype(int).tolist(),
        "chance": chance.astype(int).tolist(),
        "recurring": again.astype(int).tolist(),
        "omega": ratio,
        "per_base": {
            "naive": _per_base(_ratio(n.sum(1), chance.sum(1)), chance),
            "nonsense": float(3 * n[:, 2].sum() / nonsense) if nonsense else float("nan"),
            "mixture": float("nan") if e is None else _per_base(e, chance),
            "mixture_sd": float(np.std(spread)) if spread else float("nan"),
        },
        "ops": {} if e is None else dict(zip(OPS, map(float, e), strict=True)),
        "selection": selection(rows, flank, others),
    }
    known = true >= 0
    if known.any():
        w, g, t = want[known], got[known], true[known]
        kinds = {"error": (t == w) & (g != w), "biology": (t == g) & (g != w), "other": (t != g) & (t != w) & (g != w)}
        split = {k: counts(w[m], g[m]) for k, m in kinds.items()}
        open_ = np.tensordot(np.bincount(w[t == w], minlength=64).astype(float), CHANCE, 1)
        rate = _ratio(split["error"].sum(1), open_.sum(1))
        out["truth"] = {
            "codons": int(known.sum()),
            "labels": {k: x.astype(int).tolist() for k, x in split.items()},
            "per_base": _per_base(rate, open_),
            "ops": dict(zip(OPS, map(float, rate), strict=True)),
            "omega": omega(split["biology"], open_),
        }
    return out


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m sequencing_error_model.pfam_locus", description="Locus consensus of reads placed on Pfam."
    )
    p.add_argument("r1", type=Path, help="`pfam reads` output")
    p.add_argument("r2", type=Path)
    p.add_argument("--pairs", type=Path, required=True, help="these reads' search, pairs.tsv")
    p.add_argument("--prep", type=Path, required=True)
    p.add_argument("--sam", type=Path, help="the reads aligned to the sample's genomes, for truth")
    p.add_argument("--reference", type=Path, help="the genomes, FASTA")
    p.add_argument(
        "--mask-variants",
        action="store_true",
        help="truth as pfam_truth judges P2: MAPQ >= 20, and not at sites where these reads carry another base at >= 20 %%",
    )
    p.add_argument("--t1", type=float, default=35.0)
    p.add_argument("--t2", type=float, default=6.0)
    p.add_argument("--output", type=Path, required=True, metavar="JSON")
    args = p.parse_args(argv)
    if (args.sam is None) != (args.reference is None):
        p.error("--sam and --reference go together")
    profiles = Profiles.load(args.prep)
    rows = table(pfam.read_pairs(args.r1, args.r2), pfam.read_hits(args.pairs, profiles), profiles, args.t1, args.t2)
    if args.sam:
        with pysam.FastaFile(str(args.reference)) as f:
            genome = {c: f.fetch(c).upper() for c in f.references}
        if args.mask_variants:
            masks = pfam_truth.sam_variant_masks(args.sam, genome, 20, 0.2)
            add_truth(rows, args.sam, genome, 20, masks)
        else:
            add_truth(rows, args.sam, genome)
    found = [report(rows, k, others) for k in range(len(FLANKS)) for others in (2, 4)]
    args.output.write_text(json.dumps({"ops": OPS, "effects": EFFECTS, "exposed": len(rows), "loci": found}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
