"""InSilicoSeq `kde` error model (.npz), as `iss.error_models.kde.KDErrorModel` loads it (InSilicoSeq 2.x).

How InSilicoSeq simulates (checked in its source, phase 9): per read it picks a mean-quality bin (4 bins of
the read's mean Q: 0-9, 10-19, 20-29, 30-39) by `mean_count_*`, draws each position's Q from that bin's CDF
over Q 0..40, then **substitutes each base with probability 10^(-Q/10)**, choosing the new base from that
position's `subst_choices_*`. Insertions after and deletions of each base come from per-position rates
(`ins_*` by inserted base, `del_*` by the base). Forward = R1, reverse = R2.

So substitutions are Q-coupled (`--q-policy` applies, by read position), and InSilicoSeq's tables hold
per-position Q given the read's mean-Q bin, the substitution target by position and base, and indel rates by
position and base. Everything else is dropped. Deletion rates are conditional on the deleted base (InSilicoSeq's
own modeller divides by all bases at the position, which undercounts by about 4x; we export what its simulator
reads). Its fixed-length reads lose read-length variation from indels: it pads or trims to `read_length`.

The npz holds pickled Python objects because that is InSilicoSeq's format; this module only writes it.
"""

from pathlib import Path
from typing import Any

import numpy as np

from sequencing_error_model.export import _D, _I, _M, Sample, columns, fidelity, remap, sample
from sequencing_error_model.fit.quality import Array
from sequencing_error_model.generate import insert_sizes
from sequencing_error_model.spec import ErrorModelSpec

Q_MAX = 40  # InSilicoSeq's CDFs span Q 0..40
_BASES = "ATCG"  # InSilicoSeq's dict order
_KEYS = ("mean_count", "quality_hist", "subst_choices", "ins", "del")


def _cdfs(q: Array, pos: Array, read_length: int) -> Array:
    h = np.zeros((read_length, Q_MAX + 1))
    np.add.at(h, (pos, q.clip(max=Q_MAX)), 1)
    empty = h.sum(axis=1) == 0
    h[empty] = h.sum(axis=0) + (not h.any())  # unseen positions: the bin's pooled histogram
    return np.cumsum(h, axis=1) / h.sum(axis=1, keepdims=True)


def tables(s: Sample, read_length: int, mate: int) -> dict[str, Any]:
    """The per-mate arrays of the npz from one sample."""
    c = columns(s)
    sel = (c["mate"] == mate) & (c["pos"] < read_length)
    op, pos, q, ref, alt = (c[k][sel] for k in ("op", "pos", "q", "ref", "alt"))
    emitted = op != _D
    reads = c["read"][sel][emitted]
    rid, inverse = np.unique(reads, return_inverse=True)
    mean_q = np.bincount(inverse, q[emitted]) / np.bincount(inverse)
    bins = np.minimum(mean_q.astype(np.int64) // 10, 3)
    read_bin = bins[inverse]
    hist: list[Any] = []
    for bin_ in range(4):
        rows = read_bin == bin_
        hist.append(list(_cdfs(q[emitted][rows], pos[emitted][rows], read_length)) if rows.any() else [])
    quality = np.empty(4, dtype=object)
    quality[:] = hist

    at = np.zeros((read_length, 256))  # template bases at each position (M and D columns)
    np.add.at(at, (pos[op != _I], ref[op != _I]), 1)
    subs = np.zeros((read_length, 256, 256))
    m = (op == _M) & (ref != alt)
    np.add.at(subs, (pos[m], ref[m], alt[m]), 1)
    ins = np.zeros((read_length, 256))
    i = op == _I
    np.add.at(ins, ((pos[i] - 1).clip(min=0), alt[i]), 1)  # inserted after the previous read base
    dels = np.zeros((read_length, 256))
    np.add.at(dels, (pos[op == _D], ref[op == _D]), 1)
    bases_at = at[:, np.frombuffer(_BASES.encode(), np.uint8)].sum(axis=1)
    choices, ins_rates, del_rates = [], [], []
    for p in range(read_length):
        choice = {}
        for b in _BASES:
            alts = [x for x in _BASES if x != b]
            counts = subs[p, ord(b), np.frombuffer("".join(alts).encode(), np.uint8)]
            choice[b] = (alts, list(counts / counts.sum()) if counts.sum() else [1 / 3] * 3)
        choices.append(choice)
        ins_rates.append({b: float(ins[p, ord(b)] / bases_at[p]) if bases_at[p] else 0.0 for b in _BASES})
        del_rates.append({b: float(dels[p, ord(b)] / at[p, ord(b)]) if at[p, ord(b)] else 0.0 for b in _BASES})
    return {
        "mean_count": [int((bins == i).sum()) for i in range(4)],
        "quality_hist": quality,
        "subst_choices": choices,
        "ins": ins_rates,
        "del": del_rates,
        "reads": len(rid),
    }


def export(
    model: ErrorModelSpec,
    out: Path,
    rng: np.random.Generator,
    *,
    q_policy: str = "preserve-quality",
    read_length: int = 150,
    reads: int = 20000,
    base_profile: list[Path] | None = None,
    insert_mean: float | None = None,
    insert_sd: float = 0.0,
    **_: Any,
) -> dict[str, Any]:
    """Write `out/iss.npz` (use `iss generate --model kde --model_file out/iss.npz`); returns the fidelity report."""
    if insert_mean is not None:
        sizes = insert_sizes(insert_mean, insert_sd)
    elif "insert_size" in model.marginals:
        sizes = model.marginals["insert_size"]
    else:
        raise ValueError("InSilicoSeq needs a fragment size: pass insert_mean, or use a spec with insert_size")
    inner = np.asarray(sizes[0], np.int64) - 2 * read_length
    pmf = np.bincount(inner.clip(min=0), weights=sizes[1])
    spec_sample = sample(model, reads, read_length, rng, base_profile=base_profile)
    emitted = remap(spec_sample, "substitution", True, Q_MAX) if q_policy == "preserve-errors" else spec_sample
    r1, r2 = tables(emitted, read_length, 1), tables(emitted, read_length, 2)
    np.savez_compressed(
        out / "iss.npz",
        model="kde",
        read_length=read_length,
        insert_size=np.cumsum(pmf) / pmf.sum(),
        **{f"{k}_{d}": v[k] for d, v in (("forward", r1), ("reverse", r2)) for k in _KEYS},
    )
    c = columns(spec_sample)
    clipped = float(np.mean(c["q"][c["op"] != _D] > Q_MAX))
    return fidelity(
        "insilicoseq",
        model,
        spec_sample,
        emitted,
        q_policy,
        "substitution",
        [
            "base context (substitution target depends on position and base only)",
            "neighbouring-Q and Q x context effects in head E",
            "Q Markov dependence (positions independent given the read's mean-Q bin)",
            "indel length (InSilicoSeq indels are 1 bp) and indel dependence on Q",
        ],
        ["--model", "kde", "--model_file", str(out / "iss.npz")],
        base_profile,
        read_length,
        q_above_40_fraction=clipped,
    )
