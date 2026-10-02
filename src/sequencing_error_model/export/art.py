"""ART / art_modern Illumina quality profiles (text, one per mate), plus indel-rate arguments.

Format (as ART's `art_profiler_illumina` writes it, and art_modern ships): for each symbol `.` (all bases), then
A, T, G, C and N (separate profiles per base, used with `-sp`), and each 0-based read position, two
tab-separated lines: `<symbol> <pos> <Q values...>` and `<symbol> <pos> <cumulative counts...>`.

ART draws each position's Q independently from that profile and **substitutes the base with probability
10^(-Q/10)**, so `--q-policy` applies, by read position. Indels come from per-mate rates (`-ir/-dr/-ir2/-dr2`), unrelated to Q or context.

Per-base blocks are keyed by the true base. They are written, but the exported arguments leave `-sp` off: ART 2.3.7 keys them by the
forward-strand reference base (measured in phase 9: reverse-strand reads drew G's profile at C), while head Q
conditions on the base in read orientation. Without `-sp` ART reproduces the `.` profile exactly.
"""

from pathlib import Path
from typing import Any

import numpy as np

from sequencing_error_model.export import _D, _I, _M, Sample, columns, fidelity, remap, sample
from sequencing_error_model.spec import ErrorModelSpec

Q_MAX = 41
_SYMBOLS = ".ATGCN"


def profile(s: Sample, read_length: int, mate: int) -> str:
    c = columns(s)
    sel = (c["mate"] == mate) & (c["pos"] < read_length) & (c["op"] != _D)
    pos, q, ref = c["pos"][sel], c["q"][sel].clip(max=Q_MAX), c["ref"][sel]
    lines = []
    for symbol in _SYMBOLS:
        keep = np.ones(len(pos), bool) if symbol in ".N" else ref == ord(symbol)
        h = np.zeros((read_length, Q_MAX + 1), np.int64)
        np.add.at(h, (pos[keep], q[keep]), 1)
        h[h.sum(axis=1) == 0] = h.sum(axis=0) + (not h.any())  # unseen positions: the pooled histogram
        for p in range(read_length):
            values = np.flatnonzero(h[p])
            lines.append("\t".join(map(str, (symbol, p, *values))))
            lines.append("\t".join(map(str, (symbol, p, *np.cumsum(h[p, values])))))
    return "\n".join(lines) + "\n"


def indel_rates(s: Sample, mate: int) -> tuple[float, float]:
    """(insertion, deletion) events per template base for one mate."""
    c = columns(s)
    sel = c["mate"] == mate
    op, read = c["op"][sel], c["read"][sel]
    start = np.r_[True, (op[1:] != op[:-1]) | (read[1:] != read[:-1])]  # first column of each run
    bases = max(int((op != _I).sum()), 1)
    return float(((op == _I) & start).sum() / bases), float(((op == _D) & start).sum() / bases)


def export(
    model: ErrorModelSpec,
    out: Path,
    rng: np.random.Generator,
    *,
    q_policy: str = "preserve-quality",
    read_length: int = 150,
    reads: int = 20000,
    base_profile: list[Path] | None = None,
    **_: Any,
) -> dict[str, Any]:
    """Write `out/art_R1.txt` and `out/art_R2.txt`; returns the fidelity report, whose `arguments` are the ART flags."""
    spec_sample = sample(model, reads, read_length, rng, base_profile=base_profile)
    emitted = remap(spec_sample, "substitution", True, Q_MAX) if q_policy == "preserve-errors" else spec_sample
    args = ["-l", str(read_length)]
    for mate, (flag_i, flag_d) in ((1, ("-ir", "-dr")), (2, ("-ir2", "-dr2"))):
        path = out / f"art_R{mate}.txt"
        path.write_text(profile(emitted, read_length, mate))
        ins, dels = indel_rates(spec_sample, mate)
        args += [f"-{mate}", str(path), flag_i, f"{ins:.3g}", flag_d, f"{dels:.3g}"]
    c = columns(spec_sample)
    return fidelity(
        "art",
        model,
        spec_sample,
        emitted,
        q_policy,
        "substitution",
        [
            "base context (substitutions are uniform over the other bases)",
            "neighbouring-Q and Q x context effects in head E",
            "Q Markov dependence and read-level Q heterogeneity (positions independent)",
            "Q dependence on the base and its context",
            "indel position, context, Q and length dependence (flat per-mate rates)",
        ],
        args,
        base_profile,
        read_length,
        q_above_max_fraction=float(np.mean(c["q"][c["op"] == _M] > Q_MAX)),
    )
