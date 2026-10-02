"""PBSIM3 QSHMM (Q-bearing) and ERRHMM (error-only) model files.

Format (PBSIM3's `set_qshmm`/`set_errhmm`): space-separated lines `<accuracy> IP <state> <p>`,
`<accuracy> EP <state> <p...>` and `<accuracy> TP <state> <p...>`, states 1-based, one HMM per accuracy level.

How PBSIM3 simulates (checked in its source, phase 9), per read, after drawing an accuracy level from a fixed
distribution over [0.75, 1.05] x `--accuracy-mean`:

- QSHMM: a Q per read base from the level's HMM (EP over Q 0..; IP, EP and TP quantised to 1/100); then an
  error with probability 10^(-Q/10), split by `--difference-ratio` (sub:ins:del). Q-coupled, so `--q-policy`
  applies (positions pooled: PBSIM3 has no read position).
- ERRHMM: one HMM step per alignment column; EP gives P(match, sub, ins, del) (quantised to 1/1000), so every
  column's op comes from the hidden state. Qualities are `!`. Level 100 is error-free.

Both are written as one HMM repeated at every accuracy level 1..100, so `--accuracy-mean` does not change the
model (keep it <= 0.95 for ERRHMM, where level 100 means no errors):

- QSHMM: states are Q values (adjacent values merged into 50 states when the alphabet is larger, PBSIM3's
  `STATE_MAX`), each emitting its own Q values; IP and TP are the first-base Q and the lag-1 Q transitions of
  the sampled reads. This is the order-1 Q Markov marginal of head Q.
- ERRHMM: four states, one per op (match, sub, ins, del), each emitting its op; IP and TP are the first
  column's op and the lag-1 op transitions, so error clustering is kept at order 1.

Dropped: base context (PBSIM3 picks substitution and inserted bases uniformly), read position, Q dependence
on sequence context, read-level heterogeneity beyond the HMM, and (QSHMM) every effect on errors except Q.
"""

from pathlib import Path
from typing import Any

import numpy as np

from sequencing_error_model.export import _D, _I, _M, Sample, columns, fidelity, remap, sample
from sequencing_error_model.fit.quality import Array
from sequencing_error_model.spec import ErrorModelSpec

STATE_MAX, LEVELS = 50, range(1, 101)


def _chain(states: list[Array], n: int) -> tuple[Array, Array]:
    """(initial, transition) probabilities of an order-1 chain over `n` states from per-read state tracks;
    unseen rows get the marginal."""
    init, trans = np.zeros(n), np.zeros((n, n))
    for s in states:
        if len(s):
            init[s[0]] += 1
            np.add.at(trans, (s[:-1], s[1:]), 1)
    marginal = trans.sum(axis=0) + init
    trans[trans.sum(axis=1) == 0] = marginal
    return init / init.sum(), trans / trans.sum(axis=1, keepdims=True)


def _lines(ip: Array, ep: Array, tp: Array) -> list[str]:
    out = []
    for level in LEVELS:
        out += [f"{level} IP {i + 1} {p:.6e}" for i, p in enumerate(ip)]
        out += [f"{level} EP {i + 1} " + " ".join(f"{p:.6e}" for p in row) for i, row in enumerate(ep)]
        out += [f"{level} TP {i + 1} " + " ".join(f"{p:.6e}" for p in row) for i, row in enumerate(tp)]
    return out


def qshmm(s: Sample) -> str:
    tracks = [np.frombuffer(r.quality.encode(), np.uint8).astype(np.int64) - 33 for r in s.reads]
    values = np.unique(np.concatenate(tracks))
    groups = np.array_split(np.arange(len(values)), min(len(values), STATE_MAX))
    state_of = np.zeros(94, np.int64)
    for g, members in enumerate(groups):
        state_of[values[members]] = g
    ip, tp = _chain([state_of[t] for t in tracks], len(groups))
    counts = np.bincount(np.concatenate(tracks), minlength=94).astype(float)
    ep = np.zeros((len(groups), int(values.max()) + 1))
    for g, members in enumerate(groups):
        q = values[members]
        ep[g, q] = counts[q] / counts[q].sum()
    return "\n".join(_lines(ip, ep, tp)) + "\n"


def _op_tracks(s: Sample) -> list[Array]:
    c = columns(s)
    code = np.where(c["op"] == _I, 2, np.where(c["op"] == _D, 3, (c["ref"] != c["alt"]).astype(np.int64)))
    return np.split(code, np.flatnonzero(np.diff(c["read"])) + 1)


def errhmm(s: Sample) -> str:
    ip, tp = _chain(_op_tracks(s), 4)
    return "\n".join(_lines(ip, np.eye(4), tp)) + "\n"


def difference_ratio(s: Sample) -> str:
    """PBSIM3's `--difference-ratio` (sub:ins:del, integers 0-1000) from the sample's error columns."""
    c = columns(s)
    counts = np.array(
        [((c["op"] == _M) & (c["ref"] != c["alt"])).sum(), (c["op"] == _I).sum(), (c["op"] == _D).sum()], float
    )
    return ":".join(str(int(round(x))) for x in 1000 * counts / max(counts.sum(), 1))


def export(
    model: ErrorModelSpec,
    out: Path,
    rng: np.random.Generator,
    *,
    q_policy: str = "preserve-quality",
    read_length: int = 5000,
    reads: int = 1000,
    base_profile: list[Path] | None = None,
    **_: Any,
) -> dict[str, Any]:
    """Write `out/pbsim3_qshmm.model` and `out/pbsim3_errhmm.model`; the report's `arguments` are for QSHMM and
    `errhmm_arguments` for ERRHMM."""
    s = sample(model, reads, read_length, rng, paired=False, base_profile=base_profile)
    emitted = remap(s, "any", False, 93) if q_policy == "preserve-errors" else s
    (out / "pbsim3_qshmm.model").write_text(qshmm(emitted))
    (out / "pbsim3_errhmm.model").write_text(errhmm(s))
    common = ["--accuracy-mean", "0.9", "--length-mean", str(read_length)]
    return fidelity(
        "pbsim3",
        model,
        s,
        emitted,
        q_policy,
        "any",
        [
            "base context (substitution and inserted bases are uniform)",
            "read position and mate effects in both heads",
            "Q dependence on sequence context; QSHMM: Q Markov dependence beyond order 1, probabilities "
            "quantised to 1/100",
            "QSHMM: every effect on errors except Q; ERRHMM: Q dependence of errors, and qualities ('!')",
            "indel length beyond the order-1 op chain",
        ],
        [
            *("--method", "qshmm", "--qshmm", str(out / "pbsim3_qshmm.model")),
            *common,
            *("--difference-ratio", difference_ratio(s)),
        ],
        base_profile,
        errhmm_arguments=["--method", "errhmm", "--errhmm", str(out / "pbsim3_errhmm.model"), *common],
    )
