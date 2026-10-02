"""Badread error model and Q-score model (text, gzip-compatible), counted as `badread error_model` and
`badread qscore_model` count them, on the spec's sampled reads and their true alignments.

How Badread simulates (checked in its v0.4 source, phase 9): a read's target identity is drawn from
`--identity mean,max,stdev`; errors are added at random k-mers, each replaced by an alternative drawn from the
error model, until the read reaches that identity. Qualities are then drawn per read base from the Q-score
model, keyed by the local alignment CIGAR around it. So errors don't depend on Q (no `--q-policy`), the error
*rate* comes from `--identity` (exported from the sample's per-read identities, edlib-style: matches over
alignment columns), and error *placement* and type come from 7-mer context.

- Error model: one line per 7-mer, `KMER,p;ALT,p;...`: the k-mer itself, then up to 25 alternatives (read
  k-mers sharing its first and last base), by frequency.
- Q-score model: an `overall` line, then `CIGAR;count;q:p,...` for every odd window up to 9 read bases (runs of
  more than 6 deletions collapsed to 6) seen at least `min_occur` times; 1-base windows are always kept, as
  Badread requires `=`, `X` and `I`.

Dropped: read position, mate, neighbouring-Q and Q-dependence of errors, and Q dependence on sequence context.
Badread's own artefacts (junk/random reads, chimeras, glitches, adapters) are not modelled by the spec; the
exported arguments switch them off.
"""

import gzip
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from sequencing_error_model.export import Sample, fidelity, sample
from sequencing_error_model.fit.quality import Array
from sequencing_error_model.generate import _CIGAR
from sequencing_error_model.spec import ErrorModelSpec

K_ERROR, MAX_ALT = 7, 25
K_QSCORE, MAX_DEL, MIN_OCCUR, MAX_OUTPUT = 9, 6, 100, 10000


def aligned(template: str, sequence: str, quality: str, cigar: str) -> tuple[str, str, str, str]:
    """(aligned template, aligned read, aligned Q, column CIGAR over = X I D), gaps as '-' and Q gaps as ' '."""
    t = r = 0
    ref, read, qual, ops = [], [], [], []
    for count, op in _CIGAR.findall(cigar):
        for _ in range(int(count)):
            if op == "M":
                ref.append(template[t])
                read.append(sequence[r])
                qual.append(quality[r])
                ops.append("=" if template[t] == sequence[r] else "X")
                t, r = t + 1, r + 1
            elif op == "I":
                ref.append("-")
                read.append(sequence[r])
                qual.append(quality[r])
                ops.append("I")
                r += 1
            else:
                ref.append(template[t])
                read.append("-")
                qual.append(" ")
                ops.append("D")
                t += 1
    return "".join(ref), "".join(read), "".join(qual), "".join(ops)


def error_model(s: Sample) -> str:
    """Badread's `error_model` counting: each window of K_ERROR template bases against the read bases aligned
    inside it, kept when both share their first and last base."""
    alternatives: dict[str, Counter[str]] = defaultdict(Counter)
    for template, read in zip(s.templates, s.reads, strict=True):
        ref, seq, _, _ = aligned(template, read.sequence, read.quality, read.cigar)
        at = [i for i, b in enumerate(ref) if b != "-"]
        # ponytail: Python loop over every template base, ~1 s per Mbp; vectorise if exports get large.
        for i in range(len(at) - K_ERROR + 1):
            start, end = at[i], at[i + K_ERROR - 1] + 1
            kmer, alt = ref[start:end].replace("-", ""), seq[start:end].replace("-", "")
            if len(alt) > 1 and kmer[0] == alt[0] and kmer[-1] == alt[-1] and set(kmer + alt) <= set("ACGT"):
                alternatives[kmer][alt] += 1
    lines = []
    for kmer in sorted(alternatives):
        counts = alternatives[kmer]
        total = sum(counts.values())
        alts = sorted(((k, n) for k, n in counts.items() if k != kmer), key=lambda x: -x[1])[:MAX_ALT]
        lines.append("".join(f"{k},{n / total:.6f};" for k, n in [(kmer, counts[kmer]), *alts]))
    return "\n".join(lines) + "\n"


def _fractions(name: str, qs: Counter[int]) -> str:
    total = sum(qs.values())
    return f"{name};{total};" + "".join(f"{q}:{qs[q] / total:.6f}".rstrip("0").rstrip(".") + "," for q in sorted(qs))


def qscore_model(s: Sample) -> str:
    """Badread's `qscore_model` counting: the Q of the middle read base of every odd window of read bases (with
    the deletions between them), keyed by that window's CIGAR."""
    overall: Counter[int] = Counter()
    per_cigar: dict[str, Counter[int]] = defaultdict(Counter)
    for template, read in zip(s.templates, s.reads, strict=True):
        _, _, qual, ops = aligned(template, read.sequence, read.quality, read.cigar)
        emit = [i for i, o in enumerate(ops) if o != "D"]
        for k in range(1, K_QSCORE + 1, 2):
            for j in range(len(emit) - k + 1):
                start, end = emit[j], emit[j + k - 1] + 1
                cigar = ops[start:end]
                while "D" * (MAX_DEL + 1) in cigar:
                    cigar = cigar.replace("D" * (MAX_DEL + 1), "D" * MAX_DEL)
                q = ord(qual[emit[j + k // 2]]) - 33
                if k == 1:
                    overall[q] += 1
                per_cigar[cigar][q] += 1
    ranked = sorted(per_cigar, key=lambda c: -sum(per_cigar[c].values()))
    kept = [c for c in ranked if len(c) == 1 or sum(per_cigar[c].values()) >= MIN_OCCUR][:MAX_OUTPUT]
    kept += [c for c in ("=", "X", "I") if c in per_cigar and c not in kept]
    return "\n".join([_fractions("overall", overall), *(_fractions(c, per_cigar[c]) for c in kept)]) + "\n"


def identities(s: Sample) -> Array:
    out = []
    for template, read in zip(s.templates, s.reads, strict=True):
        ops = aligned(template, read.sequence, read.quality, read.cigar)[3]
        out.append(ops.count("=") / max(len(ops), 1))
    return np.array(out)


def export(
    model: ErrorModelSpec,
    out: Path,
    rng: np.random.Generator,
    *,
    read_length: int = 5000,
    reads: int = 1000,
    base_profile: list[Path] | None = None,
    **_: Any,
) -> dict[str, Any]:
    """Write `out/badread_error_model.gz` and `out/badread_qscore_model.gz`; returns the fidelity report."""
    s = sample(model, reads, read_length, rng, paired=False, base_profile=base_profile)
    for name, text in (("error", error_model(s)), ("qscore", qscore_model(s))):
        with gzip.open(out / f"badread_{name}_model.gz", "wt") as handle:
            handle.write(text)
    ident = 100 * identities(s)
    mean, sd = float(ident.mean()), float(ident.std())
    top = min(100.0, max(float(ident.max()), mean + 1e-3))
    args = [
        *("--error_model", str(out / "badread_error_model.gz"), "--qscore_model", str(out / "badread_qscore_model.gz")),
        *("--identity", f"{mean:.3f},{top:.3f},{max(sd, 1e-3):.3f}"),
        *("--junk_reads", "0", "--random_reads", "0", "--chimeras", "0", "--glitches", "0,0,0"),
        *("--start_adapter", "0,0", "--end_adapter", "0,0"),
    ]
    return fidelity(
        "badread",
        model,
        s,
        s,
        None,
        "any",
        [
            "read position and mate effects in both heads",
            "Q dependence of errors (Badread places errors by 7-mer, then draws Q from the local CIGAR)",
            "Q dependence on sequence context and Q Markov dependence beyond the 9-base CIGAR window",
            "error rate per read beyond a beta-distributed identity",
        ],
        args,
        base_profile,
        identity_mean=mean,
        identity_sd=sd,
    )
