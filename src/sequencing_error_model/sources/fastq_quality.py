"""Quality-process statistics from raw FASTQ(.gz), for head Q and the default-mode exposure.

A feature/output source (plan §5.2, §6.1). Reported qualities are counted as what the
generator must reproduce and as the exposure for marginal matching (§5.6). Nothing here
is error evidence, and no count may be turned into an error rate. Bases are *observed*
bases, so Q | context is conditioned on reads, not on the true template (§6.2).

Profile each mate separately; pass all files (lanes) for one mate in a single call.
"""

import gzip
from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass, field
from itertools import chain, islice
from pathlib import Path

PHRED_OFFSET = 33
MAX_Q = 93
PAD = "."  # context positions beyond either read end
_DECODE = bytes((c - PHRED_OFFSET) % 256 for c in range(256))


class FastqFormatError(ValueError):
    """A FASTQ record is truncated or malformed."""


@dataclass
class QualityProfile:
    order: int
    flank: tuple[int, int]  # (L, R): observed bases before and after the centre base
    n_reads: int = 0
    lengths: Counter[int] = field(default_factory=Counter)
    # (0-based position from read start, Q)
    position: Counter[tuple[int, int]] = field(default_factory=Counter)
    # (q_{t-order}, ..., q_{t-1}, q_t) for t >= order; the first positions come from `position`
    transitions: Counter[tuple[int, ...]] = field(default_factory=Counter)
    # (observed x_{t-L..t+R}, q_t): Q | observed context, and the centre-Q x context exposure
    context: Counter[tuple[str, int]] = field(default_factory=Counter)

    @property
    def alphabet(self) -> list[int]:
        return sorted({q for _, q in self.position})


def read_fastq(path: Path) -> Iterator[tuple[str, bytes]]:
    """Yield (bases, Phred scores) per record; iterating the scores gives ints."""
    with path.open("rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    with gzip.open(path, "rb") if gz else path.open("rb") as fh:
        line = 0
        while header := fh.readline():
            line += 1
            if not header.strip():
                continue
            seq, plus, qual = (fh.readline().rstrip(b"\r\n") for _ in range(3))
            if not header.startswith(b"@") or not plus.startswith(b"+") or len(seq) != len(qual):
                raise FastqFormatError(f"{path}:{line}: malformed or truncated record")
            if qual and not PHRED_OFFSET <= min(qual) <= max(qual) <= PHRED_OFFSET + MAX_Q:
                raise FastqFormatError(f"{path}:{line + 3}: quality outside Phred+33 0-{MAX_Q}")
            line += 3
            yield seq.decode("ascii"), qual.translate(_DECODE)


def profile_fastq(
    *paths: Path, order: int = 1, flank: tuple[int, int] = (2, 2), max_reads: int | None = None
) -> QualityProfile:
    """Count quality-process statistics over the first `max_reads` reads of `paths`, in order."""
    left, right = flank
    if order < 0 or left < 0 or right < 0:
        raise ValueError(f"order and flank must be non-negative, got order={order}, flank={flank}")
    p = QualityProfile(order, flank)
    for seq, q in islice(chain.from_iterable(map(read_fastq, paths)), max_reads):
        p.n_reads += 1
        p.lengths[len(q)] += 1
        p.position.update(enumerate(q))
        p.transitions.update(zip(*(q[i:] for i in range(order + 1)), strict=False))
        padded = PAD * left + seq + PAD * right
        p.context.update((padded[t : t + left + right + 1], qt) for t, qt in enumerate(q))
    return p
