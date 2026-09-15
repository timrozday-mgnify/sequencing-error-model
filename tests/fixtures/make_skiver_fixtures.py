"""Regenerate the skiver analyze fixtures from a seeded synthetic genome and reads.

Usage: python tests/fixtures/make_skiver_fixtures.py <skiver-binary> <out-dir>

Writes <out-dir>/reads.fastq (not committed) and runs `skiver analyze` on it, so
<out-dir>/analyze.*.csv can be diffed against the committed fixtures. Stdlib only.
"""

import random
import subprocess
import sys
from pathlib import Path

SEED = 7
GENOME_LEN = 4000
READ_LEN = 150
N_READS = 1600  # ~60x coverage
QUALS = (2, 12, 23, 37)  # binned NovaSeq alphabet
# Injected per-base substitution probability by quality bin. Deliberately not
# 10^(-Q/10): the fixtures must never agree with a Q-derived error rate by accident.
SUB_P = {2: 0.05, 12: 0.02, 23: 0.008, 37: 0.002}
INDEL_P = 0.0005


def revcomp(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def make_reads(rng: random.Random) -> str:
    genome = "".join(rng.choice("ACGT") for _ in range(GENOME_LEN))
    records: list[str] = []
    for i in range(N_READS):
        start = rng.randrange(GENOME_LEN - READ_LEN - 5)
        template = genome[start : start + READ_LEN + 5]
        if rng.random() < 0.5:
            template = revcomp(template)
        seq: list[str] = []
        qual: list[int] = []
        t = 0
        while len(seq) < READ_LEN:
            base = template[t]
            # quality decays along the read; a G in the previous template base lowers it
            p_low = 0.05 + 0.25 * len(seq) / READ_LEN + (0.1 if t and template[t - 1] == "G" else 0.0)
            q = rng.choices(QUALS, weights=(p_low / 2, p_low / 2, 0.3, 1.0 - p_low - 0.3))[0]
            r = rng.random()
            if r < INDEL_P:  # deletion: skip template base
                t += 1
                continue
            if r < 2 * INDEL_P:  # insertion before this base
                seq.append(rng.choice("ACGT"))
                qual.append(q)
                continue
            if rng.random() < SUB_P[q]:
                base = rng.choice([b for b in "ACGT" if b != base])
            seq.append(base)
            qual.append(q)
            t += 1
        records.append(f"@r{i}\n{''.join(seq)}\n+\n{''.join(chr(x + 33) for x in qual)}\n")
    return "".join(records)


def main() -> None:
    skiver, out = sys.argv[1], Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    reads = out / "reads.fastq"
    reads.write_text(make_reads(random.Random(SEED)))
    cmd = [skiver, "analyze", str(reads), "-k", "11", "-v", "13", "-c", "8", "-t", "1", "-o", str(out / "analyze")]
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
