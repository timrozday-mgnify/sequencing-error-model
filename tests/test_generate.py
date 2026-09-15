import re
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from sequencing_error_model import generate as gen
from sequencing_error_model import recovery


def templates(n: int, seed: int) -> tuple[list[str], list[int]]:
    rng = np.random.default_rng(seed)
    seqs = ["".join(rng.choice(list("ACGT"), size=rng.integers(30, 41))) for _ in range(n)]
    return seqs, [int(m) for m in rng.integers(1, 3, size=n)]


def test_reads_are_alignment_consistent() -> None:
    model = recovery.example_spec()
    seqs, mates = templates(300, seed=1)
    seqs[0] = "ACGTNNNNACGT"
    reads = gen.generate(model, seqs, mates, np.random.default_rng(2), max_ins_run=3, error_rate_scale=8)
    ops: Counter[str] = Counter()
    for template, read in zip(seqs, reads, strict=True):
        n = Counter[str]()
        for count, op in re.findall(r"(\d+)([MID])", read.cigar):
            n[op] += int(count)
        ops += n
        assert n["M"] + n["D"] == len(template) == len(read.q_track)
        assert n["M"] + n["I"] == len(read.sequence) == len(read.quality)
        assert {ord(c) - 33 for c in read.quality} <= set(model.quality_alphabet)
        _, tq = gen.align(template, read)
        assert all(q == read.q_track[t] for t, q in enumerate(tq) if q is not None)
        assert sum(q is None for q in tq) == n["D"]
    assert "NNNN" in reads[0].sequence  # non-ACGT bases pass through as matches, with no insertions between
    assert ops["I"] and ops["D"], ops


def test_cli_matches_contract(tmp_path: Path) -> None:
    recovery.example_spec().save(tmp_path / "spec")
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">0/1\nACGTACGTACGTACGTACGT\n>1/2\nGGGCCCAAATTTGGGCCC\n>2/1\n\n")
    outputs: list[Any] = []
    for i in range(2):
        out = tmp_path / f"out{i}.fastq"
        args = ["--model", str(tmp_path / "spec"), "--input", str(fasta), "--output", str(out), "--paired"]
        assert gen.main([*args, "--seed", "7"]) == 0
        outputs.append(out.read_text())
    assert outputs[0] == outputs[1]
    lines = outputs[0].splitlines()
    assert len(lines) == 12 and [lines[i].split()[0] for i in (0, 4, 8)] == ["@0/1", "@1/2", "@2/1"]
    assert all(re.fullmatch(r"@\S+ cigar:(\d+[MID])*", lines[i]) for i in (0, 4, 8))
