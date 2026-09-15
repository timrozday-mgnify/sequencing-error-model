import re
from collections import Counter
from dataclasses import replace
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


def test_fragments_place_mates_on_both_strands() -> None:
    rng = np.random.default_rng(3)
    contigs = {"chr": "".join(rng.choice(list("ACGT"), size=500)), "tiny": "acgt"}
    pairs = gen.fragments(list(contigs.items()), gen.insert_sizes(40, 15), 3000, 30, rng)
    sizes = []
    for name, r1, r2 in pairs:
        match = re.fullmatch(r"(\w+):(\d+)-(\d+)#\d+", name)
        assert match
        frag = contigs[match[1]][int(match[2]) - 1 : int(match[3])].upper()
        assert r1 == (frag + gen.ADAPTERS[0] + "N" * 30)[:30]
        assert r2 == (frag[::-1].translate(str.maketrans("ACGT", "TGCA")) + gen.ADAPTERS[1] + "N" * 30)[:30]
        sizes.append(len(frag))
    assert abs(np.mean(sizes) - 40) < 1 and abs(np.std(sizes) - 15) < 1
    assert any(s < 30 for s in sizes)  # read-through into the adapter
    try:
        gen.fragments([("x", "ACG")], gen.insert_sizes(10, 0), 1, 5, rng)
        raise AssertionError("an insert longer than every contig must fail")
    except ValueError:
        pass


def test_cli_samples_standalone_pairs(tmp_path: Path) -> None:
    model = recovery.example_spec()
    replace(model, marginals={"insert_size": gen.insert_sizes(50, 10)}).save(tmp_path / "spec")
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">chr desc\n" + "".join(np.random.default_rng(0).choice(list("ACGT"), size=300)) + "\n")
    out = tmp_path / "out.fastq"
    args = ["--model", str(tmp_path / "spec"), "--input", str(fasta), "--output", str(out), "--seed", "1"]
    assert gen.main([*args, "--pairs", "5", "--read-length", "35", "--batch", "3"]) == 0
    lines = out.read_text().splitlines()
    names = [lines[i].split()[0] for i in range(0, len(lines), 4)]
    assert len(lines) == 40 and len(set(names)) == 10
    assert all(re.fullmatch(rf"@chr:\d+-\d+#{i // 2}/{i % 2 + 1}", n) for i, n in enumerate(names))


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
