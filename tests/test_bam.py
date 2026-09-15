from pathlib import Path

import numpy as np
import pysam

from sequencing_error_model import generate as gen
from sequencing_error_model import recovery
from sequencing_error_model import spec as spec_io
from sequencing_error_model.sources import bam


def _write(tmp_path: Path, n: int) -> tuple[Path, Path, list[str], list[gen.Read], list[int]]:
    """Generated reads written as a BAM with their true CIGARs, plus records the source must skip."""
    rng = np.random.default_rng(3)
    genome = "".join(rng.choice(list("ACGT"), size=3000))
    ref = tmp_path / "ref.fa"
    ref.write_text(f">chr\n{genome}\n")
    starts, reverse = rng.integers(0, 2950, size=n), rng.random(n) < 0.5
    mates = [int(m) for m in rng.integers(1, 3, size=n)]
    fwd = [genome[s : s + int(rng.integers(30, 45))] for s in starts]
    templates = [gen._revcomp(t) if r else t for t, r in zip(fwd, reverse, strict=True)]
    reads = gen.generate(recovery.example_spec(), templates, mates, rng)
    path = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr", "LN": len(genome)}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (start, rev, mate, read) in enumerate(zip(starts, reverse, mates, reads, strict=True)):
            seq, qual, cigar = read.sequence, [ord(c) - 33 for c in read.quality], gen._CIGAR.findall(read.cigar)
            if rev:
                seq, qual, cigar = gen._revcomp(seq), qual[::-1], cigar[::-1]
            else:  # a soft clip must not add rows
                seq, qual, cigar = "TTT" + seq, [30] * 3 + qual, [("3", "S"), *cigar]
            for flag, mapq in ((0, 60), (0x100, 60), (0, 5), (0x400, 60)):  # kept, secondary, low MAPQ, duplicate
                a = pysam.AlignedSegment(out.header)
                a.query_name, a.reference_id, a.reference_start, a.mapping_quality = f"r{i}", 0, int(start), mapq
                a.flag = flag | 0x1 | (0x10 if rev else 0) | (0x80 if mate == 2 else 0x40)
                a.cigarstring = "".join(n + o for n, o in cigar)
                a.query_sequence = seq
                a.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in qual))
                out.write(a)
    return path, ref, templates, reads, mates


def test_records_reproduce_generator_tuples(tmp_path: Path) -> None:
    path, ref, templates, reads, mates = _write(tmp_path, 300)
    flank, m = (2, 2), 1
    expected = gen.observations(zip(templates, reads, mates, strict=True), flank, m)
    got = gen.observations(bam.records(path, ref), flank, m)
    assert any("-" in str(k[-1]) for k in expected.counts)  # indels exercised
    assert got.counts == expected.counts


def test_cli_fits_spec(tmp_path: Path) -> None:
    path, ref, *_ = _write(tmp_path, 200)
    out = tmp_path / "spec"
    args = [str(path), str(ref), "--output", str(out), "--error-tokens", "QualityWindow(1)", "Context(1,1)"]
    assert bam.main([*args, "--quality-tokens", "QualityMarkov(1)", "Position(3)"]) == 0
    spec = spec_io.load(out)
    assert spec.provenance["mode"] == "reference" and spec.quality_alphabet == (2, 12, 23, 37)
