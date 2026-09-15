import gzip
from pathlib import Path

import pytest

from sequencing_error_model.sources import fastq_quality as fq

# Q = 0 10 20 30, then 20 20
FASTQ = "@a\nACGT\n+\n!+5?\n\n@b\nGG\n+b\n55\n"


def test_profile(tmp_path: Path) -> None:
    path = tmp_path / "r1.fastq"
    path.write_text(FASTQ)
    p = fq.profile_fastq(path, order=1, flank=(1, 1))
    assert p.n_reads == 2
    assert p.lengths == {4: 1, 2: 1}
    assert p.alphabet == [0, 10, 20, 30]
    assert p.position == {(0, 0): 1, (1, 10): 1, (2, 20): 1, (3, 30): 1, (0, 20): 1, (1, 20): 1}
    assert p.transitions == {(0, 10): 1, (10, 20): 1, (20, 30): 1, (20, 20): 1}
    assert p.context == {(".AC", 0): 1, ("ACG", 10): 1, ("CGT", 20): 1, ("GT.", 30): 1, (".GG", 20): 1, ("GG.", 20): 1}


def test_order_zero_is_marginal(tmp_path: Path) -> None:
    path = tmp_path / "r1.fastq"
    path.write_text(FASTQ)
    assert fq.profile_fastq(path, order=0).transitions == {(0,): 1, (10,): 1, (20,): 3, (30,): 1}


def test_gzip_multiple_files_and_max_reads(tmp_path: Path) -> None:
    plain, gz = tmp_path / "r1.fastq", tmp_path / "r1.fq.gz"
    plain.write_text(FASTQ)
    gz.write_bytes(gzip.compress(FASTQ.encode()))
    assert fq.profile_fastq(gz) == fq.profile_fastq(plain)
    assert fq.profile_fastq(plain, gz).n_reads == 4
    assert fq.profile_fastq(plain, gz, max_reads=3).lengths == {4: 2, 2: 1}


@pytest.mark.parametrize(
    ("text", "match"),
    [
        ("@a\nACGT\n+\n!!!\n", ":1: malformed"),
        ("@a\nACGT\n+\n!!!!\n@b\nAC\n", ":5: malformed"),
        ("@a\nAC\n+\n! \n", ":4: quality outside"),
    ],
)
def test_malformed(tmp_path: Path, text: str, match: str) -> None:
    path = tmp_path / "bad.fastq"
    path.write_text(text)
    with pytest.raises(fq.FastqFormatError, match=match):
        fq.profile_fastq(path)
