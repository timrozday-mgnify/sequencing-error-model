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
    assert p.position == {(1, 0): 1, (2, 10): 1, (3, 20): 1, (4, 30): 1, (1, 20): 1, (2, 20): 1}
    assert p.transitions == {(0, 10): 1, (10, 20): 1, (20, 30): 1, (20, 20): 1}
    assert p.context == {(".AC", 0): 1, ("ACG", 10): 1, ("CGT", 20): 1, ("GT.", 30): 1, (".GG", 20): 1, ("GG.", 20): 1}


def test_order_zero_is_marginal(tmp_path: Path) -> None:
    path = tmp_path / "r1.fastq"
    path.write_text(FASTQ)
    assert fq.profile_fastq(path, order=0).transitions == {(0,): 1, (10,): 1, (20,): 3, (30,): 1}


def test_tables(tmp_path: Path) -> None:
    path = tmp_path / "r2.fastq"
    path.write_text(FASTQ)
    t = {x.source: x for x in fq.tables(fq.profile_fastq(path, order=2), mate=2)}
    assert all(not x.truth and "op" not in x.fields and x.fields[-1] == "mate" for x in t.values())
    assert t["fastq_quality:transitions"].fields == ("q-2", "q-1", "q", "mate")
    q_hist = {(0, 2): 1, (10, 2): 1, (20, 2): 3, (30, 2): 1}
    assert t["fastq_quality:position"].marginal("q", "mate").counts == q_hist
    assert t["fastq_quality:context"].marginal("q", "mate").counts == q_hist
    assert t["fastq_quality:length"].counts == {(4, 2): 1, (2, 2): 1}


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


def test_every_n_strides_within_each_file(tmp_path: Path) -> None:
    # Six reads of distinct lengths, so the profile says which were kept.
    text = "".join(f"@r{i}\n{'A' * (i + 1)}\n+\n{'I' * (i + 1)}\n" for i in range(6))
    a, b = tmp_path / "a.fastq", tmp_path / "b.fastq"
    a.write_text(text)
    b.write_text(text)
    assert fq.profile_fastq(a, every_n=2).lengths == {1: 1, 3: 1, 5: 1}
    # Striding is per file, so a second mate profiled the same way keeps the same records.
    assert fq.profile_fastq(a, b, every_n=2).lengths == {1: 2, 3: 2, 5: 2}
    assert fq.quality_table(a, m=1, every_n=2).marginal().counts == {(): 1 + 3 + 5}


def test_write_subset_round_trips(tmp_path: Path) -> None:
    text = "".join(f"@r{i} some comment\nACGT\n+\nIIII\n" for i in range(10))
    src, out = tmp_path / "r1.fastq", tmp_path / "sub.fastq.gz"
    src.write_text(text)
    assert fq.write_subset(src, out, every_n=4) == 3
    kept = gzip.decompress(out.read_bytes()).decode().splitlines()[::4]
    assert kept == ["@r0 some comment", "@r4 some comment", "@r8 some comment"]
    # The written subset is what the fitters would have strided to in process.
    assert fq.profile_fastq(out) == fq.profile_fastq(src, every_n=4)


def test_write_subset_rejects_truncated(tmp_path: Path) -> None:
    src = tmp_path / "bad.fastq"
    src.write_text("@a\nACGT\n+\nIIII\n@b\nACGT\n")
    with pytest.raises(fq.FastqFormatError, match="truncated record after 1 kept"):
        fq.write_subset(src, tmp_path / "out.fastq", every_n=1)
