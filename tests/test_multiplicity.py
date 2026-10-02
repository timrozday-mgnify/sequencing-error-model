"""Key multiplicity of an assembly, and skiver's kvmer.csv joined to it (phase 6b)."""

import json
import random
import sys
from pathlib import Path

import numpy as np

from sequencing_error_model import multiplicity as mult
from sequencing_error_model.sources import skiver_analyze

sys.path.insert(0, str(Path(__file__).parent / "fixtures"))
from make_skiver_fixtures import GENOME_LEN, SEED  # noqa: E402

FIXTURES = Path(__file__).parent / "fixtures" / "skiver-v0.3.2"


def _random(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_one_edit_finds_a_substitution_deletion_or_insertion_and_nothing_further() -> None:
    b = "ACGTTGCAAGTCA"
    cases = {
        "ACGTAGCAAGTCA": True,  # substitution
        "ACGTGCAAGTCAG": True,  # b's base 4 deleted, next base shifted in
        "ACGTCTGCAAGTC": True,  # extra base at 4
        "TCGTTGCAAGTCT": False,  # two substitutions
        "GTTGCAAGTCAAC": False,  # two deletions
    }
    got = mult.one_edit(np.array([mult.encode(a) for a in cases]), np.full(len(cases), mult.encode(b)), len(b))
    assert got.tolist() == list(cases.values())


def test_a_repeat_makes_multi_locus_keys_and_a_unique_genome_has_none() -> None:
    v = 13
    unique = [("u", _random(4000, 1), 1.0)]
    assert mult.summary(*mult.occurrences(unique, 21, v), v)["multi_locus"] == 0.0

    x = _random(40, 2)  # a 40 bp repeat whose two copies are followed by different sequence
    repeat = [("r", _random(500, 3) + x + _random(500, 4) + x + _random(500, 5), 1.0)]
    s = mult.summary(*mult.occurrences(repeat, 21, v), v)
    # 40 - 21 + 1 = 20 key positions per copy sit inside the repeat, of 1580 - 21 - 13 + 1 per strand
    assert s["multi_locus"] == 40 / 1547
    assert 0 < s["off_majority"] < s["multi_value"] <= s["multi_locus"]
    # a weight on the contig moves no share
    heavy = [(n, seq, 30.0) for n, seq, _ in repeat]
    assert mult.summary(*mult.occurrences(heavy, 21, v), v) == s


def test_kvmer_fixture_joined_to_its_own_genome_is_single_locus() -> None:
    genome = _random(GENOME_LEN, SEED)  # make_skiver_fixtures draws the genome first from Random(SEED)
    a = skiver_analyze.read_analyze(FIXTURES / "analyze")
    key, val, _ = mult.occurrences([("g", genome, 1.0)], a.k, a.v)
    out = mult.by_loci(a.kvmer, key, val)
    one = out["classes"]["1"]
    assert one["keys"] >= 0.95 * out["keys"]
    assert one["consensus_true_share"] > 0.99
    assert out["filter_removes_single_edit_mass"] == 0.0  # the clonal fixture's filter removes nothing


def test_cli_reads_cov_headers_and_writes_a_report(tmp_path: Path) -> None:
    fasta = tmp_path / "asm.fa"
    fasta.write_text(f">c1 len=4000 cov=23.7\n{_random(GENOME_LEN, SEED)}\n>c2\n{_random(300, 9)}\n")
    assert [w for _, _, w in mult.contigs(fasta)] == [23.7, 1.0]
    out = tmp_path / "report.json"
    assert mult.main([str(fasta), "-k", "11", "21", "--kvmer", str(FIXTURES / "analyze"), str(FIXTURES / "analyze"),
                      "--output", str(out)]) == 0  # fmt: skip
    report = json.loads(out.read_text())
    assert set(report) == {"11", "21"}
    assert report["11"]["kvmer"]["classes"]["1"]["keys"] > 600
    assert report["21"]["assembly"]["keys"] > report["11"]["assembly"]["keys"] * 0.9
