"""`kmer` default mode composed from real skiver v0.3.2 outputs, and the synthetic-recovery harness.

No skiver binary is needed: the composition runs against the committed fixture and the reads it was made from
(regenerated, not committed), and the harness runs against a stub binary that hands back that same fixture. What
the stub checks is the loop's plumbing, not its numbers - the real numbers need the CI job with a released binary.
"""

import json
import random
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from sequencing_error_model import generate as gen
from sequencing_error_model import recovery
from sequencing_error_model import spec as spec_io
from sequencing_error_model.sources import fastq_quality, skiver_analyze
from sequencing_error_model.sources import kmer as kmer_source

sys.path.insert(0, str(Path(__file__).parent / "fixtures"))
from make_skiver_fixtures import SEED, make_reads  # noqa: E402

FIXTURES = Path(__file__).parent / "fixtures" / "skiver-v0.3.2"
QUALITY_TOKENS = ["QualityMarkov(1)", "Position(8)", "Context(1,1)"]


@pytest.fixture(scope="module")
def reads_fastq(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The reads the committed fixture was produced from, regenerated from the same seed."""
    path = tmp_path_factory.mktemp("kmer") / "reads.fastq"
    path.write_text(make_reads(random.Random(SEED)))
    return path


def test_filter_stats_report_what_the_outlier_filter_removed() -> None:
    a = skiver_analyze.read_analyze(FIXTURES / "analyze")
    stats = kmer_source.filter_stats(a)
    assert stats["keys"] == len(a.kvmer) and stats["error_mass"] > 0
    assert stats["keys_dropped"] == sum(not r.passes_filter for r in a.kvmer)
    assert 0.0 <= stats["error_mass_dropped_fraction"] <= 1.0
    # Clonal synthetic reads: the filter has no variation to remove, and on this run it removes nothing.
    assert stats["keys_dropped"] == 0 and stats["error_mass_dropped"] == 0.0


def test_refuses_a_value_too_short_to_observe_errors() -> None:
    a = skiver_analyze.read_analyze(FIXTURES / "analyze")
    empty = fastq_quality.QualityProfile(1, (2, 2))

    def fit(v: int) -> None:
        kmer_source.fit(replace(a, v=v), empty, None, ["QualityWindow(0)"], QUALITY_TOKENS, {})  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="too short to observe errors"):
        fit(2)
    # v passes the guard, then the empty profile stops the fit for want of an alphabet
    with pytest.warns(UserWarning, match="biased low"), pytest.raises(ValueError, match="no qualities"):
        fit(8)


def test_rejects_tokens_default_mode_cannot_identify() -> None:
    a = skiver_analyze.read_analyze(FIXTURES / "analyze")
    profile = fastq_quality.QualityProfile(1, (2, 2))
    with pytest.raises(ValueError, match="starts with QualityWindow"):
        kmer_source.fit(a, profile, None, ["Context(1,1)"], QUALITY_TOKENS, {})  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="not identifiable"):
        kmer_source.fit(a, profile, None, ["QualityWindow(0)", "Mate"], QUALITY_TOKENS, {})  # type: ignore[arg-type]


def test_cli_composes_a_spec_the_generator_draws_from(reads_fastq: Path, tmp_path: Path) -> None:
    out = tmp_path / "spec"
    argv = [str(FIXTURES / "analyze"), str(reads_fastq), "--output", str(out), "--quality-tokens", *QUALITY_TOKENS]
    assert kmer_source.main(argv) == 0
    spec = spec_io.load(out)
    assert spec.provenance["mode"] == "kmer" and spec.provenance["v"] == 13
    assert spec.quality_alphabet == (2, 12, 23, 37)
    assert [c.name for c in spec.error_head] == ["QualityWindow", "Context", "Position", "Strand", "GC"]
    assert "skiver_weibull" in spec.marginals
    rng = np.random.default_rng(5)
    templates = ["".join(rng.choice(list("ACGT"), size=150)) for _ in range(20)]
    reads = gen.generate(spec, templates, [1] * len(templates), rng)
    assert len(reads) == 20 and all(len(r.sequence) == len(r.quality) for r in reads)


def _stub_skiver(tmp_path: Path) -> str:
    """A `skiver analyze` stand-in: ignores the reads and copies the committed fixture to `-o <prefix>`."""
    path = tmp_path / "skiver"
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import shutil, sys\n"
        "from pathlib import Path\n"
        f"src = Path({str(FIXTURES)!r})\n"
        'prefix = Path(sys.argv[sys.argv.index("-o") + 1])\n'
        'for f in src.glob("analyze.*.csv"):\n'
        '    shutil.copy(f, prefix.with_name(prefix.name + f.name.removeprefix("analyze")))\n'
    )
    path.chmod(0o755)
    return str(path)


def test_synthetic_recovery_loop_runs_end_to_end(tmp_path: Path) -> None:
    """Generate, "run skiver", profile, fit, compare: the loop's plumbing, with a stub for the binary."""
    stub = _stub_skiver(tmp_path)
    truth = recovery.example_spec()
    fitted, doc = recovery.skiver_recovery(
        truth, stub, tmp_path / "run", n_reads=200, read_length=60, genome_length=3000, q_reads=400
    )
    assert doc["v"] == kmer_source.GOOD_V and doc["reads"] == 200
    assert [c.name for c in fitted.error_head] == ["QualityWindow", "Context", "Position", "Strand", "GC"]
    assert doc["outlier_filter"]["keys"] > 0 and set(doc["scalars"]) >= {"rate_ratio", "op_tv"}
    assert isinstance(doc["failures"], list)  # the stub's CSVs describe other reads, so the numbers will not pass
    assert json.dumps(doc)  # the report is serialisable, as the CLI writes it


def test_outlier_filter_cost_compares_both_runs_per_component(tmp_path: Path) -> None:
    stub = _stub_skiver(tmp_path)
    doc = recovery.outlier_filter_cost(
        recovery.example_spec(), stub, tmp_path / "cost", n_reads=200, read_length=60, genome_length=3000, q_reads=400
    )
    assert set(doc["runs"]) == {"filtered", "use_all"}
    assert doc["removed"]["keys_dropped"] == 0  # nothing to remove on clonal reads, so no cost on this fixture
    # Same kvmer keys either way, so the two specs agree; a filter that dropped keys would move these.
    assert doc["components"]["op"]["tv"] < 1e-9, doc["components"]["op"]
    assert abs(doc["components"]["context"]["log_odds_slope"] - 1) < 1e-6, doc["components"]["context"]


def test_stub_matches_the_real_binarys_interface(tmp_path: Path) -> None:
    """The stub is only useful if `run_skiver` calls it the way it calls skiver."""
    prefix = tmp_path / "out" / "analyze"
    prefix.parent.mkdir()
    (tmp_path / "reads.fastq").write_text("@r\nACGT\n+\nIIII\n")
    recovery.run_skiver(_stub_skiver(tmp_path), tmp_path / "reads.fastq", prefix, 11, 13, 8)
    assert skiver_analyze.read_analyze(prefix).v == 13
    with pytest.raises(FileNotFoundError, match="no skiver binary"):
        recovery.run_skiver("/nonexistent/skiver", tmp_path / "reads.fastq", prefix, 11, 13, 8)
