import math
import os
import shutil
from pathlib import Path

import pytest

from sequencing_error_model.sources import skiver_analyze as sa

# The skiver-compat workflow points this at freshly regenerated outputs.
FIXTURES = Path(os.environ.get("SKIVER_ANALYZE_DIR", Path(__file__).parent / "fixtures" / "skiver-v0.3.2"))


@pytest.fixture(scope="module")
def analyze() -> sa.SkiverAnalyze:
    return sa.read_analyze(FIXTURES / "analyze")


def test_kvmer(analyze: sa.SkiverAnalyze) -> None:
    assert (analyze.k, analyze.v) == (11, 13)  # make_skiver_fixtures.py
    assert any(r.passes_filter for r in analyze.kvmer)
    for r in analyze.kvmer:
        assert (len(r.key), len(r.consensus_value)) == (analyze.k, analyze.v)
        assert list(r.op_counts) == list(sa.OPS)
        assert len(r.consensus_count_up_to_v) == analyze.v
        assert r.consensus_count_up_to_v[-1] == r.consensus_count


def test_t_axes_follow_default_trims(analyze: sa.SkiverAnalyze) -> None:
    k, v = analyze.k, analyze.v
    ts = list(range(k + 3, k + v - 1))  # k + [1 + ignore_smallest_t, v - ignore_largest_t]
    assert [h.t for h in analyze.hazard] == ts
    for r in analyze.spectrum_by_t:
        assert list(r.freq_at_t) == ts
        assert r.total == sum(r.freq_at_t.values())


def test_survival_is_fitted_weibull(analyze: sa.SkiverAnalyze) -> None:
    er = analyze.error_rate
    for t in range(1, 11):
        assert analyze.survival[t] == pytest.approx(math.exp(-er.lambda_ * t**er.beta), abs=1e-5)


def test_marginals(analyze: sa.SkiverAnalyze) -> None:
    assert {r.op for r in analyze.spectrum} <= set(sa.OPS)
    assert all(r.forward <= r.total for r in analyze.spectrum)
    assert [b.lo for b in analyze.phred] == [2, 12, 23, 37]  # injected quality alphabet
    assert all(b.lo < b.hi for b in analyze.gc_content)
    assert {r.from_start for r in analyze.read_position} == {True, False}
    er = analyze.error_rate
    proportions = er.substitution_error_proportion + er.insertion_error_proportion + er.deletion_error_proportion
    assert proportions == pytest.approx(1, abs=1e-5)


def test_hazard_na(tmp_path: Path) -> None:
    path = tmp_path / "x.hazard_rate.csv"
    path.write_text("t,num_candidates,num_survival,hazard_ratio,5th_percentile,95th_percentile\n3,0,0,NA,0,0\n")
    assert sa.read_hazard_rate(path)[0].hazard_ratio is None


def test_unexpected_header(tmp_path: Path) -> None:
    path = tmp_path / "x.summary_phred.csv"
    path.write_text("qscore,num_correct\n2,5\n")
    with pytest.raises(sa.SkiverFormatError, match="unexpected header"):
        sa.read_phred(path)


def test_short_row(tmp_path: Path) -> None:
    path = tmp_path / "x.survival_rate.csv"
    path.write_text("t,survival_rate\n1\n")
    with pytest.raises(sa.SkiverFormatError, match=":2: expected 2 fields"):
        sa.read_survival_rate(path)


def test_pre_0_3_outputs_rejected(tmp_path: Path) -> None:
    for f in FIXTURES.glob("analyze.*.csv"):
        if "gc_content" not in f.name:
            shutil.copy(f, tmp_path / f.name)
    with pytest.raises(sa.SkiverFormatError, match="older than skiver 0.3"):
        sa.read_analyze(tmp_path / "analyze")
