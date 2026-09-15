import json
from pathlib import Path

from sequencing_error_model import recovery


def test_cigar_mode_recovers_example_spec() -> None:
    truth = recovery.example_spec()
    fitted, report = recovery.recover(truth, n_reads=3000, seed=3)
    assert [c.token for c in fitted.error_head] == [c.token for c in truth.error_head]
    assert report.failures() == [], report.scalars
    assert abs(report.scalars["q_lag1_fit"] - report.scalars["q_lag1_true"]) < 0.05, report.scalars


def test_cli_writes_report(tmp_path: Path) -> None:
    out = tmp_path / "report.json"
    code = recovery.main(["--reads", "400", "--seed", "1", "--output", str(out)])
    doc = json.loads(out.read_text())
    assert code == (1 if doc["failures"] else 0)
    assert {"op_tv", "rate_ratio", "q_position_tv_max"} <= doc["scalars"].keys()
    assert len(doc["curves"]["rate_by_q_true"]) == len(doc["curves"]["reported_q"])
