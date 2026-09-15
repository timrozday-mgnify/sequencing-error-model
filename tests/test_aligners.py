import os
import shutil
from pathlib import Path

import pytest

from sequencing_error_model import recovery


def _require(tool: str) -> None:
    """CI sets REQUIRE_ALIGNERS=1 and installs both aligners; elsewhere a missing aligner skips."""
    if shutil.which(tool) is None:
        if os.environ.get("REQUIRE_ALIGNERS") == "1":
            pytest.fail(f"{tool} is not on PATH")
        pytest.skip(f"{tool} is not on PATH")


@pytest.mark.parametrize(
    ("aligner", "scale", "unclip"), [("minibwa", 0.1, False), ("minibwa", 0.1, True), ("minimap2", 1.0, False)]
)
def test_aligner_bias(tmp_path: Path, aligner: str, scale: float, unclip: bool) -> None:
    _require(aligner)
    truth = recovery.example_spec()
    report = recovery.aligner_bias(truth, aligner, tmp_path, n_reads=200, error_rate_scale=scale, unclip=unclip)
    assert report["mapped_fraction"] > 0.9
    assert 0.9 < report["op_rate_ratio"]["substitution"] < 1.1
    assert report["edits_hidden"] > 0
    assert report["failures_observable"] == [], {k: v for k, v in report.items() if k.endswith("_observable")}
