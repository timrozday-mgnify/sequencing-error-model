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


@pytest.mark.parametrize(("aligner", "scale"), [("minibwa", 0.1), ("minimap2", 1.0)])
def test_aligner_bias(tmp_path: Path, aligner: str, scale: float) -> None:
    _require(aligner)
    report = recovery.aligner_bias(recovery.example_spec(), aligner, tmp_path, n_reads=200, error_rate_scale=scale)
    assert report["mapped_fraction"] > 0.9
    assert 0.9 < report["op_rate_ratio"]["substitution"] < 1.1
