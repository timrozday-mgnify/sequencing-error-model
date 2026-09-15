import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from sequencing_error_model.spec import Component, ErrorModelSpec, SpecError, load


def make(**overrides: Any) -> ErrorModelSpec:
    kwargs: dict[str, Any] = {
        "quality_alphabet": (2, 12, 23, 37),
        "provenance": {"mode": "default", "sources": ["skiver_analyze", "fastq_quality"], "k": 31, "v": 13},
        "quality_head": (
            Component("QualityMarkov(1)", {"transitions": np.arange(16, dtype=np.int64).reshape(4, 4)}),
            Component("Position(4)", {"knots": np.linspace(0, 1, 4)}, meta={"from": ["start", "end"]}),
        ),
        "error_head": (
            Component("Context(2,2)", {"weights": np.random.default_rng(0).normal(size=(5, 4, 10))}),
            Component("QualityWindow(1)", {"weights": np.zeros((3, 10))}, identified=False),
            Component("FragmentOverdispersion", {"phi": np.array(50.0)}, generative=False),
        ),
        "latent": Component("Latent(2)", {"transitions": np.eye(2)}),
        "marginals": {"rate": np.array([0.001])},
    }
    return ErrorModelSpec(**(kwargs | overrides))


def test_round_trip(tmp_path: Path) -> None:
    spec = make()
    spec.save(tmp_path / "model")
    got = load(tmp_path / "model")
    assert (got.quality_alphabet, got.provenance, got.latent is not None) == (
        spec.quality_alphabet,
        spec.provenance,
        True,
    )
    assert got._arrays().keys() == spec._arrays().keys()
    for key, arr in spec._arrays().items():
        np.testing.assert_array_equal(got._arrays()[key], arr)
        assert got._arrays()[key].dtype == arr.dtype
    for a, b in zip(got._components(), spec._components(), strict=True):
        assert (a[1].token, a[1].generative, a[1].identified, a[1].meta) == (
            b[1].token,
            b[1].generative,
            b[1].identified,
            b[1].meta,
        )
    assert spec.error_head[0].name == "Context" and spec.error_head[0].args == (2, 2)
    assert spec.error_head[2].args == ()


@pytest.mark.parametrize(
    ("overrides", "match"),
    [
        ({"quality_alphabet": (37, 2)}, "quality alphabet"),
        ({"quality_alphabet": (2, 94)}, "quality alphabet"),
        ({"provenance": {}}, "mode"),
        ({"error_head": (Component("QualityMarkov(1)"),)}, "not a head E"),
        ({"quality_head": (Component("Contxt(2,2)"),)}, "not a head Q"),
        ({"error_head": (Component("Context(1,1)"), Component("Context(2,2)"))}, "repeats"),
        ({"latent": Component("Latent")}, "Latent"),
        ({"marginals": {"x": np.array([object()])}}, "non-object"),
    ],
)
def test_invalid(overrides: dict[str, Any], match: str) -> None:
    with pytest.raises(SpecError, match=match):
        make(**overrides)


def test_load_rejects_other_schema_and_stray_arrays(tmp_path: Path) -> None:
    make().save(tmp_path)
    doc = json.loads((tmp_path / "spec.json").read_text())
    (tmp_path / "spec.json").write_text(json.dumps(doc | {"schema_version": 999}))
    with pytest.raises(SpecError, match="schema version"):
        load(tmp_path)
    (tmp_path / "spec.json").write_text(json.dumps(doc | {"marginals": []}))
    with pytest.raises(SpecError, match="not referenced"):
        load(tmp_path)
    (tmp_path / "spec.json").write_text(json.dumps(doc | {"marginals": ["missing"]}))
    with pytest.raises(SpecError, match="malformed"):
        load(tmp_path)
