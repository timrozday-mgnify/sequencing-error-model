from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from sequencing_error_model import generate as gen
from sequencing_error_model import recovery
from sequencing_error_model.spec import Component, ErrorModelSpec, SpecError, load


def latent_spec() -> ErrorModelSpec:
    """The example spec plus a 30 % class of reads with lower qualities and more errors at the same qualities."""
    base = recovery.example_spec()
    e_w, q_w = np.zeros((2, 10)), np.zeros((2, 4))
    e_w[1, 1:] = 1.5
    q_w[1] = [1.0, 0.5, 0.0, -1.0]
    return replace(
        base,
        quality_head=(*base.quality_head, Component("Latent(2)", {"weights": q_w})),
        error_head=(*base.error_head, Component("Latent(2)", {"weights": e_w})),
        latent=Component("Latent(2)", {"prior": np.array([0.7, 0.3])}),
    )


def _draw(rng: np.random.Generator, n: int) -> tuple[list[str], list[int]]:
    templates = ["".join(rng.choice(list("ACGT"), size=rng.integers(30, 41))) for _ in range(n)]
    return templates, [int(m) for m in rng.integers(1, 3, size=n)]


def _per_read(spec: ErrorModelSpec, seed: int) -> tuple[float, float]:
    """Variance/mean of edits per read and the sd of per-read mean Q, on 5000 generated reads."""
    rng = np.random.default_rng(seed)
    templates, mates = _draw(rng, 5000)
    reads = gen.generate(spec, templates, mates, rng)
    edits = np.array([gen._edits(t, r) for t, r in zip(templates, reads, strict=True)])
    mean_q = [np.mean([ord(c) - 33 for c in r.quality]) for r in reads]
    return float(edits.var() / edits.mean()), float(np.std(mean_q))


def test_em_recovers_latent_classes(tmp_path: Path) -> None:
    truth = latent_spec()
    rng = np.random.default_rng(0)
    templates, mates = _draw(rng, 1500)
    reads = gen.generate(truth, templates, mates, rng)
    fitted = recovery.cigar_mode(templates, reads, mates, truth)
    flat = recovery.cigar_mode(templates, reads, mates, recovery.example_spec())
    fitted.save(tmp_path)
    fitted = load(tmp_path)

    assert fitted.latent is not None
    np.testing.assert_allclose(fitted.latent.params["prior"], [0.7, 0.3], atol=0.03)
    e, q = fitted.error_head[-1].params["weights"], fitted.quality_head[-1].params["weights"]
    shift = e[1] - e[0]
    np.testing.assert_allclose(shift[1:] - shift[0], 1.5, atol=0.35)
    shift = q[1] - q[0]
    np.testing.assert_allclose(shift - shift.mean(), [0.875, 0.375, -0.125, -1.125], atol=0.2)

    # Generated reads carry the truth's per-read heterogeneity; without the class they don't.
    (d_true, s_true), (d_fit, s_fit), (d_flat, s_flat) = (_per_read(s, 1) for s in (truth, fitted, flat))
    assert abs(d_fit / d_true - 1) < 0.1 and abs(s_fit / s_true - 1) < 0.05
    assert d_flat < 0.6 * d_true and s_flat < 0.9 * s_true


def test_heads_must_carry_the_latent_token() -> None:
    truth = latent_spec()
    with pytest.raises(SpecError, match="latent layer's token"):
        replace(truth, latent=None)
    with pytest.raises(SpecError, match="latent layer's token"):
        replace(truth, latent=Component("Latent(3)", {"prior": np.ones(3) / 3}))
