import numpy as np

from sequencing_error_model import baseline, recovery
from sequencing_error_model import generate as gen

GENOME = "".join(np.random.default_rng(0).choice(list("ACGT"), size=3000))


def _profile(seed: int, error_rate_scale: float = 1.0) -> baseline.Profile:
    """Reads from the example spec at their true placements, as a simulator's aligned reads."""
    rng = np.random.default_rng(seed)
    starts = rng.integers(0, len(GENOME) - 40, size=1500)
    templates = [GENOME[s : s + 40] for s in starts]
    mates = [int(m) for m in rng.integers(1, 3, size=len(templates))]
    reads = gen.generate(recovery.example_spec(), templates, mates, rng, error_rate_scale=error_rate_scale)
    qualities = ((r.sequence, [ord(c) - 33 for c in r.quality]) for r in reads)
    ref = np.unique(baseline.kmers([GENOME]))
    return baseline.profile(zip(templates, reads, mates, strict=True), qualities, ref)


def test_kmers_are_canonical_and_skip_n() -> None:
    seq = GENOME[:60]
    assert np.array_equal(baseline.kmers([seq]), baseline.kmers([gen._revcomp(seq)])[::-1])
    assert len(baseline.kmers([seq[:30] + "N" + seq[31:]])) == 10 + 9  # windows before and after the N
    assert np.array_equal(np.unique(baseline.kmers([GENOME], chunk=500)), np.unique(baseline.kmers([GENOME])))


def test_distances_rank_simulators() -> None:
    real, same, noisy = _profile(1), _profile(2), _profile(3, error_rate_scale=3.0)
    assert noisy.summary()["error_rate"] > 1.5 * real.summary()["error_rate"]
    near, far = baseline.distances(real, same), baseline.distances(real, noisy)
    for metric in ("rate", "substitution_rate", "rate_by_q", "rate_by_cycle", "edits_per_read_tv", "kmer_absent"):
        assert near[metric] < far[metric] / 2, metric
        assert baseline.verdict(near[metric], far[metric]) == "beats"
    # Error scaling leaves the Q process alone.
    assert near["q_cycle_tv"] < 0.05 and far["q_cycle_tv"] < 0.05
    assert baseline.verdict(0.30, 0.31) == "matches" and baseline.verdict(0.5, 0.2) == "trails"
