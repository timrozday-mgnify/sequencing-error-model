import gzip
import json
from pathlib import Path

import numpy as np
import pytest

from sequencing_error_model import export
from sequencing_error_model.export import art, badread, iss, pbsim3, roundtrip
from sequencing_error_model.generate import Read


def test_remap_moves_errors_into_q() -> None:
    s = export.sample(roundtrip.short_read_spec(), 2000, 50, np.random.default_rng(0))
    spec_rate = export.error_rate(export.columns(s), "substitution")
    # The example's Q is miscalibrated, so the simulators' 10^(-Q/10) is off before the remap and on after it.
    assert export.implied_rate(export.columns(s), "substitution") / spec_rate > 2
    remapped = export.columns(export.remap(s, "substitution", True, 40))
    assert export.implied_rate(remapped, "substitution") == pytest.approx(spec_rate, rel=0.05)


def test_art_profile_and_base_profile_fallback(tmp_path: Path) -> None:
    model = roundtrip.short_read_spec()
    report = art.export(model, tmp_path, np.random.default_rng(0), read_length=40, reads=1000)
    assert report["q_tv"] == 0 and "-sp" not in report["arguments"]
    qs, probs = export.read_art_profile(tmp_path / "art_R1.txt")
    assert set(qs) <= set(model.quality_alphabet) and probs.shape == (40, len(qs))
    np.testing.assert_allclose(probs.sum(axis=1), 1)
    lines = (tmp_path / "art_R1.txt").read_text().splitlines()
    assert [ln[0] for ln in lines[:: 2 * 40]] == list(".ATGCN")

    # A spec whose Q came from nowhere: draw it from that profile instead, snapped to the alphabet.
    fallback = iss.export(
        model,
        tmp_path,
        np.random.default_rng(1),
        read_length=40,
        reads=1000,
        insert_mean=200.0,
        base_profile=[tmp_path / "art_R1.txt", tmp_path / "art_R2.txt"],
    )
    assert fallback["base_profile"] == [str(tmp_path / f"art_R{m}.txt") for m in (1, 2)]
    assert fallback["spec_mean_q"] == pytest.approx(report["spec_mean_q"], rel=0.05)


def test_iss_npz(tmp_path: Path) -> None:
    report = iss.export(
        roundtrip.short_read_spec(),
        tmp_path,
        np.random.default_rng(0),
        read_length=40,
        reads=1000,
        insert_mean=200.0,
        insert_sd=10.0,
        q_policy="preserve-errors",
    )
    assert report["error_rate_ratio"] == pytest.approx(1, abs=0.05)
    with np.load(tmp_path / "iss.npz", allow_pickle=True) as npz:
        assert npz["model"] == "kde" and npz["read_length"] == 40
        assert npz["insert_size"][-1] == pytest.approx(1)
        for bin_, n in zip(npz["quality_hist_forward"], npz["mean_count_forward"], strict=True):
            assert bool(len(bin_)) == bool(n)
            assert all(cdf[-1] == pytest.approx(1) for cdf in bin_)
        choice = npz["subst_choices_reverse"][3]["G"]
        assert choice[0] == ["A", "T", "C"] and sum(choice[1]) == pytest.approx(1)
    kde = pytest.importorskip("iss.error_models.kde")  # InSilicoSeq's own loader, where installed
    model = kde.KDErrorModel(str(tmp_path / "iss.npz"))
    assert len(model.gen_phred_scores(model.quality_forward, "forward")) == 40


def test_badread_models(tmp_path: Path) -> None:
    report = badread.export(roundtrip.long_read_spec(), tmp_path, np.random.default_rng(0), read_length=400, reads=30)
    assert report["q_policy"] is None and "--identity" in report["arguments"]
    with gzip.open(tmp_path / "badread_error_model.gz", "rt") as handle:
        for line in handle:
            entries = [e.split(",") for e in line.strip().split(";") if e]
            assert len(entries[0][0]) == badread.K_ERROR and len(entries) <= badread.MAX_ALT + 1
            assert sum(float(p) for _, p in entries) <= 1 + 1e-5
    with gzip.open(tmp_path / "badread_qscore_model.gz", "rt") as handle:
        cigars = [line.split(";")[0] for line in handle]
    assert cigars[0] == "overall" and {"=", "X", "I"} <= set(cigars)
    pytest.importorskip("edlib")
    models = pytest.importorskip("badread.error_model"), pytest.importorskip("badread.qscore_model")
    assert models[0].ErrorModel(str(tmp_path / "badread_error_model.gz")).kmer_size == badread.K_ERROR
    assert models[1].QScoreModel(str(tmp_path / "badread_qscore_model.gz")).kmer_size == badread.K_QSCORE


def _hmm(path: Path) -> dict[int, dict[str, dict[int, list[float]]]]:
    out: dict[int, dict[str, dict[int, list[float]]]] = {}
    for line in path.read_text().splitlines():
        level, kind, state, *p = line.split()
        out.setdefault(int(level), {}).setdefault(kind, {})[int(state)] = [float(x) for x in p]
    return out


def test_pbsim3_models(tmp_path: Path) -> None:
    report = pbsim3.export(roundtrip.long_read_spec(), tmp_path, np.random.default_rng(0), read_length=400, reads=30)
    for name, states in (("qshmm", pbsim3.STATE_MAX), ("errhmm", 4)):
        hmm = _hmm(tmp_path / f"pbsim3_{name}.model")
        assert sorted(hmm) == list(pbsim3.LEVELS) and hmm[1] == hmm[99]
        assert len(hmm[50]["TP"]) <= states
        assert sum(hmm[50]["IP"][s][0] for s in hmm[50]["IP"]) == pytest.approx(1, abs=1e-4)
        assert all(sum(row) == pytest.approx(1, abs=1e-4) for row in hmm[50]["TP"].values())
    ratio = report["arguments"][report["arguments"].index("--difference-ratio") + 1]
    assert sum(map(int, ratio.split(":"))) == pytest.approx(1000, abs=2)


def test_qshmm_merges_q_into_50_states() -> None:
    quality = "".join(chr(33 + q) for q in range(61))
    s = export.Sample(["A" * 61], [Read("A" * 61, quality, "61M", np.zeros(0))], [1])
    hmm = {}
    for line in pbsim3.qshmm(s).splitlines():
        level, kind, state, *p = line.split()
        if level == "1":
            hmm[kind, int(state)] = [float(x) for x in p]
    assert max(state for _, state in hmm) == pbsim3.STATE_MAX
    emitted = [np.flatnonzero(hmm["EP", i]) for i in range(1, pbsim3.STATE_MAX + 1)]
    assert np.array_equal(np.concatenate(emitted), np.arange(61))  # adjacent Q values share a state


def test_cli(tmp_path: Path) -> None:
    spec_dir = tmp_path / "spec"
    roundtrip.long_read_spec().save(spec_dir)
    out = tmp_path / "out"
    argv = ["pbsim3", "--model", str(spec_dir), "--output", str(out), "--reads", "10", "--read-length", "300"]
    assert export.main([*argv, "--q-policy", "preserve-errors"]) == 0
    report = json.loads((out / "fidelity.json").read_text())
    assert report["preserves"] == "errors" and report["dropped"] and report["errhmm_arguments"]
