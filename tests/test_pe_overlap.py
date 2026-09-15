from dataclasses import replace
from pathlib import Path

import numpy as np

from sequencing_error_model import generate as gen
from sequencing_error_model import recovery
from sequencing_error_model import spec as spec_io
from sequencing_error_model.sources import pe_overlap


def test_place_finds_offset_read_through_and_indels() -> None:
    rng = np.random.default_rng(0)
    for size in (30, 45, 60, 75):
        frag = "".join(rng.choice(list("ACGT"), size=size))
        r1, r2 = (
            (s + a)[:50].ljust(50, "N") for s, a in ((frag, gen.ADAPTERS[0]), (gen._revcomp(frag), gen.ADAPTERS[1]))
        )
        assert pe_overlap.place(r1, r2) == (size - 50, False)
        mid = 49 - (max(0, size - 50) + min(50, size)) // 2 + (size - 50)
        mutated = r2[:mid] + r2[mid + 1 :] + "A"  # delete a base mid-overlap
        s, indel = pe_overlap.place(r1, mutated)
        assert s is None or indel, size
    assert pe_overlap.place("ACGT" * 10, "") == (None, False)


def test_cli_fits_spec_from_fastq(tmp_path: Path) -> None:
    model = recovery.example_spec()
    rng = np.random.default_rng(1)
    templates, mates = recovery.paired_draw(40, 45, 8, genome_length=5000)(rng, 400)
    reads = gen.generate(model, templates, mates, rng, error_rate_scale=0.2)
    for mate in (1, 2):
        records = (f"@p{i}/{mate}\n{r.sequence}\n+\n{r.quality}\n" for i, r in enumerate(reads[mate - 1 :: 2]))
        (tmp_path / f"r{mate}.fastq").write_text("".join(records))
    out = tmp_path / "spec"
    args = [str(tmp_path / "r1.fastq"), str(tmp_path / "r2.fastq"), "--output", str(out), "--iterations", "2"]
    assert pe_overlap.main(args) == 0
    spec = spec_io.load(out)
    assert spec.provenance["stats"]["pairs"] == 200 and spec.provenance["identified_ops"] == ["substitution"]
    assert np.isneginf(spec.error_head[0].params["bias"][5:]).all()


def test_pe_overlap_recovers_substitution_head() -> None:
    example = recovery.example_spec()
    head0 = example.error_head[0]
    # ~3.5% errors, not the example's 14%: overlaps that noisy fail placement. (At ~1% the per-Q rates need
    # ~10k pairs to settle within 10%, too slow for CI.)
    bias = head0.params["bias"] - np.r_[0, np.full(9, 1.5)]
    quiet = (replace(head0, params={**head0.params, "bias": bias}), *example.error_head[1:])
    truth = replace(example, error_head=quiet)
    fitted, report = recovery.recover(
        truth,
        n_reads=6000,
        seed=4,
        mode=pe_overlap.mode,
        draw=recovery.paired_draw(40, 45, 8),
        substitutions_only=True,
    )
    stats = fitted.provenance["stats"]
    assert stats["disagree"] and stats["indel"] and stats["no_overlap"] < stats["pairs"] / 10, stats
    assert fitted.provenance["mode"] == "pe-overlap"
    assert report.failures() == [], (report.scalars, report.curves["rate_by_q_true"], report.curves["rate_by_q_fit"])
