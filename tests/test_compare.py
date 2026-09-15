import json
import re
from dataclasses import replace
from pathlib import Path

import numpy as np
import pysam

from sequencing_error_model import compare, recovery
from sequencing_error_model import generate as gen
from sequencing_error_model.fit import error
from sequencing_error_model.sources import bam, pe_overlap

PCR = 0.01  # library substitutions per fragment base, shared by both mates


def test_pe_overlap_vs_reference_estimates_library_errors() -> None:
    example, rng = recovery.example_spec(), np.random.default_rng(7)
    head0 = example.error_head[0]
    bias = head0.params["bias"] - np.r_[0, np.full(9, 1.5)]  # ~3.5% errors, as in the pe-overlap recovery test
    truth = replace(
        example, error_head=(replace(head0, params={**head0.params, "bias": bias}), *example.error_head[1:])
    )
    genome, length = "".join(rng.choice(list("ACGT"), size=20_000)), 40

    def ends(frag: str) -> list[str]:
        return [
            (f + a + "N" * length)[:length] for f, a in ((frag, gen.ADAPTERS[0]), (gen._revcomp(frag), gen.ADAPTERS[1]))
        ]

    true_templates, read_templates = [], []
    for name, t1, t2 in gen.fragments([("g", genome)], gen.insert_sizes(45, 8), 3000, length, rng):
        match = re.fullmatch(r"g:(\d+)-(\d+)#\d+", name)
        assert match
        frag = genome[int(match[1]) - 1 : int(match[2])]
        assert ends(frag) == [t1, t2]
        hit = rng.random(len(frag)) < PCR
        library = "".join(
            str(rng.choice([x for x in "ACGT" if x != c])) if h else c for c, h in zip(frag, hit, strict=True)
        )
        true_templates += [t1, t2]
        read_templates += ends(library)
    mates = [1, 2] * (len(true_templates) // 2)
    reads = gen.generate(truth, read_templates, mates, rng)

    # Both mates read the library's substitutions, so they cancel in the overlap and show against the genome.
    pairs = (
        (a.sequence, [ord(c) - 33 for c in a.quality], b.sequence, [ord(c) - 33 for c in b.quality])
        for a, b in zip(reads[::2], reads[1::2], strict=True)
    )
    ev = pe_overlap.collect(pairs, 1, (2, 2))
    tokens = [c.token for c in truth.error_head]
    overlap_head = pe_overlap.fit(ev, tokens, ["QualityMarkov(1)"], iterations=5).error_head
    overlap = pe_overlap.table(ev, overlap_head)
    reference = gen.observations(zip(true_templates, reads, mates, strict=True), (2, 2), 1, "reference")

    found = compare.evidence(overlap, reference, by=("q", "mate"))
    assert abs(found["excess"] - PCR) < 0.2 * PCR, found

    # Model level: the reference-fitted head E carries the library errors on top of the overlap's.
    spec_a = replace(truth, error_head=overlap_head)
    spec_b = replace(truth, error_head=error.fit(reference, tokens, truth.quality_alphabet))
    fitted = compare.models(spec_a, spec_b, reference)
    assert abs(fitted["excess"] - PCR) < 0.3 * PCR, fitted
    assert fitted["log_likelihood_per_row"][1] > fitted["log_likelihood_per_row"][0]


def test_cli_reports_both_levels(tmp_path: Path) -> None:
    rng, length = np.random.default_rng(2), 50
    genome = "".join(rng.choice(list("ACGT"), size=5000))
    ref = tmp_path / "ref.fa"
    ref.write_text(f">g\n{genome}\n")
    placed = []  # (reference start, reverse, template) per mate, pairs interleaved
    for name, t1, t2 in gen.fragments([("g", genome)], gen.insert_sizes(70, 4), 300, length, rng):
        match = re.fullmatch(r"g:(\d+)-(\d+)#\d+", name)
        assert match
        placed += [(int(match[1]) - 1, False, t1), (int(match[2]) - length, True, t2)]
    mates = [1, 2] * (len(placed) // 2)
    reads = gen.generate(recovery.example_spec(), [t for *_, t in placed], mates, rng, error_rate_scale=0.2)
    for mate in (1, 2):
        records = (f"@p{i}\n{r.sequence}\n+\n{r.quality}\n" for i, r in enumerate(reads[mate - 1 :: 2]))
        (tmp_path / f"r{mate}.fastq").write_text("".join(records))
    path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(path), "wb", header={"SQ": [{"SN": "g", "LN": len(genome)}]}) as out:
        for i, ((start, rev, _), mate, read) in enumerate(zip(placed, mates, reads, strict=True)):
            seq, qual, cigar = read.sequence, read.quality, gen._CIGAR.findall(read.cigar)
            if rev:
                seq, qual, cigar = gen._revcomp(seq), qual[::-1], cigar[::-1]
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.reference_id, a.reference_start, a.mapping_quality = f"p{i // 2}", 0, start, 60
            a.flag = 0x1 | (0x10 if rev else 0) | (0x80 if mate == 2 else 0x40)
            a.cigarstring, a.query_sequence = "".join(n + o for n, o in cigar), seq
            a.query_qualities = pysam.qualitystring_to_array(qual)
            out.write(a)

    fastqs = [str(tmp_path / f"r{m}.fastq") for m in (1, 2)]
    spec_a, spec_b, report = tmp_path / "a", tmp_path / "b", tmp_path / "report.json"
    quality = ["--quality-tokens", "QualityMarkov(1)", "Position(3)"]
    assert pe_overlap.main([*fastqs, "--output", str(spec_a), "--iterations", "2", *quality]) == 0
    assert bam.main([str(path), str(ref), "--output", str(spec_b), *quality]) == 0
    args = [*fastqs, str(path), str(ref), "--overlap-spec", str(spec_a), "--reference-spec", str(spec_b)]
    assert compare.main([*args, "--output", str(report), "--min-exposure", "10"]) == 0
    doc = json.loads(report.read_text())
    assert doc["overlap_stats"]["pairs"] == 300 and doc["evidence"]["by"] == ["q", "mate"]
    # No library errors: the excess is small next to the rates themselves.
    assert abs(doc["evidence"]["excess"]) < 0.1 * doc["evidence"]["rate_b"], doc["evidence"]
    assert set(doc["models"]) == {"reference_rows", "overlap_rows"}
