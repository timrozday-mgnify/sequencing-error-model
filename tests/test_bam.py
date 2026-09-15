from dataclasses import replace
from pathlib import Path

import numpy as np
import pysam

from sequencing_error_model import generate as gen
from sequencing_error_model import recovery
from sequencing_error_model import spec as spec_io
from sequencing_error_model.observations import CountTable
from sequencing_error_model.sources import bam


def _write(tmp_path: Path, n: int) -> tuple[Path, Path, list[str], list[gen.Read], list[int]]:
    """Generated reads written as a BAM with their true CIGARs, plus records the source must skip."""
    rng = np.random.default_rng(3)
    genome = "".join(rng.choice(list("ACGT"), size=3000))
    ref = tmp_path / "ref.fa"
    ref.write_text(f">chr\n{genome}\n")
    starts, reverse = rng.integers(0, 2950, size=n), rng.random(n) < 0.5
    mates = [int(m) for m in rng.integers(1, 3, size=n)]
    fwd = [genome[s : s + int(rng.integers(30, 45))] for s in starts]
    templates = [gen._revcomp(t) if r else t for t, r in zip(fwd, reverse, strict=True)]
    reads = gen.generate(recovery.example_spec(), templates, mates, rng)
    path = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr", "LN": len(genome)}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (start, rev, mate, read) in enumerate(zip(starts, reverse, mates, reads, strict=True)):
            seq, qual, cigar = read.sequence, [ord(c) - 33 for c in read.quality], gen._CIGAR.findall(read.cigar)
            if rev:  # BAM-leading clip: after the aligned part in read orientation
                seq, qual, cigar = "GG" + gen._revcomp(seq), [20, 21] + qual[::-1], [("2", "S"), *cigar[::-1]]
            else:
                seq, qual, cigar = "TTT" + seq, [30] * 3 + qual, [("3", "S"), *cigar]
            for flag, mapq in ((0, 60), (0x100, 60), (0, 5), (0x400, 60)):  # kept, secondary, low MAPQ, duplicate
                a = pysam.AlignedSegment(out.header)
                a.query_name, a.reference_id, a.reference_start, a.mapping_quality = f"r{i}", 0, int(start), mapq
                a.flag = flag | 0x1 | (0x10 if rev else 0) | (0x80 if mate == 2 else 0x40)
                a.cigarstring = "".join(n + o for n, o in cigar)
                a.query_sequence = seq
                a.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in qual))
                out.write(a)
    return path, ref, templates, reads, mates


def test_records_reproduce_generator_tuples(tmp_path: Path) -> None:
    path, ref, templates, reads, mates = _write(tmp_path, 300)
    got_reads = [r for _, r, _ in bam.records(path, ref)]
    # Clipped bases shift positions and fill Q windows; strand comes from the alignment, not the mate.
    clipped = [
        replace(r, clipped=("", "65"), strand="-") if g.strand == "-" else replace(r, clipped=("???", ""), strand="+")
        for r, g in zip(reads, got_reads, strict=True)
    ]
    assert {g.strand for g in got_reads} == {"+", "-"}
    flank, m = (2, 2), 1
    expected = gen.observations(zip(templates, clipped, mates, strict=True), flank, m)
    got = gen.observations(bam.records(path, ref), flank, m)
    assert any("-" in str(k[-1]) for k in expected.counts)  # indels exercised
    fields = expected.fields
    assert any(k[fields.index("pos_start")] == 4 for k in got.counts)  # position counts the 3-base clip
    assert got.counts == expected.counts


def _subs(table: CountTable) -> int:
    return sum(n for k, n in table.counts.items() if ">" in str(k[-1]) and "-" not in str(k[-1]))


def test_pileup_and_masks(tmp_path: Path) -> None:
    path, ref, *_ = _write(tmp_path, 300)
    counts = bam.pileup(path, ref)
    table = gen.observations(bam.records(path, ref), (2, 2), 1)
    # The pileup counts the same draws as the tuples, in reference orientation.
    assert counts["chr"].sum() == sum(table.counts.values())
    genome = ref.read_text().split("\n")[1]
    matches = counts["chr"][np.arange(len(genome)), ["ACGT".index(b) for b in genome]].sum()
    assert matches == table.marginal("op").counts[("=",)]
    depth = counts["chr"][:, :5].sum(axis=1)
    p = int(depth.argmax())
    bad = tmp_path / "bad.fa"  # a reference error at the deepest site
    bad.write_text(f">chr\n{genome[:p]}{'C' if genome[p] == 'A' else 'A'}{genome[p + 1 :]}\n")
    masked, report = bam.masks(bam.pileup(path, bad), bad, 0.5)
    # The fixture has ~14% errors at depth ~4, so chance sites are masked too; the clonal cost is a later item.
    assert masked["chr"][p] and report[0]["masked_sites"] == masked["chr"].sum()
    before = gen.observations(bam.records(path, bad), (2, 2), 1)
    after = gen.observations(bam.records(path, bad, masked=masked), (2, 2), 1)
    assert _subs(before) - _subs(after) >= 0.8 * depth[p]
    assert sum(before.counts.values()) - sum(after.counts.values()) >= depth[p]
    assert bam.masks(counts, ref, None, min_contig_depth=1e6)[0] == {}
    assert list(bam.records(path, ref, masked={})) == []


def test_unclip_recovers_clipped_errors(tmp_path: Path) -> None:
    rng = np.random.default_rng(11)
    genome = "".join(rng.choice(list("ACGT"), size=2000))
    ref = tmp_path / "ref.fa"
    ref.write_text(f">chr\n{genome}\n")
    errors = list(genome[500:560])
    for p in (50, 55):  # two substitutions in the 15 bases the aligner clipped
        errors[p] = "A" if errors[p] != "A" else "C"
    adapter = genome[800:845] + "".join(rng.choice(list("ACGT"), size=15))  # 15 bases that aren't the genome
    short = list(genome[1200:1248])
    for p in (45, 47):  # two substitutions in a 3-base clip: too short to tell from an adapter, so realigned
        short[p] = "A" if short[p] != "A" else "C"
    path = tmp_path / "clipped.bam"
    with pysam.AlignmentFile(str(path), "wb", header={"SQ": [{"SN": "chr", "LN": len(genome)}]}) as out:
        reads = ((500, "".join(errors), "45M15S"), (800, adapter, "45M15S"), (1200, "".join(short), "45M3S"))
        for i, (start, seq, cigar) in enumerate(reads):
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.reference_id, a.reference_start, a.mapping_quality = f"r{i}", 0, start, 60
            a.cigarstring, a.query_sequence = cigar, seq
            a.query_qualities = pysam.qualitystring_to_array("?" * len(seq))
            out.write(a)
    scores = (2, 8, 12, 2)
    assert [r.cigar for _, r, _ in bam.records(path, ref)] == ["45M", "45M", "45M"]
    # The error-rich clips are realigned end to end; the adapter clip is kept.
    assert [r.cigar for _, r, _ in bam.records(path, ref, unclip=scores)] == ["60M", "45M", "48M"]
    assert _subs(gen.observations(bam.records(path, ref), (2, 2), 1)) == 0
    assert _subs(gen.observations(bam.records(path, ref, unclip=scores), (2, 2), 1)) == 4


def test_cli_fits_spec(tmp_path: Path) -> None:
    path, ref, *_ = _write(tmp_path, 200)
    out, report = tmp_path / "spec", tmp_path / "contigs.tsv"
    args = [str(path), str(ref), "--output", str(out), "--error-tokens", "QualityWindow(1)", "Context(1,1)"]
    args += ["--mask-alt-freq", "0.5", "--contig-report", str(report)]
    assert bam.main([*args, "--quality-tokens", "QualityMarkov(1)", "Position(3)"]) == 0
    spec = spec_io.load(out)
    # The Q±1 window reaches the clipped base beside each aligned end: Q30 (forward) and Q21 (reverse).
    assert spec.provenance["mode"] == "reference" and spec.quality_alphabet == (2, 12, 21, 23, 30, 37)
    assert spec.provenance["contigs_kept"] == 1 and report.read_text().startswith("contig\tlength\tmean_depth")
