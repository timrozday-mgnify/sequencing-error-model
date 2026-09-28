import json
import re
from dataclasses import replace
from pathlib import Path
from typing import cast

import numpy as np
import pysam

from sequencing_error_model import compare, recovery
from sequencing_error_model import generate as gen
from sequencing_error_model.fit import error
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.sources import bam, pe_overlap
from sequencing_error_model.sources import skiver_analyze as sa

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
    skiver = Path(__file__).parent / "fixtures" / "skiver-v0.3.2" / "analyze"
    # No `kmer` spec on these reads (no skiver run), so the overlap spec stands in: this checks the plumbing.
    extra = ["--skiver", str(skiver), "--kmer-spec", str(spec_a), "--generated-reads", "100", "--read-length", "40"]
    assert compare.main([*args, "--output", str(report), "--min-exposure", "10", *extra]) == 0
    doc = json.loads(report.read_text())
    assert set(doc["kmer"]) == {"overlap", "reference", "generated"}
    assert doc["kmer"]["overlap"]["components"]["op"]["tv"] < 1e-12  # the stand-in spec against itself
    assert doc["kmer"]["generated"]["reference"]["reads"] == 100
    # Different reads, so only the shape of the skiver section is checked here, not agreement.
    for mode in ("overlap", "reference"):
        assert set(doc["skiver"][mode]) == {
            "phred", "gc_content", "read_position_start", "read_position_end", "spectrum", "hazard",
        }  # fmt: skip
        assert "skipped" not in doc["skiver"][mode]["phred"], doc["skiver"][mode]["phred"]
    assert doc["overlap_stats"]["pairs"] == 300 and doc["evidence"]["by"] == ["q", "mate"]
    # No library errors: the excess is small next to the rates themselves.
    assert abs(doc["evidence"]["excess"]) < 0.1 * doc["evidence"]["rate_b"], doc["evidence"]
    assert set(doc["models"]) == {"reference_rows", "overlap_rows"}


def _rate_bins(table: CountTable, field: str, half: bool = False) -> list[sa.RateBin]:
    """`summary_phred`/`summary_gc_content` rows at the table's own rate per value of `field`.

    With `half`, each 10 % GC bin is written as skiver's two 5 % bins, at rates 0.8x and 1.2x the true one on
    equal exposure, so the comparison has to pool them by exposure to get back the rate.
    """
    m = table.marginal(field, "op")
    rows: dict[Key, list[float]] = {}
    for (value, op), n in m.counts.items():
        r = rows.setdefault((value,), [0.0, 0.0])
        r[0] += n
        r[1] += n * (op != "=")
    bins = []
    for (value,), (total, errors) in sorted(rows.items(), key=str):
        lo, hi = value if isinstance(value, tuple) else (int(cast("int", value)), int(cast("int", value)))
        rate = errors / total
        halves = [(lo, (lo + hi) // 2, 0.8), ((lo + hi) // 2, hi, 1.2)] if half else [(lo, hi, 1.0)]
        for a, b, scale in halves:
            e = total / len(halves)
            bins.append(
                sa.RateBin(a, b, rate * scale, (0.0, 0.0), round(e * (1 - rate * scale)), round(e * rate * scale))
            )
    return bins


def _position_rows(table: CountTable) -> list[sa.ReadPositionRow]:
    rows = []
    for field, start in (("pos_start", True), ("pos_end", False)):
        for b in _rate_bins(table, field):
            rows.append(sa.ReadPositionRow(b.lo, start, b.num_correct, b.num_error))
    return rows


def _spectrum_rows(pairs: list[tuple[str, gen.Read]]) -> list[sa.SpectrumRow]:
    """`summary_error_spectrum` rows walked straight off the alignments, independently of `compare`."""
    counts: dict[tuple[str, str, str], int] = {}
    for template, read in pairs:
        padded = "." + template.upper() + "."
        for t, op in gen.align(template, read)[0]:
            if op == "=":
                continue
            prev, nxt = (padded[t], padded[t + 1]) if op.startswith("-") else (padded[t], padded[t + 2])
            counts[(op, prev, nxt)] = counts.get((op, prev, nxt), 0) + 1
    return [sa.SpectrumRow(op, p, n, total, total) for (op, p, n), total in counts.items()]


def test_skiver_marginals_recomputed_from_reference_truth() -> None:
    """Every marginal skiver reports, rebuilt from the same reads' truth, lands back on skiver's own numbers."""
    rng = np.random.default_rng(11)
    genome = "".join(rng.choice(list("ACGT"), size=4000))
    templates = [genome[s : s + 45] for s in rng.integers(0, len(genome) - 45, 800)]
    mates = [1, 2] * (len(templates) // 2)
    reads = gen.generate(recovery.example_spec(), templates, mates, rng)
    pairs = list(zip(templates, reads, strict=True))
    table = gen.observations(zip(templates, reads, mates, strict=True), (2, 2), 1, "reference")

    analyze = sa.SkiverAnalyze(
        k=11,
        v=13,
        error_rate=sa.ErrorRate(
            0.01,
            (0.0, 0.0),
            0.01,
            (0.0, 0.0),
            0.01,
            (0.0, 0.0),
            0.9,
            (0.0, 0.0),
            20,
            (0.0, 0.0),
            20.0,
            (0.0, 0.0),
            0.9,
            0.05,
            0.05,
        ),  # fmt: skip
        hazard=[],
        survival={},
        spectrum=_spectrum_rows(pairs),
        spectrum_by_t=[],
        phred=_rate_bins(table, "q"),
        gc_content=_rate_bins(table, "gc", half=True),
        read_position=_position_rows(table),
        kvmer=[],
    )

    found = compare.skiver_evidence(analyze, table, min_exposure=5)
    for name in ("phred", "gc_content", "read_position_start", "read_position_end"):
        m = found[name]
        assert "skipped" not in m, (name, m)
        assert abs(m["rate_ratio"] - 1) < 1e-4, (name, m)  # rounded bin counts, not a convention gap
        assert m["coverage"][1] > 0.9, (name, m["coverage"])
    s = found["spectrum"]
    assert {"->A", "A>-"} <= set(s["ops"]), s["ops"]  # the insertion and deletion context conventions are exercised
    assert s["tv"] < 1e-9 and s["r"] > 0.9999, s
    assert found["hazard"]["table_op_proportions"]["substitution"] > 0.5, found["hazard"]


def test_components_localise_a_moved_context_effect(tmp_path: Path) -> None:
    """Identical specs agree everywhere; halving head E's Context term shrinks the context log-odds slope."""
    truth, rng = recovery.example_spec(), np.random.default_rng(3)
    genome = "".join(rng.choice(list("ACGT"), size=3000))
    templates = [genome[s : s + 40] for s in rng.integers(0, len(genome) - 40, 600)]
    mates = [1, 2] * (len(templates) // 2)
    reads = gen.generate(truth, templates, mates, rng)
    table = gen.observations(zip(templates, reads, mates, strict=True), (2, 2), 1, "reference")

    same = compare.components(truth, truth, table, min_exposure=20)
    assert same["op"]["tv"] < 1e-12 and same["context"]["r"] > 0.9999
    assert abs(same["context"]["log_odds_slope"] - 1) < 1e-9, same["context"]
    for field in ("q", "pos_start", "mate"):
        assert max(abs(r - 1) for r in same[field]["ratio"]) < 1e-9, (field, same[field])

    ctx = next(c for c in truth.error_head if c.name == "Context")
    half = replace(
        truth,
        error_head=tuple(
            replace(c, params={**c.params, "weights": c.params["weights"] * 0.5}) if c is ctx else c
            for c in truth.error_head
        ),
    )
    moved = compare.components(truth, half, table, min_exposure=20)
    # A shrunk effect, but not halved: the per-context marginal carries every context-varying component (here
    # `Homopolymer` too, unchanged), so the slope reads the whole context shape, not one component's parameters.
    assert 0.5 < moved["context"]["log_odds_slope"] < 0.95, moved["context"]
    assert moved["context"]["r"] > 0.9 and moved["context"]["mean_abs_log_ratio"] > 0.02, moved["context"]

    ref = tmp_path / "ref.fa"
    ref.write_text(f">g\n{genome}\n")
    doc = compare.generated(truth, truth, ref, 200, 40, q_reads=400)
    assert abs(doc["scalars"]["rate_ratio"] - 1) < 1e-9 and doc["scalars"]["op_tv"] < 1e-12, doc["scalars"]
