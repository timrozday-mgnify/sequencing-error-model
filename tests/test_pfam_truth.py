import json
from collections import Counter
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from test_pfam import ERROR_TOKENS, QUALITY_TOKENS, World, substitutions_only

from sequencing_error_model import pfam_truth
from sequencing_error_model.generate import _revcomp
from sequencing_error_model.sources import pfam


class Mock(World):
    """A `World` whose genes are the mock's genomes; remembers each fragment's gene."""

    def fragments(self, n: int, insert: int) -> list[tuple[str, tuple[int, int, int, bool]]]:
        out = []
        self.source: list[tuple[str, int, bool]] = []
        for _ in range(n):
            f = int(self.rng.integers(len(self.genes)))
            g = int(self.rng.integers(len(self.genes[f])))
            start = int(self.rng.integers(len(self.genes[f][g]) - insert + 1))
            strand = bool(self.rng.integers(2))
            frag = self.genes[f][g][start : start + insert]
            self.source.append((f"g{f}_{g}", start, strand))
            out.append((_revcomp(frag) if strand else frag, (f, start, insert, strand)))
        return out

    def sam(self, path: Path, pairs: list[pfam.ReadPair], insert: int, drop: set[int], clip: int = 0) -> None:
        with path.open("w") as out:
            for f, genes in enumerate(self.genes):
                out.writelines(f"@SQ\tSN:g{f}_{i}\tLN:{len(g)}\n" for i, g in enumerate(genes))
            for number, ((s1, _, s2, _), (contig, start, strand)) in enumerate(zip(pairs, self.source, strict=True)):
                for mate, s in enumerate((s1, s2)):
                    if 2 * number + mate in drop:
                        out.write(f"{number}\t{69 + 64 * mate}\t*\t0\t0\t*\t*\t0\t0\t{s}\t*\n")
                        continue
                    reverse = (mate == 1) != strand
                    lo = start + insert - len(s) if reverse else start
                    flag = 1 + 64 * (mate + 1) + 16 * reverse
                    seq = _revcomp(s) if reverse else s
                    cigar = f"{clip}S{len(s) - clip}M" if clip else f"{len(s)}M"
                    out.write(f"{number}\t{flag}\t{contig}\t{lo + clip + 1}\t60\t{cigar}\t*\t0\t0\t{seq}\t*\n")

    def genome(self) -> dict[str, str]:
        return {f"g{f}_{i}": g for f, genes in enumerate(self.genes) for i, g in enumerate(genes)}


Scored = tuple[dict[str, Any], dict[str, int], Path]


def with_stop(gene: str, codon: int) -> str:
    return gene[: 3 * codon] + "TAA" + gene[3 * codon + 3 :]


@pytest.fixture(scope="module")
def scored(tmp_path_factory: pytest.TempPathFactory) -> Scored:
    """A library with planted genome stops, a reference error at a stop, and unaligned reads."""
    tmp = tmp_path_factory.mktemp("truth")
    world = Mock(3)
    reference = world.genome()
    for g in range(4):  # real stops, in the sample and in the reference
        world.genes[0][g] = reference[f"g0_{g}"] = with_stop(world.genes[0][g], 17)
    for g in range(4):  # in the sample only: the reference differs from the strain sequenced
        world.genes[1][g] = with_stop(world.genes[1][g], 17)
    truth = substitutions_only(0.05)
    pairs, hits, templates, _ = world.pairs(truth, 6000, 210, 100)
    drop = set(range(0, 12000, 10))
    world.sam(tmp / "aln.sam", pairs, 210, drop)
    ev = pfam.collect(pairs, hits, world.profiles, t1=30, quality_every=10)
    found = pfam.drop_replicated(ev, pairs, hits)
    assert found["stops"] == ev.stats.stops and found["replicated"] + len(ev.stops) == ev.stats.stops
    sourced = [s for s in ev.stops if s.source]
    assert len(sourced) == found["sourced"] > 0 and all(s.source in s.sources for s in sourced)
    spec, _ = pfam.fit(ev, world.profiles, ERROR_TOKENS, QUALITY_TOKENS)
    report = pfam_truth.score(
        pairs,
        hits,
        world.profiles,
        tmp / "aln.sam",
        reference,
        {"both": spec, "true": truth, "stops": spec, "dose0.002": spec},
        t1=30,
        max_alt_freq=0.5,  # 15 reads deep deep, where two errors alike happen
    )
    # what each read's own template says about its exposed stops
    want = dict.fromkeys(pfam_truth.CLASSES, 0)
    for number, (s1, _, s2, _) in enumerate(pairs):
        for mate, s in enumerate((s1, s2)):
            rid = 2 * number + mate
            contig, _, _ = world.source[number]
            for c in pfam.exposed(s, hits.get(rid, []), world.profiles, 30, 6):
                if c.bases not in pfam.STOPS:
                    continue
                sequenced = "".join(templates[rid][x] for x in c.positions)
                sequenced = _revcomp(sequenced)[::-1] if c.reverse else sequenced
                if rid in drop:
                    want["unaligned"] += 1
                elif sequenced in pfam.STOPS:
                    want["genome" if contig.startswith("g0_") else "variant"] += 1
                else:
                    want["error"] += 1
    return report, want, tmp


def test_stops_are_classified_by_the_genome_under_them(scored: Scored) -> None:
    report, want, _ = scored
    left, gone = report["stops"], report["replicated_stops"]
    got = {k: left[k] + gone[k] for k in want}
    moved = got["variant"] - want["variant"]  # an error on a shallow site where most reads erred alike
    assert 0 <= moved <= 3 and got == want | {"variant": want["variant"] + moved, "error": want["error"] - moved}
    assert want["error"] > 50 and want["genome"] > 20 and want["variant"] > 20 and want["unaligned"] > 5
    # planted stops recur at their locus; errors do not, so the filter leaves the floor and takes little else
    assert left["genome"] + left["variant"] <= 5 and gone["error"] <= 0.03 * want["error"]
    assert report["floor"]["rate"] < 0.1 * report["floor_before"]["rate"]
    assert {"PF00000", "PF00001"}.isdisjoint(report["fits"]["both"]["downweighted_families"])
    assert report["fits"]["both"]["replication"]["replicated"] == sum(gone.values())


def test_recurrence_is_judged_against_what_errors_could_do() -> None:
    assert not pfam.recurs(1, 0) and not pfam.recurs(1, 50)
    assert pfam.recurs(2, 0) and pfam.recurs(3, 7)
    assert not pfam.recurs(2, 30) and not pfam.recurs(2, 98)
    assert pfam.recurs(20, 80)


def test_replicates_count_fragments_between_the_same_flanks() -> None:
    left, right = "ACGTTGCAACGT", "GGATCCGATTAC"
    locus = (7, "TAA", (left, right))
    stop, sense, other = (f"CC{left}{c}{right}GG" for c in ("TAA", "CAA", "TAA"))
    pairs = [
        (stop, b"", _revcomp(stop), b""),  # its own fragment, both mates: counted once
        (sense, b"", "", b""),
        (_revcomp(stop), b"", "", b""),
        (other.replace(right, right[::-1]), b"", "", b""),  # another locus
        (stop, b"", "", b""),  # no hit to the family: not searched
    ]
    hit = pfam.Hit(7, False, 0, 0, ((9, "M"),))
    hits = {2 * i + m: [hit] for i in range(4) for m in (0, 1)}
    assert pfam.replicates([locus, (7, "TAA", None)], pairs, hits) == [(2, 1), (0, 0)]
    # the other fragment's codon is kept: it names the stop's source
    assert pfam.at_loci([locus, (7, "TAA", None)], pairs, hits) == [(2, Counter({"CAA": 1})), (0, Counter())]


def test_a_locus_source_takes_the_stops_weight() -> None:
    class Flat:
        def weight(self, column: object, codon: str) -> float:
            return 2.0

    stop = pfam.Stop(7, (7, 3), "TAA", ("CAA", "AAA", "GAA"), (), (), (), 40, None)
    assert np.allclose(pfam.codon_weight(stop, Flat()), 2.0)  # type: ignore[arg-type]
    w = pfam.codon_weight(replace(stop, source="AAA"), Flat())  # type: ignore[arg-type]
    assert w[1] == pfam.LOCUS_RIGHT and np.isclose(w.sum(), 1.0) and w[0] == w[2] > 0


def test_observed_errors_match_the_true_head_per_cell(scored: Scored) -> None:
    report, _, _ = scored
    cells = report["cells"]
    assert {c.split(":")[0] for c in cells} == {"op", "q", "mate", "position"}
    assert sum(c["observed"] for n, c in cells.items() if n.startswith("op:") and n[3:] in pfam_truth.TEN) == sum(
        c["observed"] for n, c in cells.items() if n.startswith("mate:")
    )
    for name, c in cells.items():
        if c["observed"] > 30:
            got, sd = c["predicted"]["true"] / c["observed"], 1 / np.sqrt(c["observed"])
            assert got == pytest.approx(1, abs=4 * sd + 0.05), name


def test_checks_and_verdict(scored: Scored) -> None:
    report, _, tmp = scored
    by = {c["check"]: c for c in pfam_truth.checks(report, "true")}
    assert by["total"]["ratio"] == pytest.approx(1, abs=0.05) and by["total"]["result"] == "within"
    assert "op:A>G" not in by and by["floor"]["result"] in ("within", "unclear", "outside")
    # a fit that is 60 % too high fails, whatever its library's power allows elsewhere
    high = json.loads(json.dumps(report))
    for c in high["cells"].values():
        c["predicted"]["both"] *= 1.6
    assert next(c for c in pfam_truth.checks(high) if c["check"] == "total")["result"] != "within"
    assert pfam_truth.passed(pfam_truth.checks(high)) is not True
    # no power: neither pass nor fail
    weak = json.loads(json.dumps(report))
    for c in weak["cells"].values():
        c["stops"] = 1
    assert pfam_truth.passed(pfam_truth.checks(weak)) is None
    assert pfam_truth.verdict({str(i): weak for i in range(7)})["verdict"] == "undecided"
    assert pfam_truth.verdict({"only": high})["verdict"] == "undecided"  # too few libraries to decide
    good = json.loads(json.dumps(report))
    good["floor_kept_families"] = {"stops": 0, "rate": 0.0, "upper95": 1e-5}
    for c in good["cells"].values():
        c["predicted"]["both"], c["stops"], c["observed"] = c["observed"], 1000, c["observed"] or 1
    v = pfam_truth.verdict({str(i): good if i < 5 else high for i in range(7)})
    assert v["verdict"] == "go" and v["passed"] == 5
    assert pfam_truth.verdict({str(i): high for i in range(7)})["verdict"] != "go"
    path = tmp / "lib.json"
    path.write_text(json.dumps(report))
    assert pfam_truth.main(["verdict", str(path), "--output", str(tmp / "v")]) == 0
    assert "lib" in (tmp / "v.md").read_text() and "Locus" not in (tmp / "v.md").read_text()
    # the locus reports are shown and change nothing
    one = {"flank": 12, "others": 2, "labels": [[1, 2, 3]], "omega": 0.1}
    for dose, rate in (("0", 1e-3), ("0.01", 1.1e-2)):
        per_base = {"mixture": rate, "mixture_sd": 1e-4, "naive": 2e-3, "nonsense": 1e-3}
        (tmp / f"lib.{dose}.json").write_text(json.dumps({"loci": [one | {"per_base": per_base}]}))
    with_loci = ["--loci", str(tmp / "lib.0.json"), str(tmp / "lib.0.01.json")]
    assert pfam_truth.main(["verdict", str(path), *with_loci, "--output", str(tmp / "w")]) == 0
    text = (tmp / "w.md").read_text()
    assert "| lib | 6 | - | nan | 1.00e-03" in text and "0.01: 1.00" in text
    plain, shown = (json.loads((tmp / f"{x}.json").read_text()) for x in "vw")
    assert plain["runs"] == shown["runs"] and plain["verdict"] == shown["verdict"]


def test_codons_off_the_alignment_are_not_labelled() -> None:
    codon = pfam.Codon((0, 0), (3, 4, 5), False, "TAA", 40.0)
    genome, masks = {"c": "ACGTCAGGG"}, {"c": np.zeros(10, bool)}
    at = np.arange(9)
    placed = pfam_truth.Placement("c", False, at)
    assert pfam_truth.true_bases(codon, placed, genome) == ("TCA", "TCA")
    assert pfam_truth.stop_class(codon, placed, genome, masks) == "error"
    clipped = replace(placed, sites=np.where(at < 4, -1, at))
    gapped = replace(placed, sites=np.where(at > 4, at + 1, at))
    for p in (clipped, gapped):
        assert pfam_truth.true_bases(codon, p, genome) is None
        assert pfam_truth.stop_class(codon, p, genome, masks) == "indel"
    assert pfam_truth.stop_class(codon, None, genome, masks) == "unaligned"
    masks["c"][4] = True
    assert pfam_truth.stop_class(codon, placed, genome, masks) == "variant"


def test_reads_are_written_for_an_aligner_and_one_anchor_is_dropped(tmp_path: Path) -> None:
    pairs = [
        ("ACGT", bytes([2, 12, 23, 37]), "GGCA", bytes([37] * 4)),
        ("TTTT", bytes([2] * 4), "AAAA", bytes([2] * 4)),
    ]
    assert pfam.write_pairs(pairs, tmp_path / "x") == 2
    back = list(pfam.read_pairs(tmp_path / "x_1.fastq.gz", tmp_path / "x_2.fastq.gz"))
    assert back == pairs
    (tmp_path / "staged_1.fastq.gz").symlink_to(tmp_path / "x_1.fastq.gz")  # as a workflow stages its input
    with pytest.raises(SystemExit):
        pfam.main(["reads", *(str(tmp_path / f"x_{i}.fastq.gz") for i in (1, 2)), "--output", str(tmp_path / "staged")])
    assert list(pfam.read_pairs(tmp_path / "x_1.fastq.gz", tmp_path / "x_2.fastq.gz")) == pairs
    ev = pfam.Evidence(1, (1, 1), 35)
    ev.certain[0][((1,), (1, 1, 1))] += 1
    ev.overlap[0][(1,)] += 1
    assert not pfam.only(ev, "overlap").certain and pfam.only(ev, "overlap").overlap
    assert not pfam.only(ev, "stops").overlap and pfam.only(ev, "stops").certain
