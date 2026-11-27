import json
from pathlib import Path

import numpy as np
import pytest
from test_pfam import CODONS, substitutions_only
from test_pfam_truth import Mock

from sequencing_error_model import pfam_locus, pfam_truth
from sequencing_error_model.sources import pfam

IX = pfam_locus.INDEX


def arrays(groups: list[list[str]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = [(g, IX[c], 100 * g + k) for g, codons in enumerate(groups, 1) for k, c in enumerate(codons)]
    group, codon, pair = np.array(rows).T
    return group, codon, pair


def test_consensus_leaves_the_codon_itself_out() -> None:
    group, codon, pair = arrays(
        [
            ["AAA", "AAA", "AAA"],  # each has two others that agree
            ["AAA", "AAA", "AAA", "AAG"],  # the odd one is judged; the rest see a disagreement
            ["AAA", "AAG"],  # one other: which erred is open
            ["AAA", "AAA", "AAG", "AAG"],
        ]
    )
    row, want = pfam_locus.consensus(group, codon, pair)
    assert sorted(row.tolist()) == [0, 1, 2, 6]
    assert (want == IX["AAA"]).all()
    # mates of one fragment are one vote
    row, _ = pfam_locus.consensus(group, codon, np.where(group == 1, 7, pair))
    assert 0 not in row and 6 in row
    with pytest.raises(ValueError, match="consensus"):
        pfam_locus.consensus(group, codon, pair, 1)


def test_recurring_alleles_are_counted_once_per_locus() -> None:
    group, codon, pair = arrays(
        [["AAA"] * 5 + ["AAG"] * 3, ["AAA"] * 5 + ["AAG"] * 2, ["CTT"] * 3 + ["CGT"] * 4, ["AAA"] * 3 + ["CCC"] * 3]
    )
    got = pfam_locus.recurring(group, codon, pair)
    want = np.zeros((12, 3))
    want[pfam_locus.OPS.index("A>G"), 0] = 1  # AAA > AAG, synonymous
    want[pfam_locus.OPS.index("G>T"), 1] = 1  # CGT > CTT: the commoner is the first allele
    assert (got == want).all()


def test_effects_and_opportunities_follow_the_genetic_code() -> None:
    assert pfam_locus.CHANCE[[i for i in range(64) if not pfam_locus.STOP[i]]].sum((1, 2)) == pytest.approx(9)
    tgg = pfam_locus.CHANCE[IX["TGG"]]
    assert tgg[:, 2].sum() == 2 and tgg[pfam_locus.OPS.index("G>A"), 2] == 2  # TAG, TGA
    assert tgg[:, 0].sum() == 0


def test_the_mixture_takes_the_biology_out() -> None:
    rng = np.random.default_rng(4)
    chance = np.tensordot(np.full(64, 4000.0) * ~pfam_locus.STOP, pfam_locus.CHANCE, 1)
    e, b, ratio = rng.uniform(2e-3, 6e-3, 12), rng.uniform(1e-2, 4e-2, 12), 0.1
    n = rng.poisson(chance * (e[:, None] + b[:, None] * np.array([1, ratio, 0])))
    again = rng.poisson(chance * 5e-3 * np.array([1, ratio, 0]))
    assert pfam_locus.omega(again, chance) == pytest.approx(ratio, rel=0.15)
    got = pfam_locus.mixture(n, chance, pfam_locus.omega(again, chance))
    assert pfam_locus._per_base(got, chance) == pytest.approx(pfam_locus._per_base(e, chance), rel=0.06)
    naive = n.sum(1) / chance.sum(1)
    assert pfam_locus._per_base(naive, chance) > 1.8 * pfam_locus._per_base(e, chance)


class Strains(Mock):
    """Six copies of each gene, and eight rare strains of it with mostly synonymous differences."""

    def __init__(self, seed: int) -> None:
        super().__init__(seed, genes=12)
        self.genes = [[*(6 * genes), *(self.strain(g) for g in genes for _ in range(8))] for genes in self.genes]

    def strain(self, gene: str) -> str:
        codons = [gene[i : i + 3] for i in range(0, len(gene), 3)]
        for k in self.rng.choice(len(codons), 6, replace=False):
            near = [c for c in pfam.CODE if sum(a != b for a, b in zip(c, codons[k], strict=True)) == 1]
            same = [c for c in near if c in CODONS[pfam.CODE[codons[k]]]]
            other = [c for c in near if pfam.CODE[c] not in ("*", pfam.CODE[codons[k]])]
            codons[k] = str(self.rng.choice(same if same and self.rng.random() < 0.9 else other))
        return "".join(codons)


def test_errors_are_told_from_a_strains_variants_end_to_end(tmp_path: Path) -> None:
    world = Strains(6)
    pairs, hits, _, _ = world.pairs(substitutions_only(0.005), 4000, 210, 100)
    world.sam(tmp_path / "a.sam", pairs, 210, set(range(0, 6000, 25)))
    rows = pfam_locus.table(pairs, hits, world.profiles, t1=30)
    blind = pfam_locus.report(rows, 1)
    assert "truth" not in blind
    pfam_locus.add_truth(rows, tmp_path / "a.sam", world.genome())
    assert 0.9 < (rows[:, pfam_locus.TRUE] >= 0).mean() < 0.99  # the reads left unaligned
    got = pfam_locus.report(rows, 1)
    assert got["labels"] == blind["labels"] and got["per_base"] == blind["per_base"]  # truth is only reported
    truth = got["truth"]
    labels = {k: np.array(x) for k, x in truth["labels"].items()}
    # the variants that are seen once are called biology by the genome under the read, the rest are errors
    assert labels["error"].sum() > 150 and labels["biology"].sum() > 300 and labels["other"].sum() <= 3
    assert labels["biology"][:, 2].sum() == 0 and np.array(got["recurring"])[:, 0].sum() > 300
    # naive labels are 2.6 times too many; the mixture leaves 7-17 % over three seeds, always high: recurring
    # alleles here have a lower omega than the singletons, which is the bias the module warns of
    true, rate = truth["per_base"], got["per_base"]
    assert rate["naive"] > 2 * true
    assert true < rate["mixture"] < 1.25 * true
    assert got["omega"] < truth["omega"]
    codons = rows[rows[:, pfam_locus.TRUE] >= 0]
    all_codons = pfam_locus.DIST[codons[:, pfam_locus.TRUE], codons[:, pfam_locus.CODON]].sum() / (3 * len(codons))
    assert true == pytest.approx(all_codons, rel=0.15)  # what clean flanks and agreement select


def test_the_command_writes_a_report(tmp_path: Path) -> None:
    rng = np.random.default_rng(0)
    profiles = pfam.Profiles(("PF00001",), np.array([0, 4]), rng.normal(size=(4, 21)).astype(np.float32))
    (tmp_path / "accs.txt").write_text("PF00001\n")
    np.array([0, 4], np.uint32).tofile(tmp_path / "offs.u32")
    profiles.scores.tofile(tmp_path / "profile.f32")
    pfam.write_pairs([("ACGT" * 10, bytes([30] * 40), "ACGT" * 10, bytes([30] * 40))], tmp_path / "x")
    (tmp_path / "p.tsv").write_text("orf\tacc\thmm_from\tali_from\tcigar\n")
    reads = [str(tmp_path / f"x_{i}.fastq.gz") for i in (1, 2)]
    args = [*reads, "--pairs", str(tmp_path / "p.tsv"), "--prep", str(tmp_path), "--output", str(tmp_path / "o.json")]
    assert pfam_locus.main(args) == 0  # nothing placed: an empty report, not a crash
    assert {x["codons"] for x in json.loads((tmp_path / "o.json").read_text())["loci"]} == {0}
    with pytest.raises(SystemExit):
        pfam_locus.main([*args, "--sam", "a.sam"])


def test_selection_counts_stops_on_keyed_codons_against_all() -> None:
    def rows(spec: list[tuple[int, str, int, str]]) -> np.ndarray:
        out = np.zeros((len(spec), pfam_locus.FIRST + len(pfam_locus.FLANKS)), np.int64)
        for i, (key, codon, pair, true) in enumerate(spec):
            out[i, [pfam_locus.CODON, pfam_locus.PAIR, pfam_locus.TRUE]] = IX[codon], pair, IX[true]
            out[i, pfam_locus.FIRST :] = key
        return out

    spec = [(11, "AAA", p, "AAA") for p in range(4)]  # keyed
    spec += [(11, "AAA", 0, "AAA")]  # the other mate of fragment 0: keyed with it
    # the error stop is keyed, its neighbours are not: they each have an other fragment that disagrees
    spec += [(44, "AAA", p, "AAA") for p in range(20, 23)] + [(44, "TAA", 23, "AAA")]
    spec += [(22, "TAA", p, "TAA") for p in range(5, 8)]  # the locus agrees on a stop: biology, in neither
    spec += [(33, "AAA", 8, "AAA"), (33, "TAA", 9, "AAA")]  # one other fragment: not keyed
    spec += [(0, "TAA" if p < 13 else "AAA", p, "AAA") for p in range(10, 20)]  # off a locus: three stops
    got = pfam_locus.selection(rows(spec), 1)
    assert (got["codons_keyed"], got["codons"], got["stops_keyed"], got["stops"]) == (6, 21, 1, 5)
    assert got["keyed_over_all"] == pytest.approx((1 / 6) / (5 / 21))
    assert got["truth"]["keyed_over_all"] == pytest.approx(got["keyed_over_all"])
    assert got["truth"]["error_stops_keyed_over_all"] == pytest.approx(got["keyed_over_all"])


def test_selection_and_masked_truth_end_to_end(tmp_path: Path) -> None:
    world = Strains(6)
    pairs, hits, _, _ = world.pairs(substitutions_only(0.005), 4000, 210, 100)
    world.sam(tmp_path / "a.sam", pairs, 210, set())
    rows = pfam_locus.table(pairs, hits, world.profiles, t1=30)
    genome = world.genome()
    masks = pfam_truth.sam_variant_masks(tmp_path / "a.sam", genome, 20, 0.2)
    assert not any(m.any() for m in masks.values())  # every read on its own strain: no variant against it
    pfam_locus.add_truth(rows, tmp_path / "a.sam", genome, 20, masks)
    got = pfam_locus.report(rows, 1)["selection"]
    # uniform errors: the keyed codons are as clean as the rest, by the stops and by the genomes
    assert abs(got["keyed_over_all"] - 1) < 2.5 * got["sd"]
    assert got["truth"]["keyed_over_all"] == pytest.approx(1, abs=0.1)
    everywhere = {c: np.ones(len(g) + 1, bool) for c, g in genome.items()}
    blank = rows.copy()
    blank[:, pfam_locus.TRUE] = -1
    pfam_locus.add_truth(blank, tmp_path / "a.sam", genome, 20, everywhere)
    assert (blank[:, pfam_locus.TRUE] < 0).all()
