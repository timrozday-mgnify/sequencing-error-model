from collections import Counter
from pathlib import Path

import numpy as np

from sequencing_error_model import generate as gen
from sequencing_error_model import recovery
from sequencing_error_model import variation as var


def genome(length: int, seed: int, name: str = "c1") -> list[tuple[str, str]]:
    rng = np.random.default_rng(seed)
    return [(name, "".join(rng.choice(list("ACGT"), size=length)))]


def population(
    n: int = 3, ani: float = 0.99, *, tree: str = "star", length: int = 20_000, seed: int = 0, **kwargs: float
) -> var.Population:
    rng = np.random.default_rng(seed)
    consensus = genome(length, seed)
    haps = var.haplotypes(consensus, n, ani, rng, tree=tree)
    minor = var.minor_alleles(consensus, float(kwargs.get("density", 0.0)), rng)
    return var.Population(tuple(consensus), haps, var.abundances(list(haps), [0.2]), minor)


def test_star_tree_hits_the_target_ani() -> None:
    # `ani` is the pairwise identity of the deepest pair, so root-to-tip divergence is half of 1 - ani.
    pop = population(n=3, ani=0.99)
    pairs = [var.ani(pop, a, b) for a in pop.haplotypes for b in pop.haplotypes if a < b]
    assert all(abs((1 - x) - 0.01) < 0.002 for x in pairs), pairs
    one = [var.ani(pop, a, a) for a in pop.haplotypes]
    assert one == [1.0, 1.0, 1.0]


def test_random_tree_shares_alleles_and_star_does_not() -> None:
    shared = {
        tree: sum(len(s.haplotypes) > 1 for s in population(n=4, tree=tree, seed=3).sites() if s.haplotypes)
        for tree in ("star", "random")
    }
    # A star tree shares an allele only by coincidence (two tips drawing the same mutation); descent is not
    # coincidence, so the random tree shares two orders of magnitude more.
    assert shared["star"] < 10 < 50 < shared["random"], shared


def test_transition_transversion_ratio_and_indels() -> None:
    pop = population(n=1, ani=0.98, length=40_000, seed=5)
    kinds = Counter[str]()
    for s in pop.sites():
        if len(s.ref) == len(s.alt) == 1:
            kinds["ti" if var._TRANSITION[s.ref] == s.alt else "tv"] += 1
        else:
            kinds["indel"] += 1
    assert 1.6 < kinds["ti"] / kinds["tv"] < 2.5, kinds
    # 10% of events are indels by default, and they are the only sites that change length.
    assert 0.05 < kinds["indel"] / sum(kinds.values()) < 0.16, kinds


def test_alleles_are_derived_from_the_sequence() -> None:
    # Consensus ACGTACGTAC: T>A at 3, GG inserted after it, consensus bases 5-6 (CG) deleted. Every row is
    # left-anchored as in VCF, so the insertion and the deletion carry the base before them.
    contig = var.Contig("h|c", "c", "ACGAGGATAC", np.array([0, 1, 2, 3, -1, -1, 4, 7, 8, 9], np.int64))
    assert var._alleles(contig, "ACGTACGTAC") == [(3, "T", "A"), (3, "T", "TGG"), (4, "ACG", "A")]


def test_truth_labels_every_read_base() -> None:
    pop = population(n=2, ani=0.98, length=30_000, seed=7, density=1e-3)
    model, rng = recovery.example_spec(), np.random.default_rng(11)
    frags = var.fragments(pop, "sample1", gen.insert_sizes(180, 20), 400, 100, rng)
    records = [(f, m) for f in frags for m in (0, 1)]
    reads = gen.generate(model, [f.templates[m] for f, m in records], [m + 1 for _, m in records], rng)
    counts = Counter[str]()
    for (frag, mate), read in zip(records, reads, strict=True):
        template, consensus = frag.templates[mate], frag.consensus[mate]
        labels = var.truth(template, consensus, read)
        assert len(labels) == len(read.sequence)
        for base, cons, label in zip(read.sequence, _read_consensus(template, consensus, read), labels, strict=True):
            # The exit criterion: a read base differing from the consensus is an error or a variant, and one
            # labelled a match agrees with the consensus.
            if cons not in (".", "-") and base != cons:
                assert label in ("error", "variant"), (base, cons, label)
            if label == "match":
                assert base == cons
        counts.update(labels)
    assert counts["variant"] / sum(counts.values()) > 0.005  # ~1% divergence plus minor alleles
    assert counts["error"] > 0 and counts["adapter"] < counts["match"]


def _read_consensus(template: str, consensus: str, read: gen.Read) -> list[str]:
    """The consensus base under each read base, so the test compares bases without reusing `truth`'s logic."""
    out = []
    for t, op in gen.align(template, read)[0]:
        if op.startswith("->"):
            out.append("-")
        elif not op.endswith("-"):
            out.append(consensus[t])
    return out


def test_minor_alleles_reach_their_frequency_in_both_mates() -> None:
    rng = np.random.default_rng(13)
    consensus = genome(5_000, 17)
    haps = var.haplotypes(consensus, 1, 1.0, rng)  # clonal: every difference is a minor allele
    minor = var.minor_alleles(consensus, 0.02, rng, frequency=(0.299, 0.301))
    pop = var.Population(tuple(consensus), haps, var.abundances(list(haps), [0.0]), minor)
    frags = var.fragments(pop, "sample1", gen.insert_sizes(150, 0), 600, 150, rng)
    carried = both = covered = 0
    for f in frags:
        r1, r2 = f.templates
        c1 = f.consensus[0]
        for o, (base, cons) in enumerate(zip(r1, c1, strict=True)):
            if cons in (".", "-") or (f.consensus[0][o] == "."):
                continue
            covered += 1
            if base != cons:
                carried += 1
                # The same molecule, so mate 2 reads the same allele at the mirrored offset.
                both += r2[len(r1) - 1 - o] == gen._revcomp(base)
    assert abs(carried / (covered * 0.02) - 0.3) < 0.06, (carried, covered)
    assert both == carried


def test_repeat_copy_is_divergent_and_appended() -> None:
    consensus = genome(4_000, 19)
    out, info = var.repeat_copy(consensus, 1_200, 0.97, np.random.default_rng(23))
    assert len(out[0][1]) == 4_000 + 1_200 and out[0][1][:4_000] == consensus[0][1]
    assert abs(info.identity - 0.97) < 0.02
    assert out[0][1][slice(*info.copy)] != consensus[0][1][slice(*info.source)]


def test_abundances_are_a_series_summing_to_one() -> None:
    a = var.abundances(["h1", "h2", "h3"], [0.5, 0.05])
    assert a["sample1"] == {"h1": 0.5, "h2": 0.25, "h3": 0.25}
    assert abs(sum(a["sample2"].values()) - 1) < 1e-12 and a["sample2"]["h1"] == 0.95


def test_cli_writes_the_population(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">c1 a contig\n{genome(6_000, 29)[0][1]}\n")
    out = tmp_path / "pop"
    assert (
        var.main(
            [
                str(fasta),
                "--output",
                str(out),
                "--haplotypes",
                "3",
                "--ani",
                "0.99",
                "--tree",
                "random",
                "--minor-fraction",
                "0.2",
                "0.05",
                "--minor-density",
                "1e-3",
                "--repeat",
                "800",
                "0.95",
            ]
        )
        == 0
    )
    sites = (out / "sites.tsv").read_text().splitlines()
    assert sites[0].split("\t") == ["contig", "position", "ref", "alt", "haplotypes", "sample1", "sample2"]
    assert len(sites) > 50
    minor = [r for r in (line.split("\t") for line in sites[1:]) if not r[4]]
    assert minor and all(0.01 <= float(r[5]) <= 0.5 for r in minor)
    header = (out / "abundance.tsv").read_text().splitlines()[0]
    assert header.split("\t") == ["haplotype", "sample1", "sample2"]
    # The repeat lengthened the consensus, so it is written beside the haplotypes.
    assert len(next(iter(gen._fasta((out / "consensus.fa").open())))[1]) == 6_800
    assert sorted(p.name for p in (out / "haplotypes").iterdir()) == ["hap1.fa", "hap2.fa", "hap3.fa"]
    repeat = (out / "repeats.tsv").read_text().splitlines()[1].split("\t")
    assert repeat[0] == "c1" and (int(repeat[3]), int(repeat[4])) == (6_000, 6_800) and float(repeat[5]) > 0.9
