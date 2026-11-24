import math
from collections import Counter
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from sequencing_error_model import compare, recovery
from sequencing_error_model.fit import error
from sequencing_error_model.generate import Read, _revcomp, generate, observations
from sequencing_error_model.observations import CountTable, Key
from sequencing_error_model.sources import pfam
from sequencing_error_model.sources.pfam import Hit, Profiles
from sequencing_error_model.spec import ErrorModelSpec

ERROR_TOKENS = ["QualityWindow(1)", "Context(0,0)", "Mate"]
QUALITY_TOKENS = ["QualityMarkov(1)", "Position(3)", "Mate", "Context(1,1)"]
COMPLEMENT = dict(zip("ACGT", "TGCA", strict=True))
LENGTH = 70  # columns per family
CODONS: dict[str, list[str]] = {a: [c for c, x in pfam.CODE.items() if x == a] for a in pfam.RESIDUES}
PENALTY = 8.0  # HMMER's domain score sits this far below the raw profile bits
Source = tuple[int, int, int, bool]  # family, gene start of the fragment, its length, strand


def substitutions_only(scale: float) -> ErrorModelSpec:
    truth = recovery.scale_error_rate(recovery.example_spec(), scale)
    head = truth.error_head[0]
    bias = head.params["bias"].copy()
    bias[5:] = -np.inf
    kept = (
        replace(head, params={**head.params, "bias": bias}),
        *(c for c in truth.error_head[1:] if c.name in ("Mate",)),
    )
    return replace(truth, error_head=kept)


class World:
    """Random profiles, and diverged genes of each family: biology that varies freely by column."""

    def __init__(self, seed: int, families: int = 10, genes: int = 40) -> None:
        self.rng = rng = np.random.default_rng(seed)
        top = rng.integers(0, 20, families * LENGTH)
        main = rng.uniform(0.35, 0.9, families * LENGTH)
        p = rng.dirichlet([0.5] * 20, families * LENGTH) * (1 - main)[:, None]
        p[np.arange(len(top)), top] += main
        p = p / p.sum(1, keepdims=True)
        self.profiles = Profiles(
            tuple(f"PF{i:05d}" for i in range(families)),
            np.arange(0, (families + 1) * LENGTH, LENGTH, dtype=np.int64),
            np.log2(np.maximum(p, 1e-4) / pfam.BACKGROUND),
        )
        self.genes = [[self.gene(p[f * LENGTH : (f + 1) * LENGTH]) for _ in range(genes)] for f in range(families)]

    def gene(self, p: np.ndarray) -> str:
        residues = [pfam.RESIDUES[self.rng.choice(20, p=row)] for row in p]
        return "".join(CODONS[a][self.rng.integers(len(CODONS[a]))] for a in residues)

    def fragments(self, n: int, insert: int) -> list[tuple[str, Source]]:
        out = []
        for _ in range(n):
            f = int(self.rng.integers(len(self.genes)))
            g = self.genes[f][int(self.rng.integers(len(self.genes[f])))]
            start = int(self.rng.integers(len(g) - insert + 1))
            strand = bool(self.rng.integers(2))
            frag = g[start : start + insert]
            out.append((_revcomp(frag) if strand else frag, (f, start, insert, strand)))
        return out

    def search(self, read: int, seq: str, source: Source, mate: int) -> list[Hit]:
        """HMMER's part: the best local stretch of every ORF on the true diagonal, cut well below t1."""
        f, start, insert, strand = source
        reverse = (mate == 2) != strand
        lo = start + insert - len(seq) if reverse else start  # the read covers the gene from here
        hits = []
        for name, protein in pfam.orfs(read, seq):
            _, sign, at = name.split("_")
            if (sign == "-") != reverse:
                continue
            pos = int(at)  # on the hit's own strand
            gene_at = lo + pos
            if gene_at % 3:
                continue
            k = gene_at // 3
            scores = np.array(
                [
                    self.profiles.scores[f * LENGTH + k + i, pfam.RESIDUES.index(a)]
                    if 0 <= k + i < LENGTH and a in pfam.RESIDUES
                    else pfam.UNALIGNED
                    for i, a in enumerate(protein)
                ]
            )
            best, span, run, first = 0.0, (0, 0), 0.0, 0
            for i, x in enumerate(scores):
                if run <= 0:
                    run, first = 0.0, i
                run += x
                if run > best:
                    best, span = run, (first, i + 1)
            if best - PENALTY >= 12:
                hits.append(Hit(f, reverse, pos + 3 * span[0], k + span[0], ((span[1] - span[0], "M"),)))
        return hits

    def pairs(
        self, truth: ErrorModelSpec, n: int, insert: int, length: int, dose: float = 0.0, seed: int = 5
    ) -> tuple[list[pfam.ReadPair], dict[int, list[Hit]], list[str], list[Read]]:
        made = self.fragments(n, insert)
        templates = [t for frag, _ in made for t in (frag[:length], _revcomp(frag)[:length])]
        rng = np.random.default_rng(seed)
        reads = generate(truth, templates, [1, 2] * n, rng)
        out: list[pfam.ReadPair] = []
        hits: dict[int, list[Hit]] = {}
        for i, (_, source) in enumerate(made):
            seqs = [pfam.add_errors(reads[2 * i + k].sequence, dose, (i, k)) for k in (0, 1)]
            q = [[ord(c) - 33 for c in reads[2 * i + k].quality] for k in (0, 1)]
            out.append((seqs[0], q[0], seqs[1], q[1]))
            for k in (0, 1):
                if found := self.search(2 * i + k, seqs[k], source, k + 1):
                    hits[2 * i + k] = found
        return out, hits, templates, reads


def op_rates(spec: ErrorModelSpec, ev: pfam.Evidence) -> dict[str, float]:
    """Mean P(op) by true base and target over the evidence's own rows, at their counts."""
    rows: Counter[Key] = Counter()
    for exposed in ev.certain.values():
        for (k, _), x in exposed.items():
            rows[pfam._plain(k)] += x
    for agreed in ev.overlap.values():
        for k, x in agreed.items():
            rows[pfam._plain(k)] += x
    names = pfam.fields(ev.m)
    t = CountTable("t", names, "base", True, rows, {"flank": ev.flank})
    p = error.probabilities(spec.error_head, spec.quality_alphabet, t)
    centre, n = np.array([str(k[names.index("context")])[ev.flank[0]] for k in rows]), np.array(list(rows.values()))
    return {
        f"{x}>{y}": float(np.average(p[centre == x, 1 + pfam.BASES.index(y)], weights=n[centre == x]))
        for x in pfam.BASES
        for y in pfam.BASES
        if x != y and (centre == x).any()
    }


def test_stops_see_ten_of_twelve_substitutions() -> None:
    seen = {f"{c[p]}>{b}" for c in pfam.CODE if c not in pfam.STOPS for p, b, to in pfam.changes(c) if to in pfam.STOPS}
    flip = {"".join(COMPLEMENT.get(x, x) for x in op) for op in seen}
    assert {f"{a}>{b}" for a in "ACGT" for b in "ACGT" if a != b} - (seen | flip) == set(pfam.UNIDENTIFIED_BY_STOPS)


def test_path_follows_the_cigar() -> None:
    # 2M 1I 2M 2D 1M from codon 1, column 10
    columns, gaps, inserted = pfam._path(Hit(0, False, 3, 10, ((2, "M"), (1, "I"), (2, "M"), (2, "D"), (1, "M"))), 8)
    assert columns.tolist() == [9, 10, 11, -1, 12, 13, 16, 17]
    assert inserted.tolist() == [False] * 3 + [True] + [False] * 4
    assert gaps.tolist() == [0, 0, 0, -6, 0, 0, -7, 0]


def test_add_errors_is_reproducible_and_at_the_dose() -> None:
    seq = "ACGT" * 25_000
    a = pfam.add_errors(seq, 0.01, (1, 2))
    assert a == pfam.add_errors(seq, 0.01, (1, 2)) != pfam.add_errors(seq, 0.01, (1, 3))
    assert sum(x != y for x, y in zip(a, seq, strict=True)) / len(seq) == pytest.approx(0.01, rel=0.1)
    assert pfam.add_errors("NNNN", 1.0, (0, 0)) == "NNNN"


def test_a_codon_cannot_decide_its_own_exposure() -> None:
    """Replace an exposed codon by a stop, or by anything, and search again: it is exposed all the same."""
    world = World(1)
    truth = substitutions_only(0.01)
    pairs, hits, _, _ = world.pairs(truth, 60, 150, 120)
    made = World(1)
    sources = [s for _, s in made.fragments(60, 150)]
    rng, checked = np.random.default_rng(0), 0
    for i, (seq, *_) in enumerate(pairs):
        for c in pfam.exposed(seq, hits.get(2 * i, []), world.profiles, 30, 6)[::5]:
            for codon in ("TAA", "TGA", "".join(rng.choice(list("ACGT"), 3))):
                bases = list(seq)
                for at, b in zip(c.positions, codon, strict=True):
                    bases[at] = COMPLEMENT[b] if c.reverse else b
                changed = "".join(bases)
                again = pfam.exposed(changed, world.search(2 * i, changed, sources[i], 1), world.profiles, 30, 6)
                assert any(x.positions == c.positions and x.column == c.column for x in again), (i, c, codon)
                checked += 1
    assert checked > 100


@pytest.fixture(scope="module")
def stop_evidence() -> tuple[World, pfam.Evidence, ErrorModelSpec]:
    world, truth = World(2), substitutions_only(0.1)
    pairs, hits, _, _ = world.pairs(truth, 6000, 210, 100)  # 10 bp short of overlapping
    return world, pfam.collect(pairs, hits, world.profiles, 1, (1, 1), t1=30, quality_every=4), truth


def test_stops_recover_the_rate_and_the_identified_ops(
    stop_evidence: tuple[World, pfam.Evidence, ErrorModelSpec],
) -> None:
    world, ev, truth = stop_evidence
    assert ev.stats.stops > 150 and not ev.overlap
    spec, table = pfam.fit(ev, world.profiles, ERROR_TOKENS, QUALITY_TOKENS)
    got, want = op_rates(spec, ev), op_rates(truth, ev)
    ops = [op for op in want if op not in pfam.UNIDENTIFIED_BY_STOPS]
    sd = 1 / math.sqrt(ev.stats.stops)
    assert sum(got[o] for o in ops) == pytest.approx(sum(want[o] for o in ops), rel=max(3 * sd, 0.1))
    assert spec.provenance["unidentified_ops"] == ["A>G", "T>C"]
    assert (spec.error_head[0].params["bias"][5:] == -np.inf).all()
    assert table.meta["flank"] == (1, 1)
    # per-op: ten rates on a few hundred stops, so each within a factor of 2.5 and unbiased together
    ratios = np.array([got[o] / want[o] for o in ops])
    assert np.all((ratios > 0.4) & (ratios < 2.5)), dict(zip(ops, ratios, strict=True))
    assert float(np.exp(np.mean(np.log(ratios)))) == pytest.approx(1, abs=0.25)


def test_standard_additions_come_back_end_to_end() -> None:
    """Added errors, with the search rerun on them, raise the crude stop rate by the dose."""
    truth = substitutions_only(0.05)
    rates = []
    for dose in (0.0, 0.02):
        world = World(3)
        pairs, hits, _, _ = world.pairs(truth, 1500, 210, 100, dose=dose)
        ev = pfam.collect(pairs, hits, world.profiles, quality_every=10_000)
        rates.append(ev.crude_rate(30))
    added = rates[1][0] - rates[0][0]
    assert added == pytest.approx(0.02, rel=3 / math.sqrt(rates[1][1]) + 0.05), rates


def test_mate_overlap_identifies_all_twelve() -> None:
    world, truth = World(4), substitutions_only(1.0)
    pairs, hits, _, _ = world.pairs(truth, 1200, 130, 100)
    ev = pfam.collect(pairs, hits, world.profiles, 1, (1, 1), min_outside=15, max_outside=0.2, quality_every=2)
    assert ev.stats.overlaps > 300 and ev.stats.disagree > 100
    spec, _ = pfam.fit(ev, world.profiles, ERROR_TOKENS, QUALITY_TOKENS)
    got, want = op_rates(spec, ev), op_rates(truth, ev)
    assert spec.provenance["unidentified_ops"] == []
    assert sum(got.values()) == pytest.approx(sum(want.values()), rel=0.2)
    for op in pfam.UNIDENTIFIED_BY_STOPS:
        assert 0.4 < got[op] / want[op] < 2.5, op


def test_the_mate_that_made_a_stop_is_the_one_that_erred() -> None:
    world, truth = World(4), substitutions_only(1.0)
    pairs, hits, _, _ = world.pairs(truth, 1200, 130, 100)
    ev = pfam.collect(pairs, hits, world.profiles, 1, (1, 1), min_outside=15, max_outside=0.2, quality_every=2)
    assert ev.stats.stop_one_mate > 0
    prior = pfam.Prior(ev, world.profiles)
    assert 0 < prior.lam <= max(pfam.LAMBDAS)
    for _, column, one, two in ev.disputed:
        if one in pfam.STOPS and two not in pfam.STOPS:
            assert prior.weight(column, one, 1) < 0.01 * prior.weight(column, two, 1)


def test_overlap_needs_the_bases_outside_to_agree() -> None:
    world, truth = World(5), substitutions_only(0.05)
    pairs, hits, _, _ = world.pairs(truth, 200, 130, 100)
    scrambled = [
        (a, q, b[: len(b) - 30] + "".join(np.random.default_rng(i).choice(list("ACGT"), 30)), r)
        for i, (a, q, b, r) in enumerate(pairs)
    ]
    assert pfam.collect(pairs, hits, world.profiles).stats.overlaps > 30
    bad = pfam.collect(scrambled, hits, world.profiles, min_outside=10).stats
    assert bad.overlap_rejected > bad.overlaps


def test_families_whose_mates_disagree_too_often_lose_their_overlaps() -> None:
    """Tandem repeats: the mates of one family read different loci, so it must not set the rate."""
    world, truth = World(8), substitutions_only(0.05)
    pairs, hits, _, _ = world.pairs(truth, 1500, 130, 100)
    ev = pfam.collect(pairs, hits, world.profiles, min_outside=15, max_outside=0.2, quality_every=5)
    ev.disputed += [(ev.disputed[0][0], (3, 10), *ev.disputed[0][2:])] * 400
    spec, _ = pfam.fit(ev, world.profiles, ERROR_TOKENS, QUALITY_TOKENS)
    assert list(spec.provenance["downweighted_overlaps"]) == ["PF00003"]


def test_a_stop_in_both_mates_is_the_floor_not_an_error() -> None:
    """Pseudogenes: a family whose genes carry a stop. Both mates read it, and its family is weighted down."""
    world, truth = World(6), substitutions_only(0.1)
    broken = world.genes[0]
    for i, g in enumerate(broken):
        broken[i] = g[: 3 * 35] + "TAA" + g[3 * 36 :]
    pairs, hits, _, _ = world.pairs(truth, 2500, 130, 100)
    ev = pfam.collect(pairs, hits, world.profiles, min_outside=15, max_outside=0.2, quality_every=5)
    assert ev.stats.stop_both_mates > 10
    spec, _ = pfam.fit(ev, world.profiles, ERROR_TOKENS, QUALITY_TOKENS)
    floor = spec.provenance["floor"]
    assert floor["rate"] > 0 and floor["upper95"] > floor["rate"]
    assert "PF00000" in spec.provenance["downweighted_families"]
    want = op_rates(truth, ev)
    assert sum(op_rates(spec, ev).values()) == pytest.approx(sum(want.values()), rel=0.4)


def test_family_weights_cost_little_when_there_is_nothing_to_find() -> None:
    rng = np.random.default_rng(0)
    size = dict(enumerate(rng.lognormal(6, 1.5, 3000)))
    for rate in (5e-5, 5e-3):  # Illumina-like, and enough stops per family for a small misfit to show
        errors = {f: n * rate for f, n in size.items()}
        stops = {f: int(rng.poisson(e * rng.gamma(25, 1 / 25))) for f, e in errors.items()}  # off by 20 % per family
        w = pfam.family_weights(stops, size, errors)
        kept = sum(n * w.get(f, 1.0) for f, n in stops.items()) / sum(stops.values())
        assert kept > 0.95, rate
    found = pfam.family_weights({**stops, 0: 60}, {**size, 0: 1000.0}, {**errors, 0: 1.0})
    assert found[0] < 0.01


def test_exposure_counts_only_codons_that_are_probably_clean() -> None:
    ev = pfam.Evidence(0, (0, 0), 30)
    row = (2, "A", 1, 9, 1, "G", "=")
    ev.certain[0].update({(row, (2, 2, 37)): 10})
    w = pfam._Weights({0: 0.5}, [], np.zeros(0), {2: 0.8}, {})
    assert pfam.expected(ev, w).counts[row] == pytest.approx(10 * 0.5 * 0.8 * 0.8)


def test_plateau_flags_chance_homology() -> None:
    ev = pfam.Evidence(1, (0, 0), 25)
    ev.by_anchor.update({(0, 27): [400, 100_000], (1, 45): [40, 100_000]})
    assert ev.check_plateau() is not None
    assert ev.check_plateau({0: 0.0}) is None  # the family behind it is already weighted down
    clean = pfam.Evidence(1, (0, 0), 25)
    clean.by_anchor.update({(0, 27): [42, 100_000], (0, 45): [40, 100_000]})
    assert clean.check_plateau() is None
    assert set(clean.plateau()) == {"30", "35", "40", "45"}


def test_orf_names_and_pairs_table_round_trip(tmp_path: Path) -> None:
    world = World(7)
    gene = world.genes[0][0]
    for seq, reverse in ((gene[6:126], False), (_revcomp(gene)[6:126], True)):
        (name, protein), *_ = sorted(pfam.orfs(3, seq), key=lambda x: -len(x[1]))
        assert name.startswith("R3_")
        table = tmp_path / "pairs.tsv"
        table.write_text(f"orf\tacc\thmm_from\tali_from\tcigar\n{name}\tPF00000\t3\t2\t{len(protein) - 1}M\n")
        hits = pfam.read_hits(table, world.profiles)
        assert hits[3][0].reverse is reverse and hits[3][0].column == 2
        assert not pfam.read_hits(table, world.profiles, exclude={"PF00000"})
        with pytest.raises(ValueError, match="another Pfam release"):
            pfam.read_hits(table, replace(world.profiles, names=("PF99999",) * len(world.profiles.names)))
    (tmp_path / "bad.tsv").write_text("orf\tacc\thmm_from\tali_from\tcigar\nR3_0\tPF00000\t3\t2\t9M\n")
    with pytest.raises(ValueError, match="not named by this module"):
        pfam.read_hits(tmp_path / "bad.tsv", world.profiles)


def test_op_rates_compare_a_pfam_table_with_exact_tuples(
    stop_evidence: tuple[World, pfam.Evidence, ErrorModelSpec],
) -> None:
    """Restricted to the ops each row can show, so the two tables' rates are of the same thing."""
    world, ev, truth = stop_evidence
    _, table = pfam.fit(ev, world.profiles, ERROR_TOKENS, QUALITY_TOKENS)
    w = World(2)
    _, _, templates, reads = w.pairs(truth, 6000, 210, 100)
    exact = observations(zip(templates, reads, [1, 2] * 6000, strict=True), (1, 1), 1, "truth")
    a, b = compare.op_rates(table), compare.op_rates(exact)
    assert set(a) == set(b) - set(pfam.UNIDENTIFIED_BY_STOPS)
    pooled = [sum(r[op]["rate"] for op in a) for r in (a, b)]
    assert pooled[0] == pytest.approx(pooled[1], rel=0.35)
