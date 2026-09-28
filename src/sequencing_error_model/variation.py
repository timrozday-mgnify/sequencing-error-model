"""Test-only biological variation simulator (plan §9 phase 7): a base genome → strain haplotypes, minor
alleles, divergent repeat copies, and per-base truth.

It exists so phase 8's separation methods have truth for every base: a read base differs from the consensus
because of a sequencing error *or* because the molecule really differs there (§6.6). Nothing here is part of the
error model; it is the test rig that says which is which.

- `haplotypes` mutates the consensus down a tree: `star` (every haplotype an independent tip) or `random` (a
  random ultrametric topology, so haplotypes share alleles along internal branches, which is what linkage
  methods need). `ani` is the *pairwise* identity of the deepest pair, so root-to-tip divergence is
  `(1 - ani) / 2`. Events are substitutions (by `ti_tv`) and short indels (`indel_fraction`, geometric length).
- `minor_alleles` adds within-population substitutions that no strain owns: each sits at a consensus position
  with a frequency drawn log-uniformly over `frequency`, which is the neutral 1/f spectrum. They are realised
  per fragment (`fragments`), so both mates of a pair carry the same allele - a real molecule - and
  `pe-overlap` stays blind to them (§6.6).
- `repeat_copy` appends a diverged copy of a window to its contig: the ERR10889147 failure (two copies of a
  repeat longer than the read, phase 4), which is genome structure, not a variant, so it goes in the consensus
  before haplotyping and is reported separately.
- `Population.sites` is the truth-site table, *derived* from the haplotype sequences rather than recorded as
  they are built, so it cannot drift from them.
- `fragments` samples molecules by abundance and `truth` labels every read base match / variant / error /
  adapter from the CIGAR and the fragment's consensus counterpart. Mate conventions, adapter read-through and
  the fragment name format follow `generate.fragments`.

Deferred: recombination (plan: added once phase 8's linkage test needs a hard case); external haplotypes from a
reference plus VCF; indel minor alleles; per-sample minor-allele frequencies.
"""

import argparse
import sys
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from sequencing_error_model.fit.quality import Array
from sequencing_error_model.generate import ADAPTERS, Read, _fasta, _open, _revcomp, align

_TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}
_BASES = ("A", "C", "G", "T")


@dataclass(frozen=True)
class Contig:
    """One haplotype contig: its sequence, plus the consensus position of every base (-1 if inserted)."""

    name: str  # "<haplotype>|<consensus contig>", so a read name carries its haplotype
    reference: str  # the consensus contig it descends from
    sequence: str
    source: Array  # int64, one consensus index per base; -1 for an inserted base


@dataclass(frozen=True)
class Site:
    """One truth row: a consensus position and an alternative allele, left-anchored as in VCF."""

    contig: str
    position: int  # 0-based in the consensus
    ref: str
    alt: str
    haplotypes: tuple[str, ...]  # carriers; empty for a minor allele, which no strain owns
    frequency: dict[str, float]  # per sample


@dataclass(frozen=True)
class Repeat:
    """A divergent repeat copy: structure a read cannot resolve, not an allele."""

    contig: str
    source: tuple[int, int]  # the window that was copied, 0-based half-open
    copy: tuple[int, int]  # where the diverged copy sits
    identity: float  # achieved, not requested


@dataclass(frozen=True)
class Population:
    consensus: tuple[tuple[str, str], ...]
    haplotypes: dict[str, tuple[Contig, ...]]
    abundance: dict[str, dict[str, float]]  # sample -> haplotype -> fraction of molecules
    minor: dict[tuple[str, int], tuple[str, float]] = field(default_factory=dict)  # (contig, pos) -> (alt, freq)
    repeats: tuple[Repeat, ...] = ()

    @property
    def samples(self) -> tuple[str, ...]:
        return tuple(self.abundance)

    @property
    def sequences(self) -> dict[str, str]:
        return dict(self.consensus)

    def sites(self) -> list[Site]:
        """The truth-site table, derived from the haplotype sequences (so it always agrees with them)."""
        carriers: dict[tuple[str, int, str, str], list[str]] = {}
        for hap, contigs in self.haplotypes.items():
            for contig in contigs:
                for pos, ref, alt in _alleles(contig, self.sequences[contig.reference]):
                    carriers.setdefault((contig.reference, pos, ref, alt), []).append(hap)
        rows = [
            Site(
                c, pos, ref, alt, tuple(haps), {s: sum(a.get(h, 0.0) for h in haps) for s, a in self.abundance.items()}
            )
            for (c, pos, ref, alt), haps in carriers.items()
        ]
        rows += [
            Site(c, pos, self.sequences[c][pos], alt, (), dict.fromkeys(self.abundance, freq))
            for (c, pos), (alt, freq) in self.minor.items()
        ]
        return sorted(rows, key=lambda s: (s.contig, s.position, s.alt))


def _substitute(base: str, ti_tv: float, rng: np.random.Generator) -> str:
    """One substitution: `ti_tv` is the ratio of transitions to *all* transversions, so P(ti) = r / (1 + r)."""
    if rng.random() < ti_tv / (ti_tv + 1):
        return _TRANSITION[base]
    return str(rng.choice([b for b in _BASES if b != base and b != _TRANSITION[base]]))


def _mutate(
    sequence: str,
    source: Array,
    rate: float,
    rng: np.random.Generator,
    *,
    ti_tv: float,
    indel_fraction: float,
    max_indel: int,
) -> tuple[str, Array]:
    """One branch: each base draws an event with probability `rate` (a substitution, or a short indel)."""
    events = rng.random(len(sequence)) < rate
    bases: list[str] = []
    src: list[int] = []
    skip = 0
    for i, base in enumerate(sequence):
        if skip:
            skip -= 1
            continue
        if not events[i] or base not in _BASES:
            bases.append(base)
            src.append(int(source[i]))
            continue
        if rng.random() < indel_fraction:
            length = min(int(rng.geometric(0.5)), max_indel)
            if rng.random() < 0.5:  # deletion of `length` bases from here
                skip = length - 1
                continue
            bases += [base, *rng.choice(_BASES, size=length)]  # insertion after this base
            src += [int(source[i]), *([-1] * length)]
            continue
        bases.append(_substitute(base, ti_tv, rng))
        src.append(int(source[i]))
    return "".join(bases), np.array(src, np.int64)


def _alleles(contig: Contig, consensus: str) -> list[tuple[int, str, str]]:
    """(position, ref, alt) per difference between a haplotype contig and its consensus, left-anchored."""
    out: list[tuple[int, str, str]] = []
    inserted: list[str] = []
    prev = -1
    for base, s in zip(contig.sequence, contig.source.tolist(), strict=True):
        if s < 0:
            inserted.append(base)
            continue
        if s > prev + 1:  # consensus bases prev+1 .. s-1 are deleted
            out.append((prev, consensus[prev:s], consensus[prev]) if prev >= 0 else (0, consensus[:s], ""))
        if inserted:
            ins = "".join(inserted)
            out.append((prev, consensus[prev], consensus[prev] + ins) if prev >= 0 else (0, "", ins))
            inserted = []
        if base != consensus[s]:
            out.append((s, consensus[s], base))
        prev = s
    if prev + 1 < len(consensus):  # trailing deletion
        out.append((prev, consensus[prev:], consensus[prev]) if prev >= 0 else (0, consensus, ""))
    if inserted:
        ins = "".join(inserted)
        out.append((prev, consensus[prev], consensus[prev] + ins) if prev >= 0 else (0, "", ins))
    return out


def haplotypes(
    consensus: Sequence[tuple[str, str]],
    n: int,
    ani: float,
    rng: np.random.Generator,
    *,
    tree: str = "star",
    ti_tv: float = 2.0,
    indel_fraction: float = 0.1,
    max_indel: int = 6,
    prefix: str = "hap",
) -> dict[str, tuple[Contig, ...]]:
    """`n` haplotypes at pairwise identity `ani` for the deepest pair (root-to-tip divergence `(1 - ani) / 2`).

    `star`: every haplotype mutated from the consensus independently. `random`: a random ultrametric binary
    topology, so haplotypes share the alleles of the internal branches they descend from.
    """
    if n < 1 or not 0 < ani <= 1 or tree not in ("star", "random"):
        raise ValueError("need n >= 1, 0 < ani <= 1 and tree in {'star', 'random'}")
    if not 0 <= indel_fraction <= 1 or max_indel < 1 or ti_tv <= 0:
        raise ValueError("need 0 <= indel_fraction <= 1, max_indel >= 1 and ti_tv > 0")
    depth = (1 - ani) / 2
    names = [f"{prefix}{i + 1}" for i in range(n)]
    out: dict[str, tuple[Contig, ...]] = {name: () for name in names}

    def branch(seq: str, src: Array, rate: float) -> tuple[str, Array]:
        return _mutate(seq, src, rate, rng, ti_tv=ti_tv, indel_fraction=indel_fraction, max_indel=max_indel)

    for ref, seq in consensus:
        root = np.arange(len(seq), dtype=np.int64)

        def descend(tips: list[str], s: str, src: Array, left: float, ref: str = ref) -> None:
            """A random binary split: this branch's mutations are shared by every tip below it."""
            if len(tips) == 1:
                s, src = branch(s, src, left)
                out[tips[0]] += (Contig(f"{tips[0]}|{ref}", ref, s, src),)
                return
            step = left * float(rng.uniform(0.2, 0.8))
            s, src = branch(s, src, step)
            order = [str(t) for t in rng.permutation(tips)]
            cut = int(rng.integers(1, len(tips)))
            for half in (order[:cut], order[cut:]):
                descend(half, s, src.copy(), left - step)

        if tree == "star":
            for name in names:
                s, src = branch(seq.upper(), root.copy(), depth)
                out[name] += (Contig(f"{name}|{ref}", ref, s, src),)
        else:
            descend(names, seq.upper(), root, depth)
    return out


def minor_alleles(
    consensus: Sequence[tuple[str, str]],
    density: float,
    rng: np.random.Generator,
    *,
    frequency: tuple[float, float] = (0.01, 0.5),
    ti_tv: float = 2.0,
) -> dict[tuple[str, int], tuple[str, float]]:
    """Within-population substitutions at `density` per base, frequencies log-uniform over `frequency` (1/f)."""
    lo, hi = frequency
    if not 0 <= density <= 1 or not 0 < lo <= hi < 1:
        raise ValueError("need 0 <= density <= 1 and 0 < frequency low <= high < 1")
    out: dict[tuple[str, int], tuple[str, float]] = {}
    for name, seq in consensus:
        seq = seq.upper()
        hits = np.flatnonzero(rng.random(len(seq)) < density)
        for pos in hits.tolist():
            if seq[pos] not in _BASES:
                continue
            f = float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
            out[(name, pos)] = (_substitute(seq[pos], ti_tv, rng), f)
    return out


def repeat_copy(
    consensus: Sequence[tuple[str, str]],
    length: int,
    identity: float,
    rng: np.random.Generator,
    *,
    ti_tv: float = 2.0,
) -> tuple[list[tuple[str, str]], Repeat]:
    """Append a diverged copy of a random window to its own contig: two copies a read cannot tell apart.

    Structure, not variation: it goes into the consensus before haplotyping, so every haplotype carries both
    copies and the truth-site table stays about alleles.
    """
    out = [(name, seq.upper()) for name, seq in consensus]
    fits = [i for i, (_, seq) in enumerate(out) if len(seq) >= length]
    if length < 1 or not 0 < identity <= 1 or not fits:
        raise ValueError(f"need 1 <= length <= the longest contig and 0 < identity <= 1, got {length}, {identity}")
    i = int(rng.choice(fits))
    name, seq = out[i]
    start = int(rng.integers(len(seq) - length + 1))
    window = seq[start : start + length]
    copy, _ = _mutate(
        window, np.arange(length, dtype=np.int64), 1 - identity, rng, ti_tv=ti_tv, indel_fraction=0.0, max_indel=1
    )
    out[i] = (name, seq + copy)
    identity = sum(a == b for a, b in zip(window, copy, strict=True)) / length
    return out, Repeat(name, (start, start + length), (len(seq), len(seq) + len(copy)), identity)


def abundances(names: Sequence[str], minor: Sequence[float]) -> dict[str, dict[str, float]]:
    """One sample per minor-strain fraction: `names[0]` takes `1 - m`, the rest split `m` equally."""
    if not names or any(not 0 <= m < 1 for m in minor) or not minor:
        raise ValueError("need at least one haplotype and minor fractions in [0, 1)")
    rest = names[1:]
    return {
        f"sample{i + 1}": {names[0]: 1.0 - m if rest else 1.0, **{n: m / len(rest) for n in rest}}
        for i, m in enumerate(minor)
    }


@dataclass(frozen=True)
class Fragment:
    """One molecule: the templates the generator reads, and the consensus counterpart of each base."""

    name: str
    templates: tuple[str, str]  # mate 1, mate 2 (reverse complement), adapter and N padded as `generate.fragments`
    consensus: tuple[str, str]  # consensus base per template base; "-" where the haplotype inserted, "." past it
    haplotype: str


def _counterpart(source: Array, consensus: str) -> str:
    return "".join(consensus[s] if s >= 0 else "-" for s in source.tolist())


def fragments(
    population: Population,
    sample: str,
    insert_size: Array,
    n: int,
    read_length: int,
    rng: np.random.Generator,
    first: int = 0,
) -> list[Fragment]:
    """`n` molecules drawn by haplotype abundance, uniform over the placements that fit, with minor alleles.

    ponytail: minor alleles are drawn per fragment, so they are linked within a molecule but not across them.
    A per-molecule haplotype draw is the upgrade path when phase 8's linkage test needs one.
    """
    sizes, p = np.asarray(insert_size[0], np.int64), np.asarray(insert_size[1], float)
    abundance = population.abundance[sample]
    keys = [(hap, c) for hap, contigs in population.haplotypes.items() for c in contigs if abundance.get(hap)]
    if read_length < 1 or not len(sizes) or sizes.min() < 1 or p.min() < 0 or not p.sum() or not keys:
        raise ValueError("need read_length >= 1, an insert_size marginal of sizes >= 1, and an abundant haplotype")
    lengths = np.array([len(c.sequence) for _, c in keys], np.int64)
    weight = np.array([abundance[hap] for hap, _ in keys], float)
    sequences = population.sequences
    out = []
    for i, size in enumerate(rng.choice(sizes, size=n, p=p / p.sum())):
        slots = np.clip(lengths - size + 1, 0, None) * weight
        if not slots.sum():
            raise ValueError(f"insert size {size} is longer than every abundant contig")
        j = int(rng.choice(len(keys), p=slots / slots.sum()))
        hap, contig = keys[j]
        start = int(rng.integers(lengths[j] - size + 1))
        source = contig.source[start : start + size]
        frag = _minor(population, contig, source, contig.sequence[start : start + size].upper(), sample, rng)
        cons = _counterpart(source, sequences[contig.reference])
        r1, r2 = (s[:read_length].ljust(read_length, "N") for s in (frag + ADAPTERS[0], _revcomp(frag) + ADAPTERS[1]))
        c1, c2 = (s[:read_length].ljust(read_length, ".") for s in (cons, _revcomp(cons)))
        name = f"{contig.name}:{start + 1}-{start + size}#{first + i}"
        out.append(Fragment(name, (r1, r2), (c1, c2), hap))
    return out


def _minor(
    population: Population, contig: Contig, source: Array, frag: str, sample: str, rng: np.random.Generator
) -> str:
    """Flip this molecule's minor-allele sites, where the haplotype still carries the consensus allele."""
    if not population.minor:
        return frag
    consensus = population.sequences[contig.reference]
    bases = list(frag)
    for o, s in enumerate(source.tolist()):
        hit = population.minor.get((contig.reference, s)) if s >= 0 else None
        if hit and bases[o] == consensus[s] and rng.random() < hit[1]:
            bases[o] = hit[0]
    return "".join(bases)


def truth(template: str, consensus: str, read: Read) -> list[str]:
    """One label per read base: "error" (the read differs from its molecule), "variant" (the molecule differs
    from the consensus), "adapter" (no consensus base) or "match".

    So every read base that differs from the consensus is an error or a variant, and the two are told apart by
    the CIGAR, not by the bases.
    """
    labels: list[str] = []
    for t, op in align(template, read)[0]:
        if op.startswith("->"):  # an inserted read base: the generator put it there
            labels.append("error")
        elif op.endswith("-"):  # a deleted template base emits no read base
            continue
        elif op != "=":
            labels.append("error")
        elif consensus[t] == ".":
            labels.append("adapter")
        else:
            labels.append("variant" if template[t].upper() != consensus[t] else "match")
    return labels


def ani(population: Population, a: str, b: str) -> float:
    """Identity of two haplotypes: 1 - differing positions per consensus base (an indel counts as one)."""
    total = sum(len(seq) for _, seq in population.consensus)
    alleles: dict[tuple[str, int], dict[str, str]] = {}
    for s in population.sites():
        if s.haplotypes:
            at = alleles.setdefault((s.contig, s.position), {})
            for hap in s.haplotypes:
                at[hap] = s.alt
    differ = sum(at.get(a, "=") != at.get(b, "=") for at in alleles.values())
    return 1 - differ / total if total else 1.0


def write(population: Population, directory: Path) -> None:
    """Haplotype FASTAs, the consensus (a repeat copy changes it), the truth sites and the abundances."""
    (directory / "haplotypes").mkdir(parents=True, exist_ok=True)
    with _open(directory / "consensus.fa", "w") as out:
        for name, seq in population.consensus:
            out.write(f">{name}\n{seq}\n")
    for hap, contigs in population.haplotypes.items():
        with _open(directory / "haplotypes" / f"{hap}.fa", "w") as out:
            for contig in contigs:
                out.write(f">{contig.name}\n{contig.sequence}\n")
    samples = population.samples
    with _open(directory / "sites.tsv", "w") as out:
        out.write("\t".join(("contig", "position", "ref", "alt", "haplotypes", *samples)) + "\n")
        for s in population.sites():
            freq = (f"{s.frequency.get(name, 0.0):.6g}" for name in samples)
            out.write("\t".join((s.contig, str(s.position + 1), s.ref, s.alt, ",".join(s.haplotypes), *freq)) + "\n")
    if population.repeats:
        with _open(directory / "repeats.tsv", "w") as out:
            out.write("\t".join(("contig", "source_start", "source_end", "copy_start", "copy_end", "identity")) + "\n")
            for r in population.repeats:
                out.write("\t".join((r.contig, *map(str, (*r.source, *r.copy)), f"{r.identity:.6g}")) + "\n")
    with _open(directory / "abundance.tsv", "w") as out:
        out.write("\t".join(("haplotype", *samples)) + "\n")
        for hap in population.haplotypes:
            out.write("\t".join((hap, *(f"{population.abundance[n].get(hap, 0.0):.6g}" for n in samples))) + "\n")


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="sem-variation", description="Simulate strain haplotypes, minor alleles and divergent repeats."
    )
    p.add_argument("genome", type=Path, help="base genome FASTA[.gz] (the consensus)")
    p.add_argument("--output", type=Path, required=True, metavar="DIR")
    p.add_argument("--haplotypes", type=int, default=2)
    p.add_argument("--ani", type=float, default=0.99, help="pairwise identity of the deepest haplotype pair")
    p.add_argument("--tree", choices=("star", "random"), default="star")
    p.add_argument("--ti-tv", type=float, default=2.0)
    p.add_argument("--indel-fraction", type=float, default=0.1, help="share of events that are short indels")
    p.add_argument("--minor-fraction", type=float, nargs="+", default=[0.2], help="one sample per fraction")
    p.add_argument("--minor-density", type=float, default=0.0, help="within-population alleles per base")
    p.add_argument("--minor-frequency", type=float, nargs=2, default=[0.01, 0.5], metavar=("LOW", "HIGH"))
    p.add_argument("--repeat", type=float, nargs=2, metavar=("LENGTH", "IDENTITY"), help="add a divergent repeat")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args(argv)
    rng = np.random.default_rng(args.seed)
    with _open(args.genome, "r") as handle:
        consensus: list[tuple[str, str]] = [(name.split()[0], seq) for name, seq in _fasta(handle)]
    if not consensus:
        p.error("the genome FASTA is empty")
    repeats: tuple[Repeat, ...] = ()
    try:
        if args.repeat:
            consensus, repeat = repeat_copy(consensus, int(args.repeat[0]), args.repeat[1], rng, ti_tv=args.ti_tv)
            repeats = (repeat,)
        haps = haplotypes(
            consensus,
            args.haplotypes,
            args.ani,
            rng,
            tree=args.tree,
            ti_tv=args.ti_tv,
            indel_fraction=args.indel_fraction,
        )
        minor = minor_alleles(
            consensus, args.minor_density, rng, frequency=tuple(args.minor_frequency), ti_tv=args.ti_tv
        )
        population = Population(tuple(consensus), haps, abundances(list(haps), args.minor_fraction), minor, repeats)
    except ValueError as e:
        p.error(str(e))
    write(population, args.output)
    sites = population.sites()
    print(
        f"{len(haps)} haplotypes, {len(sites)} truth sites ({len(minor)} minor), "
        f"{len(population.samples)} samples -> {args.output}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
