"""Variation vs error at the site level (plan §6.6, phase 8): the conservative mask.

`reference` mode counts every mismatch as an error, so a co-existing strain's alleles inflate head E (phase 7:
1.11x at ANI 99%, 1.55x at 95%). The conservative mask keeps only sites whose reads can be believed clonal,
using **Q-free covariates that are not the outcome** - coverage, coverage consistent with a single copy, and
allele linkage - so head E stays unbiased on what is left. It is the first of §6.6's two layers; the joint
latent-site model is the second and is still to come.

- **Coverage.** A site below `min_depth` has no power to bound anything, so it goes.
- **Single copy.** Two collapsed repeat copies pile up at roughly twice the contig's median depth and disagree
  at their differences, so sites outside a factor `copy_ratio` of that median go, in both directions.
- **Linkage** (§6.6 signal 3). Two non-consensus alleles that really sit on one haplotype are read together by
  every read spanning both; two errors coincide at a rate set by the error rate squared. So the same *pair of
  alleles* recurring on `min_pairs` reads is a variant pair, and both its sites go. This is the one criterion
  that looks at mismatches, and it looks at their co-occurrence rather than their frequency: an error hotspot
  (which frequency-based masking drops, §6.6) has no partner allele to recur with.

Every mask reports what it dropped and what it kept, by observed op and by reference base, so a mask that eats
one op class or one context shows up before head E is fitted from what survives.

ponytail: pair counts are held in one in-memory `Counter` over (contig, site, allele, site, allele), which is
fine for a genome-scale pileup (reads x edits^2 entries) and wrong for a large metagenome; count per contig
window if that arrives. Linkage is within a single record, not across a pair's two mates, so an insert longer
than the read links nothing across its gap; mate linkage is the upgrade path when paired `reference` runs need it.
"""

from collections import Counter, defaultdict
from collections.abc import Callable, Iterable
from typing import Any

import numpy as np

from sequencing_error_model.fit.quality import Array
from sequencing_error_model.generate import Read, align
from sequencing_error_model.sources.bam import _BASE, _COMPLEMENT, Alignment

_OP_NAMES = ("alt_A", "alt_C", "alt_G", "alt_T", "deletion", "insertion")


def _alleles(contig_start: int, n: int, reverse: bool, template: str, read: Read) -> list[tuple[int, str]]:
    """A record's substituted (reference site, read base) in reference orientation. Indels are left out: their
    site is ambiguous by a homopolymer's length, which linkage cannot afford to get wrong."""
    out = []
    for t, op in align(template, read)[0]:
        if op == "=" or op[0] == "-" or op[-1] == "-":
            continue
        base = op[-1]
        if reverse:
            base = _COMPLEMENT.get(base, base)
        out.append((contig_start + n - 1 - t if reverse else contig_start + t, base))
    return out


def linked(alignments: Iterable[Alignment], min_pairs: int = 3) -> dict[str, set[int]]:
    """Per contig, the sites carrying a non-consensus allele that recurs with another one on `min_pairs` reads."""
    pairs: Counter[tuple[str, int, str, int, str]] = Counter()
    for contig, start, end, reverse, template, read, _ in alignments:
        alleles = sorted(_alleles(start, end - start, reverse, template, read))
        for i, (s1, a1) in enumerate(alleles):
            for s2, a2 in alleles[i + 1 :]:
                pairs[(contig, s1, a1, s2, a2)] += 1
    out: dict[str, set[int]] = defaultdict(set)
    for (contig, s1, _, s2, _), n in pairs.items():
        if n >= min_pairs:
            out[contig].update((s1, s2))
    return dict(out)


def _composition(alt: Array, ref: Array, which: Array) -> dict[str, Any]:
    """The op and reference-base composition of a set of sites, as shares."""
    ops = alt[which].sum(axis=0)
    bases = np.bincount(ref[which][ref[which] >= 0], minlength=4)
    return {
        "sites": int(which.sum()),
        "ops": {name: round(float(v) / max(1, int(ops.sum())), 4) for name, v in zip(_OP_NAMES, ops, strict=True)},
        "context": {b: round(float(v) / max(1, int(bases.sum())), 4) for b, v in zip("ACGT", bases, strict=True)},
    }


def conservative(
    counts: dict[str, Array],
    linked_sites: dict[str, set[int]],
    sequence: Callable[[str], str],
    *,
    min_depth: float = 5.0,
    copy_ratio: float = 2.0,
) -> tuple[dict[str, Array], list[dict[str, Any]]]:
    """Boolean masks per contig (True: drop the site, as `sources.bam.apply_masks` reads them), and one report
    row per contig: how many sites each criterion dropped, and the composition of the dropped and kept sites.

    `counts` is a `sources.bam.count_alleles` pileup, `linked_sites` a `linked` result, `sequence(contig)` the
    reference. ponytail: the median depth is the whole contig's, so a contig that is mostly one repeat family
    calls the repeat single-copy; per-window medians are the upgrade path.
    """
    masks: dict[str, Array] = {}
    report = []
    for contig, table in counts.items():
        ref = np.array([_BASE.get(b, -1) for b in sequence(contig).upper()] + [-1])
        depth = table[: len(ref), :5].sum(axis=1)
        alt = table[: len(ref)].copy()
        has_ref = ref >= 0
        alt[np.flatnonzero(has_ref), ref[has_ref]] = 0
        median = float(np.median(depth[:-1])) if len(ref) > 1 else 0.0
        low = depth < min_depth
        copies = np.zeros(len(ref), bool)
        if median > 0:
            copies = (depth > copy_ratio * median) | (depth * copy_ratio < median)
        link = np.zeros(len(ref), bool)
        sites = [s for s in linked_sites.get(contig, ()) if 0 <= s < len(ref)]
        link[sites] = True
        mask = low | copies | link
        masks[contig] = mask
        report.append(
            {
                "contig": contig,
                "length": len(ref) - 1,
                "median_depth": median,
                "dropped": int(mask.sum()),
                "dropped_low_depth": int(low.sum()),
                "dropped_multi_copy": int((copies & ~low).sum()),
                "dropped_linked": int((link & ~low & ~copies).sum()),
                "dropped_composition": _composition(alt, ref, mask),
                "kept_composition": _composition(alt, ref, ~mask),
            }
        )
    return masks, report
