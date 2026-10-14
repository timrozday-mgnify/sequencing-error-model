"""Variation vs error at the site level (plan §6.6, phase 8): the conservative mask and the joint
latent-site model.

`reference` mode counts every mismatch as an error, so a co-existing strain's alleles inflate head E (phase 7:
1.11x at ANI 99%, 1.55x at 95%). The conservative mask keeps only sites whose reads can be believed clonal,
using **Q-free covariates that are not the outcome** - coverage, coverage consistent with a single copy, and
allele linkage - so head E stays unbiased on what is left. It is the first of §6.6's two layers.

- **Coverage.** A site below `min_depth` has no power to bound anything, so it goes.
- **Single copy.** Two collapsed repeat copies pile up at roughly twice the contig's median depth and disagree
  at their differences, so sites outside a factor `copy_ratio` of that median go, in both directions.
- **Linkage** (§6.6 signal 3). Two non-consensus alleles that really sit on one haplotype are read together by
  every read spanning both; two errors coincide at a rate set by the error rate squared. So the same *pair of
  alleles* recurring on `min_pairs` reads is a variant pair, and both its sites go. This is the one criterion
  that looks at mismatches, and it looks at their co-occurrence rather than their frequency: an error hotspot
  (which frequency-based masking drops, §6.6) has no partner allele to recur with.

`joint` is the second layer. It keeps the two coverage criteria, drops linkage as a hard criterion, and asks
of every site whether its own alleles are more counts than head E predicts *for the reads that cover it* -
head E fitted, in the same loop, on the rows the site posteriors leave clonal. Linkage becomes a term in that
posterior rather than a verdict.

Every mask reports what it dropped and what it kept, by observed op and by reference base, so a mask that eats
one op class or one context shows up before head E is fitted from what survives.

ponytail: pair counts are held in one in-memory `Counter` over (contig, site, allele, site, allele), which is
fine for a genome-scale pileup (reads x edits^2 entries) and wrong for a large metagenome; count per contig
window if that arrives. Linkage is within a single record, not across a pair's two mates, so an insert longer
than the read links nothing across its gap; mate linkage is the upgrade path when paired `reference` runs need it.
"""

from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Sequence
from dataclasses import replace
from typing import Any

import numpy as np
from scipy import special

from sequencing_error_model.fit import error, kmer
from sequencing_error_model.fit.quality import Array
from sequencing_error_model.generate import Read, align, observations
from sequencing_error_model.sources import skiver_analyze
from sequencing_error_model.sources.bam import _BASE, _COMPLEMENT, Alignment, apply_masks
from sequencing_error_model.sources.skiver_analyze import SkiverAnalyze
from sequencing_error_model.spec import Component

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


def _coverage(
    table: Array, contig_sequence: str, min_depth: float, copy_ratio: float
) -> tuple[Array, Array, Array, float, Array, Array]:
    """The two Q-free coverage criteria: reference base indices, non-reference allele counts, depth, the
    contig's median depth, and the low-depth and off-median (collapsed copy or coverage hole) site masks."""
    ref = np.array([_BASE.get(b, -1) for b in contig_sequence.upper()] + [-1])
    depth = table[: len(ref), :5].sum(axis=1)
    alt = table[: len(ref)].copy()
    has_ref = ref >= 0
    alt[np.flatnonzero(has_ref), ref[has_ref]] = 0
    median = float(np.median(depth[:-1])) if len(ref) > 1 else 0.0
    low = depth < min_depth
    copies = np.zeros(len(ref), bool)
    if median > 0:
        copies = (depth > copy_ratio * median) | (depth * copy_ratio < median)
    return ref, alt, depth, median, low, copies


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
        ref, alt, _, median, low, copies = _coverage(table, sequence(contig), min_depth, copy_ratio)
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


def joint(
    alignments: Sequence[Alignment],
    counts: dict[str, Array],
    sequence: Callable[[str], str],
    error_tokens: Sequence[str],
    *,
    flank: tuple[int, int],
    m: int,
    min_depth: float = 5.0,
    copy_ratio: float = 2.0,
    min_pairs: int = 3,
    prior: float = 1e-3,
    linkage_log_odds: float = 5.0,
    iterations: int = 3,
    threshold: float = 0.5,
    tol: float = 1e-3,
) -> tuple[dict[str, Array], list[dict[str, Any]]]:
    """The joint latent-site model (§6.6, layer two): site posteriors scored against head E's own expected
    counts, alternating with head E fitted on the rows those posteriors leave clonal. Returns the same
    `(masks, report)` as `conservative`, so it drops into the same `sources.bam.apply_masks` path.

    The two coverage criteria are kept as hard drops (a site below `min_depth` has no power, one off the
    contig's median depth is a collapsed copy), and linkage stops being a hard drop and becomes a term in the
    posterior. Per site and per error category, the count test is inStrain's with head E's expected count in
    place of an assumed Q: `log BF` between the category's observed count under head E's own predicted rate at
    *those* reads' contexts, qualities and positions, and under the MLE rate that the observation implies.
    Q enters only through fitted head E; no rate is ever read off a reported Q.

    `error_tokens` are head E's tokens, `flank` and `m` its context and Q windows (`sources.bam.window`); the
    quality alphabet is read off the rows, as `sources.bam.main` does. The loop starts from the conservative
    mask (linked sites at posterior 1) and stops when no posterior moves by `tol`.

    ponytail: hard weights on the way out (a site is dropped or kept), soft weights only inside the loop; the
    multi-sample term of §6.6 is left out until several samples of one community are on hand; and the
    posterior is the best single category's, not a sum over alleles, which costs a little power at a site with
    two real alternative alleles.
    """
    alignments = [a for a in alignments if a[0] in counts]
    ref_of: dict[str, Array] = {}
    alt_of: dict[str, Array] = {}
    parts: dict[str, tuple[float, Array, Array]] = {}
    hard: dict[str, Array] = {}
    offset: dict[str, int] = {}
    total = 0
    for contig, pileup in counts.items():
        ref, alt, _, median, low, copies = _coverage(pileup, sequence(contig), min_depth, copy_ratio)
        ref_of[contig], alt_of[contig], parts[contig] = ref, alt, (median, low, copies)
        hard[contig] = low | copies
        offset[contig], total = total, total + len(ref)

    table = observations(apply_masks(alignments, hard), flank, m, "bam", read_ids=True, template_index=True)
    if not table.counts:
        raise ValueError("the coverage criteria left no row to fit head E on")
    keys = list(table.counts)
    n = np.array(list(table.counts.values()), float)
    ri, ti, oi, ci = (table.fields.index(f) for f in ("read", "t", "op", "context"))
    rec = np.array([k[ri] for k in keys], np.int64)
    tpl = np.array([k[ti] for k in keys], np.int64)
    starts = np.array([a[1] for a in alignments], np.int64)
    spans = np.array([a[2] - a[1] for a in alignments], np.int64)
    backwards = np.array([a[3] for a in alignments], bool)
    bases = np.array([offset[a[0]] for a in alignments], np.int64)
    site = bases[rec] + np.where(backwards[rec], starts[rec] + spans[rec] - 1 - tpl, starts[rec] + tpl)
    cat = np.array([error._category(k[oi], str(k[ci])[flank[0]]) for k in keys], np.int64)
    qs = {v for k in keys for f, v in zip(table.fields, k, strict=True) if f.startswith("q") and v is not None}
    alphabet = sorted(int(v) for v in qs)  # type: ignore[call-overload]

    link = np.zeros(total, bool)
    for contig, found in linked(alignments, min_pairs).items():
        link[[offset[contig] + s for s in found if 0 <= s < len(ref_of[contig])]] = True
    depth = np.bincount(site, n, minlength=total)
    k_cat = len(error.CATEGORIES)
    gamma = link.astype(float)
    head: tuple[Component, ...] | None = None
    step = 0
    while step < iterations:
        step += 1
        weighted = replace(table, counts=Counter(dict(zip(keys, n * (1 - gamma[site]), strict=True))))
        head = error.fit(weighted, error_tokens, alphabet, init=head)
        predicted = error.probabilities(head, alphabet, table)
        expected = np.stack([np.bincount(site, n * predicted[:, c], minlength=total) for c in range(1, k_cat)], 1)
        observed = np.stack([np.bincount(site, n * (cat == c), minlength=total) for c in range(1, k_cat)], 1)
        d = np.maximum(depth, 1.0)[:, None]
        p0 = np.clip(expected / d, 1e-9, 0.5)
        p1 = np.clip(observed / d, p0, 1 - 1e-9)  # the MLE rate under "this site carries the allele"
        log_bf = observed * np.log(p1 / p0) + (d - observed) * (np.log1p(-p1) - np.log1p(-p0))
        logit = np.log(prior / (1 - prior)) + log_bf.max(axis=1) + linkage_log_odds * link
        moved = special.expit(np.where(depth > 0, logit, -np.inf))
        gap = float(np.abs(moved - gamma).max())
        gamma = moved
        if gap < tol:
            break

    masks: dict[str, Array] = {}
    report = []
    for contig, ref in ref_of.items():
        median, low, copies = parts[contig]
        g = gamma[offset[contig] : offset[contig] + len(ref)]
        variant = (g >= threshold) & ~hard[contig]
        mask = hard[contig] | variant
        masks[contig] = mask
        report.append(
            {
                "contig": contig,
                "length": len(ref) - 1,
                "median_depth": median,
                "iterations": step,
                "dropped": int(mask.sum()),
                "dropped_low_depth": int(low.sum()),
                "dropped_multi_copy": int((copies & ~low).sum()),
                "dropped_variant": int(variant.sum()),
                "dropped_variant_linked": int((variant & link[offset[contig] : offset[contig] + len(ref)]).sum()),
                "dropped_composition": _composition(alt_of[contig], ref, mask),
                "kept_composition": _composition(alt_of[contig], ref, ~mask),
            }
        )
    return masks, report


def keys(
    a: SkiverAnalyze,
    context_tokens: Sequence[str],
    *,
    flank: tuple[int, int] | None = None,
    copy_ratio: float = 2.0,
    far_fraction: float = 0.10,
    prior: float = 1e-3,
    iterations: int = 3,
    threshold: float = 0.5,
    tol: float = 1e-3,
) -> tuple[SkiverAnalyze, dict[str, Any]]:
    """`joint`'s count test for `kmer` default mode (§6.6): each `kvmer.csv` key's op counts against head E's
    expected counts at that key's coverage, alternating with head E refitted on the keys left clonal.

    The same `log BF` as `joint`, per key and op: the observed count under `fit.kmer.predicted`'s share for the
    key's own contexts, against the MLE share the count implies. A key is the locus here, so coverage is its
    observations, and the one coverage criterion kept is single copy (a key off the median coverage by
    `copy_ratio` is a collapsed copy); skiver's own `-c` sets the depth floor. Linkage needs read ids, which only
    enhanced mode has. The loop starts from skiver's outlier filter (its rejected keys at posterior 1), never
    from Q, and `context_tokens` are `fit.kmer.fit`'s.

    The conservative layer has a second, reference-free member: `far_fraction` drops a key whose observations
    sit that often at a value **two or more edits** from the consensus
    (`total_count - consensus_count - neighbor_count`, over `total_count`). Two independent errors cost the
    per-base rate squared, so a frequent far value is another locus, not noise - measured at 4.5 % of
    observations on single-locus keys against 72 % on keys at five or more loci (plan phase 6b item 1). It is a
    function of the *other* values, not of the single-edit outcome head E is fitted on, so it selects on
    multiplicity rather than on the rate, which is what lets it sit beside `copy_ratio` instead of inside the
    count test. Set it to 1.0 to switch it off.

    Returns `a` with each key's `passes_filter` replaced by the test's verdict, so `sources.kmer.fit` takes it
    unchanged, and a report. `summary_phred.csv`, which sets default mode's level, was counted by skiver over
    *its* filter's keys, so the level is carried across by `level`: the ratio of the single-edit share over the
    test's kept keys to that over skiver's (`sources.kmer.fit(level=...)`).

    ponytail: that level ratio is one number for every Q, while a variant's excess sits at the Q of a correct
    base; default mode has no per-key Q to do better (§6.2), enhanced mode's per-read output would.
    """
    full = next(t for t in skiver_analyze.tables(a, use_all=True) if t.source.endswith(":kvmer"))
    loci = list(dict.fromkeys(str(locus) for locus, _ in full.counts))
    index = {locus: i for i, locus in enumerate(loci)}
    row_of = np.array([index[r.key + r.consensus_value] for r in a.kvmer])
    skiver_drop = np.zeros(len(loci), bool)
    skiver_drop[row_of] = [not r.passes_filter for r in a.kvmer]

    total = np.zeros(len(loci))
    for (locus, _), n in full.counts.items():
        total[index[str(locus)]] += n
    median = float(np.median(total))
    copies = (total > copy_ratio * median) | (total * copy_ratio < median)

    far_obs = np.zeros(len(loci))
    far = np.zeros(len(loci))
    for r, i in zip(a.kvmer, row_of, strict=True):
        far_obs[i] = r.total_count - r.consensus_count - r.neighbor_count
        far[i] = far_obs[i] / max(r.total_count, 1)
    multi_locus = far > far_fraction
    conservative = copies | multi_locus

    gamma = skiver_drop.astype(float)
    step = 0
    while step < iterations:
        step += 1
        keep = 1.0 - np.maximum(gamma, conservative)
        weighted = Counter(
            {k: n * keep[index[str(k[0])]] for k, n in full.counts.items() if keep[index[str(k[0])]] > 0}
        )
        head = kmer.fit(replace(full, counts=weighted), context_tokens, flank=flank)
        _, _, _, observed, share = kmer.predicted(full, head, flank)
        d = np.maximum(total, 1.0)[:, None]
        p0 = np.clip(share, 1e-9, 0.5)
        p1 = np.clip(observed / d, p0, 1 - 1e-9)
        log_bf = observed * np.log(p1 / p0) + (d - observed) * (np.log1p(-p1) - np.log1p(-p0))
        moved = special.expit(np.log(prior / (1 - prior)) + log_bf.max(axis=1))
        gap = float(np.abs(moved - gamma).max())
        gamma = moved
        if gap < tol:
            break

    variant = (gamma >= threshold) & ~conservative
    drop = variant | conservative

    def edits(which: Array) -> float:
        """The single-edit share of the observations at the keys `which`."""
        errors = sum(n for (locus, op), n in full.counts.items() if op != "=" and which[index[str(locus)]])
        return float(errors / max(total[which].sum(), 1.0))

    report = {
        "keys": len(loci),
        "median_coverage": median,
        "iterations": step,
        "dropped": int(drop.sum()),
        "dropped_multi_copy": int(copies.sum()),
        "dropped_multi_locus": int(multi_locus.sum()),
        "dropped_variant": int(variant.sum()),
        "far_fraction": far_fraction,
        "far_fraction_median": float(np.median(far)),
        # Of every observation sitting two or more edits off a consensus, the share in keys the flag dropped.
        "far_observations_dropped": float(far_obs[multi_locus].sum() / max(far_obs.sum(), 1e-12)),
        "skiver_dropped": int(skiver_drop.sum()),
        "both_dropped": int((drop & skiver_drop).sum()),
        "level": edits(~drop) / max(edits(~skiver_drop), 1e-12),
    }
    tested = [replace(r, passes_filter=not drop[i]) for r, i in zip(a.kvmer, row_of, strict=True)]
    return replace(a, kvmer=tested), report


def censor(
    a: SkiverAnalyze,
    context_tokens: Sequence[str],
    *,
    flank: tuple[int, int] | None = None,
    use_all: bool = True,
    far_fraction: float = 0.10,
    prior: float = 1e-3,
    iterations: int = 3,
    threshold: float = 0.5,
    tol: float = 1e-3,
) -> tuple[dict[str, list[int]], dict[str, Any]]:
    """Variant *positions* inside `kvmer.csv` keys, to right-censor rather than drop the key (plan phase 8).

    `keys` and skiver's own outlier filter both work at key granularity - a key is kept whole or dropped whole -
    so both inherit the same bias: a key is k + v bases of observation, and dropping it for one polymorphic
    site removes the other positions' exposure too. Because polymorphic keys sit in variable sequence, the
    survivors are biased toward **conserved** sequence, whose context and GC composition differ from the
    genome's. That is a sampling bias on the context distribution, which is the thing this project exists to
    estimate, and unlike the level no external source can repair it.

    `consensus_count_up_to_v{t}` is a running survival count, so `N_{K,t-1} - N_{K,t}` is the number of
    observations whose *first* mismatch is at position `t`, over a risk set of `N_{K,t-1}`. This tests that
    count per position against head E's own predicted hazard for that key's context there, with the same
    log Bayes factor and posterior as `joint` and `keys`, alternating with head E refitted on what is left
    uncensored. `fit.kmer.fit(survival=, censored=)` consumes the result.

    What is **not** available in the default output is the per-position alt identity: `kvmer.csv`'s edit-type
    columns are aggregated over the whole value, so a sharp single-position deficit cannot be told from one
    dominant alt base (a variant) against mass spread over the ~11 neighbour types (error). A sharp deficit
    is therefore suggestive, not conclusive, and it is why `fit.kmer` splits a censored position's removal
    across the key's ops in proportion. That concentration test is the one place a fork addition would pay
    for itself on substitutions (phase 12).

    Starts from **unfiltered** keys (`use_all`), since skiver's own `p0` = 1e-12 has already discarded the
    keys most likely to carry the variation this is looking for (plan phase 6b item 2), with the
    reference-free multiplicity flag (`far_fraction`, as in `keys`) excluding the collapsed copies instead.
    """
    full = next(t for t in skiver_analyze.tables(a, use_all=use_all) if t.source.endswith(":kvmer"))
    survival = next(t for t in skiver_analyze.tables(a, use_all=use_all) if t.source.endswith(":kvmer_survival"))
    window = kmer._flank(context_tokens) if flank is None else flank
    loci = kmer._prepare(full, window, survival)
    order = list(dict.fromkeys(str(locus) for locus, _ in full.counts))
    v = loci.ctx.shape[1]

    risk = loci.exposure  # N_{K,t-1}
    failed = np.maximum(
        -np.diff(np.concatenate([risk, risk[:, -1:] - _last_failures(a, order, use_all)], axis=1), axis=1), 0.0
    )
    index = {locus: i for i, locus in enumerate(order)}
    multi = np.zeros(len(order), bool)
    for r in a.kvmer:
        row = index.get(r.key + r.consensus_value)
        if row is not None:
            far = (r.total_count - r.consensus_count - r.neighbor_count) / max(r.total_count, 1)
            multi[row] = far > far_fraction

    gamma = np.zeros((len(order), v))
    step = 0
    while step < iterations:
        step += 1
        censored = {
            order[i]: [int(t) + 1 for t in np.flatnonzero(gamma[i] >= threshold)]
            for i in np.flatnonzero((gamma >= threshold).any(axis=1))
        }
        head = kmer.fit(full, context_tokens, flank=window, survival=survival, censored=censored)
        p = error.probabilities(head, (kmer.DUMMY_Q,), kmer._contexts(full, loci, window))
        hazard = np.clip(1.0 - p[loci.ctx, 0], 1e-9, 0.5)
        d = np.maximum(risk, 1.0)
        p1 = np.clip(failed / d, hazard, 1 - 1e-9)
        log_bf = failed * np.log(p1 / hazard) + (d - failed) * (np.log1p(-p1) - np.log1p(-hazard))
        moved = special.expit(np.log(prior / (1 - prior)) + log_bf)
        moved[multi] = 0.0  # a collapsed copy is the conservative layer's business, not a position's
        gap = float(np.abs(moved - gamma).max())
        gamma = moved
        if gap < tol:
            break

    called = gamma >= threshold
    out = {order[i]: [int(t) + 1 for t in np.flatnonzero(called[i])] for i in np.flatnonzero(called.any(axis=1))}
    report = {
        "keys": len(order),
        "use_all": use_all,
        "iterations": step,
        "positions": int(called.sum()),
        "position_exposure": float(risk.sum()),
        "censored_exposure_fraction": float(risk[called].sum() / max(risk.sum(), 1e-12)),
        "keys_touched": int(called.any(axis=1).sum()),
        "dropped_multi_locus": int(multi.sum()),
        # What whole-key discard would have cost instead: every position of every key with a called one.
        "whole_key_exposure_fraction": float(risk[called.any(axis=1)].sum() / max(risk.sum(), 1e-12)),
        "composition": _key_composition(loci, order, called, window),
    }
    return out, report


def _last_failures(a: SkiverAnalyze, order: Sequence[str], use_all: bool) -> Array:
    """Observations failing at the final value position, which no later running count reveals."""
    index = {locus: i for i, locus in enumerate(order)}
    out = np.zeros((len(order), 1))
    for r in a.kvmer:
        if (use_all or r.passes_filter) and (i := index.get(r.key + r.consensus_value)) is not None:
            out[i, 0] = max(r.consensus_count_up_to_v[-1] - r.consensus_count, 0)
    return out


def _key_composition(loci: "kmer._Loci", order: Sequence[str], called: Array, flank: tuple[int, int]) -> dict[str, Any]:
    """Context and GC composition of censored positions against the keys whole-key discard would remove.

    The point of censoring is that these two differ: if the censored positions' contexts look like the genome's
    while the discarded keys' do not, whole-key discard was biasing head E's context distribution and censoring
    is not (plan phase 8, which asks for exactly this comparison, as `conservative` already reports for sites).
    """
    gc = np.array([(s.count("G") + s.count("C")) / len(s) for s in order])
    keys_touched = called.any(axis=1)

    def shares(sel: Array) -> dict[str, float]:
        if not sel.any():
            return {}
        counts = np.bincount(loci.ctx[sel].ravel(), minlength=len(loci.contexts)).astype(float)
        total = counts.sum()
        return {loci.contexts[i]: float(counts[i] / total) for i in np.argsort(-counts)[:8]}

    everything = np.ones_like(called, bool)
    return {
        "flank": list(flank),
        "top_contexts_all": shares(everything),
        "top_contexts_censored": shares(called),
        "top_contexts_whole_keys": shares(np.repeat(keys_touched[:, None], called.shape[1], axis=1)),
        "gc_mean_all": float(gc.mean()),
        "gc_mean_touched_keys": float(gc[keys_touched].mean()) if keys_touched.any() else float("nan"),
    }


def q_diagnostic(alignments: Iterable[Alignment], masks: dict[str, Array] | None) -> dict[str, dict[int, Any]]:
    """Mismatch rate by reported Q at the sites a mask dropped and at the sites it kept (§6.6, phase 8).

    A report, never an input: no site is called and no rate is fitted from it. At kept sites the curve is the
    sequencer's calibration as the reads show it; at dropped sites that are real variants a mismatch is a
    correctly read allele, so it carries a correct base's Q and the curve flattens, or rises, at high Q. A dropped
    curve that looks like the kept one says the mask took error hotspots, not variants. Substitutions only, since
    an indel has no one base to read a Q off. `masks` is `True` = dropped, as `sources.bam.apply_masks` reads it.
    """
    tally: dict[str, Counter[tuple[int, bool]]] = {"kept": Counter(), "dropped": Counter()}
    for contig, start, end, reverse, template, read, _ in alignments:
        mask = (masks or {}).get(contig)
        rows, tq = align(template, read)
        for t, op in rows:
            if op[0] == "-" or op[-1] == "-":
                continue
            q, site = tq[t], end - 1 - t if reverse else start + t
            assert q is not None  # only deletions lack a read Q
            side = "dropped" if mask is not None and site < len(mask) and mask[site] else "kept"
            tally[side][(q, op != "=")] += 1
    out: dict[str, dict[int, Any]] = {}
    for side, c in tally.items():
        qs = sorted({q for q, _ in c})
        out[side] = {
            q: {"bases": c[(q, False)] + c[(q, True)], "mismatch_rate": c[(q, True)] / (c[(q, False)] + c[(q, True)])}
            for q in qs
        }
    return out
