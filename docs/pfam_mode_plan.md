# Plan: a `pfam` evidence mode (minimal reference)

Status: **explored 2026-11-24; P0 and P1 done (results in §7, P1); P2 is ready to run on HPC (2026-11-25), and not run.** P0 is [pfam_proto.py](pfam_proto.py);
P1 is [sources/pfam.py](../src/sequencing_error_model/sources/pfam.py), started before the gate closed as a user exception.
It sits behind the §6.5.1 decision gate of [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) like the rest of
phase 11, and is not an input to it. §3's numbers are P0's; where P1 corrects one, §7 says so.

**Idea.** Stack a metagenome's reads on Pfam profiles with [kfpfam](https://github.com/timrozday-mgnify/kfpfam).
Biological variation differs from column to column; sequencing error is the same everywhere. Use that to
separate them, with the profiles as the only reference.

## 1. Verdict

**Feasible for the substitution rate and most of the spectrum, on one premise: anchor on something biology cannot
explain.** The shared-component idea does not identify the error on its own (§2). Two anchors do: **in-frame stops**
(any read, ten ops) and **mate overlap** (paired reads that overlap, all twelve ops; §3.8). They are complementary.

What it would give that the other modes cannot: an error rate from reads at **1x coverage, with no assembly, no
matching genome and no mate overlap**. ERR10889147, the development library, is such a case: only 3 % of its
pairs place in `pe-overlap`, and those are contaminated by repeats (plan phase 4).

What it will not give: `A>G` / `T>C` from stops at all, a full `Context(L,R)`, or a fine-grained model from a small
library. It should be scoped as a **rate, op-composition and Q / position / mate source**, and as a check on the
other modes. It is not a fourth full mode.

**The decisive test is not done.** Every number here is on one library without ground truth, or with errors
injected into it. P2 (mock communities, on HPC) decides go or no-go.

## 2. Why the shared component needs an anchor

For a codon of a read at column k of family f, with true codon c and observed codon y:

    P(y | f, k, z) = Σ_c π_fk(c) · E_z(y | c)

π is the biology, free per column. E is the error, shared, and a function of covariates z (Q window, cycle, mate).
Knowing only that E is shared is not enough. Any error rate can be traded for diversity in every π, so a plain
matrix factorisation (NMF, topic model) has no unique answer. It becomes unique if some columns are known to
have π = 0 at some codons. This is the "anchor word" / separability condition of topic models.

| Anchor | Where π is ~0 | Identifies | Measured biological floor |
|---|---|---|---|
| **Mate overlap on one molecule** (added 2026-11-24, user's question) | not needed: both mates carry the same biology | all 12 ops, every context and Q, at overlapped cycles. Not PCR errors | none from biology; PCR errors cancel, and misplaced pairs are the risk (§3.8) |
| **In-frame stop inside a domain** | every column | the rate; 10 of 12 ops; Q, position, mate effects | ≤ 2 × 10⁻⁴ per base (§3) |
| **Frameshift inside a domain** (indel length not a multiple of 3) | every column | indel rate, and length mod 3 | not measured; too few at 125 bp |
| **Non-consensus codon at an invariant column** | a few hundred columns | all 12 ops, in principle | 2–7 × 10⁻³ non-consensus codons (§3.6) |
| **Within-column contrast** (same column, different Q, cycle or mate) | not needed: biology does not depend on cycle | the *shape* over Q, cycle, mate. Never the level | none; it has no level |

The stop anchor is the nonsense-mutation-target method, used in virology to estimate mutation rates
([Cuevas et al. 2009](https://www.uv.es/rsanjuan/Cuevas%20JVirol%2009.pdf);
[Pauly et al. 2017](https://elifesciences.org/articles/26437)). There, the stops that cannot survive count
replication errors, and sequencing error is the nuisance. Here sequencing error is the target. A search found no
earlier use for metagenomic sequencing error. The nearest relatives are DRISEE (artifactual duplicates,
[PMC3369934](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3369934/)) and shadow regression
([PMC3495688](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3495688/)), neither of which uses coding constraint.
That search was short, and P5 repeats it before anything is claimed as new.

**Two ops are invisible to stops.** No `A>G` or `T>C` change, on either strand, turns a sense codon into a stop
(`pfam_proto.py check` enumerates them). The other ten do. Genes lie on either strand of a read, so all ten are seen in
read orientation.

## 3. What was measured

ERR10889147 (HiSeq 2 × 125, human gut, 1.07 M reads), Pfam 38.2, kfpfam Mode A, ungapped HMMER hits.

**3.1 Exposure.** A codon is *exposed* when the rest of the read still places it, without that codon: one side scores
≥ T1 = 30 bits and the other ≥ 6. A codon's own errors therefore cannot decide whether it is counted, which is the
answer to the "don't select reads by similarity" warning in phase 11. 2.8 M exposed codons, 148,000 reads,
3,689 families. The top 500 families hold 67 % of it, and 1,000 hold 86 %.

**3.2 Standard additions.** Uniform substitutions were added to 400,000 reads, then **ORF calling, kfpfam and
HMMER were all rerun**. The recovered increment over the undosed run, divided by the dose, should be 1:

| | dose 0.002 | dose 0.01 |
|---|---|---|
| T1 25, hits cut at 25 bits (first rule) | 0.51 | **0.40** |
| T1 25, hits not cut | 1.06 | 0.91 |
| **T1 30** | 1.16 ± 0.10 | **1.02 ± 0.03** |
| T1 35 | 1.20 ± 0.11 | 1.01 ± 0.04 |
| T1 40 | 1.16 ± 0.13 | 1.00 ± 0.05 |
| T1 45 | 1.24 ± 0.15 | 0.97 ± 0.06 |

- **The first rule loses half the errors, and a weaker test hid it.** Injecting errors *after* the hits were found gave
  1.00–1.05. The rule compared raw profile bits with HMMER's domain score, which is several bits lower. A
  stop-truncated ORF then passed the exposure test and still missed the hit cut.
- **Any GA-only alignment table has this problem.** See P1.
- At T1 25 the undosed run has 483 stops against 77 at T1 30: chance homology in wrong frames. **The plateau over T1
  is itself the diagnostic**, and every report prints it.
- The 0.002 column is 1.5–2 standard errors high, with a shared baseline, so its rows are not independent. P2 settles it.

**3.3 Spectrum.** An EM assigns each stop to the codon it came from. At dose 0.01 (3.3 × 10⁻³ per op) all ten ops
come back at 2.6–3.8 × 10⁻³, end to end (0.77–1.15 of the truth). Attribution needs the **sample's own stack** as prior:
with the HMM's emissions alone `C>T` / `G>A` came out 1.45x too high and `A>C`, `C>G`, `T>G` collapsed to 0. This
is the part of the original idea that pays: the column's own reads say what the codon was. In a second test with a skewed
spectrum (rates 2 × 10⁻⁴ to 3 × 10⁻³, errors added after the hits, first rule), the stack gave the eight better-exposed ops 0.83–1.24 of
the truth, but `A>C` 0.68 and `T>G` 0.36: those have the fewest events and are still biased low. The prior weight (1) was set after
seeing that table, so P1 chooses it on held-out data.

**3.4 It behaves like error.** Five families excluded (3.5), T1 35, stops per 10⁶ codons:

| lowest Q in the codon | < 15 | 15–19 | 20–29 | 30–35 | ≥ 36 | | R1 | R2 |
|---|---|---|---|---|---|---|---|---|
| | 160 (5 stops) | 550 (18) | 100 (10) | 16 (5) | 21 (33) | | 29 (30) | 40 (41) |

A 25-fold range between Q 15–19 and Q ≥ 30, with Q entering only as a covariate of counted stops, never as a label. Position is flat
within noise. Estimated base rate 4.5–7 × 10⁻⁴ (T1 40 to 30), floor included.

**3.5 The floor is biology, and partly not bacteria.** Five families hold 37 of 108 stops. The worst, `PF29843`
`ORF2p_C`, is by its name **human LINE-1** (dead retroelement copies in host reads; not checked by mapping): 19 stops in 978 codons, and `PUB_ORF2p` adds 2 in
65. The others are a transposase (`Y1_Tnp`), `Zincin_1` (8 stops in 115 codons, at one column) and `NAD-GH`. The five are 0.3 % of the codons. Stops at Q ≥ 36 stay flat at
21 per 10⁶ codons, i.e. ≤ 2 × 10⁻⁴ per base, as an upper bound on the remaining floor: it may be pseudogenes, recoded
genomes or real errors. Reported Q cannot split these, and this plan does not use it to.

**3.6 Invariant columns are weak for Illumina.** At columns that never varied in half the reads (≥ 30 codons), 0.7 % of
held-out codons differ (0.2 % at ≥ 100 codons and 99 % conserved), flat in Q, so it is biology. As a substitution rate that is 1–3 × 10⁻³, 2–6x the
error we want, and the ≥ 100 case is only 6 columns. Usable at ONT error rates, or for `A>G` / `T>C` with the mixture model and far more data.

**3.7 Indels.** About a dozen candidate frameshift pairs in 1.07 M reads. Too few at 125 bp, as expected (§6).

**3.8 Mate overlap as a second anchor.** Two facts.

*It makes the original idea work.* With π_fk free per column and E shared, one read per base does not determine E:
four starts of an EM on exact counts give 0.017–0.031 for a true 0.0168. Both mates over one base do: every start
returns 0.0168, for every op, **because columns differ**. The same two-mate data with every column alike gives
0.15–0.18. So diverse biological variation is what identifies the shared error, once each molecule is read twice
(`pfam_proto.py identify`; infinite data, no Q or context, so a statement about identifiability, not precision). The stack also tells
which mate erred, which `pe-overlap` decides from base composition and its own symmetry.

*It gives a second, independent estimate here.* Pairs whose mates sit on opposite strands of the same family and columns
propose an overlap; the bases outside the counted codons confirm it (≥ 20 of them, ≤ 2 % mismatched, so the
check never sees the counted bases). Of 1,018 pairs, 801 match outside and 161 disagree there by > 20 %: **false
overlaps from tandem-repeat families** (`SH3_5`, `Choline_bind`, `S1`, `YadA_head`, `CshA_repeat`, ...), which a panel can
exclude by name. The 751 pairs that pass at ≤ 2 % hold 20,871 bases and 31 mismatches: 7 × 10⁻⁴ per base per mate, or **2 × 10⁻⁴** without one
family's two pairs, which hold 22 of the 31. (That exclusion looks at the outcome, so it is a range, not an estimate.) The
stop anchor gives 4.5–7 × 10⁻⁴, with floor and PCR error inside and over all cycles. The two agree to within what 9–31
events allow.

*What overlap cannot do alone.* It exists only for short inserts, only at overlapped (late) cycles, never for long or single-end reads, and it is
blind to PCR errors. So the anchors serve different purposes:

| | Stops | Overlap |
|---|---|---|
| Libraries | any reads, any depth | short inserts |
| Ops, contexts | ten ops, stop-adjacent | all |
| PCR/library error | included | excluded |
| Pseudogene floor | inside the rate | **measures it**: a stop in both mates is biology or PCR, in one mate is the sequencer |
| Needs the other? | overlap splits its floor | stops check cycles and reads that do not overlap |

The last row is why both belong in one fit. A stop read by both mates is not a sequencing error, whatever its Q, so the overlap gives a
**Q-free floor**, which 3.5 could not. Here it has no power: 1.4 stops expected in the shared codons, one found, in one mate.

## 4. Inference

Tiers, each kept only if it beats the one before on held-out data and on mocks.

- **Tier 1: stop anchor (P1).** The EM of 3.3. For each stop, a fractional count per candidate source op; for
  each sense codon, the ops that *would* have made a stop. This is a partially observed categorical: a base
  exposes some ops only. It maps on to head E through one change: `fit.error._mask` already drops the impossible
  self-substitution, so let a row also carry the set of ops it can show, and condition the softmax on `=` plus that set.
  Output is the existing `CountTable` with fractional counts, `truth=True`, fields `op, q, q±1, pos_start, pos_end, mate`.
  Family random effects (a Gamma or Beta-Binomial on the stop rate) replace the Poisson cut in 3.5, so a family with
  excess stops loses weight instead of being dropped by hand.
- **Tier 1b: overlap anchor (P1).** Two-view EM over (family, column) with the stack as prior, sharing E with Tier 1, so stops
  constrain the rate and overlap the twelve ops. Reuse `pe_overlap.place` and its rival guard, and require its offset to equal
  the one the HMM columns imply: **two independent placements**, which should rescue part of what failed on ERR10889147. Both
  mates' rows then go through the existing soft-EM, with π_fk in place of base composition.
- **Tier 2: whole-codon mixture (P4, only if needed).** The model of §2 with π_fk a Dirichlet around the stack,
  fitted by EM over all codons with a point mass at "invariant". It adds `A>G` / `T>C`, and more events. It is expected to fail for Illumina by 3.6 and
  to be useful for ONT.
- **Tier 3: covariate contrast (P4; §13.2).** A logit with a fixed effect per column. It adds the shape of `A>G` / `T>C` over Q and
  cycle, with the level left `identified = false` (§5.6's convention for flagged components).
- **Not planned:** plain NMF, a VAE or a Bayesian hierarchical sampler. They re-create the identifiability problem of
  §2 and add nothing the EM lacks.

**Level.** Report the stop rate with the floor **inside** it, and a separate interval for the floor from P2. Never subtract a
floor inferred from reported Q.

## 5. Where it fits the rules and the existing modes

- Q is never a label: stops are labelled by the profile and the genetic code. 3.4 is a diagnostic.
- `kmer` default mode and its pinned skiver are untouched.
- Truth failures are unlike the other three modes': it needs coding sequence and the standard code (table 11), and it counts PCR
  errors. `pe-overlap` minus `pfam` is therefore another PCR/library estimate, where a library supports both.
- Bias to expect: coding DNA only, conserved families, and stop-adjacent contexts. On the coding strand `C>T` is only seen
  at `CAA`/`CAG`/`CGA` and `G>A` only at `TGG`. So a rate per op is an average over these contexts, and it differs from the
  genome-wide one if error is context-dependent. `Context(L,R)` is fixed at 0 and flagged.

| Component | `pfam` mode |
|---|---|
| Substitution rate, 10 ops | yes |
| `A>G`, `T>C` | where mates overlap (Tier 1b); otherwise no (Tier 2 or 3) |
| centre Q, Q window, position, mate, tile | yes, given enough stops |
| `Context(L,R)`, homopolymer, clustering | no |
| indels | long reads; short reads by locus if §13.1 pays |
| head Q | from the FASTQ, as in `kmer` |

## 6. Yield: what reads buy

At 125 bp and T1 30, one million reads give about 270 stops per 10⁻³ of error (**at P1's T1 35: about 200**, so
multiply the reads below by 1.3). The opportunity is 0.8 M per million
reads, summed over ops; rare ops get 4 %.

| Target | Stops | Reads at 10⁻³ |
|---|---|---|
| total rate to 5 % | 400 | 1.5 M |
| each of the 10 ops to 10 % | ~100 for the rarest | ~10 M |
| Q × position × mate as in `pe-overlap` | ~10,000 | ~37 M |
| `QualityWindow(m)` | ~100,000 | ~370 M: not a real option |

Illumina libraries of 10–50 M pairs fit this. Reads under 90 bp will lose most exposure, since T1 needs about 25 codons
of domain, and that was not measured. Long reads gain twice: one anchor covers several domains, and frameshifts
become usable, but **HMMER's ungapped ORF search is the wrong aligner there**. It needs a frameshift-aware
profile aligner ([BATH](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11223822/), HMMER-derived), and at ONT error rates ORFs
break too often for kfpfam's six-frame input. P4 tests this.

## 7. Phases

Each is gated on the one before. Local means a laptop; **HPC means prepared here, run by the user; agents never access HPC**.

### P0: prototype ✓ (this document, `pfam_proto.py`)

### P1: `sources/pfam.py` (local) ✓

Both anchors, gapped hits, `orfs` and `fit` commands, `pfam` as a spec mode; per-row op sets in `fit.error` (a full set reproduces the old
fit to the bit) and `compare.op_rates`. 16 tests in `tests/test_pfam.py`, 3 in `tests/test_fit_error.py`.

**What P1 changed from the plan, each found by a test or by the real run:**
- **T1 is 35, not 30.** On ERR10889147 the undosed crude rate falls 6.4, 4.1, 3.4, 2.1 × 10⁻⁴ over T1 30, 35, 40, 45 before family weights and
  3.0, 2.8, 2.1 over 35, 40, 45 after. `fit` warns when the plateau fails. P0's 4.5–7 × 10⁻⁴ (T1 30) included chance homology.
- **Exposure is weighted by P(codon clean).** Errors cluster in low-Q runs, so a codon counted as clean often is not. Without the weight the rate
  reads −4 % (0.2–0.6 % error), −9 % (1.1 %), −18 % (2.9 %). With it, six seeds each: 0.98 ± 0.03 (0.24 %), 1.01 ± 0.02 (0.6 %),
  1.02 ± 0.02 (1.1 %), 1.09 ± 0.03 (2.9 %). **Valid to about 1 % error per base.**
- **Family weights** are a spike-and-slab (negative binomial null) on each family's expected stops under head E, over its own codons. A flat rate per
  codon flags families for their codon usage. They do the panel's job, so there is no panel file (`read_hits` takes one if cost needs it).
- **Overlap needs `pe_overlap.place` to agree** with the column offset, plus ≥ 20 outside bases at ≤ 5 % mismatch, plus a family weight. Without the
  agreement, repeat families (`CshA_repeat`, `Rib`, `DivIVA`, ...) put `A>G` at 2.6 × 10⁻³.
- A codon the overlap uses is not offered to the stop anchor as well.

**ERR10889147, every 2nd pair, 535,458 reads, T1 35, errors added before ORF calling** (9 minutes and 3.2 GB per fit):

| | Undosed | +0.002 | +0.01 |
|---|---|---|---|
| stops | 52 | 313 | 1,137 |
| overlapped bases / disagreements | 9,186 / 15 | 8,775 / 43 | 7,194 / 151 |
| recovered, ten ops, full fit | | 1.10 (± 0.07 counting error) | **0.97** (± 0.03) |
| recovered, all twelve | | 1.22 | 1.05 |
| crude stop rate, T1 35 / 40 / 45 | | 1.05 / 1.18 / 1.06 | 0.95 / 0.97 / 0.95 |
| per-op range, ten ops | | 0.83–1.43 | 0.79–1.17 |
| `A>G`, `T>C` (overlap only) | 2.6, 1.2 × 10⁻⁴ | 1.42, 2.31 | 1.28, 1.58 |

- Before the overlap fix (and at T1 30) the same fits recovered 0.90 at both doses.
- The `A>G` / `T>C` excess is sampling: at 0.01 the errors that fell on the overlapped bases were 17 each, **1.49x** expected, and the fit
  follows them. Over all ops, injected over expected was 1.03 at both doses, so the overlap does not select against errors. A library with 3 % placed pairs cannot pin these two ops.
- Undosed: 3.7 × 10⁻⁴ per base (R1 2.8, R2 4.6), falling from 1.9 × 10⁻² at Q < 20 to 1.3 × 10⁻⁴ at Q ≥ 36. 52 stops, so ± 14 % before the floor. Downweighted:
  `PF29843` (LINE-1 ORF2p), `PF06262`, `PF26139`, `PF23714`; overlaps `PF08481` (`GBS_Bsp-like`), `PF05103` (`DivIVA`).
- The floor has no power here: 0 stops in both mates of 3,062 shared codons (95 % bound 9.8 × 10⁻⁴ per codon).

**Exit criteria.**
- Standard additions at 0.97–1.03 end to end: **met at 0.01 (0.97, on the edge); not testable at 0.002**, whose counting error is ± 7 %. The gap
  at 0.01 is within counting error, but so is a 3 % loss.
- Ten ops within phase 2's tolerances on generated reads: **not shown.** The total is within its counting error; per op, the test allows 0.4–2.5 and a
  geometric mean of 0.75–1.25, because 20,000 reads give too few stops per op. P2's read numbers are what can test it.

**Not done:** the RawR / hmmalign parquet adapter (question 1 below; `pairs.tsv` only), panel file, the `Tile` term, indels.

### P2: mock communities and an isolate, with truth (HPC) - **go / no-go**; prepared, not run

**Libraries** (`workflows/pfam/runs.csv`): six ZymoBIOMICS runs (ERR2984773, SRR12324253, SRR12829160, ERR2935805, SRR23926920, SRR32629731; two are
about 100 bp, where yield is uncertain) and SRR24523812, the isolate with all three modes fitted (phases 5–6), each as 3.6–4.8 M pairs by stride.

**What runs** ([workflows/pfam](../workflows/pfam/README.md)): undosed, +0.002 and +0.01 copies, searched against Pfam 39.0; a fit of both anchors at each dose and of
each anchor alone undosed; the undosed reads aligned to the genomes (hit tables and alignments kept in `results/keep/` on HPC, §13). `pfam_truth score` then, per library:
- classifies every counted stop by the genome under it: `error` (a sense codon there), `genome` (a stop there: pseudogene, gene end, misplaced read),
  `variant` (a site where these reads carry another base at ≥ 20 %), `indel` or `unaligned`;
- labels **every base of every exposed codon**, not only the stop-adjacent ones a fit saw, and compares each fit's predicted errors with the observed ones per
  op, Q, mate and read position.

**Go rule** (`pfam_truth.verdict`; **agreed 2026-11-25, before any P2 run**, so not to be changed after seeing its results). A library passes if, for the fit with both anchors:
- the ten ops' total is within 15 % of the genomes', and its floor's 95 % bound is ≤ 2 × 10⁻⁴ per base (`genome` + `variant` stops the fit kept);
- no single op, Q group, mate or position quarter is `outside` 25 %, that is off by more than 25 % + 2 sd.

A check counts only where its sd, from the stops and errors counted, is at most half its tolerance. A library whose total has no power neither passes nor fails. **Go** if five of seven
pass, **no-go** once too many have failed for five to pass. This replaces the earlier "curves inside `reference`'s intervals", which a truth measured
on millions of bases would have made impossible to meet. Reported but not judged: standard additions, each anchor alone, `A>G` and `T>C` (the anchor test).

**The laptop rehearsal: what it changed before any HPC time was spent.** Reads were simulated from the ten genomes of ZymoBIOMICS' standard with ERR10889147's fitted heads
and searched with all of Pfam 38.2.
- **The floor criterion would have failed.** Error-free reads (250,000 pairs, 3x on the bacteria) give 109 stops in 3.9 M exposed codons at T1 35: 2.3 × 10⁻⁴ per base. T1 does not remove it (2.2 at 40, 2.0 at
  45), nor does T2 (2.0 × 10⁻⁴ at 10, 2.1 at 15, rising beyond), and family weights cannot, as the stops are spread over about 100 families. They concentrate in *L. fermentum*, *Salmonella* and *E. coli*:
  real pseudogenes.
- **A stop that recurs at its locus is biology** (`pfam.replicates`, `drop_replicated`): same codon between the same 12-base flanks in other fragments, more often than errors at 1 %
  could explain (p < 0.01). No Q, no genome, 2 seconds. It drops 80 of 114 at 3x, and nearly all that stay have no second read. **So it needs coverage: at 1x most of the floor remains.**
- **Rehearsal with errors** (150,000 pairs, 5x on three bacteria, errors at 1.2 × 10⁻³; kit code throughout):

  | | |
  |---|---|
  | stops | 279 found: 149 errors (none dropped), 128 in the genomes (113 dropped), 2 unaligned (both dropped) |
  | floor per base | 4.4 × 10⁻⁴ → 5.2 × 10⁻⁵ (95 % bound 8 × 10⁻⁵) |
  | ten-op total, fit over genomes | 1.00 ± 0.08 |
  | standard additions, +0.002 / +0.01 | 1.02 / 0.99 |
  | anchors alone (ten ops) | stops 1.10, overlap 1.09 |
  | `A>G`, `T>C` | both 0.81, 1.25; stops alone 0.70, 1.34; overlap alone 0.82, 1.21 |

  149 error stops are too few to test single ops, and the stop anchor's `A>G` / `T>C` are extrapolations. 4–5 M pairs should give 1,000–4,500, depending on the platform.
- **On ERR10889147** the same filter drops 4 of 52 stops; 23 are contradicted by another read and 25 have none to compare with. Of the added errors 26 % have none, so
  the undosed excess of unsupported stops (25, where 26 % of about 31 would be 8) is probably biology in rare genomes, and the crude rate an overestimate there. To be turned into an estimate in P3, after P2's
  truth shows whether the argument holds: real errors also cluster in reads whose flanks are bad.

**Not tested:** Slurm, the real images (the `singularity` profile ran with a stand-in runtime), Pfam 39.0, memory at 10 M reads (64 GB asked for a fit and 128 GB for scoring; a fit took 1.8 GB at 300,000 reads), real truth.
- **No-go** outcome is still useful: report the floor and keep the mode as a lower-effort check, or drop it.

### P3: real libraries (HPC)
- The three biome runs long enough (freshwater 145 bp, human gut 112, marine 151), the two 250 bp runs (DRR160704 soil, DRR092732 deep sea), and
  ERR10889147, each 5 M pairs by stride (`observations.subset`; file order is biased). Soil ERR5866877 (75 bp) and bird gut ERR6785304
  (72 bp, DNBSEQ) are too short (§6), and are run only to show the mode refuses them.
- Report: rate, ops, Q / position / mate curves, the floor interval, T1 plateau, standard additions at two doses, and a stop-code check
  (per-stop `TAA`/`TAG`/`TGA` excess by family, to catch recoding: code 4, 11, 15, and selenocysteine).
- Compare with `pe-overlap` where it works. Skip long reads here.
- **Exit:** the plateau and the dose-response hold in every biome, or the biomes where they do not are named.

### P4: tiers 2 and 3, and long reads (HPC, only after P2's go)
- Tier 2 against `A>G` / `T>C` truth in the mocks. Keep it only if it gets within 25 % at Illumina rates.
- Long reads: one ONT R10 and one HiFi mock the user chooses, with BATH or frameshift-joined kfpfam hits, for indel rate and length mod 3 against `reference`.

### P5: integration
- Update §1.1, §5.6, §8 of the main plan, the README, and the user guide. Repeat the novelty search. No exporter work: exporters read a spec, and
  the mode adds a source, not a component. It earns a place in joint fitting (phase 11) only if P2 is a go.

## 8. HPC preparation

Done for P2: [workflows/pfam](../workflows/pfam/README.md) has the workflow, images, inputs, layout, commands, cost and what to bring back. The `test` profile runs it on a laptop
without containers. It reuses `sem-tools.sif` (phase 6b), and the kfpfam image, Pfam 39.0 prep and index of B11. P3 uses the same workflow with other rows in `runs.csv`, but `ALIGN`, `SCORE` and
`VERDICT` need genomes, so it will need a switch to skip them.

## 9. Risks

| Risk | Where it bites | Guard |
|---|---|---|
| **Selection on error**: recruitment drops error-rich reads | P1. Already cost a factor of 2 (3.2) | the standard-addition run, end to end, in every report |
| **False homology** at low T1 | T1 25 gives 6x the stops | the T1 plateau |
| **Biological floor** (pseudogenes, retroelements, recoding, phage): measured at 2.3 × 10⁻⁴ on error-free Zymo reads | Illumina below ~3 × 10⁻⁴ | locus replication where a locus has ≥ 2 reads; family weights; P2 truth |
| **Locus replication drops real errors** that recur (a hotspot, PCR copies) | deep libraries | none dropped of 149 in the rehearsal; P2's truth keeps them and will show it |
| **Host reads** | human LINE-1 in a gut sample | host depletion first, and the panel |
| **Stop-adjacent context only** | a rate that is not genome-wide | reported as such; compared on shared contexts |
| **Reads < 90 bp, DNBSEQ** | no exposure | yield table; refuse under a threshold |
| **Wrong source-codon prior** | per-op split (3.3) | stack prior; per-op dose-response |
| **Scooped** | novelty not checked deeply | P5 |

## 10. Open questions

1. **The hmmalign parquet:** where is it, and does it have GA-only or all hits, nucleotides or residues?
2. ~~Go / no-go rule in P2~~: agreed as proposed on 2026-11-25 (§7).
3. **The gate.** The main plan says nothing in phases 9–11 starts before §6.5.1 closes. The user made an exception for P1. P2 needs its own;
   the argument is that it does not touch `kmer` mode and can share an HPC window with the gate run.

## 11. Reproduce

Needs kfpfam's `data/` and ERR10889147 (both local, untracked). Local times: 2 minutes per dose for 400,000 reads.

```bash
python docs/pfam_proto.py check
python docs/pfam_proto.py orfs R1.fastq.gz R2.fastq.gz out --n 400000 --dose 0.002
FILTER_MIN_GAPPED=21.25 filter join data/prep data/index/k4T6.5.tidx 17 16.5 400000 out.faa out.filter.tsv
python recruit/pairs.py data/pfam/Pfam-A.h3m data/prep/models.tsv out.faa out.filter.tsv 21.25 6 out.pairs.tsv
python docs/pfam_proto.py exposure R1 R2 data/prep out.pairs.tsv out.coords.tsv out.npy --n 400000 --dose 0.002 --t1 30 --hit 0
python docs/pfam_proto.py fit data/prep out.npy baseline.npy
```

## 12. Explored: a consensus of the stack, and MGnify members (2026-11-25)

Two questions: can the dominant allele(s) of the stacked reads be a second basis for the error model, and would assigning the reads to MGnify
proteins (kmer-functional-profiler's `pfam-anchor`) make the stacks clean enough for that? Code: [pfam_locus.py](../src/sequencing_error_model/pfam_locus.py) (the locus rows; the column and MGnify rows were a prototype, not kept). Nothing here
is in the `pfam` fit. Every grouping is judged the same way: each exposed codon (T1 35) against the **unanimous codon of ≥ 2 other fragments** of its group, disagreements split by effect. Pure error gives
equal rates per opportunity for synonymous, missense and nonsense; biology gives synonymous ≫ missense ≫ nonsense.

| Grouping | Mock¹: labels that are errors | Real²: syn / mis / non × 10⁻⁴ | Codons usable (real) |
|---|---|---|---|
| Pfam column, dominant ≥ 50 % | 0.5 % | 1,850 / 170 / 1.5 | 63 % |
| Pfam column, dominant ≥ 99 %, ≥ 20 deep | 16 % | 483 / 13 / 0 | < 1 % |
| Pfam column, top alleles (any seen twice) | 4 % | not split | 78 % |
| same 4 residues each side (protein) | 50 % | 165 / 7.4 / 2.1 | 52 % |
| same 10 residues each side (protein) | 84 % | 73 / 4.0 / 2.3 | 27 % |
| **locus: same 12 bases each side** | **98 %** | 13.8 / 3.0 / 1.6 | 53 % |
| locus: same 18 bases each side | 99 % | 8.3 / 2.9 / 1.7 | 46 % |

¹ The three-genome rehearsal of §7 P2 (no strain variation). ² ERR10889147, 535,458 reads; nonsense counts are 7–13, so ± 30 %.

**1. A consensus per column does not work, and more alleles or depth do not save it.** A column holds many genomes, and a rare allele one change from a common one is nearly always another genome. At 10⁻² added error, 29–51 % of
top-allele labels were errors; at the mock's 5 × 10⁻⁴, 4 %.

**2. A consensus per locus does, as a mixture and not as hard labels.** It is the stacking idea at the resolution where a molecule's neighbours are the same gene in the same strain.
- What leaks in on real reads is strain variation, and it is mostly synonymous: 12 bases each side leave about 38 % of the labels errors, but about 70 % of the missense and nonsense ones.
- **The codon effect separates the two.** Model per substitution: rate = E + b × (1, ω, 0) for synonymous, missense, nonsense. ω is the sample's own missense-to-synonymous ratio among **alleles that recur** (≥ 3 fragments), which are biology. No
  Q, no reference, and all twelve substitutions are identified, since `A>G` and `T>C` have synonymous and missense forms.

| | Stop anchor | Locus mixture | Truth |
|---|---|---|---|
| mock, per base × 10⁻⁴ (149 stops; 2,264 labels) | ± 8 % | 5.02 ± 0.13 (12 bases), 4.87 (18), 5.18 (9) | 5.03, 4.87, 5.28 |
| mock, per substitution | 0.8–1.5 of truth (with overlaps) | 0.84–1.05 of truth | |
| ERR10889147, undosed × 10⁻⁴ | 3.8 crude; 2.4 if unsupported stops are biology | 2.15 ± 0.5 (12), 2.06 (18), 1.82 (9) | unknown |
| +0.01 recovered | 0.97 (§7 P1) | 0.98 (9 and 12), 0.99 (18) | 1 |

- **Yield:** about 15 times the errors (2,170 against 149 in the mock, about 420 against 30 in ERR10889147), from the same reads.
- On ERR10889147 the total is not more precise than the stops' (± 24 % with ω from 30 recurring alleles): the biology has to be subtracted. The gain is the spectrum, `A>G` / `T>C`, and in deeper or less diverse samples the rate.
- The naive count falls from 5.6 to 3.8 × 10⁻⁴ between 12 and 18 bases; the mixture does not move. That stability, and agreement with the stops, are the checks a real library allows.
- A cheap use needs none of the mixture: the reads that **contradict a stop give its source codon**, where the stack gives a prior (119 of 140 stops in the mock, 23 of 48 in ERR10889147).

**3. MGnify members would not make the stack cleaner, and could help elsewhere.**
- *As a grouping they are too coarse.* A unit is a protein at 90 % identity, and the error model is nucleotide: species sharing a protein differ at every synonymous site. The protein rows above are the best case (identical residues) and leave synonymous rates 5–12 times the locus'. Grouping the reads' own
  bases does better with no database. Assignment is also ambiguous: a median of 8 tied units per codon, 244 at the 90th percentile (real reads, the local slice of 40 complete families).
- *As an outside consensus where no second read exists, nearly but not enough.* 80 % of exposed codons have a member that matches every other residue of the read, and of those 2,889 of 2,894 match it at the codon too. Five
  mismatches, 1.7 × 10⁻³ per codon (0.6–4 × 10⁻³), against about 0.5 × 10⁻³ expected of errors. So most are novel biology, there is no synonymous class, and no way to estimate ω, hence nothing to subtract with. 212 reads only.
- *Where they can help:* the **source residue of a stop** (a member's residue at that column, right to 99.8 % when any tie counts, 98.6 % by majority), for the stops without a second read; and those stops' gene ends and chance placements, **not** strain pseudogenes, which have intact relatives in MGnify.
  Only 40 families could be tested here. Worth a test only if P2 shows per-substitution error from the stack prior to matter.

**Limits of what was measured.** Real strain variation has no truth here: the mock has none, and the dose only shows that errors are recovered in its presence. ω comes from recurring alleles, and singleton variants are younger and less purged, so the true ω is higher and **E is biased up**; P2's D6331 (five *E. coli* strains) tests that.
That ω depends on divergence is known from genomes: within-species pN/pS is near 1 between very close strains and falls to about 0.1 with
divergence (Rocha et al. 2006, bacterial genome pairs; Garud et al. 2019, gut strains; both cited from memory, to be checked in P5). It can be
measured three ways: (1) P2's truth, already reported (`truth.omega`, the biology labels, beside the recurring-allele ω), with power in D6331;
(2) a per-family table from public genomes, pN/pS on the same Pfam domains between conspecific pairs at the locus' divergence (identical 12-base
flanks, a few differences per kb), built once from MGnify's species catalogues or GTDB; (3) a sample's own assembly where it has one (variants
≥ 20 %, as P2's `variant` class). **Order:** (1) measures the bias; build (2) only if the bias is large, to replace or correct the sample's ω.
Requiring clean flanks drops reads where errors cluster: −2 % at 12 bases and −5 % at 18 on simulated heads, about even across Q; real reads may be worse. An error hotspot that recurs counts as biology, a PCR error as an error. It needs ≥ 3 fragments at a locus. It is `kmer` mode's idea (a key
and the consensus of what follows) with two things `kmer` lacks: flanks on both sides and a reading frame, hence the effect classes. The same split could be tried on skiver's output where a frame is known.

**Recommendation.**
1. **Add a `locus` anchor (Tier 1c)**: rows for every base of a codon with ≥ 2 agreeing fragments, weighted by the mixture's P(error | substitution, effect), in the same fit as the stops and the overlap. Nonsense rows overlap the stop anchor's and are counted once.
2. **Before P2 runs on HPC**, add its scoring only: the two tables above (counts and opportunities by substitution and effect, and recurring alleles) and the share of labels that the genomes call errors. They are a few hundred numbers per library, and D6331 and D6310 are where its one unproven part, ω, has truth. **Reported, not judged**: the agreed go rule (§7) stays as it is. **Done 2026-11-25:** [pfam_locus.py](../src/sequencing_error_model/pfam_locus.py) and a `LOCUS` step in the workflow, shown in `verdict.md`. On ERR10889147 the module reproduces the table above (2.15, 2.06, 1.82 × 10⁻⁴ undosed) and recovers +0.002 at 1.00 and +0.01 at 0.97–0.98.
3. Do now, whatever else: use contradicting reads for a stop's source codon. **Prototyped 2026-11-25** ([pfam_source_proto.py](pfam_source_proto.py)): the reads of
   a stop's locus that agree on a sense codon name its source right 99.3-100 % of the time (mock, against the genome), where the stack prior's best
   is right 88-89 %. But the spectrum hardly needs it: with equal rates per op, the stack prior puts 0.9 % of stops on the wrong op at +0.01
   (2,648 error stops) and 5.0 % undosed (149), the locus where it exists 0.5 % and 2.6 %, both under the counting error of an op at P2's depth.
   The prior's errors mostly cancel across ops. Coverage: 60 % of error stops at +0.01 and 87 % undosed (mock, 5x); 31 % and 43 % in ERR10889147,
   where the locus agrees with the prior's best for 15 of 16 and 445 of 483. **So:** a cheap hard assignment where a locus exists, worth taking
   when the fit is next touched, not a reason to; its value is where the prior is weakest (58 % right on the mock's stops with no locus) and in
   samples with more strain variation per column, which the mock lacks. **In the fit since 2026-11-25** (`Stop.source`, set by `drop_replicated`; `codon_weight`
   puts `LOCUS_RIGHT` = 0.99 on it, the rest by the stack). Mock refit: 119 of 140 kept stops sourced undosed, 1,856 of 2,345 at +0.01; against the
   genomes, the mean |log ratio| over the twelve ops 0.139 → 0.138 and +0.01 recovered 0.986 → 0.985, as the prototype predicted.
4. Not the column consensus. Not MGnify grouping; revisit members for source residues after P2.

## 13. Proposed: indels, and the codons no anchor labels (2026-11-25)

Two ideas from the user's questions, not yet tested. Both run on the hit tables and reads that P2 already makes, so neither needs a new search.

**13.1 Indels by locus (short reads).** Today the mode is substitution only. Gaps in a hit are whole codons, which are biology, and a 1–2 base
sequencing indel breaks the frame. The planned anchor is the frameshift pair: two hits of one family in one read, in different frames, at adjacent
columns. That gives about a dozen per million 125 bp reads (§3.7), so it pays only on long reads with a frameshift-aware aligner (P4).
The locus gives a second route. A read whose left 12-base flank matches a locus, and whose right flank matches the others at an offset of ±1 or ±2
bases, has an indel between them. Strain indels in coding sequence are mostly multiples of 3, so the frame splits error from biology, as the codon
effect does for substitutions. Indel rate and length mod 3 come out per locus; homopolymer length can be read from the flanks.
- Test: P2's mocks, against indels in the aligned reads. Yield is the open question: Illumina indels are about 10⁻⁵ to 10⁻⁶, so it may only show a
  rate, not covariates.
- Mate overlap could also show an indel as a change in offset along the overlap. Check whether `pe_overlap.place` allows one before using it.
- **Prototyped 2026-11-25** ([pfam_indel_proto.py](pfam_indel_proto.py)). Each read's next 24 bases after a locus' 12-base key are compared with
  the consensus of the locus' other reads; at the first difference the next 10 bases must match it exactly as a substitution, or after a shift of
  1-3 (every compared consensus base confident: letting unclear ones match anything made 3 × 10⁻³ false calls per base). Mock (150,000 pairs, three
  genomes, none with indels as simulated), 5 × 10⁻⁴ indels added (lengths 1 / 2 / 3 at 60 / 20 / 20 %, half insertions):

  | Locus membership | Indels at the positions it covers | Found | Calls that are real | False calls, mock without indels |
  |---|---|---|---|---|
  | own exposed codon, search rerun | 9.6 × 10⁻⁵, 58 % of them in frame | 90 % | 100 % | 0 in 3.5 M positions |
  | 12-base key, search rerun | 4.2 × 10⁻⁴ | 64 % | 99 % | 21 in 7.8 M (2.7 × 10⁻⁶) |
  | **12-base key, loci from the undosed reads, no second search** | **5.15 × 10⁻⁴** | **64 %** (64, 65, 63 % by length) | **98 %** | as above |

  - **Exposure selects on the error counted**: a frameshift ends the hit, so the codons before it stop being exposed, and frameshifts are seen at
    1/10 of their rate (in-frame at 1/2). Reads mostly keep a hit (87 % with a frameshift, 91 % in frame), so the loss is local. A key that joins
    reads by the 12 bases before the window does not look at the window, which is the rule of §3.1.
  - **Sensitivity is 64 % at every length**, so standard additions calibrate it, and they need no second search: loci are placed by the undosed
    reads and a dosed copy is scanned for their keys. P2's `keep/` (undosed reads and hit tables) is therefore enough.
  - **ERR10889147** (535,458 reads, 4.9 M positions covered): 45 frameshift calls (38 one-base deletions) and 15 in frame. 9.3 × 10⁻⁶ per base raw,
    1.0-1.5 × 10⁻⁵ after sensitivity and the mock's false-call rate, in Illumina's expected range. At P2's 9 M reads, about 750 a library: a rate
    and a length split, not covariates.
  - **Not shown:** real indels' truth (P2: the CIGARs of `keep/`'s alignments), the false-call rate on real reads, strain indels (D6331), and
    homopolymer context. Scanning every read rather than those with a hit removes the last 13 % of selection.

**13.2 The codons no anchor labels: anchors set the level, the rest the shape.** Every exposed codon is already used: as a row that showed no
error (the denominator, and the Q / position / mate covariates), as the stack prior for a stop's source codon, and in the overlap and the locus.
What is not used is a plain disagreement with the column consensus, which biology swamps at Illumina rates (§3.6). Three ways to use more:
1. **Within-locus contrast (Tier 3; prototyped 2026-11-25, [pfam_contrast_proto.py](pfam_contrast_proto.py)).** Biology at a site does not depend
   on the cycle, Q or mate it was read at; error does. So `mismatch with the group's consensus = b(group) + e(covariates)`, a fixed effect per group,
   fitted by alternating means. e is identified up to a constant: **the contrast gives the shape, the anchors the level** (stops for ten ops,
   overlap and locus for `A>G` / `T>C`). Q stays a covariate, never a label.
   - **Per column it fails.** 36 % of bases differ from their column's consensus (many genomes per column), and that variance swamps the error:
     on the mock the cycle curve was ± 1 × 10⁻³ against a signal of 7 × 10⁻⁴, and on ERR10889147 the dose came back at +2 × 10⁻² in one Q bin.
   - **Per locus it works** (column, codon base, strand and 12 bases each side; ≥ 3 reads). Mock (150,000 pairs, three genomes, truth from the
     alignment), level from the truth: Q curve from 6.5 × 10⁻² (Q 15–19) to 3.90 × 10⁻⁴ (Q ≥ 38) against 6.6 × 10⁻² and 3.85 × 10⁻⁴; cycle from
     8.9 to 2.3 × 10⁻⁴ against 8.8 to 1.9; **all twelve ops, `A>G` / `T>C` included**, within their sd in three Q bins (`A>G` 3.6 × 10⁻², 2.4 × 10⁻³,
     2.1 × 10⁻⁴ against 3.7 × 10⁻², 2.4 × 10⁻³, 2.0 × 10⁻⁴). The locus keeps 41 % of exposed bases, the same share in every Q bin, and their error
     is that of all exposed bases (5.1 against 5.3 × 10⁻⁴).
   - **ERR10889147** (535,458 reads, 1.4 M bases in 332,874 loci): +0.01 comes back flat over Q within ± 1.5 × 10⁻³, the most in the sparse
     low-Q bins. R1 and R2 agree to 3 × 10⁻⁶ ± 2 × 10⁻⁵; the stops' R1 2.8 / R2 4.6 × 10⁻⁴ (52 stops) is within their counting error. With the
     stops' 3.7 × 10⁻⁴ as level, the Q curve runs 9.6 × 10⁻³ (Q < 15), 2.1 × 10⁻³ (20–24), 4.3 × 10⁻⁴ (30–34), 2.3 × 10⁻⁴ (≥ 38).
   - **Not shown:** that the fixed effect removes real strain variation (the mock has none: there, the plain locus mismatch did as well), that
     errors clustered in real reads do not bias it through the clean flanks (simulated errors may cluster less), and mate on the mock (98 % of
     its loci hold one mate only, as simulated). D6331 in P2's `keep/` tests the first two.
2. **The locus mixture as Tier 1c in the same fit** (§12, recommendation 1). It is checked by the stops it overlaps: the two agree on the mock
   (5.02 against 5.03 × 10⁻⁴).
3. **Contradicting reads give a stop's source codon** (§12, recommendation 3): 119 of 140 stops in the mock, 23 of 48 in ERR10889147.

Context is a fourth use: stops see only stop-adjacent contexts, and the contrast and the locus see all, so together they test whether the stops'
rate is genome-wide (§5).

**Where they land.** P3: 13.1, and 13.2 (1) scored on P2's kept hit tables against P2's truth (D6331 for strain variation). P4: 13.2 (1) as Tier 3 in the fit, if
it beats the anchors alone on held-out data.
