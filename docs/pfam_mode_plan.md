# Plan: a `pfam` evidence mode (minimal reference)

Status: **explored 2026-11-24; P0 and P1 done (results in §7, P1), P2 onwards not started.** P0 is [pfam_proto.py](pfam_proto.py);
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
- **Tier 3: covariate contrast (P4).** A logit with a fixed effect per column. It adds the shape of `A>G` / `T>C` over Q and
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
| indels | long reads only |
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

### P2: mock communities and an isolate, with truth (HPC, §8) - **go / no-go**
- Six ZymoBIOMICS runs already have kfpfam results (`nar-review-response/hpc/B11_pfam_arms/runs.tsv`: ERR2984773, SRR12324253, SRR12829160, ERR2935805,
  SRR23926920, SRR32629731; two are about 100 bp, where yield is uncertain), plus SRR24523812, the isolate with all three modes fitted (phases 5–6).
- `reference` mode on the same reads labels **each exposed codon's truth**. So each counted stop is classified: sequencing error, a stop
  in the genome (the real floor), or misplaced read (the real false-homology rate).
- **Anchor test:** with both anchors in one fit against each alone, on the runs whose inserts overlap. Keep the overlap anchor only if it
  measures the floor (stops in both mates) or the two missing ops better than the mocks' truth allows by chance.
- **Go** if, on at least five of seven runs: total rate within 15 % of `reference` on the same codons; each identified op within 25 %;
  Q and mate curves inside `reference`'s intervals; floor ≤ 2 × 10⁻⁴ per base.
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

Everything below is to be **written, tested on a laptop (`-profile test`) and handed over** under `workflows/pfam/`. It follows `workflows/hpc/`: one Nextflow
workflow, a Singularity image built on the laptop, `site.config` copied and edited, `run.sh` submitted from a run directory on Slurm.

| Process | Work | Image | Resources |
|---|---|---|---|
| `SUBSAMPLE` | every Nth pair of each raw FASTQ (no fastp, no merging) | seqkit | 1 cpu, 4 GB |
| `DOSE` | standard-addition copies at 0, 0.002 and 0.01 (undosed first) | this repo | 1 cpu, 4 GB |
| `ORFS` | six-frame ORFs, `pfam_proto.py orfs` | this repo | 1 cpu, 8 GB |
| `FILTER` | `kfpfam filter join`, panel index, `FILTER_MIN_GAPPED=21.25` (or r999's bar) | kfpfam-filter | 8 cpu, 8 GB |
| `PAIRS` | `recruit/pairs.py` | pyhmmer | 8 cpu, 16 GB |
| `EXPOSURE` | `pfam_proto.py exposure` (later `sources/pfam.py`); output is a small table | this repo | 1 cpu, 8 GB |
| `REFERENCE` | minimap2 / minibwa to the mock's genomes, truth per exposed codon (P2) | sem-tools | 8 cpu, 16 GB |
| `FIT` | EM, curves, report JSON | this repo | 1 cpu, 8 GB |

- **Cost, from a laptop run:** 1.2 core-ms per read for the filter plus HMMER on all of Pfam (480 CPU-s for 400,000 reads, ORF calling and exposure included).
  Three doses of 5 M pairs (10 M reads) is about 10 core-hours per run on an M4, 140 for 14 runs. **kfpfam's README warns that Mode A cost about 13x its
  laptop figure per read on Codon's x86 nodes**, so budget about 2,000 core-hours, or less with a panel index (§7, P1). Use `r999` or Mode A's bar, not `top`.
- **Inputs to supply on HPC:** raw FASTQs (ENA paths are in B11's `samplesheet.csv`; fetch on a login node if compute nodes cannot reach ENA), the Pfam
  index, `prep/` and `gapped.txt` already built for B11 (39.0), the mocks' genome FASTAs (B9 has the truth references), and the shovill assembly.
- **Pfam release:** the prototype used 38.2; the HPC indexes are 39.0. Rebuild or pin one, and record it in the spec provenance.
- **Returned:** aggregated count tables and report JSONs, megabytes at most (the per-codon table is ~50 MB per 400,000 reads and stays on HPC), for local fitting and plots. Nothing else.
- **Not prepared yet:** the workflow. It is P2's first task and needs P1, the parquet question below, and the site settings.

## 9. Risks

| Risk | Where it bites | Guard |
|---|---|---|
| **Selection on error**: recruitment drops error-rich reads | P1. Already cost a factor of 2 (3.2) | the standard-addition run, end to end, in every report |
| **False homology** at low T1 | T1 25 gives 6x the stops | the T1 plateau |
| **Biological floor** (pseudogenes, retroelements, recoding, phage) | Illumina below ~3 × 10⁻⁴ | panel; family random effects; P2 truth; per-stop excess test |
| **Host reads** | human LINE-1 in a gut sample | host depletion first, and the panel |
| **Stop-adjacent context only** | a rate that is not genome-wide | reported as such; compared on shared contexts |
| **Reads < 90 bp, DNBSEQ** | no exposure | yield table; refuse under a threshold |
| **Wrong source-codon prior** | per-op split (3.3) | stack prior; per-op dose-response |
| **Scooped** | novelty not checked deeply | P5 |

## 10. Open questions

1. **The hmmalign parquet:** where is it, and does it have GA-only or all hits, nucleotides or residues?
2. **Go / no-go thresholds** in P2 are proposals. They should be agreed before the run, not after.
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
