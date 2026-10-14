# Skiver paper: summary, review and critique

Review of Gu, Sharma, Wong & Nagarajan, *"Skiver: Reference-free quality control of
metagenomic sequencing datasets using (k, v)-mer sketches"*, bioRxiv
[10.64898/2026.02.12.705514](https://doi.org/10.64898/2026.02.12.705514) **v2** (posted
2026-06-04). Converted text: [skiver_paper.md](skiver_paper.md).

Written for this project, so the emphasis is: **to what extent does the paper show that
long-tail biological variation in metagenomes can be separated from sequencing error?**
Section 3 answers that; sections 4–6 cover the rest of the review and what it means for
`docs/IMPLEMENTATION_PLAN.md`.

> Caveat on the numbers in §3.2 and §4: the paper does not report its own detection limits,
> so the tables below are computed from the paper's stated model (Algorithm 1, Note S1) at
> the paper's own estimated error rates. They are implications of the method as described,
> not results the authors claim or dispute.

---

## Verdict up front

**Skiver does not separate long-tail biological variation from sequencing error, and the
paper does not claim to have measured whether it does.** What it does is estimate the
sequencing error rate on the subset of the data where the long tail has been *removed* —
keys with coverage ≥ 10 that carry no detectable second value — and the paper gives good
evidence that this estimate is not polluted by the *high-frequency* variation inside that
subset (diploid het sites, two-strain mixtures, repeats).

Three things follow:

1. **For abundant variation the separation is real and well evidenced.** The strongest single
   piece of evidence is one the authors never point at: on HG002 PacBio HiFi, skiver's
   substitution rate (≈ 6 × 10⁻⁵/base) is **8× below human heterozygosity** (≈ 5 × 10⁻⁴/base
   of alt-allele observations), and its Ti/Tv ratio is 0.33 — transversion-dominated, the
   *opposite* of biological variation's ≈ 2. Het sites are not leaking in. See §3.3.
2. **For rare variation it fails silently, and this is structural, not incidental.** The
   outlier filter is a one-sided binomial test at p₀ = 10⁻¹². At realistic key coverage it
   cannot see a minor value below roughly 10–20% allele fraction; everything below that is
   counted as sequencing error. See §3.2 for the detection-limit table.
3. **On real metagenomes the reported error rate and plausible leaked diversity are the same
   order of magnitude**, and the paper has no ground truth for any of its three real
   metagenomes. Soil Illumina: skiver reports 0.446%/base substitutions; within-population
   nucleotide diversity in soil bacteria is commonly 0.1–1%. Nothing in the paper
   distinguishes these.

The honest summary of the method's scope is narrower than the title: it is a reliable
**run-level QC** estimator computed from the high-coverage head of a sample. Whether the
head's error profile transfers to the tail is a separate question — and a live risk for us,
because error depends on context and GC, and the tail differs from the head in both.

---

## 1 What the method actually does

Four mechanisms, worth stating plainly because the separation question turns entirely on the
third and fourth.

1. **(k, v)-mer sketch.** Take every (k + v)-mer (default k = v = 17) from the reads and
   their reverse complements, subsample by `FracMinHash` on a hash of the *key* (first k
   bases), and bucket by key. Within a bucket, the most frequent *value* (last v bases) is
   the **consensus**. The key is a positional identifier; the values are observations of the
   same genomic position.
2. **Survival framing.** T = bases from a random start until the first disagreement with the
   consensus. From per-key counts N_{K,t} (observations matching the consensus up to base t),
   ĥ(t) = 1 − ΣN_{K,t}/ΣN_{K,t−1}. This is observable *only* for t ∈ (k, k + v] = (17, 34].
   Fit a discrete Weibull S(t) = exp(−λt^β) by complementary-log-log regression on those
   points, then report ĥ(1) = 1 − exp(−λ̂) as the per-base error rate.
3. **Coverage gate.** Keys with fewer than `N_key = 10` observations are discarded (0 when a
   reference is supplied).
4. **Outlier filter (Algorithm 1).** Iteratively: fit λ̂, β̂ on surviving keys; drop any key
   whose N_{K,t} is implausibly low under Binomial(N_{K,t−1}, 1 − ĥ(t)) at p₀ = 10⁻¹²; refit;
   repeat to convergence.

Gates 3 and 4 are the *entire* mechanism by which biological variation is kept out of the
error estimate. There is no variant model, no allele-frequency estimate, no explicit
partition of the sketch into "variant" and "error" — a key is either used whole or discarded
whole.

The error **spectrum** is estimated separately: enumerate the ≤ 11v distance-1 neighbours of
each consensus value, count how often each appears, assume the composition π_e is independent
of t, and normalise. Neighbours reachable by more than one edit are marked **Ambiguous and
discarded**.

---

## 2 What the paper establishes

Credit where due; these are solid and several are genuinely useful to us.

- **Mapping-based error profiling is badly biased without the exact strain.** Two clean
  demonstrations: on the Zymo Log mock community, mapping to a single isolate inflates the
  error rate because foreign reads mis-map to it (§3.3, Fig 3.A); and on simulated *E. coli*
  reads, the estimated error rate **nearly doubles** as reference ANI falls to 96%
  (Fig 3.B). The circularity argument is sharp — tightening the alignment score threshold to
  suppress spurious mappings requires knowing the error rate you are trying to measure.
- **Quality scores alone are unusable as an error-rate estimator.** seqtk's `ErrQ` is wrong
  by orders of magnitude in both directions (Table S2: 5 × 10⁻⁹ on Zymo Gut HiFi against
  skiver's 1.1 × 10⁻³; 0.0010 on Human R10.4 against 0.015). This directly supports the
  project's no-Q-as-truth rule in `CLAUDE.md`.
- **The hazard rate is the right object, and β < 1 is a real finding.** Error rate alone
  cannot distinguish 4 consecutive errors per 16 bases from 1 error per 4 bases (Fig 1), and
  the fitted β̂ is well below 1 on Nanopore and Illumina (Fig S1: 0.595 for R10.4, 0.717 for
  *B. subtilis* Illumina). Errors cluster; read-level error rate is heterogeneous. The
  constant-hazard ablation is consistently worse on Nanopore (Table 3).
- **The consensus-from-reads strategy works where it is testable.** ĥ(t) with and without a
  reference are nearly identical (Fig 5.A), and the spectra agree (Fig S6). This is a real
  result: the most frequent value recovers the reference base.
- **Alignment-free spectra avoid aligner scoring bias.** On Badread data with a known uniform
  spectrum, skiver recovers ≈ 1/3 each for ins/del/sub while Minimap2 + BEST over-reports
  substitutions (Fig 4). This is the paper's cleanest single result, because the ground truth
  is the simulator's own parameters rather than another aligner.
- **Phred calibration without circularity.** §2.5 estimates h(t | Q) as a conditional hazard
  from consensus disagreement — Q is a covariate, never a source of error labels. Finding:
  low Q (< 10) is pessimistic, high Q is over-confident.
- **Cheap.** ~10× faster than Minimap2 + BEST, ~13× less memory than KMC + GenomeScope
  (Tables S8–S9), single-threaded.

---

## 3 The focal question: long-tail biological variation vs sequencing error

### 3.1 Where the separation happens, and what it is made of

Variation is excluded by three gates, each with a different and mostly unexamined failure
mode:

| Gate | Removes | Fails on |
|---|---|---|
| `N_key ≥ 10` coverage threshold | every key from a genome too rare to reach 10 observations | nothing *leaks* — but the entire long tail is **excluded from the estimate by construction**, so no statement about it is being made |
| `k = 17` key uniqueness | — | keys shared between related genomes collapse non-homologous positions into one bucket (the paper's own explanation for why unfiltered skiver over-estimates even on a single isolate, §3.4) |
| Algorithm 1 outlier filter | keys with a *statistically detectable* second high-count value | any minor value below the detection limit — see §3.2 |

Note what gate 1 means for the title's claim. In a soil metagenome at 33 Gb, the organisms
reaching ≥ 10× per-key coverage are the most abundant handful. (FracMinHash subsamples
*distinct keys* by key hash, so per-key multiplicity is preserved — subsampling costs
power in the number of keys, not coverage per key. The threshold still binds on abundance.)
The number skiver reports is therefore the error rate of the sample's **head**, presented as
the error rate of the dataset. The paper's own stated applicability criterion says as much —
"at least one genome in the read set has coverage exceeding 20×" (§3.4) — but the Results
and the title do not carry that qualification forward.

### 3.2 The filter's detection limit

Algorithm 1 flags a key when `Pr[X ≤ N_{K,t}] < 10⁻¹²` under `Binomial(N_{K,t−1}, 1 − ĥ(t))`.
A second value at allele fraction *f* pushes N_{K,t} down to ≈ (1 − f)·N. Solving for the
minimum key coverage N at which that deficit reaches p₀ = 10⁻¹²:

| ĥ(t) (platform) | f = 0.5 | f = 0.2 | f = 0.1 | f = 0.05 | f = 0.02 | f = 0.01 |
|---|---|---|---|---|---|---|
| 0.0015 (PacBio HiFi) | 9 | 33 | 76 | 211 | 776 | 2451 |
| 0.0040 (Illumina) | 13 | 43 | 116 | 310 | 1525 | 7950 |
| 0.0150 (Nanopore R10.4) | 17 | 73 | 236 | 971 | 32476 | > 2 × 10⁶ |
| 0.0900 (Nanopore, Zymo) | 45 | 443 | 41996 | > 2 × 10⁶ | > 2 × 10⁶ | > 2 × 10⁶ |

Read this against the default `N_key = 10` and against realistic per-key coverage for a
non-dominant organism (tens):

- **A 50/50 split is always caught.** Diploid het sites, a two-strain 1:1 mixture, a
  two-copy repeat — all comfortably inside the detection region. This is exactly the regime
  the paper tests, and it passes.
- **A 10% minor allele needs ~116× key coverage on Illumina and ~42,000× on high-error
  Nanopore.** At metagenomic coverage this is unreachable for anything but the dominant
  organism.
- **Below ~5% allele fraction the filter is effectively blind at any realistic coverage**, and
  every such observation is counted as a sequencing error.

The long tail of biological variation in a metagenome *is* the low-allele-fraction regime:
rare strains, minor variants within a species population, recent mutations. The filter is
structurally unable to see it. The paper's framing — the filter "effectively suppresses
confounding signals arising from repetitive and shared regions" (§3.4) — is accurate and
also narrower than it sounds: repeats and shared orthologous regions between two
near-equal-abundance strains are high-*f* phenomena.

A second consequence, relevant to us: the filter discards the **whole key**, i.e. up to 34
bases of observation, not just the variant position. Regions of elevated diversity lose keys
wholesale, so the surviving data is biased toward conserved sequence — which differs from
variable sequence in context and GC composition. For a context-conditioned error model this
is a sampling bias on the context distribution itself, not just a loss of power.

Third: the filter conflates three different causes of a low N_{K,t} — polymorphism, repeats,
and **consensus failure** (§4, C1) — and removes all three without distinguishing them. So
`passes_filter` in `analyze.kvmer.csv` is not a variant indicator, and the filtered fraction
cannot be read as a diversity measure.

### 3.3 The strongest evidence — which the authors do not make

Two internal consistency checks on HG002 support the claim that high-frequency variation is
genuinely excluded. Neither appears in the paper.

**(a) The HiFi substitution rate is below human heterozygosity.** From Table S2 and Fig S8,
HG002 PacBio HiFi: ĥ(1) = 1.487 × 10⁻³ with substitutions only 4.0% of the spectrum, so the
per-base substitution rate is ≈ 6 × 10⁻⁵. Human heterozygosity is ≈ 1 per 1000 bp, and at a
het site roughly half the reads carry the alt allele, so unfiltered het sites would contribute
≈ 5 × 10⁻⁴/base — **8× more than the total substitution signal reported.** The reported rate
being an order of magnitude *below* the genome's own heterozygosity is a non-trivial check
that het sites are being removed rather than absorbed. It also agrees with §3.2: at f = 0.5
and HiFi coverage, detection needs only N ≥ 9.

**(b) Ti/Tv discriminates error from variation, per platform.** Computing transition/
transversion ratios from the Fig S8 matrices (per-type, accounting for 4 transitions vs 8
transversions):

| HG002 platform | sub / ins / del | Ti/Tv (per-type) |
|---|---|---|
| PacBio HiFi | 4.2% / 39.5% / 56.3% | **0.33** |
| Nanopore R10.4 | 31.0% / 31.9% / 37.2% | 2.54 |
| Nanopore R9.4 | 30.0% / 30.2% / 39.8% | 9.76 |

Biological variation sits at Ti/Tv ≈ 2.0–2.1. HiFi's 0.33 is transversion-dominated — clearly
not residual heterozygosity. R9.4's 9.76 is far *above* the biological value, so nanopore has
a genuine intrinsic transition bias (A↔G, C↔T), which means R10.4's 2.54 — superficially
suspicious because it is so close to the germline value — most likely sits on that platform
trend rather than indicating leakage. **Ti/Tv is therefore a usable leak diagnostic, but only
after calibrating the platform's intrinsic ratio on a clonal isolate.** That is a cheap check
we should adopt (§5).

### 3.4 What is not shown

The gap between what is demonstrated and what the abstract claims ("reliably handles complex
datasets containing multiple strains, alleles, and repetitive regions") is wide.

- **The only controlled test of strain variation is a single two-strain simulation.**
  Fig 3.E: *E. coli* K-12 MG1655 + O157:H7, 98% ANI, **fixed 64× total coverage**, mixing
  proportion swept 0→1. That is one ANI, one pair, one total coverage, two strains. A
  metagenome's long tail is many strains, a power-law abundance distribution, a continuum of
  ANI, and coverage straddling the threshold. Also, because *total* coverage is fixed, the
  minor strain's own coverage falls toward zero at the extremes — so the sweep conflates
  "filter removed it" with "it fell below `N_key`", and the paper does not separate the two.
- **A 10-species mock community does not test strain variation at all.** Different species do
  not share 17-mers, so there is nothing for the filter to do. Confirmed by the paper's own
  numbers: on Zymo Log, ĥ(1) is 0.088681 with the filter and 0.089319 without — a 0.7%
  difference. The filter's real tests are only Zymo Gut (multiple co-existing *E. coli*
  strains, at unknown frequencies) and HG002 (diploid + repeats).
- **The filter makes things worse on the most metagenome-like dataset with ground truth, and
  this goes unremarked.** Table 3, Zymo Log Nanopore: MSE_S is 1.47 × 10⁻⁴ with the filter
  and **2.51 × 10⁻⁵ without** — 5.9× worse with it. Table S7 shows the same for MSE_h
  (1.75 × 10⁻⁵ vs 9.97 × 10⁻⁶). The paper's own table *bolds the no-filter column* in that
  row while the text says the filter "substantially improves robustness". Both are true of
  different datasets, but the narrative is stated unconditionally.
- **No ground truth for any real metagenome.** The three actual microbiome datasets (human
  gut, soil, water) appear only in Table S5, as skiver vs skiver-no-filter. The soil number
  moves from 6.565% to 0.448% — a **14.7× swing** attributable entirely to the filter — and
  the caption asserts the unfiltered value is "biased heavily upwards" without evidence for
  where the truth lies. 0.448% is plausible-ish for NovaSeq but at the top of the usual range;
  the within-population nucleotide diversity it would have to be distinguished from is the
  same order of magnitude. "The filtered number looks reasonable" is not validation.
- **A weak internal signal worth noting in skiver's favour.** All three Illumina estimates
  cluster tightly: *B. subtilis* isolate 0.400%, soil 0.448%, water 0.462%. If long-tail
  diversity were leaking heavily into the metagenomes you would expect them well above the
  clonal isolate; they are 12% and 16% above. That is mild evidence the filter is working —
  or evidence that ~0.4% is a characteristic scale of the method on Illumina rather than a
  property of the data. The two readings are not distinguished, and the second is made more
  plausible by C7 below.
- **No sensitivity analysis on p₀.** 10⁻¹² is asserted once and never swept, despite being the
  single parameter that sets the variation/error boundary (§3.2). No convergence proof for
  Algorithm 1's fixed point either.
- **The decisive experiment is absent**: spike known variation into a clonal isolate at
  controlled allele frequencies and measure how much of it lands in ĥ(1) and in the spectrum.

### 3.5 Verdict

| Regime | Separated from error? | Evidence |
|---|---|---|
| Diploid het sites, high coverage | **Yes** | HG002 substitution rate 8× below heterozygosity; Ti/Tv wrong sign for variation (§3.3) |
| Two strains, near-equal abundance, ≥ 64× | **Yes** | Fig 3.E |
| Repeats / shared orthologous regions | **Yes**, within the head | ablation Tables 3, S6, S7 on human data |
| Minor strains below ~10% allele fraction | **No** — counted as error | §3.2; no experiment in the paper |
| Rare variants in a species population | **No** — counted as error | §3.2; no experiment in the paper |
| Genomes below the coverage threshold | **Excluded, not separated** | `N_key = 10`; §3.4 criterion |

The paper shows that *abundant* biological variation can be separated from sequencing error
in high-coverage data. It provides no evidence about the long tail, and the method's own
structure — a one-sided binomial test at p₀ = 10⁻¹² plus a hard coverage cut — implies the
long tail is absorbed into the error estimate. For a metagenome, "long-tail variation is
separated from sequencing error" is better stated as "long-tail variation is excluded from
the data the error rate is estimated on, or silently counted as error."

---

## 4 Other methodological criticisms

**C1 — Supplementary Note S1's union bound is stated incorrectly, and `N_key = 10` is below
the paper's own bound anyway.** The note concludes
`Pr[C* > C_w for all w] ≥ 11v · Pr[C* > C_w] = 11v · 0.99/(11v) = 0.99`. A union bound gives
`Pr[∃w : C* ≤ C_w] ≤ 11v · Pr[C* ≤ C_w]`, hence
`Pr[all] ≥ 1 − 11v·exp(−N(p₀−q)²/2)`; for 99% confidence the requirement is
`exp(...) ≤ 0.01/(11v)`, i.e. `log(11v/0.01)` in the numerator, not `log(11v/0.99)`. As
written the chain also derives ≥ 0.99 from a quantity that is ≤ 1, and the preceding line
drops a sign (`Pr[S_N − E[S_N] ≤ −N(p₀−q)]`). The D_j are claimed to take values in {−1, 1}
but actually take values in {−1, 0, 1}; Hoeffding with range 2 still gives the stated
exponent, so that part survives. The asymptotic `N = Ω((1−ε)^(−2v) log v)` is unaffected —
only the constant is wrong, roughly a factor of 2.

The bigger issue is that the bound, evaluated at the paper's own estimated error rates, does
not support the default threshold:

| Dataset (ĥ(1) from Table S2) | N required (as printed) | N required (corrected) | P(≥2 errors in a 17-base value) |
|---|---|---|---|
| Zymo Gut HiFi, 0.15% | 11 | 21 | 0.0% |
| *B. subtilis* Illumina, 0.40% | 12 | 23 | 0.2% |
| Human R10.4, 1.5% | 18 | 34 | 2.6% |
| Zymo Log Nanopore, 8.87% | 303 | 568 | 45.3% |
| *B. subtilis* RSII, 15.0% | 3931 | 7378 | 74.9% |

The default is `N_key = 10` — below the requirement for **every** dataset in the paper, and
below it by 1–3 orders of magnitude on the two high-error datasets. The paper's defence is
empirical ("we show that this is sufficient for the evaluated datasets"), which is fine as
far as it goes, but Note S1 is presented as theoretical justification for a threshold it does
not justify. The practical consequence is a **wrong consensus** on low-coverage keys at high
error rates: when the consensus is itself error-bearing, every correct observation at that key
is counted as an error. Note this interacts with the filter in skiver's favour — a wrong
consensus produces a large N_{K,t} deficit, so the filter tends to remove those keys too.
Part of what the filter is doing is cleaning up consensus failure, not polymorphism.

**C2 — ĥ(1) is a 17-position extrapolation, and the model-form choice moves it ~2×.** h(t) is
observable only on t ∈ (17, 34]; the headline per-base error rate ĥ(1) = 1 − exp(−λ̂) lies 17
positions below the lower edge of the fitting window, reached by a two-parameter Weibull fitted
by *ridge regression with Huber loss on a log-log transform* — not a maximum-likelihood fit of
the survival model, and the ridge penalty strength is never reported. Sensitivity is severe:
fixing β = 1 roughly halves ĥ(1) across the board (Table S6: Zymo Log 0.0887 → 0.0459, R10.4
0.0150 → 0.0066, RSII 0.1503 → 0.0905). The bootstrap intervals in Fig 5.A are sampling
intervals within the window and do not propagate this model-form uncertainty. The released
v0.3.2 also defaults `--ignore-smallest-t 2 --ignore-largest-t 2`, narrowing the already-small
window to 13 points — a detail absent from the paper.

**C3 — The ground truth is Minimap2, the method under criticism.** §3.2 defines ground-truth
S(t) by counting t-mers that align error-free in the Minimap2 BAM, while §§3.3–3.4 argue at
length that Minimap2 alignment is biased (scoring favours mismatches over indels, inflated by
divergence, soft clipping). The benchmark therefore cannot demonstrate accuracy — only
*consistency with alignment given a perfect reference*. Where skiver and Minimap2 disagree,
the disagreement is scored against skiver. The spectrum results escape this for simulated data
(Fig 4's truth is the simulator), which is why Fig 4 is the paper's strongest evidence; for
real data the spectrum claim rests on agreement with "previously reported error profiles",
which is much weaker.

**C4 — The MSE_S comparison handicaps the baselines by construction.** §3.2 forces
`Ŝ(t) = exp(−ε̂t)` on any method reporting a single number. skiver gets two parameters fitted
to the curve's shape; seqtk and GenomeScope get one number forced into an exponential they
never claimed. The resulting 10–100× MSE_S margins overstate the methodological gap. The fair
comparison for the *rate* is Table S2, where the picture is much tighter and sometimes
unflattering: on *B. subtilis* RSII, skiver = 0.1503, skiver-with-reference = 0.1141,
BEST = 0.1033 — skiver's own two modes disagree by 32% on a clonal isolate with a perfect
reference. On Human R10.4: 0.0150 / 0.0172 / 0.0285.

**C5 — The spectrum is estimated from the cleanest minority of the data, and discards
Ambiguous neighbours non-uniformly.** Only distance-1 neighbours of the consensus are counted,
so at RSII rates **75% of observed values carry ≥ 2 errors and are invisible** to the spectrum
estimator (45% on Zymo Nanopore — see the table in C1). The visible 25% is a biased subsample
conditioned on being nearly clean, which correlates with read quality — the very read-level
heterogeneity that β < 1 reveals. The `π_e` independent-of-t assumption is checked (Fig S3);
independence of *read quality* is not, and under heterogeneity the clean subsample's spectrum
need not equal the population's.

Separately, Ambiguous neighbours are **discarded** (20% of the neighbourhood in the Table S1
example), and ambiguity is a function of the consensus's own composition — homopolymers and
short repeats generate far more of it, because a deletion is modelled as delete-plus-append and
so collides with insertions and substitutions. So the method discards disproportionately many
events in exactly the contexts where PacBio and Nanopore indels concentrate, biasing π̂ toward
**understating indels in homopolymers**. Any context-conditioned indel model built on skiver's
spectrum inherits this.

**C6 — The covariate model is marginal and proportional-hazards, so its θ's cannot be
composed.** §2.5 estimates each θᵢ as the *mean of the difference* between
log(−log(1 − ĥ(t | Xᵢ))) and log(−log(1 − ĥ(t))) — a per-covariate offset against a shared
baseline shape, fitted one covariate at a time. Two problems for anyone wanting a joint model:
correlated covariates double-count (Q and read position are strongly correlated; so are Q and
GC), and the proportional-hazards assumption says a covariate shifts the hazard multiplicatively
without changing β. For Illumina the latter is visibly wrong — Fig S2 shows a *bathtub* hazard
(decreasing then increasing), so position changes the shape, not just the level, and Q is
confounded with position. **The marginal θ's in `summary_phred.csv` / `summary_gc_content.csv` /
`summary_read_position.csv` must not be multiplied together to build a joint model.**

**C7 — Illumina is skiver's weakest case, and the paper's own tables show it without saying
so.** *B. subtilis* Illumina is the one row in Table S3 where skiver is **not** best on MSE_h,
and where all four methods land within 12% of each other (seqtk 2.51 × 10⁻⁴, KMC
2.24 × 10⁻⁴ — the winner, skiver 2.33 × 10⁻⁴, skiver+ref 2.32 × 10⁻⁴). It is also the dataset
where the filter *hurts* MSE_S (Table 3: 6.80 × 10⁻⁴ default vs 5.59 × 10⁻⁴ no-filter). And
skiver's rate (0.0040) is 45% *below* BEST's (0.0072) on a clonal isolate with an exact
reference, where alignment should be trustworthy and the paper's aligner-bias explanation (an
indel-scoring artefact) is least applicable to a substitution-dominated 100 bp Illumina run.
The under-estimate's direction is what over-filtering, or a monotone Weibull extrapolated
backwards through a bathtub hazard (C2), would both produce. The text mentions only that
Illumina hazard rates are "slightly biased upward". **For a paired-end Illumina metagenome
project this is the most important finding in the paper, and it is the least foregrounded** —
but see C7a: the 0.0040 does **not** reproduce on the pinned v0.3.2 binary, which returns
0.0081 on the same dataset. The MSE_h and filter results above stand; the rate bias appears
to be specific to the benchmarked v0.2.2.

### C7a — Why Illumina reads low: ridge shrinkage of β removes the frailty lift

Tested against the paper's own Illumina dataset (ERR2935851, *B. subtilis*, HiSeq 1500,
101 bp) with the **pinned v0.3.2** binary. Three findings, the last of which explains the
paper's result.

**1. The paper's Illumina number reproduces. My first attempt did not — my sampling was at
fault.** Running the full R1 (23.16 M reads, 557×) gives **ĥ(1) = 0.004854, β̂ = 0.938**,
against the paper's v0.2.2 value of 0.0040 — a 21% gap, not the 2× I reported earlier from a
contiguous head slice. The rate is essentially **flat in depth**:

| subsample | depth | key median cov | ĥ(1) | β̂ |
|---|---|---|---|---|
| every 9th read | 62× | 38 | 0.005100 | 0.934 |
| every 4th read | 139× | 86 | 0.005189 | 0.929 |
| every 2nd read | 279× | 172 | 0.005104 | 0.931 |
| all reads | 557× | 343 | 0.004854 | 0.938 |

So no version fix is involved, and upstream's changelog agrees: between the benchmarked v0.2.2
and v0.3.2 there is **no change to the hazard fit or the per-base rate**. v0.2.1 added
`--ignore-{smallest,largest}-t` (so v0.2.2 already had the window trim), v0.2.2 changed
*Q-score* calibration to "clog-log regression with a common beta", v0.3.0 added output files,
v0.3.1 added FASTA input, and **v0.3.2 is a pure performance release** (`seq_io_parallel`,
`rapidgzip`, ~7× sketching speedup). **The Illumina bias was not accounted for upstream; it is
still present in the pinned binary.**

**2. Within-run read-block heterogeneity spans 2.5×, which is a trap for any harness.** At
*identical* read counts (2,573,803 reads, ~62×, key coverage 36–40), the answer depends
entirely on *which* reads:

| selection | ĥ(1) | β̂ |
|---|---|---|
| first 2.57 M reads (contiguous head) | **0.008081** | 0.918 |
| every 9th read across the run | 0.005100 | 0.934 |
| last 2.57 M reads (contiguous tail) | **0.003208** | 0.954 |

A 2.5× spread from flowcell position alone, and β̂ moves inversely with it (the dirtier block
has the lower β̂) exactly as the frailty story predicts. The per-base rate is therefore not a
property of "the run" but of the reads you fed in. **Any depth sweep, test fixture, or
cross-mode comparison built by slicing a FASTQ contiguously is biased by up to 2.5×**, and
cross-mode comparisons must consume the *same* reads, not different slices. My earlier 0.0081
was precisely this mistake.

**3. The paper's "skiver reads low on Illumina" is the ridge penalty discarding the frailty
lift.** `h(t)` is observable only on t ∈ (k, k+v], narrowed by v0.3.2's undocumented
`--ignore-{smallest,largest}-t 2` to 13 points (t = 20…32 at k = v = 17), and ĥ(1) is an
extrapolation ~19 positions below the window's lower edge whose value is set almost entirely
by β. On the full-R1 window (hazard 0.003979 → 0.003620):

| β | reported ĥ(1) |
|---|---|
| 1.000 | 0.00375 |
| **0.938** (skiver's ridge + Huber fit) | **0.00489** |
| 0.900 | 0.00576 |
| **0.831** (plain OLS on the same 13 points) | **0.00779** |
| 0.850 | 0.00717 |
| 0.800 | 0.00896 |
| 0.700 | 0.01413 |

**Plain OLS on skiver's own window returns 0.00779 — essentially BEST's 0.0072 — while
skiver's regularised fit returns 0.00489.** The ratio, 0.63, is the whole of the "skiver reads
low" effect, and it comes from a ridge penalty whose strength §2.3 never states. Ridge shrinks
the slope toward 0, i.e. β toward 1, and β = 1 applies *no* lift; since the window is measured
on a population already enriched for clean reads (that is what β < 1 means), removing the lift
under-reports. The headline number is therefore not reproducible from the published window
without an unreported hyperparameter, and the direction of its error is predictable.

**Why β controls the level.** `h(t)` is the chance base t is wrong *given the previous t−1 were
right*. With a homogeneous error rate that would not depend on t at all and β would be 1. It
depends on t because reads differ: conditioning on a long clean run progressively selects the
*cleaner* reads, so the surviving population gets cleaner as t grows and the hazard falls.
β < 1 is that selection. Consequently the hazard skiver can see, at t = 20…32, sits **below**
the unconditioned per-base rate; ĥ(1) is that unconditioned rate, and the extrapolation back
to t = 1 is the step meant to undo the enrichment. β is the only thing setting how much undo
you get.

**The aliasing mechanism is real but second-order.** A (k, v)-mer spans 34 bases — a third of
a 101 bp read — so larger t sits systematically later in the read, where the rate is higher,
pushing β̂ up against the frailty-driven decline. The measured per-cycle ramp is **6.4× from
read start to end** (from `summary_read_position.csv`, i.e. consensus disagreement, never
reported quality). Simulating the quantity skiver fits
([skiver_illumina_sim.py](skiver_illumina_sim.py)):

| scenario | true h(1) | skiver ĥ(1) | ratio | β̂ |
|---|---|---|---|---|
| 10 kb read, flat rate, no frailty | 0.00600 | 0.00600 | 1.00× | 1.000 |
| 10 kb read, frailty, no ramp | 0.00510 | 0.00707 | 1.39× | 0.899 |
| 101 bp, frailty, no ramp | 0.00510 | 0.00707 | 1.39× | 0.899 |
| 101 bp, ramp at ¼ strength | 0.00510 | 0.00480 | 0.94× | 1.001 |
| 101 bp, **measured** ramp | 0.00510 | 0.00181 | 0.35× | 1.263 |
| 151 bp, measured ramp | 0.00510 | 0.00292 | 0.57× | 1.126 |
| 301 bp, measured ramp | 0.00510 | 0.00476 | 0.93× | 0.992 |
| 10 kb read, measured ramp | 0.00510 | 0.00849 | 1.66× | 0.842 |

Direction confirmed — a ramp inflates β̂, read length monotonically rescues it (301 bp already
recovers 0.93×), and on 10 kb reads the ramp dilutes across the window and the bias reverses to
over-estimation. But full measured strength (0.35×, β̂ = 1.26) is far worse than reality, where
v0.3.2 returns β̂ = 0.938; the real data behaves like the ¼-strength row, because requiring a
clean 17-base key preferentially selects early cycles and low-error reads, so the ramp *as
sampled inside a (k, v)-mer* is much flatter than the cycle-level ramp. The simulation brackets
the mechanism; it does not pin it.

Also innocent, contrary to my first guess: Note S2's Taylor approximation recovers λ and β to
within 3–5% even at β = 0.6.

**Mitigations for us**, cheapest first: never slice a FASTQ contiguously (finding 2); treat β̂
rather than the filter as the main Illumina uncertainty and carry a sensitivity band; refit the
window without ridge shrinkage, or take the level from `kvmer.csv` single-edit counts directly,
to recover the lift skiver discards; use `--hazard-model constant` only as a floor (0.00375
here). The fit-free per-cycle mean in `summary_read_position.csv` is a useful anchor but a
lower bound, being conditioned on a clean neighbouring key.

**C8 — Reverse-complement double counting.** Step 1 sketches reads *and* their reverse
complements, so every read position is observed twice under two different keys. Equation 1
pools these as independent, and Algorithm 1's binomial test treats N_{K,t} as a binomial
sample. This roughly doubles the effective N, making the nominal p₀ = 10⁻¹² test more
aggressive than stated and tightening the detection limits in §3.2 by a corresponding factor.
`--forward-only` exists; the paper does not quantify the difference.

**C9 — The released tool exceeds the paper, and the paper does not document the release.**
v0.3.2 emits a per-key `passes_filter` flag, `homopolymer_length`, and a spectrum conditioned
on `prev_base`/`next_base` — i.e. the context-conditioned spectrum the paper marginalises away,
computed from a sketch that always carried the context. It also exposes `--outlier-threshold`,
`--use-all`, `--hazard-model` and the undocumented `--ignore-{smallest,largest}-t`. Meanwhile
the paper benchmarks **v0.2.2**, so its numbers do not characterise the v0.3.2 binary this
project pins.

---

## 5 Implications for this project

Tying to `docs/IMPLEMENTATION_PLAN.md` and the three evidence modes.

1. **Do not treat skiver's kmer mode as covering the long tail.** It characterises the
   high-coverage head. The plan already notes skiver needs high coverage; the sharper
   statement is that `N_key = 10` plus a p₀ = 10⁻¹² filter means the model is trained on
   *conserved sequence from abundant organisms*, which is a biased sample of contexts (§3.2),
   not merely a smaller one. Worth stating explicitly wherever the kmer mode's scope is
   documented.
2. **The cross-mode check we already have is the experiment the paper is missing.** pe-overlap
   needs no reference and no consensus, so it is not subject to the filter at all. Running
   both modes on the same Illumina metagenome and comparing ĥ(1) and the substitution spectrum
   gives a direct measurement of how much variation skiver's filter is absorbing — the
   spike-in experiment of §3.4, for free. If skiver reads systematically *higher* than
   pe-overlap on a diverse sample and they agree on a clonal isolate, that difference is the
   leak.
3. **Adopt Ti/Tv as a leak diagnostic, calibrated per platform.** Per §3.3 it separates
   biological variation (≈ 2.0) from HiFi error (0.33) cleanly, but nanopore's intrinsic ratio
   is itself 2.5–9.8, so the baseline must come from a clonal isolate on the same chemistry.
   Cheap to compute from `summary_error_spectrum.csv` and a genuinely informative check on any
   fitted model. This connects to the clonal point already landed on `phase6b-per-mate`.
4. **Expect ~30% disagreement between error-rate sources and set round-trip tolerances
   accordingly.** skiver's own two modes differ by 32% on a clonal isolate with a perfect
   reference, and by 2× against BEST on Human R10.4 (C4). A round-trip test demanding tighter
   agreement than that is testing the harness, not the model.
5. **Never compose skiver's marginal covariate effects.** Per C6, the per-covariate θ's are
   fitted marginally against a shared baseline shape and the Q/position confound is real for
   Illumina. Build the joint base × Q × context model from counts, not from skiver's summary
   offsets.
6. **Expect indels to be under-called in homopolymers** if any exporter's indel model is
   derived from skiver's spectrum (C5), and prefer v0.3.2's `homopolymer_length` and
   `prev_base`/`next_base` columns over the paper's marginal π_e (C9).
7. **On Illumina, β̂ is the uncertainty that matters, not the filter** (C7a). The pinned
   v0.3.2 returns a defensible rate on the paper's own Illumina dataset (0.0081 vs BEST's
   0.0072), so the paper's 2× under-estimate is a v0.2.2 artefact and not a reason to distrust
   the mode. But ĥ(1) swings 6× across a 0.3 range in β, and a 1.9× slice of that comes from
   an unreported ridge penalty — so carry a β-sensitivity band, don't quote a point estimate,
   and keep reference and pe-overlap primary for Illumina with skiver as the cross-check.
8. **`passes_filter` is not a variant flag.** It conflates polymorphism, repeats and consensus
   failure (§3.2, C1), so don't use the filtered fraction as a diversity estimate.

---

## 6 Experiments that would settle the focal question

In rough order of cost, all runnable with what the project already has:

1. **Spike-in titration.** Take a clonal isolate at high coverage, introduce known variants at
   controlled allele fractions (50%, 20%, 10%, 5%, 2%, 1%), and measure what fraction of each
   lands in ĥ(1) and in the spectrum. This directly measures the §3.2 detection curve on real
   data and is the single experiment the paper most needs.
2. **pe-overlap vs kmer on the same sample**, clonal isolate and diverse metagenome (item 2
   above). No new data required.
3. **Sweep `-e/--outlier-threshold`** across several orders of magnitude on a diverse sample
   and plot ĥ(1). A plateau would suggest the variation/error boundary is well separated; a
   smooth slope would show the reported rate is an artefact of where p₀ was placed.
4. **Sweep `-l/--lower-bound`** (the `N_key` threshold) and check whether ĥ(1) drifts with it.
   Drift means the estimate depends on which organisms clear the gate, i.e. the head/tail
   extrapolation is unsafe.
5. **Realistic multi-strain mixture**, many strains on a power-law abundance distribution
   across a range of ANI, with the minor strains' *absolute* coverage swept independently of
   total — separating "filtered" from "below threshold", which Fig 3.E conflates.
6. **GC-stratified error rate on a metagenome** using skiver's own GC covariate, to test
   whether the head's error profile transfers to the tail's sequence composition (§3.1).
