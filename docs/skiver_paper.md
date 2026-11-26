# Skiver: Reference-free quality control of metagenomic sequencing datasets using (k, v)-mer sketches

Zhenhao Gu, Puru Sharma, Limsoon Wong, Niranjan Nagarajan

National University of Singapore / Genome Institute of Singapore, A*STAR.

bioRxiv preprint doi: [10.64898/2026.02.12.705514](https://doi.org/10.64898/2026.02.12.705514),
**version 2, posted June 4, 2026**. CC-BY-NC-ND 4.0.

> Converted from the bioRxiv PDF (main text + supplementary material) for offline
> reference. Figures are not reproduced; their captions are retained inline so the
> figure references in the text resolve. The companion critique is
> [skiver_paper_review.md](skiver_paper_review.md).
>
> **Version note.** `docs/IMPLEMENTATION_PLAN.md` cites both v1 and v2. v1 was titled
> "Alignment-free Estimation of Sequencing Error Rates and Spectra using (k, v)-mer
> Sketches"; v2 is retitled to foreground metagenomic QC and adds the Phred-calibration
> results. This file is **v2**. The benchmarks in it were run with **skiver v0.2.2**, not
> the v0.3.2 release this project pins.

---

## Abstract

**Background.** Quality control of sequencing datasets is an important first step in
numerous bioinformatics pipelines such as mapping, variant calling, and assembly.
Existing methods typically rely on alignment results or quality scores. However, the
reference genome is not always available for mapping, and uncalibrated quality scores may
yield biased estimates of error rates.

**Results.** We present skiver, a reference-free and alignment-free algorithm that
estimates sequencing error rates and calibrates Phred quality scores using (k, v)-mer
sketches. By identifying the consensus from the sketched (k, v)-mers, skiver estimates
survival and hazard rates that capture positional information of sequencing errors. Across
simulated and real datasets from various sequencing platforms, skiver accurately recovers
sequencing error rates and the proportion of different error types. We further demonstrate
its ability to calibrate Phred scores. It also reliably handles complex datasets containing
multiple strains, alleles, and repetitive regions through an iterative outlier filtering
strategy. Skiver is computationally efficient and supports tools that need accurate
sequencing error rate estimates or quality scores as prior knowledge.

**Availability and Implementation.** An implementation of skiver is available at
<https://github.com/GZHoffie/skiver>, and dataset and scripts for reproducibility are
available at <https://github.com/GZHoffie/skiver-test>.

**Keywords:** sequencing error, survival analysis, alignment-free methods, sequence quality

---

## 1 Introduction

Given a set of sequenced reads, accurately estimating the overall sequencing error rate and
the spectrum of error types (substitutions, insertions, and deletions) is a fundamental
task in computational biology. These statistics are central to routine quality control
(QC). They are used to decide whether a run is usable, to compare runs across different
flow cells, chemistries, or library preparations, and to detect systematic failure modes
[1]. They also directly affect downstream inference. In particular, error profiles influence
alignment scoring and chaining heuristics in read mapping [2], and serve as key priors or
likelihood components in sensitive variant calling [3] and strain-aware analysis such as
phasing [4]. Because many pipelines implicitly assume a particular error model, biased
error rate estimates can propagate and lead to avoidable miscalls, unstable parameter
tuning, and misleading biological conclusions.

The challenge is amplified in metagenomic settings. Real samples often contain a mixture of
organisms with highly uneven abundance, including low-coverage genomes and closely related
coexisting strains. These properties can confound error rate estimation in two ways: (i)
low-coverage genomes and strain heterogeneity lead to many k-mers with low coverage, which
can be mistaken for sequencing errors, and (ii) reference genomes may be incomplete or
unavailable because organisms are missing from databases or differ substantially from
available references. Consequently, methods that rely on an incomplete reference genome
database can produce biased estimates precisely in the scenarios where robust error
profiling is most needed.

### 1.1 Previous work

Existing approaches for estimating sequencing error rates and spectra can be broadly divided
into reference-based and reference-free methods. Reference-based methods infer errors from
alignments to an external reference, whereas reference-free methods avoid mapping and
instead rely on internal read statistics such as k-mer frequency patterns or basecalling
quality scores.

**Reference-based methods.** A common and often effective strategy is to map the reads to a
reference genome and compute error statistics from the resulting alignments (e.g., via CIGAR
operations) [2, 5]. When a high-quality, closely matching reference is available, this
approach can provide accurate estimates. However, it has two major drawbacks. First, mapping
and alignment can be computationally expensive for large read sets and for long references,
making it costly as a routine QC step. Second, and more critically for metagenomics,
suitable references may be missing, incomplete, or diverge from the sequenced genomes
because organisms are not represented in databases or because strains differ substantially
from available references. In such cases, true biological differences are conflated with
sequencing errors, inflating apparent mismatch and indel rates and biasing the estimated
sequencing error rates upward [6]. One possible mitigation is to assemble the reads to
obtain a sample-specific consensus and then compare the reads to this consensus [7], but
assembly can be computationally intensive and may be unstable for low-coverage components
and complex mixtures.

**Reference-free methods.** To avoid dependence on an external reference, many approaches
estimate error rates directly from read-derived statistics. For example, shadow regression
[8] and SequencErr [9] utilize overlapping (paired-end) short reads to quantify
disagreements between reads. Another prominent family of methods leverages k-mer frequency
information to infer error-related quantities without explicit alignment [10–17]. While
scalable, these methods are often tailored to specific regimes (e.g., single-genome
assumptions or specific sequencing platforms) and can struggle to recover the full error
spectrum (i.e., distribution of substitutions vs. insertions vs. deletions). Moreover, the
severe coverage imbalance typical of metagenomic samples introduces a key pitfall:
low-coverage genomes yield sparse and highly variable k-mer counts, making frequency-based
inference unstable and causing rare error-free k-mers to be difficult to distinguish from
erroneous k-mers induced by sequencing noise.

Another class of reference-free methods uses Phred quality scores, which theoretically
satisfy Q = −10 log₁₀ Pr[error] and are used by many downstream tools [3, 18]. However, in
practice, quality scores can be miscalibrated and vary with sequencing technology [19] and
library preparation [20]. Some calibration procedures do exist [21], but require additional
assumptions and data. Consistent with these observations, our experiments show that
uncalibrated quality scores can substantially underestimate or overestimate true error
rates. Moreover, they do not directly provide reliable estimates of the full
substitution/insertion/deletion spectrum.

### 1.2 Our contribution

In this work, we propose skiver, which introduces a novel (k, v)-mer sketching approach for
reference-free and alignment-free profiling of sequencing error rates and spectra. The core
idea is to construct sketches that group (k, v)-mers sharing the same key and then find the
error patterns in the associated values. Using concepts in survival analysis, we use the
hazard rate and the survival rate to describe the frequencies and the distribution of
sequencing errors. Across simulated and real datasets from multiple sequencing platforms,
our method yields accurate hazard rate and survival rate estimates while remaining robust in
settings where reference genomes are unavailable or incomplete, including metagenomic
regimes with uneven coverage and strain heterogeneity.

---

## 2 Methods

### 2.1 Problem formulation

The definition of sequencing error rate varies across the literature. It is usually defined
as the number of errors divided by the alignment length [22]. In this paper, we borrow
concepts from survival analysis. Let T (T ≥ 1) be the random variable that is the number of
bases from a random starting position until the first failure (disagreement between the
sequenced base and the true underlying base).

> **Figure 1.** The alignment of two reads (bottom sequence) to the reference (top
> sequence). Both reads have an error rate of 1/4 but different hazard rates h(t) (right).

> **Definition 1.** The hazard rate at the t-th base from the starting position, denoted by
> h(t) := Pr[T = t | T ≥ t], is the conditional probability of the t-th base disagreeing
> with the underlying genome given that the previous (t − 1) bases agree.

The benefits of using the hazard rate instead of the error rate inferred from alignment
include the following. Firstly, alignment can be slow and ambiguous. Different alignment
scoring schemes may yield different numbers of aligned bases, resulting in varying estimates
using different mappers [23]. Secondly, the hazard rate offers a more comprehensive view of
the position-dependent error process. For example, imagine two reads in Figure 1, where read
A has 4 consecutive errors every 16 bases, and read B has 1 error every 4 bases. The error
rates of the two reads are the same and are equal to 1/4, but their hazard rates are
drastically different. Specifically, read A has a much lower hazard rate at small t,
indicating a higher chance of observing a long stretch of consecutive matches.

> **Definition 2.** The survival rate at the t-th base, denoted S(t) := Pr[T > t], is the
> probability of the first t bases from the point of observation being free of sequencing
> errors.

Note that the survival rate

$$S(k) = \prod_{t=1}^{k}\bigl(1 - h(t)\bigr)$$

also indicates the probability of a k-mer in the read being free of sequencing error. The
per-base error rate can be calculated by ε_per-base := h(1) = 1 − S(1), which is the
probability that a random base in the read is erroneous.

In this work, we aim to propose an efficient and accurate estimator for hazard rate and
survival rate using (k, v)-mer sketches, and show its ability to estimate the frequency of
each error type (substitutions, insertions, and deletions) at the same time.

### 2.2 (k, v)-mer sketches

> **Figure 2.** The workflow of skiver. **A.** The structure of a (k, v)-mer. **B.**
> Extracting (k, v)-mer sketches from the read set. (k, v)-mers with the same keys are
> gathered together to infer the consensus value and identify sequencing errors. **C.** The
> sketched (k, v)-mers go through a hazard rate outlier filter to exclude keys that are
> associated with multiple values with high counts. **D.** The hazard rate is calculated for
> each key and used to fit the parameters of a Discrete Weibull distribution. **E.** We find
> all the distance-1 neighbors of the consensus and infer the frequency of each type of
> sequencing error (insertion, deletion, and substitution).

> **Definition 3.** A (k, v)-mer (as shown in Figure 2.A) is a segment of DNA of length
> k + v, with the first k bases being the key and the last v bases being the value.

The structure of a (k, v)-mer is similar to the idea of anchor-target k-mer pairs in SPLASH
[24, 25], but different in that our key and value must be adjacent for accurate error
profiling, whereas the anchor-target pair can be separated in the reads. In this paper, we
focus on using this structure for sequencing error rate and spectra estimation.

A (k, v)-mer sketch of a set of sequenced reads is created with the following steps
(Figure 2.B):

- **Step 1.** We extract all the (k, v)-mers from the reads and their reverse complement.
  Optionally, only the forward strand of the reads is used.
- **Step 2.** We subsample roughly 1/c (k, v)-mers with FracMinHash [26]. In particular,
  given a random hash function H : Σᵏ → (0, 1), we subsample a (k, v)-mer if the hash value
  of its key is less than 1/c.
- **Step 3.** All the subsampled (k, v)-mers are stored in a hash table that maps a key to
  the set of associated values along with the number of times they appear in the read set.
  The most frequently appearing value is identified as the consensus.

Here, k is chosen to be large such that the keys are mostly unique in the set of sequenced
genomes. The keys are used as positional identifiers in the genomes. We then identify
variation in the associated values to determine possible sequencing errors or mutations. By
default, both k and v are chosen to be 17, which generally works well across different real
sequencing reads with various error rates (Supplementary Figure S9).

It is possible to show that given a per-base error rate of ε, when the number of times a key
appears in the read set is N_key = Ω((1 − ε)^−2v log v), the value with the highest count
(consensus) matches the true value from the sequenced genome with high probability
(**Supplementary Note S1**). In practice, we choose a default threshold N_key = 10, and only
use keys above this coverage threshold for subsequent tasks, and show that this is sufficient
for the evaluated datasets.

If a reference genome is provided, we also extract the set of (k, v)-mers from the genome
using the same hash function. If a key K is associated with multiple different values in the
reference genome, which indicates a repeat sequence, the key is discarded. Otherwise, the
unique value associated with the key is regarded as the consensus. Since the consensus is
obtained from the reference, a lower threshold for key multiplicity is set (N_key = 1) to
enable profiling of lower coverage read sets.

### 2.3 Estimating hazard rate and survival rate

In this work, we assume that T follows a discrete Weibull distribution with parameters λ and
β, which is often used for discrete-time survival analysis. This assumption comes from the
empirical observation that hazard rates h(t) in real sequencing datasets decrease with t
(Figure 5.A), and the survival rate S(t) fits well to the curve S(t) = exp(−λt^β)
(Supplementary Figure S1).

Our model differs from existing error models which typically assume a constant error rate
[15], essentially assuming β = 1. In real datasets, especially Nanopore and Illumina, the
best fit parameter β̂ is much smaller than 1. This indicates a decreasing hazard rate, and is
often interpreted as evidence for heterogeneity in the failure process [27], or heterogeneous
sequencing error rate across reads in our case.

Given a (k, v)-mer sketch of the sequenced reads, we count N_{K,t}, the number of (k, v)-mers
that have key K and match with the consensus up to the t-th base. For example, in the example
of Figure 2.B, 4 (k, v)-mers share the key K = TTACATTGGCAG. The consensus value is taken as
the value with the highest count, in this case, AGCG. Two of the (k, v)-mers differ from the
consensus in the 14th and 15th base, respectively. We therefore have N_{K,12} = N_{K,13} = 4,
N_{K,14} = 3, N_{K,15} = N_{K,16} = 2.

This process is repeated for all keys. Assuming the hazard rate at the t-th base is h(t), we
should have N_{K,t} ~ Binomial(N_{K,t−1}, 1 − h(t)). We take the maximum likelihood estimate,

$$\hat h(t) = 1 - \frac{\sum_K N_{K,t}}{\sum_K N_{K,t-1}} \quad \text{for } k < t \le k+v. \tag{1}$$

This allows us to estimate h(t) in a small interval between k and k + v. Under the assumption
that T follows a discrete Weibull distribution, we have for large t (**Supplementary Note
S2**),

$$h(t) = 1 - \exp\bigl(-\lambda\,(t^\beta - (t-1)^\beta)\bigr)$$
$$\Rightarrow \log(1 - h(t)) = -\lambda\,(t^\beta - (t-1)^\beta) \approx -\lambda\beta t^{\beta-1}$$
$$\Rightarrow \log(-\log(1 - h(t))) \approx \log(\lambda\beta) + (\beta - 1)\log t.$$

This transformation (complementary log-log) is widely used in discrete survival analysis
[28]. We can then perform a **ridge regression with Huber loss** of log(−log(1 − h(t)))
vs. log t in the range [k + 1, k + v]. Let a and b be the estimated slope and intercept.
Then, we have

$$\hat\beta = a + 1, \qquad \hat\lambda = \exp(b)/\hat\beta.$$

The survival rate can then be estimated by directly plugging the estimated parameters,
Ŝ(t) = exp(−λ̂t^β̂), and the sequencing error rate is estimated to be
ĥ(1) = 1 − Ŝ(1) = 1 − exp(−λ̂). This process takes time and space that are linear in v and
the number of keys in the sketch.

### 2.4 Estimating error spectra

To estimate the composition of errors (types of substitutions, insertions, and deletions), we
make an additional assumption that the error composition is roughly constant and is
independent of t. In other words, for a given single base error type e (such as substitution
C → A), the hazard rate of this type of error occurring can be expressed as

$$h_e(t) = \pi_e h(t),$$

where π_e is a constant value. The probability of the type of error e happening exactly once
in the value, while no other error is present, is

$$\Pr[\text{error } e \text{ happens exactly once in the value}] = \sum_{t=k+1}^{k+v} S(t-1)\,h_e(t)\,S(k+v-t) \propto \pi_e,$$

given that k and v are fixed. This assumption is based on the observation that the
proportions of each error type happening exactly once at position t are roughly constant
across t (**Supplementary Figure S3**). Under this assumption, we can estimate π_e by the
frequency of e happening exactly once in the value.

After identifying the consensus value for each key, we find the set of all pairs
(value, edit type), where the value can be obtained from the consensus via one edit of edit
type (such as A → C). An example of the neighbor set can be found in **Supplementary Table
S1**. If a neighbor can be reached via multiple types of edit, we mark the edit type as
**Ambiguous**.

We then count the number of times each distance-1 neighbor appears in the (k, v)-mer sketch.
In the example in Figure 2.E, the neighbors corresponding to the single base substitutions
C → A and G → C appear once respectively. This process is repeated for all keys, which gives
us the total number of times each edit type has appeared. The frequency of an edit type is
estimated to be its count divided by the sum of counts of all the profiled error types. **The
Ambiguous neighbors are not counted.**

There are at most 11v distance-1 neighbors given a consensus value. As a result, the time
needed to count the error frequencies is only O(v · #keys + #values), where #keys and #values
are the total number of keys and values in the (k, v)-mer sketch. This is much more efficient
than the case in which we align all values to the consensus, which can take time
O(v² · #values).

### 2.5 Estimating the dependence of error rate on covariates

We can also estimate the conditional hazard rates h(t | Xᵢ) assuming a Cox proportional
hazard rate model with the complementary-log-log link [28],

$$\log(-\log(1 - h(t \mid X_i))) = \log(-\log(1 - h_0(t))) + \theta_i X_i,$$

where Xᵢ represents the covariates such as Phred quality scores, GC-content of the read, or
the pore/well from which the read was generated. We choose the baseline hazard h₀(t) to be
the hazard rate ĥ(t) = 1 − exp(−λ̂(t^β̂ − (t − 1)^β̂)) estimated using all the data points. We
then calculate N_{K,t|Xᵢ} to be the number of bases in key K that satisfies Xᵢ and agrees
with the consensus up to the t-th base, and estimate

$$\hat h(t \mid X_i) = 1 - \frac{\sum_K N_{K,t\mid X_i}}{\sum_K N_{K,t-1\mid X_i}} \quad \text{for } k < t \le k+v.$$

The individual scale parameters θᵢ are then estimated using the **mean of the difference**
between log(−log(1 − ĥ(t | Xᵢ))) and log(−log(1 − ĥ(t))).

### 2.6 Dealing with repeats, multiple alleles or strains

In real datasets, it is common that the sequenced genomes contain long repetitive regions and
heterozygous sites in the case of the human genome, or multiple co-existing strains of the
same species in the case of metagenomic samples. In these cases, a key can be associated with
multiple values with high count in the reads, causing overestimation of the hazard rate
(Figure 3.E).

A key observation is that if a key K is associated with multiple values that have high
counts, N_{K,t} is going to be significantly smaller for some t ∈ [k + 1, k + v]. We
therefore use a simple outlier filter as shown in Algorithm 1 that filters all K that have a
significantly smaller N_{K,t}, assuming that N_{K,t} ~ Binomial(N_{K,t−1}, 1 − ĥ(t)).

If the reference genome is provided, this filter is disabled by default, as keys that are
associated with multiple values are already discarded.

**Algorithm 1: Outlier filter in hazard rate estimation**

```
Input:  N_{K,t} for all keys K and k < t <= k+v; threshold eps to indicate
        convergence of parameters; P-value threshold p0 (1e-12 by default).
Output: A set of keys that pass the filter.

 1  keep_key <- {}
 2  foreach key K do
 3      keep_key[K] <- true
 4  end
    /* Iterative outlier filtering */
 5  repeat
        /* estimate lambda_hat, beta_hat using the clog-log regression
           described in the previous section */
 6      lambda_hat, beta_hat <- estimate_weibull({N_{K,t} : k+1 <= t <= k+v,
                                                  keep_key[K] is true})
 7      for t <- k+1 to k+v do
 8          h_t <- 1 - exp(-lambda_hat * (t^beta_hat - (t-1)^beta_hat))
 9          foreach key K do
10              if Pr_{X ~ Binomial(N_{K,t-1}, 1-h_t)}[X <= N_{K,t}] < p0 then
11                  keep_key[K] <- false
12              end
13          end
14      end
15      lambda_new, beta_new <- estimate_weibull({N_{K,t} : k+1 <= t <= k+v,
                                                  keep_key[K] is true})
16  until |lambda_new - lambda_hat| < eps and |beta_new - beta_hat| < eps
17  return {K : keep_key[K] is true}
```

---

## 3 Results

### 3.1 Baselines and datasets

To test the ability of our algorithm to estimate h(t) and S(t), we benchmarked **skiver
v0.2.2** against popular state-of-the-art baselines.

For mapping-based tools, we used BEST [5] to analyze the BAM output from Minimap2 [2], with
ground truth reference genomes. By default, BEST uses only the primary alignment for each
read for read error profiling. For reference-free methods, we used GenomeScope2.0 [16] to fit
the k-mer histogram output of KMC 3.2.4 [33], and used the `Read Error Rate` field in the
summary file. For quality-score-based methods, we used seqtk [34] and the `ErrQ` field to
compute sequencing error rates, ε̂_seqtk = 10^(−ErrQ/10).

We selected a wide range of sequencing datasets from multiple sources and sequencing
platforms (Table 1) [29–32, 35], including mock bacterial communities, bacterial isolates,
human reads, and different microbiomes. All of the chosen datasets (except the microbiomes)
have well-defined reference genomes for fair comparisons.

All experiments were run on a server with 16 Intel(R) Xeon(R) Platinum 8259CL CPUs @
2.50GHz. All baseline tools were run with their default parameters.

**Table 1** Sequencing datasets used for various experiments in this study.

| Sample Type | Details | Platform | # of bases | Source |
|---|---|---|---|---|
| Mock Community | Zymo Log [29] | Nanopore GridION | 16.5 Gb | NCBI: ERR3152366 |
| Mock Community | Zymo Gut Microbiome | PacBio HiFi | 18 Gb | NCBI: SRR13128014 |
| Microbial Isolates | *Bacillus subtilis* [29] | Illumina | 4.7 Gb | NCBI: ERR2935851 |
| Microbial Isolates | *Bacillus subtilis* [30] | PacBio RSII | 1.6 Gb | NCBI: SRR7498042 |
| Human DNA | HG002 [31] | Nanopore R10.4 | 60.6 Gb | github.com/marbl/HG002 |
| Human DNA | HG002 [31] | Nanopore R9.4 | 64.2 Gb | github.com/marbl/HG002 |
| Human DNA | HG002 [31] | PacBio HiFi | 109.06 Gb | github.com/marbl/HG002 |
| Microbiome | Human Gut [32] | Nanopore PromethION | 5.4 Gb | NCBI: ERR7625321 |
| Microbiome | Soil Microbiome | Illumina NovaSeq | 33.1 Gb | NCBI: SRR11742949 |
| Microbiome | Water Microbiome | Illumina NovaSeq | 6.4 Gb | NCBI: SRR14560391 |

### 3.2 Evaluation metrics

In all of the experiments, we calculate the **ground truth** survival rate S(t) by counting
the proportion of t-mers in the reads that align with the reference without any error in the
BAM output of Minimap2. The ground truth hazard rate was calculated as h(t) = 1 − S(t)/S(t−1)
for t ≥ 2 and h(1) = 1 − S(1).

The main metric used to evaluate the tools was the mean squared error (MSE) of the estimated
survival rate Ŝ(t) vs. the ground truth S(t), calculated by

$$\mathrm{MSE}_S := \frac{1}{n}\sum_{t=1}^{n}\bigl(S(t) - \hat S(t)\bigr)^2,$$

where the upper bound n is chosen to be 100. This is similar to the definition of Brier
scores, which are used extensively in the evaluation of survival curve estimation [36].
Intuitively, this metric measures how well the tools predict the proportion of the t-mers in
the read set that are error-free for various t. Similarly, we define the MSE for hazard rate
to be MSE_h := (1/n) Σ (h(t) − ĥ(t))². Here, ĥ(t) = 1 − exp(−λ̂(t^β̂ − (t−1)^β̂)) and
Ŝ(t) = exp(−λ̂t^β̂) for skiver, while **ĥ(t) = ε̂, Ŝ(t) = exp(−ε̂t) for the other methods that
only report one error rate estimate ε̂**.

> **Figure 3.** Performance of various tools under simulated conditions. **A.** Profiled
> error rate for BEST + Minimap2 on the Zymo Log mock community dataset, if only one
> bacterial isolate is given as the reference (on the x-axis). The lines representing the
> abundance and the % reads mapped use the axis on the right of the figure. **B.** Profiled
> error rate for BEST + Minimap2 for a simulated read from an *E. coli* K-12 MG1655 strain,
> using different *E. coli* strains as references. **C.** Estimation of survival rates S(t)
> for simulated datasets of *E. coli* with 100× coverage. The dashed line represents the
> ground truth S(t) at different read identities (%), the solid line represents the mean of
> the estimates Ŝ(t), and the shaded area covers the range between the min and max of Ŝ(t)
> in the 20 simulated experiments. **D.** MSE_S of skiver estimates at different coverage
> values with simulated reads that have 96% read identity. **E.** Estimated error rates for
> skiver when mixing two *E. coli* strains together at different proportions, with and
> without the outlier filter.

### 3.3 Results for mapping-based methods may be biased by wrong or missing references

In this section, we evaluate the behavior of mapping-based methods under simulated conditions
in which the reference genome is missing or incorrect.

We first consider the scenario of an incomplete reference database, a common situation in
metagenomic analyses. Using the Zymo Log, Nanopore GridION dataset, we applied Minimap2 to
map the reads under two reference configurations: (i) using all 10 known genomes in the mock
community as references, and (ii) using each individual genome alone as the reference
(Figure 3.A).

Interestingly, the estimated error rate obtained using the complete reference set is the
lowest (approximately 12%), whereas estimates derived from individual reference genomes are
consistently biased upwards. Examination of the mapping results reveals that Minimap2 assigns
a disproportionately large number of reads to certain isolate genomes relative to their
expected abundances in the mock community. This indicates that many reads originating from
other species are incorrectly mapped to these genomes. As a result, the apparent sequencing
error rate is substantially inflated. These observations suggest that read mapping with an
incomplete reference database is inherently unreliable for error-rate profiling. While it may
be possible to reduce spurious mappings by applying more stringent alignment score
thresholds, the optimal threshold itself depends on the unknown sequencing error rate,
creating a **circular dependency** that limits the robustness of mapping-based approaches.

Next, we examine the case in which the reference genome is present but does not exactly match
the sequenced strain. Reads were simulated from *E. coli* K-12 MG1655 at 96% read identity
using Badread [37]. Error rates were then estimated using Minimap2 and BEST, with reference
genomes chosen from other *E. coli* strains sampled from GTDB [38]. The average nucleotide
identity (ANI) between each reference genome and *E. coli* K-12 MG1655 was computed using
skani [39], and all selected references were verified to have ANI values greater than 95%.
Despite this high sequence similarity, **the estimated error rate almost doubled as the ANI
between the reference genome and the true sequenced genome decreased to 96%** (Figure 3.B).
These results demonstrate that without an additional step of variant calling, even modest
strain-level divergence can severely bias error-rate estimates derived from read mapping.
Consequently, accurate profiling using mapping-based methods requires the exact strain, or a
near-identical genome, to be present in the reference database, a condition that is rarely
met in real metagenomic samples.

### 3.4 Skiver accurately estimates hazard rates and survival rates in simulated reads

Next, we tested the behavior of skiver on simulated reads.

We first evaluated skiver under a simplified setting in which the read set contains only a
single genome. Reads from the *E. coli* K-12 MG1655 strain were simulated across a range of
read identities (90% to 100%) and sequencing coverages (10× to 100×). We performed 20
experiments on each set of parameters.

Across these conditions, skiver accurately recovered both the survival rate and the hazard
rate. For example, at 100× coverage (Figure 3.C), the estimated survival curve Ŝ(t) (solid
line) closely matches the ground-truth S(t) (dashed line), exhibiting minimal bias and low
variance. When the read identity is fixed (Figure 3.D), the mean squared error of the
survival rate estimate (MSE_S) remains close to zero as coverage decreases, until coverage
drops to 10×. At this depth, the default skiver configuration shows increased variance,
whereas providing a reference genome restores accuracy and reduces variance.

Based on the complete set of simulations (Supplementary Figures S4-5), **skiver produces
generally reliable survival and hazard rate estimates when at least one genome in the read
set has coverage exceeding 20×, and the read error rate is below 6%.** For lower coverage or
higher error rates, a reference genome is required for an accurate estimate.

In addition to accurately estimating survival and hazard functions, skiver is also able to
recover the underlying sequencing error spectrum used by the simulator. In this experiment,
Badread was configured with a **random error model**, in which both the type of edit
operation (substitution, insertion, or deletion) and the resulting base after each edit are
selected uniformly at random. Consistent with this design, the error spectrum reported by
skiver (Figure 4) exhibits approximately equal frequencies across all single-base
substitution (SBS) types, while the overall frequencies of substitutions, insertions, and
deletions are each close to one third. These results closely match the parameters specified
in the simulator.

In contrast, when BEST is applied to Minimap2 alignments generated from the same dataset, a
higher proportion of substitutions is reported relative to insertions and deletions. This
bias is likely attributable to the alignment scoring scheme of Minimap2, which favors
mismatches over indels during alignment optimization. By comparison, skiver is an
alignment-free method and therefore remains unaffected by alignment heuristics or scoring
parameters, enabling a more faithful recovery of the true sequencing error spectra.

> **Figure 4.** Estimated sequencing error spectrum of a 128× coverage *E. coli* dataset
> with 96% read identity simulated using Badread, where all error types are equally likely.
> (k, v)-mer SBS frequencies are all ≈0.08–0.09; the insertion / deletion / substitution
> composition from (k, v)-mer sits at ≈1/3 each, matching ground truth, while
> BEST + Minimap2 over-reports substitutions.

Next, we evaluated the effect of the outlier filter in datasets containing multiple closely
related strains (Figure 3.E). We simulated reads with 64× coverage from *E. coli* K-12 MG1655
and *E. coli* O157:H7 at the same read identity. These two strains share an average
nucleotide identity (ANI) of **98%** according to skani [39]. The two read sets were then
subsampled and mixed to generate datasets with a **fixed total coverage of 64×** but varying
proportions of the two strains. For each mixture, we estimated the sequencing error rate
using skiver, defined as the first-base hazard ĥ(1), both with and without the outlier filter
enabled.

When the proportion of K-12 MG1655 is 1.0 or 0.0, the dataset effectively represents a single
isolate. In these cases, skiver without filtering reports overestimated error rates, largely
due to similar genomic regions that introduce spurious (k, v)-mer groupings. As the mixture
proportion approaches 0.5, the estimated error rate from skiver without filtering increases
sharply, reflecting the growing ambiguity between orthologous regions shared by the two
strains. In contrast, skiver with the outlier filter remains relatively robust across all
mixture proportions and consistently reports error rates close to the ground truth. These
results demonstrate that the outlier filter effectively suppresses confounding signals
arising from repetitive and shared regions in strain genomes, enabling accurate error-rate
estimation even in multi-strain datasets.

> **Figure 5.** Performance of skiver on real metagenomic datasets. **A.** Estimated hazard
> rate using skiver, with and without reference genomes. The shaded areas are the 90%
> confidence intervals for the hazard rate estimated via bootstrapping. **B.** Estimated
> survival rate using baselines and skiver. **C.** Estimated error spectra for the datasets.
> The last row represents insertions, and the last column represents deletions, while the
> rest of the entries are substitutions. **D.** Estimated true error rate of the bases with
> respect to their associated Phred scores.

### 3.5 Skiver generalizes well to reads from various sequencing platforms

We next evaluated the performance of skiver and several baseline methods on real sequencing
datasets.

Across all datasets, skiver produces accurate hazard rate estimates (Figure 5.A). For
Illumina reads, the estimated hazard rates are slightly biased upward, which is likely
attributable to known positional biases in Illumina sequencing, where errors are more
prevalent at the ends of the reads (Supplementary Figure S2). Notably, the estimated ĥ(t)
curves are nearly identical for skiver with and without a reference genome, indicating that
selecting the most frequently observed value as the consensus is a reliable strategy.

For survival rate estimation, skiver consistently outperforms the baseline methods, achieving
the lowest MSE_S and MSE_h across almost all datasets (Figure 5.B; Table 2; Supplementary
Table S3). Skiver also produces a reasonable estimation of error rates across multiple real
complex metagenomic samples that agrees with previous error profiles. Among the baselines,
KMC + GenomeScope 2.0 is not designed for metagenomic samples and therefore struggles to
identify an appropriate cutoff between erroneous and solid k-mers in the k-mer count
histogram. Its performance is nevertheless substantially better on bacterial isolates and
human sequencing data. Seqtk, on the other hand, may either under- or overestimate error
rates, reflecting the fact that base quality scores in many read sets are not well calibrated.

**Table 2** MSE_S for benchmarked tools on real sequencing datasets. Lowest in **bold**,
second lowest <u>underlined</u>.

| Dataset | seqtk | KMC + GenomeScope2.0 | skiver | skiver (with reference) |
|---|---|---|---|---|
| Zymo Log, Nanopore GridION | 6.38 × 10⁻² | 1.05 × 10⁻² | **1.47 × 10⁻⁴** | <u>1.78 × 10⁻⁴</u> |
| Zymo Gut Microbiome, PacBio HiFi | 4.36 × 10⁻³ | 7.95 × 10⁻² | <u>1.74 × 10⁻⁴</u> | **1.56 × 10⁻⁶** |
| *B. subtilis*, Illumina | 1.83 × 10⁻² | 1.24 × 10⁻² | <u>6.80 × 10⁻⁴</u> | **5.65 × 10⁻⁴** |
| *B. subtilis*, PacBio RSII | 8.53 × 10⁻⁴ | <u>3.61 × 10⁻⁴</u> | 4.71 × 10⁻⁴ | **6.16 × 10⁻⁵** |
| Human, Nanopore R10.4 | 7.11 × 10⁻² | 1.47 × 10⁻² | <u>1.63 × 10⁻³</u> | **6.33 × 10⁻⁴** |
| Human, Nanopore R9.4 | 3.99 × 10⁻² | 5.44 × 10⁻² | **1.58 × 10⁻³** | <u>3.06 × 10⁻³</u> |
| Human, PacBio HiFi | 2.58 × 10⁻³ | <u>1.70 × 10⁻⁴</u> | 1.78 × 10⁻⁴ | **5.51 × 10⁻⁶** |

Since skiver accurately estimates the survival rate S(k), i.e., the probability that a k-mer
is free of sequencing errors, its output can serve as prior knowledge for downstream k-mer
-based methods. For example, skiver's error rate estimation can be used as input to sylph
[40], a fast metagenomic profiling tool, to accurately estimate the percentage of unclassified
reads (Supplementary Table S4).

In addition to estimating error rates, skiver enables profiling of sequencing error spectra
(Figure 5.C). For nanopore datasets, substitutions are dominated by **G → A and A → G
transitions**; PacBio datasets are characterized primarily by insertions and deletions; and
Illumina datasets are dominated by substitution errors. These patterns are consistent with
previously reported error profiles for the respective sequencing platforms [19, 23, 41]. The
estimated spectra produced by skiver are also highly consistent regardless of whether a
reference genome is provided (Supplementary Figure S6), further supporting the idea that
the consensus derived from the most frequently observed value recovers the true reference base
with high probability. On the other hand, the error spectra estimated using BEST in
combination with Minimap2 differ noticeably from those obtained by skiver (Supplementary
Figure S7). In particular, substitution rates are frequently inflated, which may be due to
the alignment scoring scheme employed by Minimap2 and soft clipping.

Furthermore, skiver is able to calibrate Phred scores reported in the FASTQ files by
estimating the conditional hazard rates, reporting similar empirical error rates as Minimap2
(Figure 5.D). In Minimap2, a base and its associated Phred score are considered as an error
if it is a substitution, an insertion, or follows a deletion. In all cases, the empirical
error rate decreases as the Phred score increases. Interestingly, across multiple sequencing
platforms, **the empirical error rate for low Phred scores (below 10) tends to be lower than
the theoretical line, while high Phred scores are usually over-confident.** The ability to
calibrate Phred scores can potentially help downstream tools that rely on accurate Phred
scores, such as LoFreq [3].

Finally, we compared the memory usage and running times of all evaluated tools (Figure 6;
Supplementary Tables S8-9). Although skiver currently operates in a single-threaded mode, it
remains among the fastest and most memory-efficient methods evaluated. This efficiency
highlights the suitability of skiver as a lightweight and scalable first step in
bioinformatics pipelines, particularly for large-scale analysis of sequencing datasets.

> **Figure 6.** Wall clock time (left) and peak memory usage (right) on the Zymo Log,
> Nanopore GridION dataset. Time (s): Minimap2 + BEST 3591.82, seqtk 310.20,
> KMC + GenomeScope 934.24, skiver 322.45. Memory (GB): KMC + GenomeScope 29.81,
> Minimap2 + BEST 3.86, seqtk 0.20, skiver 2.25.

### 3.6 Ablation Studies

Lastly, we assessed whether modeling the time to first sequencing error, T, using a discrete
Weibull distribution is appropriate. To this end, we constructed a variant of skiver in which
the shape parameter is fixed at β̂ = 1. Under this constraint, the survival function simplifies
to Ŝ(t) = exp(−λ̂t), corresponding to a constant hazard rate ĥ(t) = 1 − exp(−λ̂).

Experiments conducted on the same datasets (Table 3; Supplementary Table S7; Supplementary
Figure S10) show that this constant-hazard model still outperforms all baseline methods in
terms of both MSE_S and MSE_h, and provides a reasonable single-number estimate of the
sequencing error rate. However, across all datasets, both MSE_S and MSE_h are consistently
similar to or worse than those obtained using the full skiver model. The performance gap is
particularly pronounced for Nanopore datasets, where sequencing errors are known to be
non-uniform and tend to cluster along reads. These results indicate that the discrete Weibull
distribution more closely reflects real sequencing error processes than a constant-hazard
assumption. The difference can also be seen from the per-base error rate estimations, where
the estimations of skiver with a constant hazard rate are much lower than the default skiver
and Minimap2 + BEST's estimations (Supplementary S2 and S6).

We further performed the same ablation experiments using skiver with the outlier filter
disabled. In this setting, the estimated survival rates are substantially lower
(Supplementary Figure S10), particularly for the Zymo Gut Microbiome mock community, where
multiple *E. coli* strains coexist, and for human datasets, which are diploid and contain
significant amounts of repetitive regions. The improved agreement between the default skiver
estimates and the ground truth in these challenging datasets further demonstrates that the
proposed outlier filter is robust to genomic heterogeneity, enabling accurate survival and
hazard rate estimation even in the presence of multiple strains, alleles, or repetitive
sequences.

**Table 3** MSE_S across different datasets for variants of the skiver algorithm.
"constant hazard rate": β̂ fixed to 1; "no filter": outlier filter disabled during hazard
rate estimation. Lowest in **bold**.

| Dataset | Default skiver | constant hazard rate | no filter |
|---|---|---|---|
| Zymo Log, Nanopore GridION | 1.47 × 10⁻⁴ | 2.99 × 10⁻³ | **2.51 × 10⁻⁵** |
| Zymo Gut Microbiome, PacBio HiFi | 1.74 × 10⁻⁴ | **1.70 × 10⁻⁴** | 1.81 × 10⁻³ |
| *B. subtilis*, Illumina | 6.80 × 10⁻⁴ | 6.42 × 10⁻⁴ | **5.59 × 10⁻⁴** |
| *B. subtilis*, PacBio RSII | **4.71 × 10⁻⁴** | 5.66 × 10⁻⁴ | 1.32 × 10⁻³ |
| Human, Nanopore R10.4 | **1.63 × 10⁻³** | 2.03 × 10⁻³ | 1.16 × 10⁻¹ |
| Human, Nanopore R9.4 | **1.58 × 10⁻³** | 2.88 × 10⁻³ | 4.73 × 10⁻² |
| Human, PacBio HiFi | 1.78 × 10⁻⁴ | **1.45 × 10⁻⁴** | 1.84 × 10⁻¹ |

---

## 4 Discussion and Conclusion

In this study, we present skiver, a reference-free and alignment-free framework for profiling
sequencing errors by using (k, v)-mer sketches. By modeling the time to first error as a
discrete Weibull distribution, skiver provides direct estimates of survival and hazard rates,
enabling a more expressive characterization of sequencing error processes than conventional
single-number error metrics.

Across both simulated and real datasets, skiver consistently produces accurate survival and
hazard rate estimates under a wide range of sequencing conditions. Skiver further enables
reliable profiling of sequencing error spectra. The inferred substitution, insertion, and
deletion patterns across Illumina, PacBio, and Nanopore datasets closely match previously
reported platform-specific error characteristics. The proposed outlier filter substantially
improves robustness in datasets containing multiple strains, alleles, or repetitive regions,
enabling accurate estimation of error rates of reads coming from complex microbial
communities and diploid genomes.

Finally, skiver is computationally efficient, achieving low runtime and memory usage, making
it suitable as a lightweight preprocessing step in sequencing analysis pipelines.

There remains room for further improvement in the skiver pipeline. In particular, performance
on datasets with high error rates or low sequencing coverage could potentially be enhanced
through adaptive selection of the parameters k and v. The error rate estimation results of
skiver can also support a future generation of k-mer-based algorithms for alignment-free
variant calling, read error correction, and other sequence analysis tasks.

Beyond applications to sequencing reads, the skiver framework could also be extended to
collections of genomes. In this setting, the hazard rate can be interpreted as the
probability that the next genomic position differs due to mutation. Under this formulation,
skiver could enable efficient estimation of mutational spectra across large genome
collections, providing a scalable approach for large-scale mutational spectra analysis.

---

## Declarations

**Funding.** Z.G. and P.S. are supported by NUS Research Scholarship.

**Data and code availability.** <https://github.com/GZHoffie/skiver>;
<https://github.com/GZHoffie/skiver-test>.

**Conflict of interest.** The authors declare no competing interests.

---

## Supplementary Material

### Table S1 — one-edit neighbor set of value `AACG`

`A → T` denotes a substitution of A to T, `− → A` an insertion of A, `A → −` a deletion of A.
**When a base is deleted, a new base is filled in at the end to maintain the length v**, so
there are at most 4v possible deletion neighbors. *Ambiguous* means the neighbor can be
reached by multiple edit operations (e.g. `AAGC` is reachable from `AACG` by deleting G or
inserting C).

| Neighbor | Edit | Neighbor | Edit |
|---|---|---|---|
| AAAC | − → A | ATCG | A → T |
| TACG | A → T | ACGG | A → − |
| AGCG | A → G | CAAC | − → C |
| ACAC | − → C | TAAC | − → T |
| AAGA | C → − | ACGT | A → − |
| ATAC | − → T | AACT | Ambiguous |
| AACA | Ambiguous | AGAC | − → G |
| AAGC | Ambiguous | ACGC | A → − |
| GAAC | − → G | AAAG | C → A |
| ACGA | A → − | AATC | − → T |
| CACG | A → C | AAGG | Ambiguous |
| ACCG | A → C | AAGT | C → − |
| GACG | A → G | AATG | C → T |
| AACC | Ambiguous | | |

(6 of 30 listed neighbors — 20% — are Ambiguous and therefore **discarded**.)

### Figure captions S1–S10

- **S1.** Least-square best fits of the observed survival rate, assuming T follows a Weibull
  distribution, S(t) = exp(−λt^β). β < 1 in most fits, indicating heterogeneous error rate
  among reads. Fitted values: Zymo Log Nanopore `S(t)=e^(−0.0886 t^0.8732)`; Zymo Gut HiFi
  `e^(−0.0016 t^0.9418)`; *B. subtilis* Illumina `e^(−0.0088 t^0.7170)`; *B. subtilis* RSII
  `e^(−0.1072 t^1.0109)`; Human R10.4 `e^(−0.0369 t^0.5952)`; Human R9.4
  `e^(−0.0512 t^0.7761)`; Human HiFi `e^(−0.0011 t^0.9903)`.
- **S2.** Observed survival and hazard rate of the *B. subtilis* Illumina dataset (read
  length 100). Because Illumina errors are more likely at read ends, the dataset shows a
  decreasing hazard for small t and an **increasing hazard for large t — a "bathtub" curve**.
- **S3.** Dependence of the proportion of each error type (insertion, deletion, substitution)
  on t, for Zymo Log Nanopore, Zymo Gut HiFi, and *B. subtilis* Illumina.
- **S4 / S5.** MSE_S and MSE_h of skiver on Badread simulated data across read identity and
  coverage.
- **S6.** Estimated error spectra with skiver on metagenomic data, with and without reference.
- **S7.** Estimated error type composition, Minimap2 + BEST vs skiver.
- **S8.** Estimated hazard rate, survival rate and error spectrum of HG002 across platforms.
- **S9.** MSE_S on three real datasets across combinations of k and v. MSE_S does not vary a
  lot as k and v change. For low-error reads a larger v is preferable; for high error a
  smaller v. **k = v = 17 is a good choice across all datasets.**
- **S10.** Ablation: predicted survival rate of default skiver, constant-hazard skiver, and
  skiver without the outlier filter.

### Table S2 — Estimated per-base error rates

Skiver values are ĥ(1) = 1 − exp(−λ̂). For Minimap2 + BEST the gap-compressed error rate is
reported.

| Dataset | Minimap2 + BEST | seqtk | KMC + GenomeScope2.0 | skiver | skiver (with ref) |
|---|---|---|---|---|---|
| Zymo Log, Nanopore GridION | 0.089030 | 0.0214 | 0.000373 | 0.0887 | 0.0743 |
| Zymo Gut Microbiome, PacBio HiFi | 0.001441 | 5.12 × 10⁻⁹ | 0.000079 | 0.0011 | 0.0015 |
| *B. subtilis*, Illumina | 0.007221 | 0.0005 | 0.000057 | 0.0040 | 0.0041 |
| *B. subtilis*, PacBio RSII | 0.103340 | 0.1259 | 0.001173 | 0.1503 | 0.1141 |
| Human, Nanopore R10.4 | 0.028530 | 0.0010 | 0.010760 | 0.0150 | 0.0172 |
| Human, Nanopore R9.4 | 0.047698 | 0.0110 | 0.051815 | 0.0517 | 0.0587 |
| Human, PacBio HiFi | 0.001175 | 0.0001 | 0.000814 | 0.0015 | 0.0021 |

### Table S3 — MSE_h across datasets

| Dataset | seqtk | KMC + GenomeScope2.0 | skiver | skiver (with ref) |
|---|---|---|---|---|
| Zymo Log, Nanopore GridION | 6.93 × 10⁻⁴ | 1.44 × 10⁻⁴ | 1.75 × 10⁻⁵ | **1.13 × 10⁻⁵** |
| Zymo Gut Microbiome, PacBio HiFi | 1.29 × 10⁻⁶ | 4.53 × 10⁻⁵ | 5.08 × 10⁻⁸ | **8.75 × 10⁻⁹** |
| *B. subtilis*, Illumina | 2.51 × 10⁻⁴ | **2.24 × 10⁻⁴** | 2.33 × 10⁻⁴ | 2.32 × 10⁻⁴ |
| *B. subtilis*, PacBio RSII | 1.30 × 10⁻³ | 7.67 × 10⁻⁴ | 7.04 × 10⁻⁵ | **4.62 × 10⁻⁵** |
| Human, Nanopore R10.4 | 3.79 × 10⁻⁵ | 3.41 × 10⁻⁵ | 4.87 × 10⁻⁶ | **3.85 × 10⁻⁶** |
| Human, Nanopore R9.4 | 6.87 × 10⁻⁵ | 1.23 × 10⁻³ | **1.36 × 10⁻⁵** | 2.00 × 10⁻⁵ |
| Human, PacBio HiFi | 7.96 × 10⁻⁷ | **5.38 × 10⁻⁸** | 8.64 × 10⁻⁸ | 3.60 × 10⁻⁸ |

### Table S4 — Percentage of unclassified reads

Ground truth = Minimap2 against the Zymo reference database.

| Dataset | Minimap2 (ground truth) | Default sylph | Skiver + sylph |
|---|---|---|---|
| Zymo Log, Nanopore GridION | 89.9126% | 36.3728% | 89.4377% |
| Zymo Gut Microbiome, PacBio HiFi | 96.9265% | 95.5556% | 96.9850% |

### Table S5 — Per-base error rate on real metagenomes (ins / del / sub)

**No ground truth is available for these three datasets.** The paper's caption states "Skiver
without outlier filter is biased heavily upwards."

| Dataset | Skiver | Skiver (no filter) |
|---|---|---|
| Human Gut, Nanopore PromethION | 7.27% (2.099% / 2.900% / 2.271%) | 9.628% (2.711% / 3.934% / 2.983%) |
| Soil sample, Illumina NovaSeq | 0.448% (0.001% / 0.001% / 0.446%) | 6.565% (0.055% / 0.064% / 6.446%) |
| Water sample, Illumina NovaSeq | 0.462% (0.001% / 0.002% / 0.459%) | 5.289% (0.036% / 0.030% / 5.223%) |

### Table S6 — Estimated error rates ĥ(1) for algorithm variants

| Dataset | Default skiver | constant hazard rate | no filter |
|---|---|---|---|
| Zymo Log, Nanopore GridION | 0.088681 | 0.045869 | 0.089319 |
| Zymo Gut Microbiome, PacBio HiFi | 0.001083 | 0.000948 | 0.004329 |
| *B. subtilis*, Illumina | 0.004018 | 0.002898 | 0.004162 |
| *B. subtilis*, PacBio RSII | 0.150313 | 0.090512 | 0.181092 |
| Human, Nanopore R10.4 | 0.014999 | 0.006597 | 0.078046 |
| Human, Nanopore R9.4 | 0.051685 | 0.023164 | 0.095097 |
| Human, PacBio HiFi | 0.001487 | 0.000833 | 0.046505 |

### Table S7 — MSE_h for algorithm variants

| Dataset | Default skiver | constant hazard rate | no filter |
|---|---|---|---|
| Zymo Log, Nanopore GridION | 1.75 × 10⁻⁵ | 5.80 × 10⁻⁵ | **9.97 × 10⁻⁶** |
| Zymo Gut Microbiome, PacBio HiFi | **5.08 × 10⁻⁸** | 5.09 × 10⁻⁸ | 6.66 × 10⁻⁷ |
| *B. subtilis*, Illumina | 2.33 × 10⁻⁴ | **2.32 × 10⁻⁴** | 2.32 × 10⁻⁴ |
| *B. subtilis*, PacBio RSII | **7.04 × 10⁻⁵** | 7.02 × 10⁻⁵ | 1.61 × 10⁻⁴ |
| Human, Nanopore R10.4 | **9.92 × 10⁻⁶** | 2.07 × 10⁻⁵ | 2.86 × 10⁻⁴ |
| Human, Nanopore R9.4 | **1.34 × 10⁻⁵** | 3.63 × 10⁻⁵ | 5.74 × 10⁻⁴ |
| Human, PacBio HiFi | 8.96 × 10⁻⁸ | **7.12 × 10⁻⁸** | 1.97 × 10⁻⁴ |

### Tables S8 / S9 — Running time (s) and memory (GB)

| Dataset | KMC+GS time | seqtk time | MM2+BEST time | Skiver time | KMC+GS mem | seqtk mem | MM2+BEST mem | Skiver mem |
|---|---|---|---|---|---|---|---|---|
| Zymo Log, Nanopore GridION | 934.24 | 310.20 | 3591.82 | 322.45 | 29.81 | 0.20 | 3.86 | 2.25 |
| Zymo Gut Microbiome, PacBio HiFi | 1164.23 | 221.79 | 4753.00 | 349.93 | 29.75 | 0.05 | 3.15 | 5.51 |
| *B. subtilis*, Illumina | 304.16 | 66.82 | 1078.41 | 93.11 | 23.29 | 0.00 | 0.87 | 0.65 |
| *B. subtilis*, PacBio RSII | 73.11 | 12.12 | 441.58 | 99.99 | 11.91 | 0.10 | 3.62 | 4.64 |
| Human, Nanopore R10.4 | 4631.00 | 1823.54 | 29867.00 | 1078.43 | 29.89 | 1.55 | 22.02 | 3.27 |
| Human, Nanopore R9.4 | 3858.00 | 1679.65 | 5942.00 | 1103.88 | 29.96 | 3.10 | 22.02 | 3.25 |
| Human, PacBio HiFi | 6331.00 | 1627.24 | 20533.00 | 1253.03 | 29.93 | 0.05 | 22.02 | 5.86 |

### Supplementary Note S1 — coverage needed for the consensus to be correct

An **informal** lower bound for the coverage of keys needed so that the value with the highest
count (the consensus) is indeed the true value in the sequenced genome with high probability.

Assumptions:

1. Each base in the reads has an **independent and equal** chance of being erroneous with
   probability ε.
2. ε is small such that the value with the highest count contains **at most one error almost
   surely**.

For a given key K, let V\* be the true value after the key in the underlying reference genome.
Suppose K appears N times, with observed values V₁, …, V_N. Then

$$p_0 := \Pr[V_j = V^*] = (1-\varepsilon)^v.$$

Fix a distance-1 neighbor w of V\*; the probability of observing it is

$$q := \Pr[V_j = w] \le (1-\varepsilon)^{v-1}\varepsilon.$$

Let C\* and C_w be the counts of V\* and w. Let D_j := I(V_j = V\*) − I(V_j = w), i.i.d. taking
values in {−1, 1}, and S_N := Σ D_j, so E[D_j] = p₀ − q and E[S_N] = N(p₀ − q). By Hoeffding's
inequality,

$$\Pr[C^* \le C_w] = \Pr[S_N \le 0] \le \exp\left(-\frac{2N^2(p_0-q)^2}{4N}\right) = \exp\left(-\frac{N(p_0-q)^2}{2}\right).$$

There are at most 11v distance-1 neighbors of V\*. If the multiplicity of the key is
sufficiently large, i.e.

$$N \ge \frac{2\log\frac{11v}{0.99}}{(p_0-q)^2} = \frac{2\log\frac{10v}{9}}{\bigl((1-\varepsilon)^{v-1}(1-2\varepsilon)\bigr)^2}, \tag{1}$$

then by union bound

$$\Pr[C^* > C_w \text{ for all } w] \ge 11v \cdot \Pr[C^* > C_w] = 11v\cdot\exp\left(-\frac{N(p_0-q)^2}{2}\right) \ge 11v\cdot\frac{0.99}{11v} = 0.99.$$

This means that if the multiplicity of keys N = Ω((1 − ε)^−2v log v) then C\* > C_w for all w
with high probability. The higher the coverage N, the higher the probability that C\* has the
highest count.

> *(Reproduced as printed. The union-bound step is stated incorrectly in the preprint — see
> [skiver_paper_review.md](skiver_paper_review.md) §C1.)*

### Supplementary Note S2 — the complementary log-log approximation

$$h(t) = 1 - \exp\bigl(-\lambda(t^\beta - (t-1)^\beta)\bigr)$$
$$\Rightarrow \log(-\log(1-h(t))) = \log\lambda + \log\bigl(t^\beta - (t-1)^\beta\bigr)$$

By Taylor expansion,

$$(t-1)^\beta = t^\beta\left(1-\frac{1}{t}\right)^\beta = t^\beta\left(1 - \frac{\beta}{t} + \frac{\beta(\beta-1)}{2t^2} + O(t^{-3})\right),$$

so

$$t^\beta - (t-1)^\beta = \beta t^{\beta-1}\left(1 - \frac{\beta-1}{2t} + O(t^{-2})\right),$$

and therefore

$$\log(-\log(1-h(t))) = \log\lambda\beta + (\beta-1)\log t + \log\left(1 - \frac{\beta-1}{2t} + O(t^{-2})\right) \approx \log\lambda\beta + (\beta-1)\log t.$$

This approximation is based on the observation that β ≈ 1 in the real dataset
(Supplementary Figure S1) and t ∈ [k + 1, k + v] is large in the estimated hazard rates
h(t). **It should be noted that this approximation is biased, but causes negligible error in
the estimation of λ and β.**

---

## References

1. Chen, S., Zhou, Y., Chen, Y., Gu, J.: fastp: an ultra-fast all-in-one fastq preprocessor. *Bioinformatics* **34**(17), 884–890 (2018)
2. Li, H.: Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* **34**(18), 3094–3100 (2018)
3. Wilm, A., Aw, P.P.K., Bertrand, D., Yeo, G.H.T., Ong, S.H., Wong, C.H., Khor, C.C., Petric, R., Hibberd, M.L., Nagarajan, N.: LoFreq: a sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. *Nucleic Acids Research* **40**(22), 11189–11201 (2012)
4. Shaw, J., Gounot, J.-S., Chen, H., Nagarajan, N., Yu, Y.W.: Floria: fast and accurate strain haplotyping in metagenomes. *Bioinformatics* **40**(Suppl 1), 30–38 (2024)
5. Liu, D., Belyaeva, A., Shafin, K., Chang, P.-C., Carroll, A., Cook, D.E.: BEST: a tool for characterizing sequencing errors. *bioRxiv* 2022–12 (2022)
6. Ross, M.G., Russ, C., Costello, M., Hollinger, A., Lennon, N.J., Hegarty, R., Nusbaum, C., Jaffe, D.B.: Characterizing and measuring bias in sequence data. *Genome Biology* **14**(5), 51 (2013)
7. Liu, Y., Li, Y., Chen, E., Xu, J., Zhang, W., Zeng, X., Luo, X.: Repeat and haplotype aware error correction in nanopore sequencing reads with dechat. *Communications Biology* **7**(1), 1678 (2024)
8. Victoria Wang, X., Blades, N., Ding, J., Sultana, R., Parmigiani, G.: Estimation of sequencing error rates in short reads. *BMC Bioinformatics* **13**(1), 185 (2012)
9. Davis, E.M., Sun, Y., Liu, Y., Kolekar, P., Shao, Y., Szlachta, K., Mulder, H.L., Ren, D., Rice, S.V., Wang, Z., et al.: SequencErr: measuring and suppressing sequencer errors in next-generation sequencing data. *Genome Biology* **22**(1), 37 (2021)
10. Bethune, J., Kleppe, A., Besenbacher, S.: A method to build extended sequence context models of point mutations and indels. *Nature Communications* **13**(1), 7884 (2022)
11. Blanca, A., Harris, R.S., Koslicki, D., Medvedev, P.: The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches. *Journal of Computational Biology* **29**(2), 155–168 (2022)
12. Greenfield, P., Duesing, K., Papanicolaou, A., Bauer, D.C.: Blue: correcting sequencing errors using consensus and context. *Bioinformatics* **30**(19), 2723–2732 (2014)
13. Hera, M.R., Medvedev, P., Koslicki, D., Blanca, A.: Estimation of substitution and indel rates via k-mer statistics. *bioRxiv* 2025–05 (2025)
14. Kaplinski, L., Möls, M., Puurand, T., Remm, M.: Docest — fast and accurate estimator of human NGS sequencing depth and error rate. *Bioinformatics Advances* **3**(1), 084 (2023)
15. Melsted, P., Halldórsson, B.V.: KmerStream: streaming algorithms for k-mer abundance estimation. *Bioinformatics* **30**(24), 3541–3547 (2014)
16. Ranallo-Benavidez, T.R., Jaron, K.S., Schatz, M.C.: GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nature Communications* **11**(1), 1432 (2020)
17. Wu, H., Blanca, A., Medvedev, P.: A k-mer-based estimator of the substitution rate between repetitive sequences. *bioRxiv* 2025–06 (2025)
18. Yu, Y.W., Yorukoglu, D., Peng, J., Berger, B.: Quality score compression improves genotyping accuracy. *Nature Biotechnology* **33**(3), 240–243 (2015)
19. Delahaye, C., Nicolas, J.: Sequencing DNA with nanopores: troubles and biases. *PLoS ONE* **16**(10), 0257521 (2021)
20. Schirmer, M., Ijaz, U.Z., D'Amore, R., Hall, N., Sloan, W.T., Quince, C.: Insight into biases and sequencing errors for amplicon sequencing with the Illumina MiSeq platform. *Nucleic Acids Research* **43**(6), 37–37 (2015)
21. Liao, P., Satten, G.A., Hu, Y.-J.: PhredEM: a Phred-score-informed genotype-calling approach for next-generation sequencing studies. *Genetic Epidemiology* **41**(5), 375–387 (2017)
22. Sahlin, K., Medvedev, P.: Error correction enables use of Oxford Nanopore technology for reference-free transcriptome analysis. *Nature Communications* **12**(1), 2 (2021)
23. Dohm, J.C., Peters, P., Stralis-Pavese, N., Himmelbauer, H.: Benchmarking of long-read correction methods. *NAR Genomics and Bioinformatics* **2**(2), 037 (2020)
24. Chaung, K., Baharav, T.Z., Henderson, G., Zheludev, I.N., Wang, P.L., Salzman, J.: SPLASH: a statistical, reference-free genomic algorithm unifies biological discovery. *Cell* **186**(25), 5440–5456 (2023)
25. Kokot, M., Dehghannasiri, R., Baharav, T., Salzman, J., Deorowicz, S.: Scalable and unsupervised discovery from raw sequencing reads using SPLASH2. *Nature Biotechnology* **43**(7), 1084–1090 (2025)
26. Irber, L., Brooks, P.T., Reiter, T., Pierce-Ward, N.T., Hera, M.R., Koslicki, D., Brown, C.T.: Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers. *bioRxiv* 2022–01 (2022)
27. Vaupel, J.W., Manton, K.G., Stallard, E.: The impact of heterogeneity in individual frailty on the dynamics of mortality. *Demography* **16**(3), 439–454 (1979)
28. Allison, P.D.: Discrete-time methods for the analysis of event histories. *Sociological Methodology* **13**, 61–98 (1982)
29. Nicholls, S.M., Quick, J.C., Tang, S., Loman, N.J.: Ultra-deep, long-read nanopore sequencing of mock microbial community standards. *GigaScience* **8**(5), 043 (2019)
30. McIntyre, A.B., Alexander, N., Grigorev, K., Bezdan, D., Sichtig, H., Chiu, C.Y., Mason, C.E.: Single-molecule sequencing detection of N6-methyladenine in microbial reference materials. *Nature Communications* **10**(1), 579 (2019)
31. Shafin, K., Pesout, T., Lorig-Roach, R., Haukness, M., Olsen, H.E., Bosworth, C., Armstrong, J., Tigyi, K., Maurer, N., Koren, S., et al.: Efficient de novo assembly of eleven human genomes using PromethION sequencing and a novel nanopore toolkit. *bioRxiv* 715722 (2019)
32. Gounot, J.-S., Chia, M., Bertrand, D., Saw, W.-Y., Ravikrishnan, A., Low, A., Ding, Y., Ng, A.H.Q., Tan, L.W.L., Teo, Y.-Y., et al.: Genome-centric analysis of short and long read metagenomes reveals uncharacterized microbiome diversity in Southeast Asians. *Nature Communications* **13**(1), 6044 (2022)
33. Kokot, M., Dlugosz, M., Deorowicz, S.: KMC 3: counting and manipulating k-mer statistics. *Bioinformatics* **33**(17), 2759–2761 (2017)
34. Li, H.: seqtk: toolkit for processing sequences in FASTA/Q formats. <https://github.com/lh3/seqtk> (2025)
35. Chen, Q., Zhang, B., Peng, C., Huang, J., Liu, Z., Shen, X., Jiang, C.: Kun-peng enables scalable and accurate pan-domain metagenomic classification. *Briefings in Bioinformatics* **27**(2), 119 (2026)
36. Suresh, K., Severn, C., Ghosh, D.: Survival prediction models: an introduction to discrete-time modeling. *BMC Medical Research Methodology* **22**(1), 207 (2022)
37. Wick, R.R.: Badread: simulation of error-prone long reads. *Journal of Open Source Software* **4**(36), 1316 (2019)
38. Parks, D.H., Chuvochina, M., Rinke, C., Mussig, A.J., Chaumeil, P.-A., Hugenholtz, P.: GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Research* **50**(D1), 785–794 (2022)
39. Shaw, J., Yu, Y.W.: Fast and robust metagenomic sequence comparison through sparse chaining with skani. *Nature Methods* **20**(11), 1661–1665 (2023)
40. Shaw, J., Yu, Y.W.: Rapid species-level metagenome profiling and containment estimation with sylph. *Nature Biotechnology* **43**(8), 1348–1359 (2025)
41. Minoche, A.E., Dohm, J.C., Himmelbauer, H.: Evaluation of genomic high-throughput sequencing data generated on Illumina HiSeq and Genome Analyzer systems. *Genome Biology* **12**(11), 112 (2011)
