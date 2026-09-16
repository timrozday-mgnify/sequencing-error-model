# Implementation plan: sequencing-error-model

Status: **draft, revised 2026-09-16**. Phases 0 (repo, CI and PR policy) and 1 (inputs) are done; phase 2 is in progress (`spec.py`, the head Q and head E fitters and per-head selection landed; the insertion-quality sub-head waits for `reference` tuples); phase 3 (generator, paired output and recovery harness) is done; phase 4 (`pe-overlap`) has met its exit criteria. The source, EM fit, recovery test, real-data placement guards and a real-run spec have landed, along with Q smoothing of head E and a finer default head Q (user guide: [pe_overlap.md](pe_overlap.md)). Its open items and the non-blocking ErrorProfiler comparison remain. Phase 5 (`reference`) is in progress: the core of `sources/bam.py`, site masks, the contig filter, the per-contig report, the `pe-overlap` vs `reference` comparison, the aligner bias check against an observable truth with a soft-clipping correction, and the real near-clonal run (SRR24523812: `pe-overlap` vs `reference` disagreements explained, `QualityWindow` beats centre Q) and the baseline harness with its Illumina (vs ReSeq) and ONT (vs Badread, PBSIM3, CycSim) runs have landed. The main open finding is per-read quality and error heterogeneity (`Latent(S)`). Everything else is planned. This revision adds the `pe-overlap` and `reference` evidence modes beside the `kmer` (skiver) mode, and builds them first so they can check the `kmer` evidence (§1.1, §9). The latest revision adds separating biological variation (strains, minor alleles, divergent repeats) from sequencing error (§6.6): every method is first shown unbiased on clonal simulations, then a variation simulator (phase 7) measures the problem, and separation methods (phase 8) are developed on it before real metagenomes. Later phases are renumbered 9–13.

## 1. Goal

A Python package and CLI that:

1. **learns** a sequencing error model from any of three independent kinds of evidence, the **evidence modes** (§1.1): paired-end overlap, alignment to a reference or self-assembly, and k-mer consensus (skiver). Two of the three need no external reference. Every mode trains the same kind of model (§5), so models from different modes can be compared directly;
2. **models bases and qualities jointly** (§5). The core quantity is the error profile of a base *given its own quality score, the surrounding bases and qualities*, and other context such as read position, mate, strand and read-level state;
3. **generates** reads with that model, outputting **both bases and quality scores**, and keeping the `FASTA in → FASTQ + CIGAR out` contract that genome-blender already uses; and
4. **exports** the model to the native formats of popular, maintained read simulators, so it plugs into existing pipelines (including the metagenome wrappers CAMISIM and MeSS).

**Quality scores are never treated as evidence of error.** A reported Q is a feature that may covary with the true error rate, and an output the generator must reproduce. Error labels only come from truth-bearing sources: mate overlap, alignment to a reference or self-assembly, skiver consensus, duplex reads, and so on.

Out of scope: signal-level simulation (squigulator, seq2squiggle), variant calling and mutation models as products (a minimal variation simulator and variant-site posteriors exist only to separate variation from error, §6.6), abundance/community modelling (CAMISIM and MeSS already do this), and error *correction*. Quality recalibration is produced as a diagnostic, not as a read-rewriting tool.

### 1.1 Evidence modes

| Mode | Truth for error labels | Input | External reference? | Platforms |
|---|---|---|---|---|
| **`pe-overlap`** | Agreement of R1 and R2 where they overlap | paired FASTQ | No | Illumina (inserts short enough to overlap) |
| **`reference`** | Alignment to a reference: an external genome, a spike-in or mock community (PhiX, lambda, ZymoBIOMICS), or a polished self-assembly of the sample | FASTQ + reference or contigs + BAM/CRAM | Yes, but the sample's own assembly counts | all |
| **`kmer`** | skiver (k,v)-mer consensus | skiver outputs + raw FASTQ | No | all |

The `kmer` mode has two skiver builds:

| Build | Input | skiver build |
|---|---|---|
| **default** | `skiver analyze` CSVs (+ raw FASTQ for the quality model) | unmodified upstream release, pinned (currently **v0.3.2**) |
| **enhanced** | `skiver dump` per-observation outputs (plus the default inputs) | modified skiver: a new, minimal fork of GZHoffie/skiver, seeded from the `timrozday-mgnify/skiver` fork (§12, decision 5) |

"Default mode" and "enhanced mode", here and in the repo, mean the `kmer` mode on these two builds.

**One model, three kinds of evidence.** All modes reduce to the same observation schema (§7) and fit the same two heads with the same component vocabulary (§5.4). They differ only in which components they identify (§5.6, §8). That makes cross-mode comparison a validation tool:

- **`pe-overlap` and `reference` are built first** (phases 4–5). They observe the exact error position together with per-observation quality windows, so fitting from them is direct: no latent edit position and no marginal matching.
- **They then check the `kmer` mode** (phase 6), at two levels:
  - *evidence*: on the same reads, the marginal tables skiver reports (spectrum, P(error | Q), position curve, GC, clustering) are tabulated from overlap and alignment observations and compared with skiver's CSVs;
  - *model*: the fitted specs are compared per component and on held-out likelihood of each other's observations.
- **Their truth assumptions fail differently**, so agreement between two modes is stronger evidence than either alone:
  - `pe-overlap` cancels errors both mates share (PCR, library and cluster errors) and sees only substitutions, only inside the overlap;
  - `reference` includes PCR errors and indels, but inherits reference errors, strain divergence and aligner scoring bias;
  - `kmer` includes PCR errors, needs consensus coverage and drops errors its outlier filter hides.

  Comparisons are restricted to the support two modes share (genomes, op classes, read positions). Gaps expected from these differences are reported as estimates, not failures. For example, the `reference` minus `pe-overlap` substitution rate estimates the PCR/library error contribution.

Reference-free training stays the project's niche (§4.4). The `reference` mode is a mode in its own right and the validation anchor for the other two.

---

## 2. What skiver gives us

### 2.1 Unmodified skiver (v0.3.2): the default-mode data contract

`skiver analyze <reads> -o prefix` always writes these files (from `OUTPUT_SUFFIXES` in upstream `src/analyze.rs`):

| File | Content | Model information it carries |
|---|---|---|
| `summary_error_rate.csv` | λ, β of the discrete Weibull, per-base and effective error rate, coverage, bootstrap CIs | Global rate; β < 1 indicates error clustering |
| `hazard_rate.csv`, `survival_rate.csv` | h(t) for t ∈ 1..v with CIs; fitted S(t) | Clustering shape |
| `summary_error_spectrum.csv` | counts per (20 ops × prev_base × next_base), `total` and `forward` | Marginal spectrum in trinucleotide context; strand asymmetry |
| `summary_error_spectrum_dependence_on_t.csv` | the same, per window position t | Check of the "composition constant in t" assumption |
| `summary_phred.csv` | correct/error counts per reported Q | Marginal P(error \| Q at the base), **centre Q only**, pooled over contexts |
| `summary_read_position.csv` | correct/error counts per index from read start and from read end | Marginal error-rate curve along the read |
| `summary_gc_content.csv` | error rate per GC bin | GC dependence |
| **`kvmer.csv`** | per key: `key`, `consensus_value`, `passes_filter`, `homopolymer_length`, `consensus_count`, `neighbor_count`, `total_count`, **20 per-op counts**, `consensus_count_up_to_v1..v` | **Per-locus (k+v)-base sequence context with per-op error counts.** The edit *position* inside the value is not recorded; no qualities |

**`kvmer.csv` is the key to default-mode feasibility for the base-context part.** It gives, for every sampled locus, the full k+v consensus sequence plus counts of each single-edit error type. For each key:

- the survival columns fix how many observations survive to each t;
- the op counts are a *mixture over the unknown position t* of the edit.

A context model P(op | surrounding true bases) can be fit by marginalising that latent position. The likelihood of op *e* at a key is Σₜ S(t−1) · h_e(ctx_t) · S(v−t), fit by EM or direct gradient on the mixture.

**Unmodified skiver carries no joint information about qualities.** It never pairs neighbouring qualities with errors, and never pairs Q with context. How default mode still produces a quality-aware model, and what it can't identify, is in §5.4 and §6.2.

### 2.2 Modified skiver: what enhanced mode adds

The current fork adds `skiver dump` (`src/dump.rs`) with `--raw`, `--base`, `--survival` and `--windows`. `ValueInfo` also gains `read_id` and `fragment_id`, so R1/R2 share an id.

| Output | Adds |
|---|---|
| `base_observations.tsv` | Per-base rows with the known edit position, Phred, `read_pos`, `dist_to_end`, strand, `fragment_id`. Joint op × context × **per-base Q** × position × strand (~32 GB for hq-illumina) |
| `raw_observations.tsv` | Key string and **`qual_str` for the whole v-window** per observation, i.e. the neighbouring-quality window within the value |
| `survival_observations.tsv` | Per-observation first-error time |
| `windows.bin` | Read-ordered 13-byte records (~430 MB, 75× smaller than the TSV). Enables read-level latent-state HMM fitting, but carries **no qualities** |

Candidate *new* outputs for the new minimal fork (§9, phase 12):

- **key-side qualities**, so the quality window extends left of t=1;
- a per-window quality summary in `windows.bin`, so a read-level latent state can drive both Q and errors;
- homopolymer run length at the edit, and indel length (> 1 bp);
- multi-edit (2+) values aligned to consensus;
- an explicit R1/R2 flag.

---

## 3. Lessons from developing the skiver fork

These come from `skiver/docs/hmm_error_model.md`, `docs/performance.md`, `examples/*synthetic_recovery*`, commit history and the report-integration work.

1. **The value length v controls recoverability.** Amplicon recovery test (source rate ≈ 3.4e-4, BAM truth 3.43e-4):
   - `v=13`: retrained 3.48e-4 ✓
   - `v=6`: 2.76e-4 (≈ −20%)
   - `v=1`: **0** (no errors observable)

   Keep v ≈ 13 and fit context from key bases when v < L. Reject inputs with tiny v.
2. **Mask impossible categories.** Destination-base substitutions (`sub_to_A` when the true base is A) must be masked in the likelihood, or fitted mass leaks into categories the generator discards and the marginal rate drops. Bake this into the core likelihood.
3. **Additive context scales and full tables don't.** `BaseContext(L)` is 4ᴸ rows, which is memory-bound past L≈6–8. AIC kept selecting longer context up to the L=8 screen cap. The same explosion applies to quality windows (Q-bins^m). Default to additive/low-rank features and materialise dense tables only at export time, where Badread uses 7-mers.
4. **Not every fitted component is generative.** In the fork, `PhredContext` and `Position` conditioned errors on quality and position but no model generated the qualities, so those components were training-only; `FragmentOverdispersion` only affects the training likelihood. **This project closes that gap:** the quality process is modelled explicitly (§5), so quality-conditioned error components become generative. Components that are still training-only are flagged in the model spec.
5. **Quality was two disconnected models.** Generation sampled Q from an empirical P(Q | error type) table, while P(error | Q) (qcal) was fitted separately and used only for validation. Nothing forced the two to agree, or tied Q to neighbouring qualities, context or position. A single joint model (§5) replaces both.
6. **`PhredContext` averaged Phred per context row.** Aggregating Q into per-context means discarded the per-observation pairing of quality and error. Quality-window features must stay per observation (§7 "Sparse tuples first").
7. **Synthetic recovery needs `--use-all -l 0`** (the outlier filter hides generated errors) and **reference-guided dump** to isolate model error from consensus error. Low-error Illumina needs large synthetic sets before metrics mean anything.
8. **Data size dominates engineering.** The 32 GB TSV needed an HDF5 cache and a single-aggregation model search. Plan for streaming aggregation from day one.
9. **Reports.** Quarto `-P` parameters only reach a `#| tags: [parameters]` cell. Embedding plotly.js per figure plus MathJax OOM'd pandoc at ~6 GB. Emit plotly once and skip MathJax.
10. **Determinism.** Consensus ties used to follow HashMap order; now broken by smallest value. Golden tests on exact outputs were what caught this.
11. **`is_forward` is strand, not R1/R2.** Paired-end R1/R2 differences need the fork's `fragment_id`/`read_id` or a new flag.
12. **Code structure.** The fork's Python grew as flat `scripts/` with a 2.5 k-line library, pickled torch artifacts and hard-coded platform reports. This repo starts as a proper package with a safe, versioned artifact format.

---

## 4. Landscape

### 4.1 Error *estimation* tools (context)

| Tool | Needs reference/alignment | Output |
|---|---|---|
| **skiver** | No (optionally reference-guided) | rate, hazard, spectrum, marginal Q calibration, position, GC |
| samtools stats / Qualimap / alignment-based (e.g. BEST) | Yes | mismatch/indel rates by cycle and quality. Biased by incomplete references and strain divergence (skiver paper) |
| GATK BQSR | Yes, plus known sites | Joint read group × Q × cycle × dinucleotide context mismatch table: the closest existing analogue of the §5 error head, but reference-bound |
| k-mer spectrum tools (GenomeScope-style) | No | a single error rate, no spectrum |
| Simulator "profilers" (`iss model`, `badread error_model`/`qscore_model`, NanoSim `read_analysis.py`, `art_profile_builder`, `neat model-seq-err`, ReSeq `stats`, CycSim training, LongISLND, GemSIM) | Yes, BAM/PAF against a reference (ART and NEAT can use FASTQ for qualities) | simulator-specific models |
| BEAR (DRISEE on artifactual duplicate reads) | No | per-position substitution rates, pooled indel rate, mean Q of erroneous bases vs position; feeds its own simulator |
| ErrorProfiler (PE overlap), DADA2 `learnErrors` | No | substitution × Q tables; no simulator |

**Our niche:** simulator-ready, quality-aware error models for samples where no good reference exists (metagenomes, novel isolates). We produce what the profilers would, without the alignment step, and in a representation richer than any single simulator's format. §4.4 states how narrow that niche is and what would make the project redundant.

### 4.2 Read simulators: export targets

Maintenance data is from the GitHub API on 2026-09-15. "Q coupling" describes how each simulator ties errors to qualities, which determines how much of the §5 joint model survives export.

| Simulator | Reads | Last push / ★ | Error model format | Q coupling | Where it's used | Export fidelity from our model |
|---|---|---|---|---|---|---|
| **ART** / **art_modern** | Illumina (+ MGI, AVITI, Onso profiles) | art_modern 2026-07 / 39 (original ART unmaintained, still ubiquitous) | Text quality profiles (per-position Q distributions, optionally per base); indel rates as CLI flags (`-ir/-dr/-ir2/-dr2`) | Q sampled per position, then substitution probability = 10^(−Q/10), i.e. assumes perfectly calibrated Q; no neighbour or context effects | MeSS, CAMISIM, countless pipelines | Medium. Per-position Q marginals map well. Our calibration can't be expressed without distorting either the Q distribution or the error rate (§6.3). Context and Q-window effects are lost |
| **InSilicoSeq** | Illumina (MiSeq/HiSeq/NextSeq/NovaSeq) | 2026-09 / 227 | `.npz` with `read_length`, `insert_size`, `mean_count_{forward,reverse}`, `quality_hist_{forward,reverse}`, `subst_choices_*` (per position × base), `ins_*`, `del_*` | Q-driven substitutions (to verify in phase 9) plus per-position substitution choices and indel rates | Metagenome benchmarking | Medium. Per-position Q histograms and substitution/indel tables are populated from marginals of both heads; same calibration caveat as ART |
| **NEAT v4** | Illumina | 2026-08 / 72 | gzip-pickle dict `{error_model1, error_model2, qual_score_model1, qual_score_model2}` (`SequencingErrorModel` + quality Markov model per mate) | Errors from Q-derived probabilities; Q from a Markov model | Variant-calling benchmarks | Medium. The Q head maps onto NEAT's quality Markov model (order 1, per mate). Needs NEAT's classes to pickle |
| **Mason2** (SeqAn) | Illumina/454/Sanger | seqan 2026-08 / 502 | CLI parameters only (mismatch/indel probabilities, begin/end ramps, quality mean/sd for correct and wrong bases) | Separate Q mean/sd for correct vs erroneous bases | Aligner benchmarks | Low (parametric), but trivial to emit. The correct/wrong Q split comes from the joint model |
| **wgsim** | Illumina | 2021 / 288 (unmaintained) | CLI base error rate | none | CAMISIM | Low. Emit a rate only, for CAMISIM compatibility |
| **Badread** | ONT/PacBio | 2026-07 / 301 | Error model: text, one line per 7-mer `KMER,p;ALT,p;...` (≤ 25 alternatives); Q model: `CIGAR;count;q:p,...` keyed by local CIGAR window | **Q conditioned on the local error pattern** (CIGAR window around the base) | Long-read tool benchmarks | **High.** 7-mer → alternatives materialises the base-context part; P(Q \| local CIGAR window) is a marginal of the joint model. Neighbour-Q and position effects are lost |
| **PBSIM3** | PacBio/ONT | 2025-04 / 117 | ERRHMM text (`IP`/`EP`(match,sub,ins,del)/`TP` per accuracy level) or QSHMM (quality HMM; errors from Q) | ERRHMM: none (qualities `!`). QSHMM: latent-state Q process, errors from Q | MeSS; top of the 2026 ONT simulator benchmark for length/Q realism | Medium. Export **QSHMM** for Q-bearing output (maps our latent state + Q head) and ERRHMM for error-only use. Sequence context is lost |
| **NanoSim** | ONT | 2026-03 / 311 | Directory of pickled KDEs + text Markov/histogram files (`_error_markov_model`, `_match_markov_model`, `*_hist`, `_error_rate.tsv`, quality models) from alignments | Quality models conditioned on match/error state | CAMISIM (`nanosim3`) | Low–medium. Only the error/quality Markov parts are estimable; length KDEs must come from a base model. Lowest priority |
| **ReSeq** | Illumina | 2021 paper / GitHub (maintenance to check) | Binary stats file from `reseq illuminaPE` statistics step (format to confirm in phase 10) | **Errors conditioned on Q**, position, errors so far in the read, reference base and recent dominant error; Q conditioned on previous Q, position, mate, tile, sequence quality and reference base; per-site systematic errors; stored as 2-D margins | Illumina tool benchmarks | Medium–high if the format is writable: head Q and centre-Q head E map closely. Preceding-sequence context is not modelled by ReSeq and is lost |
| **CycSim** | ONT/HiFi/Cyclone | 2025 preprint, GigaScience 2026 | k-mer sliding-window error model + error-state transition matrix, learned from BAM (format to confirm) | Q assigned after errors are placed | Long-read mapper parameter tuning | Medium: base context and error clustering map; Q process thin. Candidate, not committed |
| **genome-blender** | all | internal | `skiver-generate` subprocess contract | full joint model | This project's current consumer | Full: it uses our native generator |

Metagenome wrappers: **MeSS** wraps ART + PBSIM3, and **CAMISIM** wraps ART, NanoSim3 and wgsim. Covering **ART, InSilicoSeq, Badread and PBSIM3** therefore reaches most short-read, long-read and metagenome-simulation users. NEAT, Mason2, wgsim and NanoSim are a second wave.

The 2026 ONT simulator benchmark (Badread, LongISLND, lrsim, NanoSim, PBSIM3, SimLoRD) found that no simulator reproduces all metrics:

- **PBSIM3** was best for length, read-level Q and per-base Q calibration.
- **Badread and LongISLND** were the only ones to capture context-dependent substitution rates, which span ~2 orders of magnitude.
- **Homopolymer errors** were poorly reproduced by most simulators.

CycSim (not in that benchmark) reports better context-dependent substitution and simple-repeat error heterogeneity than Badread, NanoSim and PBSIM3, but its qualities are assigned after errors.

No existing long-read simulator combines context-dependent errors with realistic, calibrated qualities. For Illumina, ReSeq comes closest (errors conditioned on Q, a rich Q process) but has no preceding-sequence context model and is BAM-bound. That combination is the gap the native generator fills. Exports are necessarily projections of it.

### 4.3 Evidence sources beyond skiver

Skiver's gaps (§6.2) are mostly about *which axes are observed together*. It sees base context and op, but only centre-Q marginals, no quality windows, only marginal position, no R1/R2, no indel length, and no read-level structure. Other tools and data observe different axes with different truth assumptions. Many of them need no external reference.

The design consequence: **every source reduces to the same observation schema**, a (true-context window, observed bases, quality window, covariates, op) tuple with counts (§7). Each source declares:

- which fields it observes;
- **what it uses as truth for the op label**. Quality is never a truth label.

Sources can be fitted alone or jointly, and disagreement between them is a first-class diagnostic.

The first two rows of the table below are now evidence modes of their own (`pe-overlap` and `reference`, §1.1); spike-in and control references are handled by the `reference` mode.

| Source | Truth used for errors | Reference needed? | Adds what skiver lacks | Limitations | Maintained |
|---|---|---|---|---|---|
| **Paired-end overlap discordance**: ErrorProfiler (yaotianran/ErrorProfiler), or our own implementation over NGmerge/fastp merges | Mate agreement in the overlap | No | Exact **cycle position × centre Q × neighbouring Q × substitution × context × R1/R2**, jointly, with both mates' qualities at the same template base. PCR errors cancel because both mates share them, so it isolates sequencer error | Illumina only; needs short inserts (fully overlapping for ErrorProfiler); substitutions only; which mate erred is ambiguous (resolve by symmetric EM over both mates' Q, **not** by trusting the higher Q); covers the read ends of both mates | ErrorProfiler: 2025 preprint, code on GitHub; NGmerge 2025-11; fastp 2026-09 |
| **Self-reference ("assemble, then align")**: metaFlye / myloasm / hifiasm-meta (long), metaSPAdes (short) → minimap2 → pileup | Polished consensus of high-coverage contigs | No external one (the sample assembles its own) | **Everything alignment gives**: full quality windows, indel lengths, homopolymer run-length errors, joint Q × position × context, R1/R2, per-read error and quality trajectories (latent-state models), insert size. The simulators' own profilers can then run directly on the BAM | Assembly errors become false "truth" (ONT homopolymers especially), so polish (medaka) and mask sites with minor-allele frequency above a threshold, like BQSR known sites. Strain mixtures and repeats need masking. Aligner scoring bias (minimap2 favours mismatches over indels; skiver paper). Compute-heavy. Only high-coverage genomes, as with skiver | Flye 2026-04, myloasm 2026-09, hifiasm-meta 2025-11, SPAdes 2026-09, minimap2 2026-05, medaka 2026-05 |
| **ONT duplex vs simplex** (dorado duplex) | Duplex read, as near-truth for its two simplex parents | No | **Long-read truth per molecule** with full simplex Q tracks: homopolymer indels, per-read error/quality trajectories, without assembly | Needs duplex-capable chemistry and runs; only the duplex-paired fraction of reads; residual duplex errors are correlated with simplex errors in homopolymers | dorado 2026-09 |
| **Raw FASTQ quality profiling** (built in; same idea as `art_profile_builder`) | **None, never evidence of error** | No | The **quality process only**: per-position, per-mate Q distributions, Q autocorrelation (Q transitions), Q given observed bases, read-length distribution. Needed to *generate* qualities in default mode | Carries no error information and must not be read as an error rate: reported Q is often miscalibrated (Q binning; ErrorProfiler found systematic underestimation of accuracy). Conditions on *observed* rather than true bases, a good approximation when errors are rare (Illumina) and a poor one for ONT | n/a |
| **UMI / duplex consensus** (fgbio) | Molecular consensus | No | Per-molecule truth with qualities for amplicon and targeted libraries; separates PCR from sequencer error | Only UMI libraries | fgbio 2026-08 |
| **Spike-in and control references**: Illumina PhiX, ONT lambda/DCS, mock communities (e.g. ZymoBIOMICS) | Known reference | Yes, but a free, exact one | Truth-grade alignment profiles with qualities for the *same run*, including Illumina InterOp per-cycle error metrics | The control's library prep differs from the sample's; PhiX reads often sit in `Undetermined` or are removed; mock-community strain drift | Illumina/interop 2026-04 |
| **Existing error tables (importers)**: GATK BQSR recalibration tables; DADA2 `learnErrors` matrices | Reference + known sites (BQSR); denoised ASVs (DADA2) | BQSR yes; DADA2 no | BQSR: joint read group × Q × cycle × context mismatch counts. DADA2: reference-free substitution × Q error matrices for amplicons | BQSR needs known sites, which rarely exist for metagenomes. DADA2 is substitutions only, amplicon only, and its matrices are smoothed fits rather than raw counts. Both condition on centre Q only | GATK 2026-09, dada2 2026-07 |
| **k-mer spectrum** (GenomeScope2) | k-mer histogram model | No | An independent global rate check | Single rate, isolate-like genomes only (poor for metagenomes) | 2025-10 |
| **Error-corrector diffs** (Lighter, Rcorrector; BFC unmaintained since 2016) | Corrected read | No | Per-read substitution calls with Q, position and mate | Miscorrection and under-correction are coverage-dependent and unmeasured, and correctors often use Q themselves, which would leak Q into the labels. **Not recommended** | Lighter/Rcorrector 2026-01 |

**Recommendation** (matches §12, decision 6):

1. **Paired-end overlap → `pe-overlap` mode.** The cheapest reference-free evidence that identifies **neighbouring-quality and quality × context effects** for Illumina, which default-mode skiver cannot.
2. **Reference / self-assembly alignment → `reference` mode.** The most general evidence, especially for long reads (homopolymers, indel lengths, read-level quality/error regimes). It also gives a real-data **check of skiver-derived models** on the same high-coverage genomes, which the fork only ever validated on synthetic data.
3. **`kmer` mode fitting** after both, so it is checked against them from the start.
4. **ONT duplex/simplex source.** The best reference-free long-read truth when duplex data exist.
5. **Raw FASTQ quality profiling.** Low as *evidence* (it carries none), but a hard dependency for *generating* qualities in every mode, so it was built first (phase 1) as a feature/output source.
6. **Importers**: BQSR tables, DADA2 matrices, InterOp metrics.

### 4.4 Closest prior art and redundancy risks

Most parts of this project exist somewhere already. Only the combination below is new, so the plan must show it adds something on real data rather than assume it.

| Tool | Reference-free training | Simulates | Base context | Quality model | Errors vs Q |
|---|---|---|---|---|---|
| **ReSeq** (Illumina) | No (BAM) | Yes | No preceding-sequence model (stated limitation); per-site systematic errors | Rich: previous Q, position, mate, tile, sequence quality, reference base (2-D margins) | Errors conditioned on Q, learned from data |
| **Badread** | No | Yes | 7-mers | Q given local CIGAR window | Errors first, then Q |
| **CycSim** (long reads) | No (BAM) | Yes | k-mer sliding window + error-state transitions | Assigned after errors | Errors first, then Q |
| **LongISLND** | No | Yes | Per-k-mer error patterns | Position scaling | — |
| **BEAR** (DRISEE) | **Yes** (artifactual duplicate reads) | Yes | None (per-position rates) | Mean Q of erroneous bases vs position | Q fitted to errors |
| **skiver** (upstream) | **Yes** | **No** (QC only) | Marginal trinucleotide spectrum | Marginal P(error \| Q) | — |
| **This project, `kmer` default mode** | Yes | Yes | Two-sided, latent-position `Context(L,R)` | FASTQ Q Markov process | Centre-Q term by marginal matching; no Q window |
| **This project, full model** | Yes (`pe-overlap`, enhanced `kmer`) or `reference`/duplex | Yes | Two-sided + `QualityxContext` | Q conditioned on true context, shared `Latent(S)` | Q first, then errors given the Q window; Q never a label |

**Clearly distinct:**

- simulator-ready models trained without alignment. BEAR is the only precedent, and it has no context and no Q coupling;
- quality as a feature, with the calibration gap learned rather than assumed;
- one spec exported to many simulators, with `--q-policy` and fidelity reports;
- multi-source fitting with disagreement diagnostics.

**Redundancy risks, each tied to a test:**

1. **The reference-free advantage is narrower than it looks.** Skiver needs ≥ ~20× coverage on some genomes, roughly the genomes that also self-assemble. On a self-reference BAM, ReSeq, Badread, CycSim and NanoSim can train directly. Skiver's remaining advantages are speed, no aligner bias and no assembly errors used as truth. *Test:* the `kmer`-vs-`reference` cross-check (phase 6, §12 decision 7).
2. **Default mode is only modestly richer than existing profiles.** It is roughly Badread-level context plus an ART/NEAT-level quality process. *Test:* comparative checks against BAM-trained baselines (phases 5–6, §12 decision 11).
3. **The richest components come from alignment-type evidence** (BAM, duplex) that existing profilers already use. PE overlap is the reference-free exception for Illumina. On those sources the contribution is the joint model and exporter hub, not the evidence.
4. **Exports discard what is new** (§6.3). Exporter users get about what the simulator's own profiler gives, minus the alignment step. Only the native generator carries the full model.
5. **Expressiveness isn't realism.** A simple Q-HMM (PBSIM3) won on quality realism in the 2026 ONT benchmark. Components must earn their place on held-out data (§12 decision 10).

---

## 5. Model formulation: bases and qualities together

### 5.1 Notation

For one read:

- **Template:** x = x₁…x_n, the true bases.
- **Read:** observed bases y and reported qualities q, related to x by an alignment a (a sequence of per-template-position ops: match, `sub_to_*`, `ins_*` before the position, deletion).
- **Covariates z:** position from the read start and end, mate (R1/R2), strand, GC of the fragment, platform/run.
- **Latent state s:** an optional read-level state s₁…s_n (e.g. a degraded-stretch regime).

Windows around template position t:

- **Base-context window:** x_{t−L..t+R}.
- **Quality window:** q_{t−m..t+m}.

### 5.2 Two heads, one joint model

The generative model of a read is **P(y, q, a | x, z) = Σ_s P(s) · P(q | x, z, s) · P(a, y | x, q, z, s)**, i.e. a *quality head* followed by an *error head*. Both are ordinary conditional models over the same feature vocabulary.

**Error head (E).** P(op_t | x_{t−L..t+R}, q_{t−m..t+m}, z_t, s_t).

This is the error profile "given the base and quality in question and the surrounding bases and qualities". Its outputs are the 10 categories: match, 4 `sub_to_*` with the true-base mask, 4 `ins_*`, and deletion. Indel length comes from an extension head when a truth source observes it.

**Quality head (Q).** P(q_t | q_{t−1..t−m}, x_{t−L..t+R}, z_t, s_t).

This is an autoregressive quality process, a Markov chain over read positions. It is conditioned on true context (e.g. quality drops after homopolymers, GGC motifs on Illumina) and on position and mate (end-of-read decay, R2 < R1).

**Why this factorisation (Q first, then errors given the whole Q window):**

- **Head E is exactly the quantity to report.** It is what truth-bearing sources observe directly (skiver dump, PE overlap, BAM and duplex all pair an op label with the read's qualities). It is also what BQSR-style calibration and quality-aware downstream tools consume.
- **Two-sided quality windows become possible at generation time.** The whole Q track is sampled first, so an error at t may depend on q_{t+1}.
- **Reported Q stays a feature, never a label.** No component interprets 10^(−Q/10) as an error probability. Calibration is *learned* by head E (§5.5).
- **Coupling in both directions is represented.** Errors cluster at low Q through E's dependence on q_t and its neighbours. Qualities cluster around troublesome sequence through Q's dependence on x. A shared latent state s drives both, for slow-varying regimes.

The alternative factorisation, errors first and then Q | local error pattern (as in Badread's Q model), is derived from the joint model when an exporter needs it, not fitted separately.

### 5.3 Alignment and indel conventions

Qualities live on *read* positions, errors on *template* positions. The conventions match the fork's `skiver dump --base` so training data and generation agree:

- **Match / substitution** at template position t emits one read base with quality q_t.
- **Deletion** at t emits no read base and consumes no quality. The quality window for neighbouring ops is taken over read positions, so a deletion is flanked by the qualities of its neighbouring read bases.
- **Insertion** before t emits an extra read base. Its quality comes from an insertion quality sub-head, P(q_ins | q neighbours, x context, z), trained on insertion rows (`true_base = '-'`).
- **Generation order per read:**
  1. sample s;
  2. sample the template-indexed Q track from head Q;
  3. sample ops from head E using the Q window;
  4. materialise y, q and the CIGAR, drawing inserted-base qualities from the sub-head and dropping deleted positions' qualities.

A consistency test checks that re-deriving head E's training tuples from generated reads reproduces the sampled (op, Q-window) pairs exactly.

### 5.4 Parameterisation

Full tables are infeasible: 4^(L+R+1) contexts × Q-bins^(2m+1) × position × mate. So both heads are log-linear with additive and low-rank terms, generalising the fork's composable components:

| Component (spec token) | Head | Features | Replaces fork component |
|---|---|---|---|
| `Context(L,R)` | E, Q | per-offset base effects, two-sided (additive); `ContextTable(L,R)` for small dense tables | `AdditiveContext(N)`, `BaseContext(N)` (preceding only) |
| `QualityWindow(m)` | E | per-offset effects of binned or splined q_{t±k}, per op class; centre-Q term always included | `PhredContext(N)` (preceding only, averaged per row) |
| `QualityMarkov(m)` | Q | transitions from the previous m qualities (binned) | — (new) |
| `QualityxContext(rank)` | E | low-rank interaction of centre Q with centre/nearby bases and homopolymer length | — (new) |
| `Position(spline)` | E, Q | splines of distance from start and end | `Position(N)` (log1p features) |
| `Mate`, `Strand` | E, Q | categorical offsets | `Strand` |
| `Homopolymer` | E, Q | run length at and around t | `Homopolymer` |
| `Latent(S)` | E, Q | shared read-level HMM state scaling both heads | `HMM(S)` (errors only) |
| `IndelLength` | E ext. | length distribution by run length and Q | — (needs truth source) |
| `FragmentOverdispersion` | E (training only) | Dirichlet-multinomial between-fragment variance | same |
| `NeuralWindow(...)` (optional, later) | E, Q | small 1-D CNN/MLP over the joint base+quality window | — |

Model search uses the same greedy, criterion-based selection as the fork, now over both heads. The heads are fitted separately because their likelihoods factorise given the data, and their criteria are reported separately.

**Quality alphabets are platform-specific**, and the spec stores the alphabet, so generation emits only legal values:

- binned 4-level Illumina qualities (e.g. NovaSeq {2, 12, 23, 37});
- full-range ONT qualities (0–50);
- HiFi qualities (capped at 93).

### 5.5 Derived quantities

These are computed from the fitted heads, not fitted separately:

- **Empirical quality** of a base given its window: −10·log₁₀ P(error | window). This is a calibration diagnostic (reported against the reported Q) and the basis of quality-dependent export policies.
- **Recognition view**, P(true base ≠ observed | observed bases and qualities in the window, z). Obtained by Bayes over the true centre base using both heads. It is the reads-only error profile a quality-aware tool would use on real data, where the true context is unknown.
- **Marginals for exporters:**
  - per-position Q histograms (ART, InSilicoSeq);
  - Q Markov transitions per mate (NEAT);
  - P(Q | local CIGAR window) (Badread);
  - a Q-state HMM (PBSIM3 QSHMM);
  - Q mean/sd for correct vs erroneous bases (Mason2).

### 5.6 Identifiability by data source

| Parameter | `kmer` default (+ FASTQ) | `kmer` enhanced | `pe-overlap` | `reference` / duplex |
|---|---|---|---|---|
| E: base context | ✓ (latent position, `kvmer.csv`) | ✓ | ✓ (short context) | ✓ |
| E: centre Q | ✓ marginal only, joint with context by marginal matching (below) | ✓ | ✓ | ✓ |
| E: neighbouring Q, Q × context | ✗ (weights fixed at 0; flagged) | ✓ within value (key-side needs fork addition) | ✓ | ✓ |
| E: position, mate | position marginal only; mate ✗ | position ✓; mate needs fork flag | ✓ | ✓ |
| Q head | ✓ from FASTQ, conditioned on observed bases | ✓ conditioned on consensus | ✓ | ✓ |
| Q–E coupling via latent state | ✗ | partial (needs Q in `windows.bin`) | ✗ (short span) | ✓ |

**Default-mode marginal matching.** Unmodified skiver gives two error marginals, P(error | context) from `kvmer.csv` and P(error | centre Q) from `summary_phred.csv`. The FASTQ gives the exposure: the joint distribution of centre Q and observed context. Head E is fitted as log-additive context + centre-Q terms so that its implied marginals *under that exposure* reproduce both skiver marginals. This is raking/IPF-style, and it avoids double-counting when low Q and difficult context co-occur. Neighbouring-Q terms are unidentifiable here and fixed at zero. The spec records this, and exports carry it in their fidelity report.

**`summary_phred.csv`'s counting convention** (pinned in phase 1 from upstream v0.3.2 `PhredScoreSummary`): each observed value is scanned from its first base against the consensus and stops at the first mismatch. Every base up to and including that mismatch contributes its own Q (matches to `num_correct`, the mismatch to `num_error`), over value positions t ∈ [1 + ignore_smallest_t, v − ignore_largest_t] only. The exposure for marginal matching is therefore the hazard-model exposure, not all value bases. The reported `per_base_error_rate` is a per-Q Weibull λ with the global β, not the raw count ratio. Details are in `sources/skiver_analyze.py`.

---

## 6. Feasibility and limitations

### 6.1 Learnable in `kmer` default mode (unmodified skiver + raw FASTQ)

| Model element | Source | Confidence |
|---|---|---|
| Per-base error rate, clustering (λ, β) | `summary_error_rate.csv` | High (the published method) |
| Error-type composition (12 subs + 4 ins + 4 del), strand asymmetry | `summary_error_spectrum.csv` | High |
| Head E base context, two-sided within k+v | `kvmer.csv` with latent edit position (§2.1) | Medium. Identifiable because position varies across keys; needs a recovery test |
| Head E centre-Q term | `summary_phred.csv` + FASTQ exposure (marginal matching, §5.6) | Medium. Correct only under the log-additive assumption |
| Marginal error rate vs read position (start/end) | `summary_read_position.csv` | High for rate; composition along the read assumed constant |
| Head Q: per-position/mate Q distributions, Q Markov transitions, Q given observed context | raw FASTQ | High as a *quality process*; carries no information about errors |
| GC dependence | `summary_gc_content.csv` | Medium |

### 6.2 Not learnable, or only weakly, in `kmer` default mode

- **Neighbouring-quality effects and quality × context interactions** in head E. Unmodified skiver never pairs errors with quality windows or with context-specific qualities. The effects are fixed at zero and flagged; closing this needs enhanced skiver or a truth-bearing source such as PE overlap.
- **Joint effects of position, strand and mate with context or Q.** Each CSV is a separate marginal. They are combined log-additively, an assumption that is untestable without per-observation data.
- **Quality head conditioned on true context and on errors.** FASTQ conditions on observed bases. That is acceptable for Illumina (error ≈ 10⁻³), but biased for ONT, where erroneous bases are frequent and correlated with low Q.
- **Indels longer than 1 bp and multi-error windows.** 2+-edit values are dropped, so the indel length distribution is unobserved. Assume 1 bp (fine for Illumina, wrong for ONT homopolymers).
- **Homopolymer effects.** `kvmer.csv` has only the longest run in the value, and the position of an indel inside a run isn't identifiable. Only a coarse run-length covariate is possible.
- **Read-level heterogeneity** jointly affecting Q and errors. Only the aggregate β is available.
- **R1 vs R2** differences in errors are unobservable; R1 vs R2 qualities come from the FASTQ.

### 6.3 Structural limitations (the first two items are `kmer`-specific; the rest apply to all modes)

- **Reference-free truth is the consensus.** Coverage ≥ ~20× on at least some genomes is needed. Strain heterogeneity, repeats and low-coverage metagenomes inflate or filter estimates, so the outlier filter matters. Document minimum-coverage guidance and surface `passes_filter` statistics in reports. Separating variation from error in every mode is §6.6.
- **FracMinHash subsampling** (`-c`) trades precision for memory. Low-error platforms need a low `c` or large inputs for stable context and quality-window tables.
- **Reported quality is not ground truth, anywhere.** Two consequences:
  - No source may derive error labels from Q, which is why error-corrector diffs are excluded.
  - Exporters to simulators that equate Q with error probability (ART, InSilicoSeq, NEAT's error sampling, PBSIM3 QSHMM) must choose a `--q-policy`:
    - `preserve-quality` (default): emit the learned Q distribution; the simulated error rate then follows the simulator's assumption and deviates by our calibration gap.
    - `preserve-errors`: remap the emitted Q per position so the simulator's Q-derived errors match head E; the Q distribution is then distorted.

  The fidelity report quantifies whichever is violated.
- **Simulator formats are projections.** ART, InSilicoSeq, NEAT and PBSIM3 lose base context; all simulators lose neighbouring-quality effects in head E; Badread keeps base context and error-conditioned Q but loses position. Only the native generator keeps the full joint model. Every exporter reports what was dropped.
- **Base profile fallback.** When no FASTQ is available, head Q starts from a *base profile* (a simulator's built-in profile or a user file) and the fallback is recorded in the provenance.
- **Simulator quirks.** PBSIM3 ERRHMM outputs `!` qualities (use QSHMM for Q-bearing exports). NEAT's model files are Python pickles of NEAT classes, a version-coupling risk. NanoSim's model directory mixes KDE pickles tied to its scikit-learn version.
- **skiver upstream drift.** The CSV schema changed between releases (v0.3.x added GC and survival files). Parsers validate headers and pin supported version ranges; a CI matrix runs released skiver binaries.

### 6.4 `pe-overlap` and `reference` modes

**`pe-overlap`** identifies, per overlapping template base, substitution × two-sided context × centre and neighbouring Q of both mates × cycle position × mate. It does not identify:

- **indels.** The overlap is placed gaplessly; pairs whose discordance implies an indel are dropped and counted;
- **errors both mates share** (PCR, library, cluster-level), which cancel;
- **the read middle** when inserts are long. Position coverage depends on the insert-size distribution, so position effects are fitted against the recorded overlap exposure;
- **which mate erred** when they disagree. This is resolved by symmetric EM with head E over both mates' windows, never by trusting the higher Q;
- **both mates erring to the same base** at one position. This is rare for Illumina and is bounded as a nuisance term.

Overlap placement must not use Q, so Q cannot leak into the labels.

**`reference`** identifies everything in the last column of §5.6, including indel lengths, homopolymer errors and per-read trajectories. Its risks:

- **Reference errors become false errors.** External genomes diverge from the sample's strain, and self-assemblies carry consensus errors (ONT homopolymers especially). Mitigations: polish; mask sites whose minor-allele frequency is above a threshold; keep only high-coverage contigs; report per-contig rates so outliers show.
- **Strain variation.** Minor strains and within-population alleles show as mismatches against any single consensus. Minor-allele masking is outcome-dependent and can remove true error hotspots (§6.6).
- **Aligner scoring bias** between mismatches and indels (minimap2). This is measured by aligning generated reads with known CIGARs.
- **Spike-ins** (PhiX, lambda) give exact truth, but from a different library prep than the sample.

### 6.5 Verdict

Feasible, in stages:

- **`pe-overlap`** (Illumina) and **`reference`** (all platforms) yield the full §5 model directly, within the scopes above. They are the quickest route to a working generator, and the yardstick for the `kmer` mode.
- **`kmer` default mode** yields a context-aware error head with a centre-quality term, plus a realistic quality process. Generation outputs full FASTQ (bases + qualities) for the native generator. Credible Illumina exports and a context-preserving Badread export follow.
- **The full base × quality-window model** from reference-free k-mer evidence needs per-observation quality windows, i.e. enhanced skiver. For Illumina, `pe-overlap` already supplies them reference-free; `reference` and duplex supply them for all platforms. Long-read realism (homopolymers, indel lengths, read-level regimes) needs the same.

### 6.6 Biological variation vs sequencing error

In a metagenome, a read can differ from its "truth" because of an error, or because the molecule really differs from that truth: a co-existing strain, a minor within-population allele, a divergent repeat copy or paralog, or an unpolished consensus. Counted as errors, these inflate head E, and not evenly. Variation concentrates in particular genes, contexts (transitions, third codon positions) and coverage levels, so it biases context and op composition as well as the rate.

**How each mode is exposed.**

| Mode | What variation looks like | Exposure |
|---|---|---|
| `pe-overlap` | Both mates read one molecule, so a strain allele is read by both mates and agrees. It cancels, like a PCR error | Near-immune. The remaining risk is placement: mates overlapping divergent copies of a repeat longer than the read (ERR10889147, phase 4) |
| `reference`, external genome | Mismatches at every site where the sample's strain differs; reads from absent relatives mis-map | High. In the skiver paper the alignment-based error rate nearly doubled as reference ANI fell to 96% |
| `reference`, self-assembly | The consensus is the majority strain; minor strains' alleles appear as mismatches at their frequency; collapsed repeats | Medium, concentrated at minor-allele sites |
| `kmer` | A key with several true values | The outlier filter (per-t hazard above median + 3·IQR, or a binomial test against the fitted hazard) drops keys with a frequent second value. Minor values below its power pass as errors. On clonal data it can also drop true error hotspots (fork lesson 7: it hides generated errors) |

**How other tools separate them.**

| Tool | Signal | Q as error probability? | What we take |
|---|---|---|---|
| skiver | Per-key hazard outliers; reference-guided mode drops keys with several reference values | No | Kept as the `kmer` default. Measure its minor-variant recall and its cost on clonal hotspots |
| inStrain | Minimum coverage, minimum variant frequency, and a null model of error counts at that coverage that assumes Q30 | Yes (fixed Q30) | The coverage-dependent count test, with head E's expected count in place of an assumed Q |
| LoFreq | Poisson-binomial test of variant counts from per-base Q; strand-bias filter | Yes | The site count test structure, not Q as probability; strand bias only as a model check |
| DADA2 `learnErrors` | Alternates error-rate estimation and sample inference (abundance p-values under the error model) until they agree, starting from "only the most abundant sequence is correct" | No: rates per Q are learned | **The core loop:** fit head E on sites believed clonal, re-score sites under head E, repeat |
| UNOISE3, Deblur | Abundance skew of a sequence over its 1-edit neighbours; a static upper-bound error profile | No | Abundance skew is the `kmer` analogue (consensus vs neighbour counts). Static profiles are what this project replaces |
| GATK BQSR | Mask known variant sites before counting mismatches | No (tabulates by Q) | Site masking, with sites found in the sample rather than a database |
| DESMAN | Variant frequencies co-vary across samples with strain abundance; errors don't | No | Multi-sample consistency, when several samples of one community exist |
| Floria, Strainy (phasing); VeChat, DeChat (haplotype-aware correction) | True variants co-occur on reads and pairs (MEC clustering, variation or de Bruijn graphs); errors are independent across reads | No | A linkage test |
| Duplex, UMI consensus | A variant is on both strands of a molecule, an error on one | No | Already an evidence source (phase 11) |

**Signals, and the rules for using them.** A site-level decision is part of the label pipeline, so the "Q is never a label" rule applies to it. A site-level decision based on a feature that head E models also biases that component.

1. **Within-molecule concordance** (PE overlap, duplex). Q-free and variant-immune; the yardstick for how much variation contaminates the other modes.
2. **Allele count vs coverage.** At high coverage a site's error count is tightly bounded by the error model, so a variant at a frequency well above the error rate stands out. Test counts against head E's expected count for that site's reads (their contexts, Q windows, positions), never against 10^(−Q/10) or a fixed Q30. Variants near the error rate are not separable this way; phase 7 measures the floor.
3. **Linkage.** Two non-consensus alleles on one read or pair co-occur far more often than independent errors would. Q-free, and independent of every head E feature. It needs two or more variant sites within a read or insert, so it finds moderately divergent strains, not isolated SNPs; long reads extend its reach.
4. **Independence from read covariates.** A variant has a similar frequency on both strands, both mates, every cycle position and many distinct read starts; errors follow strand, mate, position and context. **Caution:** these are head E components. Filtering sites on strand or position bias removes true strand- or position-dependent errors and shrinks those components. They may be used only inside the joint site model below, where head E supplies the expected bias, and their effect is tested on clonal simulations.
5. **Reported quality.** Mismatches spread evenly over reported Q look like a variant; mismatches concentrated at low Q look like errors. Q may enter only through fitted head E in the joint site model, with its calibration learned, never as a threshold, a fixed rate or a filter of its own. Q is also a **diagnostic**: after separation, the mismatch rate by reported Q should be flat at called variant sites and follow head E at retained sites.
6. **Cross-mode agreement.** On the same reads, `pe-overlap` is variant-immune. The `reference` or `kmer` excess over `pe-overlap` on shared support estimates PCR/library error plus residual variation; separation should shrink it to the PCR/library excess measured without variation.
7. **Multi-sample consistency** (DESMAN-style). Optional, when several samples of one community exist.

**Approach.** Two layers, each tested on clonal truth before it runs on variation:

- **Unambiguous sites first (conservative mask).** Keep sites with high coverage, coverage consistent with a single copy, no linked variants, and (`kmer`) keys passing the outlier filter. Head E conditions on context, Q window and position, so selecting sites on Q-free covariates that aren't the outcome leaves it unbiased. Selecting on the observed mismatches (minor-allele frequency, hazard outliers) is outcome-dependent and removes error hotspots. Every mask reports the op and context composition of what it drops, and clonal simulations measure its bias.
- **Joint latent-site model** (DADA2-style self-consistency). Each site has latent alleles with frequencies; a read's base is drawn from the site's alleles and then passes through head E. Alternate: (i) fit head E with each read row weighted by the posterior that the read carries the consensus allele; (ii) update site posteriors from allele counts under head E, linkage and (optionally) multi-sample terms. Start from the conservative mask, never from Q. This reuses the `pe-overlap` soft-EM machinery (fractional counts, warm starts). Variant-site calls are a by-product and a diagnostic, not a product.

**Order.** Each step is gated on the one before:

1. **Clonal simulation** (phases 3–6). Every mode's fitters, filters and masks are shown unbiased on reads from single genomes, and each filter's cost when there is nothing to filter is measured.
2. **Simulated variation** (phase 7). A variation simulator with per-base truth (error or variant) measures how much strains, minor alleles and divergent repeats bias each unmodified mode.
3. **Separation** (phase 8). Masks and the joint site model are developed on the phase 7 scenarios, then validated on strain-resolved mock communities (ZymoBIOMICS D6331 has five E. coli strains at equal abundance; skiver's K-12/O157:H7 mixtures), and only then run on real metagenomes.

---

## 7. Architecture

```
src/sequencing_error_model/
  sources/skiver_analyze.py  # parse + validate v0.3.x analyze CSVs → marginal count tables (default mode)
  sources/skiver_dump.py     # streaming reader for fork dump TSV / windows.bin (enhanced mode)
  sources/fastq_quality.py   # quality process: per-position/mate Q, Q transitions, Q | observed context, lengths
  sources/pe_overlap.py      # pe-overlap mode: mate discordance with both mates' quality windows
  sources/bam.py             # reference mode: external / spike-in / self-assembly alignments → tuples (site masking)
  sources/ont_duplex.py      # simplex-vs-duplex alignments
  sources/importers.py       # GATK BQSR tables, DADA2 error matrices, InterOp metrics
  sites.py                   # variation vs error (§6.6): conservative masks, joint latent-site model
  variation.py               # test-only variation simulator: strain haplotypes, minor alleles, repeats + truth sites
  observations.py            # sparse (context window, Q window, covariates, op) → count tuples; truth label per source
  spec.py                    # ErrorModelSpec: versioned JSON manifest + .npz arrays (no pickle)
  fit/                       # quality head, error head, marginal matching, latent state, selection criteria
  select.py                  # greedy criterion-based component search over both heads
  compare.py                 # cross-mode comparison on shared support: count tables and fitted specs
  generate.py                # native generator: FASTA → FASTQ (bases + qualities) + CIGAR, single or paired
  export/{art,iss,badread,pbsim3,neat,mason,wgsim,nanosim}.py
  cli.py                     # sem fit | select | compare | generate | export <sim> | report | inspect
```

Principles:

- **One simulator-neutral model spec** (`spec.py`), the only thing exporters and the generator read. It holds:
  - provenance: evidence mode(s), sources, input hashes, base-profile fallback; for `kmer`, the skiver build and version, k, v, c;
  - the **quality alphabet** and binning;
  - **head Q** and **head E**: component tokens, parameters, per-component `generative` and `identified` flags (e.g. neighbouring-Q weights fixed in default mode);
  - the optional latent-state layer shared by both heads;
  - summary marginals cached for exporters and reports (rate, λ/β, op composition, position curves, calibration curve).

  Stored as JSON + `.npz`: safe to load, diffable, language-agnostic. It is not a torch pickle.
- **Sparse tuples first.** Every source reduces to deduplicated (context window, quality window, covariate bins, op) tuples with counts. Binned qualities make this compress hard for Illumina. When unique tuples exceed memory (long windows, ONT's full Q range), fitting falls back to streaming minibatches over the same schema. Marginal-only sources (skiver analyze, importers) enter the joint likelihood through their marginal constraints (§5.6).
- **Joint fitting across sources.** The per-source log-likelihood of the model is marginalised onto each source's observed fields and summed, with optional per-source nuisance terms (skiver consensus misses, assembly errors, mate ambiguity). The report shows per-source disagreement before anything is pooled.
- **Modes are compared on shared support.** `compare.py` restricts two sources' count tables to the fields, op classes, read positions and genomes both observe. It then compares the tables (evidence level) and the specs fitted from them (per component and held-out likelihood). Each mode pair declares its expected gaps (e.g. PCR errors are invisible to `pe-overlap`), and those gaps are reported as estimates.
- **Fitting stack.** numpy/scipy (L-BFGS) for the default mode's small, marginal problems. torch (and optionally pyro for VI uncertainty) as an `enhanced` extra for per-observation windows, latent states and the optional `NeuralWindow`. Revisit if one stack proves sufficient.
- **The generator is the reference implementation of §5.3.** It is vectorised per read (sample the Q track, then batched categorical ops) and emits bases, qualities and CIGAR. Exporters are validated against it.
- **Exporters are thin, tested adapters.** Each declares:
  - required spec fields;
  - its `--q-policy` behaviour;
  - an optional base profile;
  - a machine-readable *fidelity report* of what was dropped or approximated.
- **Optional simulator dependencies** (NEAT classes, InSilicoSeq's loader) are imported lazily and only needed for export or round-trip tests.

---

## 8. Mode capability matrix

| Component | `pe-overlap` | `reference` | `kmer` default | `kmer` enhanced | Generative | Exporters that can use it |
|---|---|---|---|---|---|---|
| Global rate, op composition | substitutions only | ✓ | ✓ | ✓ | ✓ | all |
| PCR/library errors | excluded (cancel) | included | included | included | – | (comparison diagnostic) |
| Biological variation (strains, minor alleles) | excluded (both mates agree), except divergent repeats | contaminates; separated in phase 8 | partly filtered (outlier filter); separated in phase 8 | as default, plus linkage via read ids | – | (diagnostic, §6.6) |
| Strand asymmetry | ✓ | ✓ | ✓ | ✓ | ✓ | InSilicoSeq (fwd/rev), native |
| E: `Context(L,R)` | ✓ (exact position) | ✓ | ✓ (latent position) | ✓ (exact position) | ✓ | Badread, native |
| E: centre Q | ✓ | ✓ | ✓ (marginal matching) | ✓ | ✓ | all Q-coupled exporters via `--q-policy`, Mason2 (correct/wrong Q), native |
| E: `QualityWindow(m)`, `QualityxContext` | ✓ | ✓ | – (fixed at 0) | ✓ (value window; key side needs fork addition) | ✓ | native only |
| Q: per-position/mate Q, `QualityMarkov(m)` | ✓ | ✓ | ✓ (FASTQ) | ✓ | ✓ | ART, InSilicoSeq, NEAT (order 1), PBSIM3 QSHMM, native |
| Q: conditioned on true context / errors | ✓ (in overlap) | ✓ | approx. (observed bases) | ✓ | ✓ | Badread Q model, NanoSim, native |
| Position | ✓ (overlap region, exposure-weighted) | ✓ | ✓ (marginal) | ✓ (joint) | ✓ | ART, InSilicoSeq, NEAT, Mason (ramps), native |
| GC | ✓ | ✓ | ✓ | ✓ | ✓ | native |
| Clustering β | within overlap | ✓ | ✓ | ✓ | approx. | PBSIM3 (2-state approximation), native |
| `Latent(S)` shared by E and Q | – (short span) | ✓ | – | ✓ (needs Q in `windows.bin`) | ✓ | PBSIM3 QSHMM/ERRHMM, native |
| Homopolymer / indel length | – | ✓ | coarse | ✓ (needs new fork outputs) | ✓ | Badread (via k-mers), native |
| Fragment overdispersion, R1/R2 | R1/R2 ✓ | ✓ | – | ✓ | training only / R1–R2 profiles | ART/ISS/NEAT R2 profiles |

---

## 9. Phases

Each phase lands as one or more PRs with green CI. Exit criteria are testable. The order follows §1.1: the model and generator (phases 2–3), then `pe-overlap` and `reference` (4–5), then `kmer`, checked against both (6). Phases 3–6 establish each mode on clonal simulations and near-clonal real data (isolates, spike-ins, amplicon mocks). Biological variation is then simulated (7) and separated from error (8) before any exit depends on a real metagenome (§6.6).

### Phase 0: repository, CI and PR policy ✓
uv/ruff/mypy/pytest, pre-commit (revs standardised with the MIMICC/ENA repos), Linting, Testing and Release workflows, PR template, and a `main` ruleset requiring PRs with passing `lint` + `test`.

### Phase 1: inputs (skiver analyze, FASTQ quality, observation schema) ✓
Landed before the evidence modes were added. The FASTQ statistics and the schema serve every mode; the skiver parsers are consumed in phase 6.
- ✓ `sources/skiver_analyze.py`: typed, header-validated parsers for all v0.3.x analyze CSVs. `summary_phred.csv`'s counting convention is documented (§5.6).
- ✓ `sources/fastq_quality.py`: streaming quality-process statistics from raw FASTQ(.gz), per mate (one call per mate, several lanes allowed). Labelled in code and spec as a feature/output source, never as error evidence. Collects sparse counters for:
  - per-position Q histograms (from the read start);
  - Q transition counts (order m, for t ≥ m);
  - Q | observed base context (`flank` = (L, R), `.`-padded at read ends), which is also the joint centre Q × observed context exposure used for marginal matching;
  - read-length distribution.
- ✓ `observations.py`: the sparse tuple schema shared by all sources. A `CountTable` is a sparse count over a subset of one field vocabulary (`op`, `context`, `locus`, `q`/`q±i`, `t`, `pos_start`/`pos_end`, `strand`, `mate`, `gc`, `length`); fewer fields means a marginal of the joint model. Tables declare their unit (base, value, read, error event) and whether they are truth-bearing, and the constructor rejects an `op` field on non-truth tables, so quality-only sources can't carry error labels. Positions and `t` are 1-based value/read positions throughout (skiver's k + t is rebased). `skiver_analyze.tables` and `fastq_quality.tables` convert both default-mode sources; the fitted-parameter files (`summary_error_rate.csv`, `survival_rate.csv`) stay on `SkiverAnalyze` as passthrough, and `kvmer.csv` keys failing skiver's filter are dropped.
- ✓ **Fixtures.** A tiny synthetic genome plus reads with known injected errors *and qualities* (`tests/fixtures/make_skiver_fixtures.py`, seeded), analysed by the skiver v0.3.2 x86-64 release binary (taken from the `skiver-compat` artifact; an arm64 build counts slightly differently). The CSVs are committed (`tests/fixtures/skiver-v0.3.2/`, ~130 KB); the reads are regenerated, not committed.
- ✓ **CI job `skiver-compat`.** Downloads skiver release binaries (matrix: v0.3.1, v0.3.2, latest), regenerates the fixtures, runs the parser tests on them and diffs the deterministic files (bootstrap CIs vary run to run) against the committed fixtures. It runs on PR only when `sources/skiver_*` changes, plus a weekly schedule to catch new releases.
- ✓ **Exit:** every analyze file and the FASTQ statistics round-trip into the schema; unsupported versions fail with a clear error.

### Phase 2: model spec and both heads, fitted from exact-position observations (in progress)
- ✓ `spec.py`: an `ErrorModelSpec` directory of `spec.json` (schema version, quality alphabet, provenance with required `mode`, component tokens, `generative`/`identified` flags, meta) and `arrays.npz` (parameters and cached marginals; object arrays are rejected and loading never unpickles). Tokens are validated against the §5.4 component names and the heads each may appear in, plus `GC`, `Weibull` and `InsertionQuality` from this phase; the shared `Latent(S)` layer sits outside both heads. Argument arity is left to each fitter. Loading rejects other schema versions, missing arrays and unreferenced arrays.
- ✓ Provenance `mode` is the evidence mode (`pe-overlap`, `reference`, `kmer`, or a list of distinct modes for a joint fit), with `skiver_build` (`default`/`enhanced`) required exactly when `kmer` is among them (`spec.MODES`, `spec.SKIVER_BUILDS`).
- ✓ Head Q fitters (`fit/quality.py`): `QualityMarkov(m)` (required; holds the bias), `Position(n)` (linear splines over log distance from start and end), `Mate` and `Context(L,R)`, as one log-linear softmax over the quality alphabet fitted by L2-penalised L-BFGS (scipy) from a joint `CountTable` with a `q` field. `sample` draws Q tracks from the same design matrix, for tests now and the generator in phase 3. Recovery test: on tuples sampled from a known head, lag, context and mate effects correlate r > 0.9 and per-position Q histograms are within TV < 0.05.
  - Still to do: the insertion-quality sub-head (needs insertion rows, so it waits for `reference` tuples), and a joint (lags × position × context) table from `fastq_quality`, whose current counters are separate marginals.
- Input is FASTQ statistics (every mode, conditioned on observed bases) or per-observation tuples (conditioned on true bases, from `pe-overlap` and `reference`).
- ✓ Head E fitters over exact-position tuples, the form `pe-overlap` and `reference` produce (`fit/error.py`): one log-linear softmax over 10 categories (match, 4 substitutions, 4 insertions, deletion) with the substitute-to-self category masked in the likelihood and in `probabilities`, fitted by L2-penalised L-BFGS. Components:
  - `QualityWindow(m)` (required; holds the bias and always the centre-Q term), per-offset alphabet effects plus "beyond the read end";
  - `Context(L,R)`; `Homopolymer` (run length of the centre base inside the table's context window);
  - `Position(n)`, `Mate`, `Strand`, `GC(n)` (linear spline over the GC bin midpoint);
  - `QualityxContext(r)`: a rank-r bilinear interaction of centre Q with x_{t−1..t+1}, sum-to-zero over bases so it can't absorb a centre-Q main effect.

  Each row is one template base with a true centre base; an insertion and a substitution at the same base can't both be counted yet (revisit with `reference` tuples). Recovery test on labels sampled from a known head: Q-window, context (substitutions) and homopolymer (indels) log-odds correlate r > 0.9, the error rate per reported-Q bin is within 10% and the marginal rate within 5%, with a truth far from 10^(−Q/10).

  The `kmer` default fitters (latent edit position, marginal matching) are in phase 6.
- ✓ Criterion-based selection per head (`select.py`): greedy forward search as in the fork, screening variants of the required first component (`QualityMarkov(m)` / `QualityWindow(m)`), then adding the candidate group that most improves the criterion until none does. Both heads expose `log_likelihood`; criteria are scored on a held-out table split by read. The default is held-out log-likelihood (`test-ll`); `aic`/`bic` are available, but with raw parameter counts (gauge and masked cells included) BIC rejected real position and neighbouring-Q effects in the recovery tests. Tests: head Q selects every true component; head E keeps the true Q window and context and rejects `Mate` and `Strand`, which have no true effect.
- **Exit:** on observation tuples sampled from a known spec:
  - recovered context and quality-window log-odds correlate with truth (r > 0.9 for dominant effects);
  - the error rate implied per reported-Q bin is within 10% of truth;
  - per-position Q histograms are within TV < 0.05;
  - the marginal error rate is within 5%.

### Phase 3: native generator (bases + qualities, single or paired) and recovery harness ✓
- ✓ `generate.py` (§5.3), vectorised across a batch's bases: sample the template-indexed Q track (head Q), draw head E outcomes from each base's Q window, then materialise bases, qualities and CIGAR. As in the fork, an insertion draws again at the same base (up to `max_ins_run`), so head E rows are per draw. Deletions drop q_t. Inserted bases reuse q_t until `InsertionQuality` exists. Non-ACGT template bases pass through as N matches. `Homopolymer` now stores its fitted flank in `meta` so generation rebuilds the same window.
  - `align` / `observations` re-derive head E tuples from template + read + CIGAR. The Q window is template-indexed from the read; a deleted base takes the next read base's Q (the fork's dump gives deletions no Q), so windows touching a deletion differ from the generator's track. **Open:** settle the deletion-Q convention with `reference` tuples (phase 5).
  - CLI `sem-generate` keeps the fork's `skiver-generate` flags and `@name cigar:CIGAR` output; `--model` is a spec directory (no presets or `.pt`). genome-blender switches via its generate command setting.
  - Tests: CIGAR/sequence/quality lengths agree and read-derived qualities equal the Q track at every non-deleted base; refitting both heads from the generator's own CIGARs recovers the marginal error rate (5%), per-Q rates (10%) and Q lag/context effects (r > 0.9); the CLI is deterministic per seed and matches the header contract.
- ✓ Recovery harness (`recovery.py`): `recover(truth, n_reads, seed, mode)` generates reads, fits them with a mode (a function templates, reads, mates, truth → spec; `cigar_mode` refits from the generator's own CIGARs), and `compare`s on held-out generated reads. Error metrics are expectations under both specs on the same tuples (op TV/KL, marginal rate, rate by read position and reported Q, empirical Q); quality metrics come from Q tracks generated by each spec (per-position TV, lag-1 autocorrelation). `TOLERANCES` holds the phase 2 exit criteria and `Report.failures()` checks them. `example_spec()` is the shared built-in truth. The small version runs in `task test`; `python -m sequencing_error_model.recovery --reads N --output report.json` runs at scale, via the manual `recovery` workflow.
- ✓ Standalone paired output: `fragments` samples placements uniformly over genome contigs with insert sizes from the spec's `insert_size` marginal (row 0 sizes, row 1 probabilities; `insert_sizes(mean, sd)` builds a discretised normal). Mate 2's template is the reverse complement, and short inserts read through into TruSeq adapters, then N. `sem-generate --pairs N --read-length L [--insert-mean M --insert-sd S]` writes interleaved `contig:start-end#i/1|2` records. Tests: mate templates match the named placement on both strands, including adapter read-through, and the sampled insert mean and sd match the marginal.
- `generate.py` implementing §5.3. Qualities are always emitted from head Q, and `--no-quality` is kept only for compatibility. It keeps the `skiver-generate` CLI contract so genome-blender can switch without code changes, and adds an alignment-consistency self-test.
- Paired-end output with an insert-size distribution stored in the spec. The generator is the truth-known fixture source for every mode: overlapping pairs for `pe-overlap`, reads plus their true CIGARs for `reference`, reads for skiver.
- A mode-agnostic recovery harness: generate from a known spec, run one mode's source and fitters, compare with the spec. Each mode plugs in as it lands. Metrics:
  - **errors:** TV/KL on op probabilities, marginal rate, position curve;
  - **qualities:** per-position Q TV, Q lag-1 autocorrelation, P(error | reported Q) calibration curve vs source, empirical-vs-reported Q curve.
- A small version runs in CI (a few seconds of reads); the full-scale version is a manual `workflow_dispatch` workflow.
- **Exit:** re-deriving tuples from the generator's own CIGARs and refitting recovers the spec for both heads within the phase 2 tolerances; the alignment-consistency self-test passes.

### Phase 4: `pe-overlap` mode (exit met; open items below)
- ✓ `sources/pe_overlap.py` (Illumina), streaming over paired FASTQ(.gz):
  - ✓ place the overlap per pair gaplessly, handling read-through into adapters when the insert is shorter than the read. Placement uses bases only, never Q (lowest mismatch fraction over ≥ `min_overlap` bases, then longest);
  - ✓ emit one head E row per mate per overlapping base, in the mate's own orientation: Q window, cycle position, mate, strand, and the consensus context (N where mates disagree, the mate's own bases outside the overlap);
  - ✓ where mates disagree, the template base is latent. Soft EM with head E over both mates' rows attributes the error, starting from a coin flip, never the higher Q. Rows carry expected (fractional) counts, and `fit.error.fit(init=)` warm-starts each refit;
  - ✓ the rows are the overlap exposure per read position, so position effects are fitted conditional on what the overlap covers;
  - ✓ drop and count pairs whose discordance implies an indel (a breakpoint plus a shift of ≤ 3 removes ≥ 3 mismatches, checked before the mismatch cap) or whose best overlap is shorter than `min_overlap` or above 20% mismatches. A tight cap (NGmerge's 10%) dropped the error-rich pairs and biased rates down 8% in recovery.
  - ✓ head Q is fitted from every base of every pair, placed or not, on observed bases. Fitting it on kept overlaps only biased it (per-position TV 0.067), because the dropped pairs are the low-Q ones.
- ✓ The spec records `identified_ops: ["substitution"]` and gets −inf indel logits, so generation from it emits no indels; PCR/library errors cancel and are absent. `recovery.compare(substitutions_only=True)` compares on that support; `recovery.paired_draw` feeds the harness pairs from `generate.fragments`.
- ✓ Recovery test on generated pairs (3000 pairs, ~3.5% errors) passes the phase 2 tolerances. Reads are regenerated, not committed. At ~1% errors the per-Q rates settle within 10% only by ~10k pairs; there, Q37 (rate 0.002) still comes out 14% high. The cause was L2 shrinkage of the rarest window weights, which the `smooth` + `window_l2` defaults below reduce.
- ✓ CLI: `python -m sequencing_error_model.sources.pe_overlap R1 R2 --output SPEC` writes a spec whose provenance carries the pair statistics.
  - Defaults: head E `QualityWindow(1) Context(1,1) Position(4) Mate` with `--smooth 30 --window-l2 0.01`; head Q `QualityMarkov(1) Position(48) Mate Context(1,1)`.
  - Usage and limits: [pe_overlap.md](pe_overlap.md).
- ✓ Real-data placement guards, from two runs:
  - **SRR5240881** (MiSeq 2×150, V3-V4/V9 mock amplicons): 97.7% of pairs placed, 0.15% indel-flagged, and a mismatch rate of 0.01% where both mates have Q ≥ 30.
  - **ERR10889147** (HiSeq 2×125, genomic, long inserts): the first version placed about 25% of pairs, and most were false: 7–15% mismatches even at Q ≥ 30. Two Q-free fixes:
    - an indel flag now needs the split alignment itself to pass the mismatch cap (unrelated mates had gained ≥ 3 mismatches from any shift, 525 of 2000 flagged);
    - a pair is ambiguous when another offset more than `max_shift` away also has a mismatch fraction below `max_rival` = 0.35, which catches repeats and low-complexity sequence. This threshold is absolute: a margin over the best offset rejected error-rich true overlaps in recovery (17% of 40 bp pairs).
  - After both fixes, 3% of ERR10889147 pairs place. Their Q ≥ 30 mismatch rate is still 1.6%, from about 0.3% of pairs that overlap a divergent copy of a repeat family longer than the read (they look like alpha-satellite). Bases alone can't flag these. They dominate only because true overlaps are so rare, so **long-insert genomic libraries are not `pe-overlap` evidence**; use amplicon or short-insert libraries, or the `reference` mode.
- Fitting cost: collecting rows takes about 1 ms per pair; the fits dominate (1000 MiSeq pairs, 27 Q values: 30 s of EM refits for head E, 23 s for head Q). **Open:** bin or subsample before fitting large runs.
- ✓ **Exit run** on SRR5240881, 20,000 pairs (12 min, dominated by fitting): 97.7% placed, 45 indel-flagged, 13,204 disagreements in 2.44 M overlap bases.
  - On 1,000 held-out pairs, the fitted marginal substitution rate is 0.29% per mate (raw mate disagreement is 0.6% per overlap base).
  - The fitted rate mostly falls as reported Q rises: about 7% at Q14–18, 0.1–0.5% at Q27–33, 0.02% at Q34–37 and 0.001% at Q38–39. Empirical Q is below reported at low Q and above it at the top, consistent with PCR errors cancelling.
  - **Open, error head position:** the fitted `Position(4)` curve swings from 0 to 6% per 25-bp bin and differs between mates, while the raw disagreement rate is flat at 0.45–0.62% across positions. Fixed-length amplicon reads make `pos_start` and `pos_end` collinear, so the paired splines are unidentified and only L2 pins them. Needs one position axis when read lengths don't vary, or selection to reject it.
  - **Open, head Q capacity:** on held-out reads the per-position Q TV median is 0.153 (max 0.60 at cycle 41) with `Position(4)`. With `Position(24)` it falls to 0.071 median, max 0.38. MiSeq's cycle-specific Q structure needs per-cycle terms or selection over n. The CLI default is now `Position(48)` (median 0.060); cycles 40–41 remain open (see the head Q table below).
- HiSeq-like simulation for a merge-rate target: 2×125 pairs from a gut genome (MGYG000175911) with the SRR5240881 spec. Insert mean 294, sd 50 makes ~10% of pairs overlap by ≥ 20 bases. Pilot on a random genome: 9.6% placed, as predicted. On MGYG000175911 (416 contigs), 20,000 pairs gave 10.5% placed and 22 indel-flagged, matching the theoretical 10.5%, so the placement guards cost nothing on this genome.
  - Fitting `pe-overlap` on these pairs takes 9.6 min, mostly head Q on all 5 M bases. On 1,000 held-out simulated pairs against the truth:
    - marginal substitution rate 1.07×, op TV 0.001, Q lag-1 autocorrelation 0.615 vs 0.614;
    - per-Q rates within 0.8–1.15× for most of Q12–30, with outliers at Q20 (0.34×), Q25 (1.74×) and Q27 (0.58×);
    - Q31–39 overestimated 1.7–11×. Only 670 disagreements were seen, so the rarest-Q window weights are shrunk toward the mean rate (the recovery Q37 effect, now stronger).
  - Smoothing `QualityWindow` weights across the Q alphabet (`fit.error.fit(smooth=)`) was tried against this bias. Mean |log2(fitted / true rate)| over Q ≥ 30, head E fitted by EM on the same 20,000 pairs:

    | smooth | first differences | second differences (kept) |
    |---|---|---|
    | 0 | 1.54 | 1.54 |
    | 3 | – | 1.66 |
    | 30 | 1.72 | 1.68 |
    | 300 | 2.32 | 1.61 |

    Q38–39 stay 8–14× too high at every setting, so the overestimate is bias, not sampling noise that borrowing strength could fix. First differences also flatten the curve, pulling well-observed bins toward sparse ones. On a synthetic table with a sparse Q40 bin, second differences do follow the trend (Q40 0.029 → 0.006, below Q30's 0.013), which the unit test pins.
  - **Cause: L2 on every weight.** The Q38–39 window weights must sit far below zero, and with few errors L2 pulls them back toward the overall rate; a curvature prior can't counter a pull on the level. `fit.error.fit(window_l2=)` sets a separate, weaker precision on the window weights. Same evidence and metric (mean |log2 ratio|, marginal rate ratio):

    | window_l2 | smooth | all Q | Q ≥ 30 | marginal | Q38 / Q39 |
    |---|---|---|---|---|---|
    | 1 (= l2) | 0 | 0.83 | 1.54 | 1.073 | 5.6× / 11.2× |
    | 0.01 | 0 | 1.44 | 1.10 | 1.017 | 2.9× / 7.3× |
    | 0.01 | 30 | **0.60** | 0.97 | 1.022 | 5.5× / 5.5× |
    | 0.01 | 300 | 0.64 | 0.96 | 1.020 | 5.8× / 6.0× |
    | 0.001 | 30 | 0.60 | 0.96 | 1.010 | 5.4× / 5.5× |

    `window_l2` 0.001 is indistinguishable from 0.01, so the `pe-overlap` default stays at the less extreme 0.01, with `smooth` 30.

    A weak window L2 alone collapses sparse bins (Q24 and Q26–29 at 0.03–0.13×). Smoothing alone can't move the level. Together they cut the error across Q by ~30% and at Q ≥ 30 by ~37%; smooth 300 over-smooths the middle (Q24–27 at 1.4–2.7×). The truth itself has a Q34 dip (rate below Q33 and Q35) that a shape prior can't follow (4.9× at smooth 30).
  - Head Q per-position TV reaches 0.35 even with the truth's own tokens, but only at the knots: 0.36 at cycle 5 and 0.12 at 25–27, against 0.01–0.05 elsewhere. The truth's knots sit at 1, 5.3, 28.4 and 151 (fitted on 150-bp reads) and the refit's at 1, 5, 25 and 125, so two coarse splines disagree at their kinks.
- ✓ **Finer default head Q.** Candidates fitted on 2,000 SRR5240881 pairs and scored by per-position Q TV on 1,000 held-out reads:

  | head Q (with `Mate`, `Context(1,1)`) | fit | TV median | 95th | max |
  |---|---|---|---|---|
  | `QualityMarkov(1)`, `Position(4)` (old default) | 48 s | 0.177 | 0.465 | 0.593 |
  | `QualityMarkov(1)`, `Position(12)` | 53 s | 0.106 | 0.404 | 0.570 |
  | `QualityMarkov(1)`, `Position(24)` | 69 s | 0.071 | 0.334 | 0.383 |
  | `QualityMarkov(1)`, `Position(48)` (**new default**) | 129 s | 0.060 | 0.233 | 0.322 |
  | `QualityMarkov(2)`, `Position(24)` | 110 s | 0.074 | 0.359 | 0.396 |
  | `QualityMarkov(1)`, `Position(24)`, mixed knots (log to cycle 10, then linear) | 77 s | 0.112 | 0.367 | 0.450 |
  | `QualityMarkov(1)`, `Position(48)`, mixed knots | 119 s | 0.082 | 0.314 | 0.340 |

  A second lag and linear knot spacing don't help; more log-spaced knots do. **Open:** cycles 40–41 still reach TV 0.32 at 48 knots, a narrow cycle-specific dip that a spline can't follow. It needs a per-cycle term or knots chosen by selection.
- ✓ **Q smoothing of head E** (`fit.error.fit(smooth=, window_l2=)`, tables above) is the `pe-overlap` default.
- Reported, not blocking: run ErrorProfiler on the same data and compare its substitution × Q tables.
- **Open items**, carried forward:
  - error head `Position` splines are unidentified on fixed-length reads;
  - Q38–39 are still ~5.5× too high at ~10% merging;
  - head Q misfits cycles 40–41 (TV 0.32);
  - fitting cost (bin or subsample large runs; `Position(48)` doubles head Q time).
- **Exit:**
  - ✓ on generated pairs, the substitution part of head E (context, `QualityWindow`, position, mate) and head Q are recovered within the phase 2 tolerances (`tests/test_pe_overlap.py`);
  - ✓ on one real Illumina run, a fitted spec and a report are produced (SRR5240881, above).

### Phase 5: `reference` mode and cross-mode comparison (in progress)
- ✓ `sources/bam.py` core: primary, mapped, non-duplicate records with MAPQ ≥ `--min-mapq` (default 20) become (template, read, mate) triples in read orientation and go through `generate.observations`. `bam.fit` fits both heads from those tuples, and `recovery.cigar_mode` now calls it, so the recovery harness and the BAM source share one fit path. Test: generated reads written to a BAM with their true CIGARs (forward reads soft-clipped; plus secondary, low-MAPQ and duplicate copies) reproduce the generator's tuples exactly, indels included. CLI: `python -m sequencing_error_model.sources.bam BAM REF --output SPEC`.
  - Soft-clipped bases (and an insertion after the last aligned base) give no rows but count toward `pos_start`/`pos_end` and fill Q windows at the aligned ends (`Read.clipped`); `strand` is the alignment's (`Read.strand`). The test checks both against reads clipped on either side. Hard clips aren't recoverable, and records with N/P ops are skipped.
  - ✓ MAPQ, secondary/supplementary, QC-fail and duplicate filters (in the core above).
- ✓ Site masks, contig filter and per-contig report (`bam.pileup`, `bam.masks`): a first pass over the same records and CIGAR walk counts A/C/G/T/deletion/insertion-before-site per reference position; the pileup's draws equal the tuple rows. A site is masked when a non-reference allele is seen ≥ 2 times at frequency ≥ `--mask-alt-freq` (catches reference errors as well as minor alleles), and masked sites' draws give no rows (`Read.masked`). Contigs below `--min-contig-depth` mean depth are dropped. `--contig-report TSV` writes length, mean depth, raw error rate, masked sites and kept per contig. Masks and filters are off by default until their clonal cost is measured. Test: a planted reference error is masked and its substitutions leave the table.
  - An insertion variant is masked only on the strand whose insertions attach to the masked base.
- ✓ Indel-length evidence (`generate.indel_events`): one `error`-unit count per indel event keyed by `indel` (I/D), `indel_length`, `run` (the template base's full homopolymer run, not window-limited) and `q`. Consecutive insertions before a base, or consecutive deleted bases, form one event. It reads the same (template, read, mate) triples as `observations`, so it serves the generator and BAM records alike, masks included. Test: hand-built CIGARs give the expected events.
- ✓ `IndelLength(n)` (`fit/indel.py`): P(length 1..n | kind, run capped at 8, Q) as a log-linear softmax fitted by L-BFGS from `indel_events`; lengths above n count as n.
  - With it in head E, the generator draws one head E outcome per indel event. An insertion takes its length from `IndelLength`, draws its other bases from that base's insertion mix, then makes one final no-insertion draw. A deletion removes the next length − 1 bases.
  - `observations(per_event=True)` gives the matching one row per event.
  - `bam.fit` fits it when an `IndelLength` token is given (the CLI then makes a second pass for events), and `recovery.cigar_mode` passes events through.
  - Tests: generated length distributions per kind and run group are within TV 0.05 of the truth; `recover` on a 4× indel-rate truth with longer indels in runs of 3+ passes the phase 2 tolerances, and refitted lengths are within weighted TV 0.05.
  - Known approximations: the post-insertion row is fitted as a full draw; adjacent deletions read back as one event; for deletions Q is the next read base's (the open deletion-Q convention).
- Deferred: per-read error *and quality* trajectories for `Latent(S)`. No phase fits `Latent(S)` yet, so they wait for the phase that adds its fitter.
- One code path for three kinds of reference: an external genome, a spike-in or mock community, and a self-assembly. Ship documented recipes (assemble with metaSPAdes / metaFlye / myloasm / hifiasm-meta, polish, align with minimap2 or bwa-mem2) rather than wrapping assemblers.
- ✓ **Aligner bias check** (`recovery.aligner_bias`, `python -m sequencing_error_model.recovery --aligner minibwa|minimap2`). It uses minibwa for short reads and minimap2 (`-x map-ont` by default) for long reads. The `Testing` workflow installs the pinned x64-linux release binaries (minibwa v0.7, minimap2 v2.31) and sets `REQUIRE_ALIGNERS=1`, so a missing aligner fails there; locally, `tests/test_aligners.py` skips without them (Homebrew has both at the same versions).
  - Single-end reads from one random genome, on both strands, with known CIGARs. The truth's error rate is scaled by shifting head E's non-match bias (what `generate(error_rate_scale=)` does).
  - Reported, true vs aligned: mapped fraction; op-class rates per head E row; per-site agreement (the share of true substitution, deletion and insertion counts that the alignments place at the same reference site); mean indel length per kind; and head E refitted on the aligned tuples against the truth.
  - Example spec, 50 kb genome, both at MAPQ ≥ 20 with every read mapped:

    | aligner | reads | errors per row | rate ratio (sub / del / ins) | site agreement (sub / del / ins) | refit rate ratio | phase 2 tolerances |
    |---|---|---|---|---|---|---|
    | minibwa | 3,000 × 100–150 bp | 1.0% (scale 0.1) | 0.98 / 0.88 / 0.88 | 0.96 / 0.53 / 0.60 | 0.984 | pass |
    | minimap2 `map-ont` | 400 × 1–2 kb | 6.8% | 0.95 / 1.17 / 1.06 | 0.92 / 0.46 / 0.53 | 0.975 | pass |

  - Substitutions are kept and placed well. Indels are moved: only about half land on their true site, so per-site indel evidence (`Homopolymer` and context effects on indels) is where aligner bias concentrates. minibwa loses ~12% of indels at low error; minimap2 at high error calls 17% more deletions and 6% more insertions, and makes them slightly longer (mean deletion length 1.09 vs 1.00). The likeliest cause of the site shifts is indels placed at a different base of a homopolymer run, or merged with a neighbouring substitution, but that isn't separated yet.
  - The example spec is Illumina-like with 30–40 bp position knots, so the long-read row checks the aligner rather than a realistic ONT spec. Tests (200 reads each): every read maps, and the substitution rate is within 10%.
- ✓ **Observable truth** (decided 2026-09-15). Generated CIGARs carry edits no alignment can show: a deletion and an insertion of one base in a homopolymer run cancel, and some indel pairs collapse into substitutions. So aligner bias is measured against `generate.realign`, each true read's optimal end-to-end alignment under the aligner's own affine scores (`recovery.ALIGNER_SCORES`: minibwa 2/8/12/2, minimap2 2/4/4/2), indels left-aligned and kept whole. The spec is recovered only up to that observable equivalence. `aligner_bias` reports `edits_hidden`, the aligned fit vs the observable fit (`*_observable`), and indel rates and lengths by kind and run (1, 2, 3+). Unit scores give the minimum edit distance, and a test checks that against a plain DP.
  - Rate ratios, fitted over expected on held-out reads (phase 2 tolerances: 5% marginal, 10% per Q). 50 kb genome; the `IndelLength` truth is the example spec with 4× indels and `IndelLength(4)`, longer in runs of 3+:

    | truth | aligner | reads | edits hidden | vs truth | vs observable | vs observable, `--unclip` | indels in runs of 3+ vs observable (D / I) |
    |---|---|---|---|---|---|---|---|
    | example, scale 0.1 | minibwa | 3,000 × 100–150 bp | 0.7% | 0.984 | 0.977 | 1.000 | 0.96 / 0.94; 1.00 / 1.00 (`--unclip`) |
    | `IndelLength`, scale 0.25 | minibwa | 4,000 × 100–150 bp | 3.8% | **0.926** | 0.956 | 0.998 | 0.95 / 0.91; 1.00 / 1.00 (`--unclip`) |
    | example | minimap2 `map-ont` | 400 × 1–2 kb | 2.3% | 0.975 | 0.993 | – | 0.98 / 1.00 |
    | `IndelLength` | minimap2 `map-ont` | 600 × 1–2 kb | 7.4% | **0.915** | 0.995 | – | 0.99 / 0.99 |

    Observable columns are against the forward-strand observable (below).

  - Against the observable truth on the `IndelLength` truth, with minimap2 and minibwa `--unclip`:
    - indel rates by kind and run are 0.98–1.00;
    - length TV is ≤ 0.014;
    - indel site agreement is 0.97–0.99, up from 0.73–0.78 against the read-orientation observable, whose tie-break placed reverse reads' gaps elsewhere.
  - Without `--unclip`, minibwa's indel site agreement is 0.85–0.94. Most of the site disagreement against the generator's CIGARs was equivalent placement.
  - **Free template ends** (2026-09-15). The observable truth first aligned each read end to end on its true template. That span isn't observable: an indel a few bases from a read end ties with, or loses to, a shifted end (`110M1I1M` vs `110M2D2M`). About 40% of minibwa's remaining 3+-run indel loss was this. The observable now realigns each whole read inside the genome around its true span (±16 bases) with free template ends.
  - A unit-cost observable was tried first and rejected. Substitutions and indels tie under unit costs, so its tie-break either split long indels (length TV 0.15–0.22) or inflated indels (deletion ratio 0.75).
- ✓ **Soft-clipping correction** (`bam --unclip MATCH MISMATCH OPEN EXTEND`, `recovery --aligner … --unclip`). minibwa's indel loss is clipping. On the `IndelLength` run, 14% of reads were clipped and those kept 55% of their edits, while unclipped reads matched the optimal alignment (NM 1.004×). A clipped read is realigned end to end inside the reference widened by its clips (`realign` with free template ends), unless a clipped end of 8+ bases differs from the reference at more than half its bases (adapters, chimeras). Test: errors in a clipped end, including a 3-base clip, come back as rows; a clipped adapter stays clipped.
  - The cap was 30% of any clip at first, which rejected 390 of 542 clipped reads: one error in a 2–3 base clip exceeds it. Realigned with no cap, clips holding only sequencing errors reach 0.42–0.67 divergence at the 90th percentile from 8 bases on, and random 8+ base ends 0.75–0.91 at the median. Below 8 bases the two overlap, so short clips are always realigned. At 0.5 the cap rejects 4–19% of error clips and 76–96% of random ends of 8+ bases.
  - **Open:** whether `--unclip` should default on for short reads (check on a real run with adapters first; short adapter remnants under 8 bases are now realigned); the unclip DP and the free-ends observable are unbanded, so long reads need a band (the minimap2 rows take 2–3 minutes).
- ✓ **Indel components** (`recovery.indel_components`, in `aligner_bias` as `indel_components*`). Both specs' head E on held-out tuples: expected deletions and insertions, fitted over reference, grouped by the centre base's homopolymer run (`Homopolymer`) and by the base at each context offset (`Context`). On random templates the other components average out of each group. The phase 2 per-Q tolerance (10%) applies. Test: removing the truth's run effect on indels moves only runs of 3+, and a deletion-after-A context weight moves only deletions.
  - Aligned fit over observable fit, same runs as the table above:

    | truth | aligner | D by run 1 / 2 / 3+ | I by run 1 / 2 / 3+ | context groups off by > 10% |
    |---|---|---|---|---|
    | example, scale 0.1 | minibwa `--unclip` | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 | 0 of 40 |
    | `IndelLength`, scale 0.25 | minibwa `--unclip` | 1.00 / 0.99 / 1.00 | 0.99 / 1.00 / 1.00 | 0 of 40 |
    | `IndelLength`, scale 0.25 | minibwa | 0.92 / 0.95 / 0.97 | **0.88** / 0.91 / 0.94 | 3 of 40 (0.89–0.90) |
    | example | minimap2 `map-ont` | 0.98 / 0.95 / 0.98 | 0.98 / 0.98 / 0.99 | 0 of 40 |
    | `IndelLength` | minimap2 `map-ont` | 1.00 / 0.98 / 0.99 | 0.99 / 1.00 / 0.99 | 0 of 40 |
    | `IndelLength`, `Context` zeroed, scale 0.25 | minibwa `--unclip` | 1.00 / 0.99 / 1.00 | 0.99 / 1.00 / 1.00 | 0 of 40 |
    | `IndelLength`, `Context` zeroed | minimap2 `map-ont` | 0.99 / 0.99 / 0.99 | 1.00 / 1.00 / 0.99 | 0 of 40 |

  - `Homopolymer` and `Context` on indels are met for minimap2 and for minibwa with `--unclip`: every group within 5%, most within 1%. Without `--unclip`, minibwa loses about 10% of insertions in every group, which is clipping (above), not a component shift.
- ✓ **Forward-strand observable** (`recovery.realign_observable`). The first indel components run failed `Context` on deletions with both aligners: deletions after G 1.14–1.41, before G 0.81–0.89. It persisted with the truth's G substitution hotspot zeroed, since G stays an error hotspot through head Q's low Q on G.
  - Scoring the aligned CIGARs under the aligner's scores showed ties. On the `IndelLength` truth with `Context` zeroed, minibwa aligned 2,791 reads as the observable did, 748 at an equal score but a different placement, and 7 worse; the whole G shift was in the tied reads.
  - 745 of the 748 tied reads were on the reverse strand (minimap2: 131 reverse, 40 forward, split evenly by direction). The aligners left-align gaps on the reference, but `realign` ran in read orientation, which right-aligns reverse reads on the reference.
  - The observable now realigns on the forward strand and flips back; a leading insertion there becomes a clipped read end, as in `sources.bam`. Tests: a reverse read's deletion in a GGG run sits at the run's first base on the reference, and a leading forward insertion comes back clipped.
- ✓ **Clonal cost of masks** (`recovery.mask_cost`, `python -m sequencing_error_model.recovery --mask-cost`). Reads come from one random genome, placed on both strands, so every mismatch is an error. Masks are computed with `bam.count_alleles` / `bam.site_masks` / `bam.apply_masks`, the BAM path's cores without file I/O. For each `max_alt_freq` it reports masked sites, removed rows and errors by op class, the most enriched trinucleotide contexts among removed errors, and a head E refit on the kept tuples against the truth on held-out reads: phase 2 rate tolerances, op TV, and the slope of fitted on true `Context` log-odds. The contig depth filter isn't measured: it selects whole contigs on depth, not on the outcome.
  - Example spec (14% row error rate, so a harsh case), 20 kb genome, depth 30, 16,901 reads (72 s):

    | `max_alt_freq` | masked sites | errors removed | removed rows' error rate | rate ratio | `Context` slope | phase 2 tolerances |
    |---|---|---|---|---|---|---|
    | off | 0 | 0 | – | 1.004 | 1.008 | pass |
    | 0.1 | 5,298 | 43.6% | 24% | 0.888 | 0.867 | **fail** (marginal, every Q bin 0.85–0.92) |
    | 0.2 | 220 | 2.5% | 36% | 0.998 | 0.993 | pass |
    | 0.5 | 1 | ~0 | – | 1.004 | 1.008 | pass |

  - The removed errors concentrate where the truth's substitution hotspot is (after G; `GGC`/`GCC` enriched 2× at 0.1 and 5.5× at 0.2), which is the outcome-dependent loss §6.6 predicts. At 0.2 it is too small to move a component past tolerance at this depth and error rate. Chance masking grows with the error rate and shrinks with depth.
  - **Default:** masks stay off. A threshold only makes sense with the run's error rate and depth in view; the report tells a user what a given threshold costs on a clonal simulation of their spec. Test: on a small genome, a 0.2 mask removes rows more than twice as error-rich as the rest.
- ✓ `compare.py` (§7): evidence-level and model-level comparison on shared support. First use: `pe-overlap` vs `reference` on the same Illumina reads, which also estimates the PCR/library error contribution.
  - `compare.evidence(a, b, by)`: two labelled head E tables restricted to substitution support (match and substitution rows, since `pe-overlap` identifies no indels), marginalised onto `by` and kept at values both observe with at least `min_exposure` rows each. Rates are standardised to the pooled exposure, so coverage differences (overlap positions, Q mix) don't read as rate differences. `excess` = rate_b − rate_a.
  - `compare.models(a, b, table)`: both specs' head E on one table, each conditioned on no indel: op TV, rates, excess, rate by a covariate, and log-likelihood per row.
  - `pe_overlap.table(ev, error_head)` exposes the overlap evidence: disagreements split by the posterior under a fitted head E (the EM step, now shared with `pe_overlap.fit`), or a coin flip without one.
  - Test: 3,000 generated pairs (40 bp, insert 45 ± 8, ~3.5% sequencing errors) from fragments carrying 1% library substitutions that both mates read. Those cancel in the overlap and show against the unmutated genome.
    - Evidence level, by (Q, mate): overlap 4.24% vs reference 5.21%, **excess 0.97%** (truth 1%; tolerance 20%).
    - Model level, on the reference table: excess 1.11% (tolerance 30%), op TV 0.014, and the reference-fitted head E has the higher log-likelihood.
  - ✓ CLI: `python -m sequencing_error_model.compare R1 R2 BAM REF --overlap-spec A --reference-spec B --output report.json` takes the specs the two source CLIs fitted from one run and writes both levels (`evidence` by `--by`, default Q and mate; `models` on each mode's rows), with optional `--mask-alt-freq`. Both tables share one context flank and Q window, wide enough for both specs' head E. Test: 300 generated pairs (50 bp, insert 70 ± 4, no library errors) written as FASTQ and a BAM go through both source CLIs and `compare`; the evidence excess is within 10% of the reference rate (observed 0.02% on 4.0%).
  - **Open:** per-component comparison (context log-odds, position curve) beyond rates; expected-gap declarations for other mode pairs (`kmer`, phase 6); standardise `evidence` by a read-position bin as well (below); `--unclip` in `compare` (its reference rows use the aligner's clips).
- ✓ **Real near-clonal run** (2026-09-16): SRR24523812, *Phocaeicola vulgatus* NMBE-5 isolate (SAMN35060538), MiSeq 2×150, reads adapter-trimmed to 1–151 bp, Q binned to 14/21/27/32/36. No assembly is deposited, so the reference is a shovill (SPAdes) assembly of all 8.0 M pairs: 184 contigs, 4.81 Mb, N50 119 kb. The comparison uses 100,000 pairs from mid-run (pairs 1,000,001–1,100,000), aligned with minibwa: 99.8% of reads mapped, 1.7% soft-clipped, NM 0.76% per aligned base, median insert 365.
  - Real reads carry N calls, which the generator never emits. `generate.observations` now skips a draw whose read base is not A/C/G/T, as it skips masked sites; head E has no category for a no-call.
  - Per-contig raw error rates are 0.73–1.1%, the high-depth contigs (up to 107×, repeats) included, so there is no divergent-repeat signal.
  - `pe-overlap` places 8.5% of pairs (inserts ≤ ~280): 8,537 overlaps, 4,562 disagreements, 71 indel-flagged. Fitting: 47 min for 20,000 pairs, 104 min for 100,000 (head Q on every base dominates). `reference` (`--unclip`) on 40,000 reads, 6.0 M rows: 95 min. Both ran concurrently with other fits.
  - **`QualityWindow` vs centre Q** (exit item), head E fitted on pairs 1–20,000 and scored on held-out pairs 80,001–100,000:

    | mode | held-out rows | ll/row, `QualityWindow(0)` | (1) | (2) | gain over (0), nats: (1) / (2) |
    |---|---|---|---|---|---|
    | `reference` (`Context(1,1) Homopolymer Mate`) | 5.86 M | −0.031787 | −0.031449 | −0.031246 | +1,977 / +3,163 |
    | `pe-overlap` (`Context(1,1) Position(4) Mate`, EM) | 162 k | −0.034692 | −0.034505 | −0.034375 | +30 / +51 |

    `pe-overlap` rows are scored with disagreements split by a coin flip, the same table for every model. Each window step adds ~140 weights, so the `reference` gain is not overfitting; the `pe-overlap` gain is small on ~1,700 held-out overlaps. `QualityWindow(2)` beats the default (1) in both modes; selection over m should decide the default.
  - **Evidence level** (`compare`, by Q and mate). `reference` shows 1.31× the overlap rate (0.757% vs 0.575%), and the ratio is Q-dependent (Q21 3.1–3.5×, Q36 3.5–4.3×), so it is not a PCR excess, which adds a constant rate. Restricting the aligned reads to proper pairs with insert ≤ 280 (16,605 of 200,049 records), the reads `pe-overlap` can see, brings it to 1.05×:

    | Q (mate 1 / 2) | overlap | reference, all | reference, insert ≤ 280 | ratio, insert ≤ 280 |
    |---|---|---|---|---|
    | 14 | 5.02% / 5.69% | 6.00% / 7.21% | 4.99% / 5.64% | 0.99 / 0.99 |
    | 21 | 0.361% / 0.355% | 1.10% / 1.25% | 0.804% / 0.692% | 2.23 / 1.95 |
    | 27 | 0.174% / 0.233% | 0.311% / 0.430% | 0.202% / 0.256% | 1.16 / 1.10 |
    | 32 | 0.053% / 0.043% | 0.099% / 0.125% | 0.069% / 0.069% | 1.29 / 1.61 |
    | 36 | 0.0087% / 0.0081% | 0.031% / 0.035% | 0.025% / 0.031% | 2.82 / 3.80 |

  - Each remaining disagreement has a cause:
    1. **Insert-size selection** (1.31× → 1.05×). Short-insert clusters have fewer errors at a given reported Q in this run. A `pe-overlap` spec from a long-insert library describes its short-insert subset.
    2. **Cycle selection within a Q bin** (Q21). On short inserts, Q21's mismatch rate by 25-cycle bin is 1.78%, 1.81%, 0.80%, 0.70%, 0.44%, 0.25%; overlaps cover mostly the late cycles (0.25%, near the overlap's 0.36%). No Q21 mismatch recurs at a site (0 of 286). Q36 is flat across cycles (0.023–0.034%).
    3. **A flat ~0.02 pp excess at Q27–36**: library/PCR errors (they cancel in the overlap) and reference differences. 16% of Q36 mismatches on short inserts (10% on all) recur at the same site and allele at ~6× depth, which is assembly consensus error or minor variants (~0.005 pp of the 0.019 pp); the rest are singletons.
  - **Model level**, both specs' head E conditioned on no indel, on each mode's rows (`reference` spec from 40,000 reads of all inserts):

    | `pe-overlap` spec | rows | op TV | rate ratio, reference / overlap spec | ll/row, overlap / reference spec |
    |---|---|---|---|---|
    | 20,000 pairs | reference, all | 0.0039 | 1.46 | −0.03726 / **−0.03510** |
    | 20,000 pairs | overlap | 0.0031 | 1.35 | **−0.02506** / −0.02614 |
    | 100,000 pairs | reference, all | 0.0035 | 1.62 | −0.03652 / **−0.03510** |
    | 100,000 pairs | reference, insert ≤ 280 | 0.0024 | 1.66 | **−0.02041** / −0.02057 |
    | 100,000 pairs | overlap | 0.0033 | 1.50 | **−0.02473** / −0.02607 |

    Each spec wins on its own mode's rows, and the 100,000-pair overlap spec also wins on the short-insert reference rows, which the reference spec (fitted on all inserts) doesn't specialise to: the selection effect again. Rate ratios at model level (1.35–1.66) exceed the evidence-level 1.31, and the two overlap specs differ (0.55% vs 0.61% on overlap rows). Likely cause, not yet checked: on reference rows the overlap spec's `Position(4)` extrapolates to early cycles it rarely saw, where Q21 is 7× more error-prone.
- **Baseline harness.** The phase 3 metrics computed on real reads vs reads from (a) the native generator, (b) ReSeq for Illumina, (c) Badread, PBSIM3 and CycSim for long reads, each simulator trained by its own profiler on the same alignment. Add k-mer spectrum concordance and context-dependent substitution rates, the metrics of the 2026 ONT benchmark.
  - ✓ `baseline.py` (`python -m sequencing_error_model.baseline REF REAL_BAM NAME=BAM … --output report.json`). Every read set is aligned the same way and profiled by the same code.
    - Error rows come from `bam.records` → `generate.observations` (by reported Q, cycle bin and trinucleotide, split into match, substitution, deletion and insertion), in batches of 200 records, with each aligned read's error rate (0.2% bins, capped at 10%) counted alongside. `--cycle-bin` sets the bin (10 for short reads, 500 for long reads).
    - Q comes from every base of every primary record (Q histogram per cycle bin, lag-1 pairs), and canonical 21-mers from the same records (multiplicity spectrum, and occurrences absent from the reference).
    - `distances` gives one number per metric, 0 when identical. Rate metrics are the real-exposure-weighted mean \|log rate ratio\|; the others are TVs, except `q_lag1` (absolute difference) and `kmer_absent` (\|log ratio\|).
    - `verdict` gives the native generator's result against each baseline: beats, matches (within 10% or 0.01), or trails.
    - Test: reads from the example spec rank a same-spec read set above a 3× error-scaled one on every rate, per-read error rate and k-mer-absent metric, and canonical k-mers ignore strand and skip N windows.
  - ✓ **Illumina run** (2026-09-16), SRR24523812 (as above).
    - Training: pairs 1–20,000 aligned with minibwa. The native spec is the `reference` spec above (`QualityWindow(1) Context(1,1) Homopolymer Mate`, fitted with `--unclip`). ReSeq 1.1 (bioconda container, amd64 under emulation, 1 min) got `--adapterFile` with the TruSeq pair, because its adapter auto-detection fails on adapter-trimmed reads. Its coverage bias fit didn't converge at 1.2× depth, so it simulated uniform coverage.
    - Evaluation: held-out pairs 20,001–100,000 against 80,000 simulated pairs each (native `--insert-mean 369 --insert-sd 66`, the training alignment's), all aligned with minibwa and profiled with `--unclip 2 8 12 2`, 2 min per run.

    | metric | real | native | `Position(4)` | `Position(8)` | ReSeq | distance: native / `Position(4)` / `Position(8)` / ReSeq | `Position(8)` vs ReSeq |
    |---|---|---|---|---|---|---|---|
    | error rate | 0.776% | 0.832% | 0.819% | 0.798% | 0.849% | 0.070 / 0.054 / **0.028** / 0.089 | beats |
    | deletion / insertion bases per row (×10⁻⁵) | 2.6 / 3.2 | 3.0 / 3.2 | 3.0 / 3.3 | 2.9 / 3.6 | 3.2 / 2.3 | D 0.137 / 0.146 / 0.093 / 0.206; I 0.006 / 0.060 / 0.143 / 0.308 | beats |
    | rate by reported Q | | | | | | 0.140 / 0.141 / **0.108** / 0.243 | beats |
    | rate by 10-cycle bin | | | | | | 0.122 / 0.074 / **0.056** / 0.094 | beats |
    | substitution rate by trinucleotide | | | | | | 0.180 / 0.176 / 0.168 / 0.245 | beats |
    | Q by 10-cycle bin (TV) | mean Q 32.36 | 32.30 | 32.30 | 32.30 | 32.13 | 0.010 / 0.010 / 0.010 / 0.018 | matches |
    | Q lag-1 autocorrelation | 0.200 | 0.175 | 0.175 | 0.175 | 0.188 | 0.026 / 0.026 / 0.026 / 0.012 | trails |
    | per-read error rate (TV) | median < 0.2%, p90 2.7% | 0.7%, 1.9% | 0.7%, 1.9% | 0.7%, 1.9% | < 0.2%, 2.7% | 0.373 / 0.368 / 0.360 / 0.076 | trails |
    | k-mer spectrum (TV) | | | | | | 0.057 / 0.054 / 0.051 / 0.045 | matches (within the 0.01 floor) |
    | reference-absent k-mers | 11.5% | 15.4% | 15.0% | 14.6% | 13.1% | 0.288 / 0.263 / 0.241 / 0.132 | trails |

    - **Beats** on per-base error structure: rates by Q and trinucleotide context, and indel rates.
    - **Rate by cycle** trailed because the first spec had no `Position` in head E (the Q21 cycle effect above). Refitting head E as `QualityWindow(1) Context(1,1) Position(4) Homopolymer Mate` on the same training records (14 min; +0.00022 nats per row, +1,290 on 5.98 M training rows), with head Q unchanged and generation repeated with the same seed, cuts that distance from 0.122 to 0.074. It now beats ReSeq (0.094).
    - `Position(8)` (19 min; +0.00035 nats per training row) goes further: rate by cycle 0.056, error rate 0.028, and rate by reported Q 0.108 (from 0.140), so part of the Q curve's miss was cycle confounded with Q.
      - Held-out head E log-likelihood on pairs 20,001–50,000 (9.0 M rows, same `--unclip` and window) agrees: +1,568 nats for `Position(4)` and +2,538 for `Position(8)` over no `Position`. The extra knots generalise.
      - The insertion distance grows (0.006 → 0.060 → 0.143), but it rests on about 750 insertion bases per read set, so without a noise floor it reads as noise; deletions improve (0.137 → 0.093).
      - Per-read heterogeneity doesn't move (0.373 → 0.360).
      - Default for Illumina `reference` specs: `Position(8)` over `Position(4)`. The `pe-overlap` default stays `Position(4)`: fixed-length amplicon reads leave its position axes collinear (phase 4).
    - **Trails** on per-read error heterogeneity, the largest gap. Real mismatches per read (NM) are overdispersed (variance/mean 4.0; 59% of reads have none, against 33% under Poisson). ReSeq reaches 3.4 through its errors-so-far conditioning. The native generator draws errors independently given Q, so it gives 1.05, and its lag-1 Q (0.175 vs 0.200) doesn't carry enough read-level clustering. The spread-out errors break more k-mers, which drives the k-mer spectrum and reference-absent k-mer gaps. That is the deferred per-read trajectory model (`Latent(S)`).
    - Not yet measured: a sampling-noise floor (real vs real halves), which matters for the rare indels (600–750 indel bases per read set).
  - ✓ **Long-read run** (2026-09-16): SRR30993108, *E. coli* K-12, PromethION R10.4.1 (kit14, 5 kHz, dorado sup v5), 62,469 reads, mean 7.9 kb. Reference NC_000913.3 (MG1655). Reads 1–5,000 train and reads 30,001–35,000 (39.7 Mb) are held out, all aligned with minimap2 `-ax map-ont`, profiled without `--unclip` (its DP is unbanded) and with `--cycle-bin 500`, 5.5 min for five read sets.
    - **Badread 0.4.2:** `error_model` and `qscore_model` on the training PAF; `simulate` with the training reads' length (8,356 ± 10,166) and primary-alignment identity (97.8, max 99.9, sd 4.2), and default adapters, junk, random reads, chimeras and glitches.
    - **PBSIM3 3.0.5:** `--method sample --sample train.fq` (lengths and Q strings from the training reads), `--difference-ratio 39:24:36` (its recommendation for ONT).
    - **CycSim 1.0.4:** built from source (not on bioconda; needs `CMAKE_POLICY_VERSION_MINIMUM=3.5` with current CMake). `train -r nanopore` on the training alignment, `sim -D` to the held-out base count. Its unaligned BAM output was converted to FASTQ and aligned like the others.
    - **Native:** head E `QualityWindow(1) Context(1,1) Homopolymer IndelLength(8)` fitted on training reads 1–300 (2.4 M rows, 3 min); head Q `QualityMarkov(1) Position(48) Context(1,1)` fitted on reads 1–100 (27 min). Head Q is the limit: long-read rows never share `pos_start`/`pos_end`, so nothing collapses, and a 49-value Q alphabet costs about 6 s per 1,000 rows. Templates have lengths resampled from the training reads, placed uniformly on both strands; generation took 6.5 min.

    | metric | real | native | Badread | PBSIM3 | CycSim | native vs Badread / PBSIM3 / CycSim |
    |---|---|---|---|---|---|---|
    | error rate | 1.97% | 1.44% (0.317) | 2.43% (0.210) | 1.58% (0.219) | 2.36% (0.179) | trails / trails / trails |
    | substitution / deletion / insertion | 0.83 / 0.57 / 0.57% | 0.72 / 0.44 / 0.27% | 0.97 / 0.83 / 0.63% | 0.74 / 0.53 / 0.31% | 0.89 / 0.88 / 0.59% | insertions trail all three (0.732 vs 0.107 / 0.604 / 0.044) |
    | rate by reported Q | | 0.454 | 2.126 | 0.532 | 0.856 | beats / beats / beats |
    | rate by 500-cycle bin | | 0.330 | 0.229 | 0.308 | 0.195 | trails / matches / trails |
    | substitution rate by trinucleotide | | 0.331 | 0.172 | 0.613 | 0.471 | trails / beats / beats |
    | Q by cycle bin (TV) | mean Q 34.2 | 34.7 (0.066) | 34.3 (0.036) | 34.7 (0.034) | 20.7 (0.895) | trails / trails / beats |
    | Q lag-1 autocorrelation | 0.947 | 0.932 (0.015) | 0.337 (0.610) | 0.947 (0.000) | 0.391 (0.556) | beats / trails / beats |
    | per-read error rate (TV) | median 0.5%, p90 5.3% | 1.5%, 2.1% (0.739) | 0.9%, 6.5% (0.221) | 0.5%, 3.9% (0.091) | 1.9%, 4.9% (0.742) | trails / trails / matches |
    | k-mer spectrum (TV) | | 0.441 | 0.329 | 0.428 | 0.248 | trails / matches / trails |
    | reference-absent k-mers | 14.2% | 13.2% (0.075) | 21.9% (0.431) | 13.9% (0.024) | 30.1% (0.751) | beats / trails / beats |

    Distances in brackets. From NM, held-out real reads have 1.98% errors, all training reads 2.20%, and training reads 1–300 2.45%.
    - **Beats** every baseline on rate by reported Q. Badread's Q barely tracks its errors, and CycSim's Q is off altogether (mean 20.7). PBSIM3 copies real Q strings, so it matches Q almost exactly and reaches its errors through Q, but with its fixed error ratio.
    - **Trails** on the error rate (1.44%, against 2.45% in its own training reads), insertions most of all. This is not the error head: on the training reads head E's expected rates equal the observed ones exactly (substitution 1.125%, deletion 0.481%, insertion 0.337% of per-event rows). On the same templates, generated reads give substitution 0.745%, deletion 0.317% and insertion 0.227%, and head E is again self-consistent on those rows. So the generated Q tracks carry the loss. Real reads have a heavier low-Q tail (Q ≤ 10 on 7.0% of bases, against 4.9% generated) and much wider per-read quality: per-read mean Q spans 26.7–40.8 (p10–p90), native only 33.3–36.0. Head Q's `QualityMarkov(1)` gets lag-1 right (0.932 vs 0.947) but has no per-read level, so whole low-quality reads, where errors concentrate, are never generated. It is the Illumina per-read heterogeneity gap again, this time on the Q side, and `Latent(S)` should close both.
    - The per-read error rate shows it directly: native 1.5–2.1% (p50–p90) against real 0.5–5.3%. PBSIM3 inherits per-read quality from the real Q strings and gets closest (TV 0.09).
  - **Open:**
    - `Latent(S)` (per-read quality and error level), then re-run both rows;
    - head Q fitting cost on long reads (bin `pos_start`/`pos_end` geometrically, or fit on a Q-binned alphabet);
    - a sampling-noise floor.
- **Exit:**
  - on generated reads aligned back to their genome, the full spec (including indel lengths and homopolymer effects) is recovered within the phase 2 tolerances, after correcting for the measured aligner bias. Redefined 2026-09-15 against the observable truth (aligner bias check above). Rates, indel lengths and indel rates by homopolymer run (within 10%) are met for minimap2, and for minibwa with `--unclip`, and so are `Homopolymer` and `Context` on indels (indel components above);
  - the clonal cost of each mask is reported, and masks that bias a head E component beyond the phase 2 tolerances on clonal reads are off by default;
  - ✓ on one real near-clonal Illumina dataset (isolate, spike-in or amplicon mock), `pe-overlap` and `reference` specs agree per component on shared support within tolerance, or the report explains each disagreement (SRR24523812 above: rates are explained, by insert-size and cycle selection plus a flat library/reference excess; per-component comparison beyond rates is still open);
  - ✓ `QualityWindow` beats the centre-Q-only head E on held-out likelihood, or the report shows it doesn't (SRR24523812: it does in both modes, and m=2 beats m=1);
  - ✓ the baseline harness runs on at least one Illumina and one long-read dataset, and the report states where the native generator beats, matches or trails each baseline (SRR24523812 vs ReSeq; SRR30993108 vs Badread, PBSIM3 and CycSim, above). The native generator beats every baseline on rate by reported Q; it trails on per-read heterogeneity on both platforms, and on the ONT error rate through it.

### Phase 6: `kmer` mode (default build), checked against `pe-overlap` and `reference`
- Head E fitters for unmodified skiver:
  - latent-position `Context(L,R)` from `kvmer.csv` with the true-base mask;
  - centre-Q term by marginal matching against `summary_phred.csv` under the FASTQ exposure;
  - position curve, strand, GC, Weibull passthrough.

  Head Q uses the phase 2 FASTQ fitters.
- **Evidence check.** On the same reads, tabulate from `pe-overlap` and `reference` tuples the marginals skiver reports: spectrum in trinucleotide context, P(error | Q), read-position curve, GC, and hazard/clustering. Use skiver's counting conventions where they can be emulated (§5.6; e.g. first-error stopping for P(error | Q)), otherwise compare rates. Compare with skiver's CSVs through `compare.py`. This tests skiver's evidence independently of the `kmer` fitters.
- **Model check.** Compare the `kmer` default spec with the `pe-overlap` and `reference` specs per component (context log-odds, centre-Q curve, position, op composition), on held-out likelihood of their tuples, and on phase 3 metrics of generated reads.
- **Synthetic recovery.** Generate, run skiver (released binary), profile the FASTQ, fit, compare with the spec. On clonal reads, run with skiver's outlier filter on (default) and off (`--use-all`), and report the keys, contexts and error mass the filter removes when there is no variation to remove.
- **Exit:**
  - the synthetic loop passes at v=13 within the phase 2 tolerances, and the documented failure at small v is reproduced as a guarded error or warning;
  - the outlier filter's clonal cost is quantified per head E component;
  - on the phase 5 near-clonal real datasets, skiver's marginals and the `kmer` spec are compared with both other modes on shared support, and the report states agreement per component, with expected gaps (PCR errors vs `pe-overlap`) given as estimates;
  - reported, not blocking: the default-mode model's held-out likelihood and marginals vs ReSeq (Illumina) and Badread (long reads) profiles trained on the aligned reads (§4.4, risk 2).

### Phase 7: biological variation simulator and the size of the problem
Test infrastructure only (§1 scope): it exists so separation methods have per-base truth.
- `variation.py`: from a base genome, make strain haplotypes over a star or random tree at a target divergence (ANI), with SNPs (transition/transversion ratio) and short indels; add within-population minor alleles from a frequency spectrum; insert divergent repeat copies (the ERR10889147 failure). Output haplotype FASTAs, abundances per sample, and a truth-site table (consensus position, alleles, allele per haplotype, frequency per sample).
  - No recombination at first. It weakens linkage, so it is added as a scenario once phase 8's linkage test needs a hard case.
- External haplotypes: accept a reference plus VCF (msprime, SimBac) or haplotype FASTAs with a truth table (CAMISIM's sgEvolver strains), so richer population structure plugs in without a dependency.
- `generate.fragments` samples haplotypes by abundance; read names carry the haplotype, and every base's truth is one of match, error or variant (the read's CIGAR against its haplotype, the haplotype against the consensus).
- **Scenario grid:** divergence (ANI 99.99, 99.9, 99, 98, 95%) × minor-strain fraction (0.5, 0.2, 0.05, 0.01) × per-genome coverage (10–200×) × minor-allele density × divergent repeats, single sample and a multi-sample abundance series; Illumina 2×150 with short and long inserts, and ONT/HiFi.
- **Problem size.** Run each mode unmodified over the grid (`pe-overlap`; `reference` against an external relative, the majority-strain consensus and a self-assembly recipe; `kmer` default with its outlier filter) and report bias against the clonal truth per head E component: marginal rate, op composition, context log-odds, rate per reported Q, position. This is the baseline phase 8 must beat.
- **Exit:**
  - per-base truth round-trips: every read mismatch against the consensus is labelled error or variant, consistently with the CIGAR and haplotype;
  - the grid's clonal point reproduces the phase 4–6 recovery results;
  - the bias table exists for all three modes, and `pe-overlap` is shown variant-immune outside the repeat scenarios.

### Phase 8: separating biological variation from sequencing error
- `sites.py` (§6.6):
  - **conservative mask**: coverage, single-copy coverage, no linked variants, Q-free; reports the composition of what it drops;
  - **joint latent-site model** over `reference` tuples: site allele counts under head E's expected counts, a linkage term from read and pair co-occurrence, an optional multi-sample term; alternates with head E fitting and warm starts, as `pe-overlap`'s EM does;
  - **`kmer` default:** a key-level test of each key's op counts against head E's expected counts at the key's coverage, beside or instead of skiver's outlier filter; enhanced mode adds linkage through read ids.
- **Q diagnostic:** mismatch rate by reported Q at called variant sites vs retained sites. The miscalibrated-Q guard (§10) is rerun with variation present.
- **Cross-mode check:** where inserts overlap, `pe-overlap` on the same reads bounds the residual variation in the other modes.
- **Real strain-resolved data**, in order: skiver's K-12/O157:H7 mixtures; ZymoBIOMICS D6331 (21 strains including five E. coli at equal abundance; strain genomes supplied); then one real metagenome.
- **Exit:**
  - on the phase 7 grid, head E is recovered within the phase 2 tolerances against the clonal truth in every scenario above a stated floor of minor-strain fraction, divergence and coverage; the floor is set from the phase 7 bias table before the methods are tuned, and failures below it are reported;
  - no method moves the clonal point outside the phase 2 tolerances;
  - with miscalibrated Q and variants, fitted rates follow the injected truth, not the reported Q;
  - on D6331, `reference` and `kmer` specs agree with `pe-overlap` or a clonal spike-in on shared support within tolerance, or the report explains each gap; variant-site precision and recall against the known strain differences are reported, not blocking.

### Phase 9: exporters, wave 1 (ART/art_modern, InSilicoSeq, Badread, PBSIM3)
- Each exporter:
  - implements `--q-policy`, where the simulator couples errors to Q;
  - supports a base-profile fallback;
  - emits a fidelity report.
- Badread gets both the error model and the Q model. PBSIM3 gets QSHMM (Q-bearing) and ERRHMM.
- **Round-trip tests.** Export, run the real simulator, re-estimate, and compare both the error rate and the Q marginals to the spec within tolerance. Re-estimation uses the `reference` mode against the simulator's source genome (cheap and exact), plus the `kmer` mode where the default path must be exercised. Simulators come from bioconda in a separate `export-roundtrip` workflow (weekly + `workflow_dispatch` + PRs touching `export/`), so the core Testing workflow stays fast.
- Unit tests validate generated files against each simulator's own loader where one is importable (InSilicoSeq `KDErrorModel`, Badread model parser). Confirm how InSilicoSeq couples substitutions to Q.
- **Exit:** all four simulators run on exported models; under the chosen `--q-policy`, the preserved quantity is within 10% of the spec, and the violated one is quantified in the fidelity report.

### Phase 10: exporters, wave 2, and wrapper recipes
- NEAT v4 (quality Markov + error model; optional dependency), Mason2 (including correct/wrong-base Q parameters) and wgsim parameter sets, and NanoSim error and quality model parts on a base model.
- Candidates, pending a format check: **ReSeq** (its Q process and Q-conditioned errors map closely onto both heads, a better Illumina target than ART if its stats file can be written) and **CycSim** (k-mer errors and error-state transitions for long reads). Drop either if its model file is not documented or stable.
- Documented CAMISIM (`art`/`nanosim3`/`wgsim` config) and MeSS recipes pointing at exported profiles.
- **Exit:** the same round-trip criterion (NanoSim: best-effort, documented gaps).

### Phase 11: further evidence sources and joint fitting
- `sources/ont_duplex.py` (dorado duplex pairs, simplex Q tracks).
- UMI/duplex consensus (fgbio) where libraries have UMIs.
- `sources/importers.py` for GATK BQSR recalibration tables, DADA2 `learnErrors` matrices and InterOp error metrics.
- Multi-source joint fit with per-source disagreement diagnostics, reusing `compare.py`.
- **Exit:**
  - each source recovers a known model on synthetic data;
  - on at least one real ONT metagenome, with phase 8 separation applied, the report quantifies agreement between `reference`, `kmer` and duplex per component.

### Phase 12: enhanced `kmer` mode
- `sources/skiver_dump.py`: a streaming aggregator from TSV and `windows.bin` into the observation schema (no 32 GB in memory). It reads `phred` from base observations and `qual_str` windows from raw observations.
- Port the fork's composable components under the §5.4 names, with the generative/training-only split now resolved by head Q.
- **New minimal fork** of GZHoffie/skiver (§12, decision 5): extend with distinct module(s), each output gated by a dump format version header:
  - key-side qualities, so quality windows span t < 1;
  - a per-window quality summary in `windows.bin`, for a `Latent(S)` shared by both heads;
  - homopolymer run length and indel length at edits;
  - 2+-edit values aligned to consensus;
  - R1/R2 flag;
  - a binary/Parquet `--base` to replace the 32 GB TSV.
- Offer the non-invasive ones upstream so more of the default mode improves over time.
- **Exit:** enhanced models beat default models on held-out likelihood for both heads on real ONT and Illumina datasets, and agree more closely with `pe-overlap` and `reference` on shared support; PBSIM3/Badread exports gain the latent-state and error-conditioned Q models.

### Phase 13: reports, packaging, release
- `sem report`: a single HTML page (plotly loaded once, no MathJax) with:
  - coverage/filter diagnostics;
  - both heads' fitted effects;
  - calibration (empirical vs reported Q, by context and position);
  - identifiability flags;
  - variation diagnostics (mask composition, variant-site calls, rate by Q at variant vs retained sites);
  - cross-mode comparison;
  - per-export fidelity reports;
  - selection traces.
- PyPI and bioconda recipes, plus a container image for Nextflow (nf-core-style module for the synthetic-metagenomic-benchmark pipeline).
- **Exit:** v0.1.0 tagged via the Release workflow; the bioconda PR is open.

---

## 10. Testing strategy

| Layer | Where | Runs |
|---|---|---|
| Unit + property tests: parsers, true-base mask, spec I/O, quality alphabet legality, §5.3 alignment/quality consistency (every generated read has `len(qual) == len(seq)`, deletions consume no quality, inserted bases get qualities) | `tests/` | every PR (Testing) |
| Golden fixtures from pinned skiver releases | `tests/fixtures/skiver-vX.Y.Z/` | every PR; regenerated by `skiver-compat` |
| Small statistical recovery, both heads (fixed seeds, tolerance bands) | `tests/recovery/` | every PR, under ~2 min |
| "Q is never a label" guard: a synthetic dataset with deliberately miscalibrated qualities, where the fitted error rate must follow the injected truth, not 10^(−Q/10) | `tests/recovery/` | every PR |
| Cross-mode agreement: the same generated paired reads (with their genome) through `pe-overlap`, `reference` and `kmer`; fitted specs agree within tolerance on shared support | `tests/recovery/` | `pe-overlap` + `reference` every PR; with skiver in `skiver-compat` |
| Clonal cost of filters and masks: on single-genome reads, each mask and filter (minor-allele masking, skiver's outlier filter, `sites.py`) keeps head E within the phase 2 tolerances | `tests/recovery/` | every PR |
| Variation guard: one small strain mixture (one divergence, one minor fraction) through `reference` and `kmer` with separation recovers the clonal head E within tolerance, `pe-overlap` is unchanged, and the miscalibrated-Q guard still holds; the full phase 7 grid runs manually | `tests/recovery/`; manual workflow | every PR (small); dispatch (grid) |
| Exporter round-trips with real simulators (errors + qualities, per `--q-policy`) | `export-roundtrip` workflow | weekly, dispatch, PRs touching `export/` |
| Full-scale recovery on real datasets | manual workflow or HPC | before releases |

---

## 11. CI/CD and PR policy (in place)

The standard matches `EBI-Metagenomics/mimicc-ena-submission-assistant`: uv + Python 3.12, ruff format/check, mypy, pytest, and pre-commit with the same pinned revs. That reference repo has no branch protection configured, so this repo goes further:

- **Linting** (`lint`): pre-commit hygiene hooks, ruff format --check, ruff check, mypy (strict).
- **Testing** (`test`): pytest.
- **Release**: a `v*` tag checks that the tag matches the version, then tests, builds and creates a GitHub release with sdist + wheel.
- **Ruleset on `main`**: PR required, required checks `lint` + `test` (strict, up to date), no force-push, no deletion. Approvals are not required while there is a single maintainer; raise to 1 when collaborators join.
- To be added in later phases: `skiver-compat` (phase 1), `export-roundtrip` (phase 9), PyPI trusted publishing (phase 13).

---

## 12. Open decisions

1. **Name.** `sequencing-error-model` is a working name. The package and import names (`sequencing_error_model`, CLI `sem`) should change together before v0.1.0.
2. **Licence.** `pyproject.toml` declares MIT to match upstream skiver (MIT); there is no LICENSE file yet. Confirm before ported code lands.
3. **Fitting stack.** numpy/scipy core with a torch extra (recommended), or torch throughout like the fork.
4. **Priority of wave-1 exporters.** Recommended order: Badread → InSilicoSeq → ART → PBSIM3 (context-rich first, then most-used Illumina tools).
5. **Fork strategy.** Make a new fork of GZHoffie/skiver for enhanced outputs, base some of the changes on the current `timrozday-mgnify/skiver` fork but try to minimize modifications and if possible just extend with distinct module(s).
6. **Evidence-source priority** (decided 2026-09-15). `pe-overlap` and `reference` modes first, then `kmer` default mode checked against them, then ONT duplex, importers and enhanced `kmer`. Phase 11 could move ahead of phase 10 if long-read realism matters more than second-wave exporters. Phases 7–8 (variation) come before exporters because every real-metagenome exit depends on them; exporters only need clonal specs and could run in parallel.
7. **Is skiver "primary"?** No longer assumed: `kmer` is one of three modes, checked against the other two in phase 6. If `reference` dominates in accuracy, position `kmer` as the fast, reference-free mode and `reference` as the high-fidelity path. The phase 6 comparison settles this with data. It is also the redundancy question (§4.4): if alignment sources dominate *and* existing profilers (ReSeq, Badread, CycSim) match the native generator on the same BAM, the project's value reduces to the joint Q model and the exporter hub, and scope should shrink accordingly.
8. **Factorisation.** Q-then-errors (§5.2, recommended: head E is directly the reported error profile and allows two-sided quality windows) vs errors-then-Q (Badread-style; simpler for long-read exports). The alternative is derivable from the joint model either way, so this only decides which one is fitted.
9. **Default `--q-policy`** for Q-coupled simulators: `preserve-quality` (recommended; realistic FASTQs for QC and aligner benchmarking) or `preserve-errors` (true error rates for assembly and variant benchmarking). It may be worth defaulting per use case.
10. **Quality windows vs neural window.** Start with additive `QualityWindow(m)` + low-rank `QualityxContext` and add `NeuralWindow` only if phase 5 shows held-out gains large enough to justify an export-unfriendly model.
11. **Baselines.** Which external simulators the baseline harness (phase 5) must include. Recommended: ReSeq for Illumina; Badread, PBSIM3 and CycSim for long reads. Also whether beating them should become a blocking exit criterion, rather than a reported one, before v0.1.0.
12. **Variation simulator.** Minimal in-package `variation.py` with external haplotypes accepted (recommended: exact per-base truth, no dependency), or CAMISIM/sgEvolver, SimBac or msprime as the generator.
13. **Q in site posteriors.** Whether the joint site model may use head E, and so Q windows, when scoring sites. Recommended: yes, since calibration is learned and not assumed, guarded by the miscalibrated-Q test with variants. If the guard fails, fall back to a Q-marginalised head E plus linkage.
14. **Minor-variant floor.** Which minor-strain fraction, divergence and coverage the project promises to separate. Set from the phase 7 bias table, before phase 8 methods are tuned.

## Sources

- Skiver preprint: https://www.biorxiv.org/content/10.64898/2026.02.12.705514v1
- Benchmarking long-read simulators on ONT data (2026): https://www.biorxiv.org/content/10.64898/2026.05.06.723380v1.full
- MeSS (ART + PBSIM3): https://academic.oup.com/bioinformatics/article/41/1/btae760/7935384
- Six short-read simulators evaluated (ART, DWGSIM, InSilicoSeq, Mason, NEAT, wgsim): https://www.nature.com/articles/s41437-022-00577-3
- PBSIM3: https://academic.oup.com/nargab/article/4/4/lqac092/6855700
- InSilicoSeq error models: https://insilicoseq.readthedocs.io/en/latest/iss/model.html
- Badread error / Q-score models: https://github.com/rrwick/Badread/wiki/Error-models, https://github.com/rrwick/Badread/wiki/QScore-models
- NanoSim: https://github.com/bcgsc/NanoSim
- art_modern: https://github.com/YU-Zhejian/art_modern
- NEAT: https://github.com/ncsa/NEAT
- CAMISIM configuration: https://github.com/CAMI-challenge/CAMISIM/wiki/Configuration-File-Options
- ErrorProfiler (paired-end overlap error profiling, 2025): https://www.biorxiv.org/content/10.1101/2025.07.15.665004v1.full
- GATK BaseRecalibrator: https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator
- Dorado duplex: https://software-docs.nanoporetech.com/dorado/latest/basecaller/duplex/
- Duplex basecalling for assembly (R. Wick, 2024): https://rrwick.github.io/2024/05/08/duplex_assemblies.html
- ReSeq (Genome Biology 2021): https://pmc.ncbi.nlm.nih.gov/articles/PMC7896392/, https://github.com/schmeing/ReSeq
- CycSim, context-aware long-read simulation (2025 preprint; GigaScience 2026): https://www.biorxiv.org/content/10.64898/2025.12.04.692264, https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giag079/8728720
- BEAR, sequence-read simulator for metagenomics (2014): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4168713/
- DRISEE: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3369934/
- LongISLND: https://doi.org/10.1093/bioinformatics/btw602
- GemSIM: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3305602/
- Skiver v2 preprint: https://www.biorxiv.org/content/10.64898/2026.02.12.705514v2
- inStrain SNV calling and null model: https://instrain.readthedocs.io/en/latest/user_manual.html
- LoFreq: https://pmc.ncbi.nlm.nih.gov/articles/PMC3526318/
- DADA2 `learnErrors` (self-consistent error learning): https://github.com/benjjneb/dada2/blob/master/man/learnErrors.Rd
- DESMAN: https://link.springer.com/article/10.1186/s13059-017-1309-9
- Floria: https://academic.oup.com/bioinformatics/article/40/Supplement_1/i30/7700908
- Strainy: https://www.nature.com/articles/s41592-024-02424-1
- VeChat: https://www.nature.com/articles/s41467-022-34381-8
- DeChat: https://www.nature.com/articles/s42003-024-07376-y
- CAMISIM strain simulation (sgEvolver): https://github.com/CAMI-challenge/CAMISIM/wiki/Strain-Simulation
- msprime 1.0: https://academic.oup.com/genetics/article/220/3/iyab229/6460344
- SimBac: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5049688/
- ZymoBIOMICS Gut Microbiome Standard (D6331): https://files.zymoresearch.com/protocols/_d6331_zymobiomics_gut_microbiome_standard.pdf
- Strain-level E. coli profiling benchmark on D6331 (2026): https://www.biorxiv.org/content/10.64898/2026.05.19.726160v1.full
