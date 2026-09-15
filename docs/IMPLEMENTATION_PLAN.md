# Implementation plan: sequencing-error-model

Status: **draft, revised 2026-09-15**. Phases 0 (repo, CI and PR policy) and 1 (default-mode inputs) are done; phase 2 is in progress (`spec.py` landed); everything else is planned.

## 1. Goal

A Python package and CLI that:

1. **learns** a sequencing error model from skiver outputs, with no reference genome or alignment required. Skiver is the primary evidence source; optional additional sources (§4.3) fill the gaps skiver can't cover;
2. **models bases and qualities jointly** (§5). The core quantity is the error profile of a base *given its own quality score, the surrounding bases and qualities*, and other context such as read position, mate, strand and read-level state;
3. **generates** reads with that model, outputting **both bases and quality scores**, and keeping the `FASTA in → FASTQ + CIGAR out` contract that genome-blender already uses; and
4. **exports** the model to the native formats of popular, maintained read simulators, so it plugs into existing pipelines (including the metagenome wrappers CAMISIM and MeSS).

**Quality scores are never treated as evidence of error.** A reported Q is a feature that may covary with the true error rate, and an output the generator must reproduce. Error labels only come from truth-bearing sources: skiver consensus, mate overlap, alignment to a self-assembled or control reference, duplex reads, and so on.

It runs in two modes:

| Mode | Input | skiver build |
|---|---|---|
| **default** | `skiver analyze` CSVs (+ raw FASTQ for the quality model) | unmodified upstream release, pinned (currently **v0.3.2**) |
| **enhanced** | `skiver dump` per-observation outputs (plus the default inputs) | modified skiver: a new, minimal fork of GZHoffie/skiver, seeded from the `timrozday-mgnify/skiver` fork (§12, decision 5) |

Out of scope: signal-level simulation (squigulator, seq2squiggle), variant/mutation models, abundance/community modelling (CAMISIM and MeSS already do this), and error *correction*. Quality recalibration is produced as a diagnostic, not as a read-rewriting tool.

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

Candidate *new* outputs for the new minimal fork (§9, phase 7):

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
| **InSilicoSeq** | Illumina (MiSeq/HiSeq/NextSeq/NovaSeq) | 2026-09 / 227 | `.npz` with `read_length`, `insert_size`, `mean_count_{forward,reverse}`, `quality_hist_{forward,reverse}`, `subst_choices_*` (per position × base), `ins_*`, `del_*` | Q-driven substitutions (to verify in phase 4) plus per-position substitution choices and indel rates | Metagenome benchmarking | Medium. Per-position Q histograms and substitution/indel tables are populated from marginals of both heads; same calibration caveat as ART |
| **NEAT v4** | Illumina | 2026-08 / 72 | gzip-pickle dict `{error_model1, error_model2, qual_score_model1, qual_score_model2}` (`SequencingErrorModel` + quality Markov model per mate) | Errors from Q-derived probabilities; Q from a Markov model | Variant-calling benchmarks | Medium. The Q head maps onto NEAT's quality Markov model (order 1, per mate). Needs NEAT's classes to pickle |
| **Mason2** (SeqAn) | Illumina/454/Sanger | seqan 2026-08 / 502 | CLI parameters only (mismatch/indel probabilities, begin/end ramps, quality mean/sd for correct and wrong bases) | Separate Q mean/sd for correct vs erroneous bases | Aligner benchmarks | Low (parametric), but trivial to emit. The correct/wrong Q split comes from the joint model |
| **wgsim** | Illumina | 2021 / 288 (unmaintained) | CLI base error rate | none | CAMISIM | Low. Emit a rate only, for CAMISIM compatibility |
| **Badread** | ONT/PacBio | 2026-07 / 301 | Error model: text, one line per 7-mer `KMER,p;ALT,p;...` (≤ 25 alternatives); Q model: `CIGAR;count;q:p,...` keyed by local CIGAR window | **Q conditioned on the local error pattern** (CIGAR window around the base) | Long-read tool benchmarks | **High.** 7-mer → alternatives materialises the base-context part; P(Q \| local CIGAR window) is a marginal of the joint model. Neighbour-Q and position effects are lost |
| **PBSIM3** | PacBio/ONT | 2025-04 / 117 | ERRHMM text (`IP`/`EP`(match,sub,ins,del)/`TP` per accuracy level) or QSHMM (quality HMM; errors from Q) | ERRHMM: none (qualities `!`). QSHMM: latent-state Q process, errors from Q | MeSS; top of the 2026 ONT simulator benchmark for length/Q realism | Medium. Export **QSHMM** for Q-bearing output (maps our latent state + Q head) and ERRHMM for error-only use. Sequence context is lost |
| **NanoSim** | ONT | 2026-03 / 311 | Directory of pickled KDEs + text Markov/histogram files (`_error_markov_model`, `_match_markov_model`, `*_hist`, `_error_rate.tsv`, quality models) from alignments | Quality models conditioned on match/error state | CAMISIM (`nanosim3`) | Low–medium. Only the error/quality Markov parts are estimable; length KDEs must come from a base model. Lowest priority |
| **ReSeq** | Illumina | 2021 paper / GitHub (maintenance to check) | Binary stats file from `reseq illuminaPE` statistics step (format to confirm in phase 5) | **Errors conditioned on Q**, position, errors so far in the read, reference base and recent dominant error; Q conditioned on previous Q, position, mate, tile, sequence quality and reference base; per-site systematic errors; stored as 2-D margins | Illumina tool benchmarks | Medium–high if the format is writable: head Q and centre-Q head E map closely. Preceding-sequence context is not modelled by ReSeq and is lost |
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

**Recommendation**, ordered by the gap closed in the joint base × quality model per unit effort (matches §12, decision 6):

1. **Paired-end overlap source.** The cheapest reference-free source that identifies **neighbouring-quality and quality × context effects** for Illumina, which default-mode skiver cannot.
2. **Self-reference BAM source.** The most general gap-filler, especially for long reads (homopolymers, indel lengths, read-level quality/error regimes). It also gives a real-data **cross-check of skiver-derived models** on the same high-coverage genomes, which the fork only ever validated on synthetic data.
3. **ONT duplex/simplex source.** The best reference-free long-read truth when duplex data exist.
4. **Raw FASTQ quality profiling.** Low as *evidence* (it carries none), but a hard dependency for *generating* qualities in default mode, so it is still built early (phase 1) as a feature/output source.
5. **Importers**: BQSR tables, DADA2 matrices, InterOp metrics, spike-in BAMs (a special case of item 2 with a known reference).

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
| **This project, default mode** | Yes | Yes | Two-sided, latent-position `Context(L,R)` | FASTQ Q Markov process | Centre-Q term by marginal matching; no Q window |
| **This project, full model** | Yes (enhanced skiver, PE overlap) or BAM/duplex | Yes | Two-sided + `QualityxContext` | Q conditioned on true context, shared `Latent(S)` | Q first, then errors given the Q window; Q never a label |

**Clearly distinct:**

- simulator-ready models trained without alignment. BEAR is the only precedent, and it has no context and no Q coupling;
- quality as a feature, with the calibration gap learned rather than assumed;
- one spec exported to many simulators, with `--q-policy` and fidelity reports;
- multi-source fitting with disagreement diagnostics.

**Redundancy risks, each tied to a test:**

1. **The reference-free advantage is narrower than it looks.** Skiver needs ≥ ~20× coverage on some genomes, roughly the genomes that also self-assemble. On a self-reference BAM, ReSeq, Badread, CycSim and NanoSim can train directly. Skiver's remaining advantages are speed, no aligner bias and no assembly errors used as truth. *Test:* the skiver-vs-BAM cross-check (phase 6, §12 decision 7).
2. **Default mode is only modestly richer than existing profiles.** It is roughly Badread-level context plus an ART/NEAT-level quality process. *Test:* comparative exit criteria against BAM-trained baselines (phases 2–3, §12 decision 11).
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

| Parameter | default skiver (+ FASTQ) | enhanced skiver dump | PE overlap | self-reference BAM / duplex / spike-in |
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

### 6.1 Learnable in default mode (unmodified skiver + raw FASTQ)

| Model element | Source | Confidence |
|---|---|---|
| Per-base error rate, clustering (λ, β) | `summary_error_rate.csv` | High (the published method) |
| Error-type composition (12 subs + 4 ins + 4 del), strand asymmetry | `summary_error_spectrum.csv` | High |
| Head E base context, two-sided within k+v | `kvmer.csv` with latent edit position (§2.1) | Medium. Identifiable because position varies across keys; needs a recovery test |
| Head E centre-Q term | `summary_phred.csv` + FASTQ exposure (marginal matching, §5.6) | Medium. Correct only under the log-additive assumption |
| Marginal error rate vs read position (start/end) | `summary_read_position.csv` | High for rate; composition along the read assumed constant |
| Head Q: per-position/mate Q distributions, Q Markov transitions, Q given observed context | raw FASTQ | High as a *quality process*; carries no information about errors |
| GC dependence | `summary_gc_content.csv` | Medium |

### 6.2 Not learnable, or only weakly, in default mode

- **Neighbouring-quality effects and quality × context interactions** in head E. Unmodified skiver never pairs errors with quality windows or with context-specific qualities. The effects are fixed at zero and flagged; closing this needs enhanced skiver or a truth-bearing source such as PE overlap.
- **Joint effects of position, strand and mate with context or Q.** Each CSV is a separate marginal. They are combined log-additively, an assumption that is untestable without per-observation data.
- **Quality head conditioned on true context and on errors.** FASTQ conditions on observed bases. That is acceptable for Illumina (error ≈ 10⁻³), but biased for ONT, where erroneous bases are frequent and correlated with low Q.
- **Indels longer than 1 bp and multi-error windows.** 2+-edit values are dropped, so the indel length distribution is unobserved. Assume 1 bp (fine for Illumina, wrong for ONT homopolymers).
- **Homopolymer effects.** `kvmer.csv` has only the longest run in the value, and the position of an indel inside a run isn't identifiable. Only a coarse run-length covariate is possible.
- **Read-level heterogeneity** jointly affecting Q and errors. Only the aggregate β is available.
- **R1 vs R2** differences in errors are unobservable; R1 vs R2 qualities come from the FASTQ.

### 6.3 Structural limitations (both modes)

- **Reference-free truth is the consensus.** Coverage ≥ ~20× on at least some genomes is needed. Strain heterogeneity, repeats and low-coverage metagenomes inflate or filter estimates, so the outlier filter matters. Document minimum-coverage guidance and surface `passes_filter` statistics in reports.
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

### 6.4 Verdict

Feasible, in stages:

- **Default mode** yields a context-aware error head with a centre-quality term, plus a realistic quality process. Generation outputs full FASTQ (bases + qualities) for the native generator. Credible Illumina exports and a context-preserving Badread export follow.
- **The full base × quality-window model**, the one described in §5, needs per-observation quality windows. For Illumina, PE overlap supplies them reference-free; enhanced skiver and self-reference/duplex supply them for all platforms. Long-read realism (homopolymers, indel lengths, read-level regimes) needs the same.

---

## 7. Architecture

```
src/sequencing_error_model/
  sources/skiver_analyze.py  # parse + validate v0.3.x analyze CSVs → marginal count tables (default mode)
  sources/skiver_dump.py     # streaming reader for fork dump TSV / windows.bin (enhanced mode)
  sources/fastq_quality.py   # quality process: per-position/mate Q, Q transitions, Q | observed context, lengths
  sources/pe_overlap.py      # paired-end overlap discordance (with both mates' quality windows)
  sources/bam.py             # self-reference / spike-in alignments → observation tuples (site masking)
  sources/ont_duplex.py      # simplex-vs-duplex alignments
  sources/importers.py       # GATK BQSR tables, DADA2 error matrices, InterOp metrics
  observations.py            # sparse (context window, Q window, covariates, op) → count tuples; truth label per source
  spec.py                    # ErrorModelSpec: versioned JSON manifest + .npz arrays (no pickle)
  fit/                       # quality head, error head, marginal matching, latent state, selection criteria
  select.py                  # greedy criterion-based component search over both heads
  generate.py                # native generator: FASTA → FASTQ (bases + qualities) + CIGAR
  export/{art,iss,badread,pbsim3,neat,mason,wgsim,nanosim}.py
  cli.py                     # sem fit | select | generate | export <sim> | report | inspect
```

Principles:

- **One simulator-neutral model spec** (`spec.py`), the only thing exporters and the generator read. It holds:
  - provenance: sources, skiver version, mode, k, v, c, input hashes, base-profile fallback;
  - the **quality alphabet** and binning;
  - **head Q** and **head E**: component tokens, parameters, per-component `generative` and `identified` flags (e.g. neighbouring-Q weights fixed in default mode);
  - the optional latent-state layer shared by both heads;
  - summary marginals cached for exporters and reports (rate, λ/β, op composition, position curves, calibration curve).

  Stored as JSON + `.npz`: safe to load, diffable, language-agnostic. It is not a torch pickle.
- **Sparse tuples first.** Every source reduces to deduplicated (context window, quality window, covariate bins, op) tuples with counts. Binned qualities make this compress hard for Illumina. When unique tuples exceed memory (long windows, ONT's full Q range), fitting falls back to streaming minibatches over the same schema. Marginal-only sources (skiver analyze, importers) enter the joint likelihood through their marginal constraints (§5.6).
- **Joint fitting across sources.** The per-source log-likelihood of the model is marginalised onto each source's observed fields and summed, with optional per-source nuisance terms (skiver consensus misses, assembly errors, mate ambiguity). The report shows per-source disagreement before anything is pooled.
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

| Component | Default | Enhanced | Generative | Exporters that can use it |
|---|---|---|---|---|
| Global rate, op composition | ✓ | ✓ | ✓ | all |
| Strand asymmetry | ✓ | ✓ | ✓ | InSilicoSeq (fwd/rev), native |
| E: `Context(L,R)` | ✓ (latent position) | ✓ (exact position) | ✓ | Badread, native |
| E: centre Q | ✓ (marginal matching) | ✓ | ✓ | all Q-coupled exporters via `--q-policy`, Mason2 (correct/wrong Q), native |
| E: `QualityWindow(m)`, `QualityxContext` | – (fixed at 0) | ✓ (value window; key side needs fork addition) | ✓ | native only |
| Q: per-position/mate Q, `QualityMarkov(m)` | ✓ (FASTQ) | ✓ | ✓ | ART, InSilicoSeq, NEAT (order 1), PBSIM3 QSHMM, native |
| Q: conditioned on true context / errors | approx. (observed bases) | ✓ | ✓ | Badread Q model, NanoSim, native |
| Position | ✓ (marginal) | ✓ (joint) | ✓ | ART, InSilicoSeq, NEAT, Mason (ramps), native |
| GC | ✓ | ✓ | ✓ | native |
| Clustering β | ✓ | ✓ | approx. | PBSIM3 (2-state approximation), native |
| `Latent(S)` shared by E and Q | – | ✓ (needs Q in `windows.bin`) | ✓ | PBSIM3 QSHMM/ERRHMM, native |
| Homopolymer / indel length | coarse | ✓ (needs new fork outputs) | ✓ | Badread (via k-mers), native |
| Fragment overdispersion, R1/R2 | – | ✓ | training only / R1–R2 profiles | ART/ISS/NEAT R2 profiles |

---

## 9. Phases

Each phase lands as one or more PRs with green CI. Exit criteria are testable.

### Phase 0: repository, CI and PR policy ✓
uv/ruff/mypy/pytest, pre-commit (revs standardised with the MIMICC/ENA repos), Linting, Testing and Release workflows, PR template, and a `main` ruleset requiring PRs with passing `lint` + `test`.

### Phase 1: default-mode inputs ✓
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

### Phase 2: model spec, both heads, default-mode fitting (in progress)
- ✓ `spec.py`: an `ErrorModelSpec` directory of `spec.json` (schema version, quality alphabet, provenance with required `mode`, component tokens, `generative`/`identified` flags, meta) and `arrays.npz` (parameters and cached marginals; object arrays are rejected and loading never unpickles). Tokens are validated against the §5.4 component names and the heads each may appear in, plus `GC`, `Weibull` and `InsertionQuality` from this phase; the shared `Latent(S)` layer sits outside both heads. Argument arity is left to each fitter. Loading rejects other schema versions, missing arrays and unreferenced arrays.
- Head Q fitters: `QualityMarkov(m)`, `Position`, `Mate`, `Context` on observed bases, insertion-quality sub-head (base-profile prior in default mode).
- Head E fitters:
  - latent-position `Context(L,R)` from `kvmer.csv` with the true-base mask;
  - centre-Q term by marginal matching against `summary_phred.csv` under the FASTQ exposure;
  - position curve, strand, GC, Weibull passthrough.
- Criterion-based selection per head (AIC/BIC on held-out keys and reads).
- **Exit:** on simulated data with known context and quality effects:
  - recovered context log-odds correlate with truth (r > 0.9 for dominant effects);
  - the error rate implied per reported-Q bin is within 10% of truth;
  - per-position Q histograms are within TV < 0.05;
  - the marginal error rate is within 5%.
- **Comparative check (not blocking, reported):** on a dataset with a reference, the default-mode model's held-out likelihood and marginals are compared with ReSeq (Illumina) and Badread (long reads) profiles trained on the aligned reads (§4.4, risk 2).

### Phase 3: native generator (bases + qualities) and recovery harness
- `generate.py` implementing §5.3. Qualities are always emitted from head Q, and `--no-quality` is kept only for compatibility. It keeps the `skiver-generate` CLI contract so genome-blender can switch without code changes, and adds an alignment-consistency self-test.
- A recovery harness: generate, run skiver (released binary), profile the generated FASTQ, fit, compare. Metrics:
  - **errors:** TV/KL on op probabilities, marginal rate, position curve;
  - **qualities:** per-position Q TV, Q lag-1 autocorrelation, P(error | reported Q) calibration curve vs source, empirical-vs-reported Q curve.
- A small version runs in CI (a few seconds of reads); the full-scale version is a manual `workflow_dispatch` workflow.
- **Baseline harness.** The same metrics computed on real reads vs reads from (a) the native generator, (b) ReSeq for Illumina, (c) Badread, PBSIM3 and CycSim for long reads, each simulator trained by its own profiler on an alignment of the same sample. Add k-mer spectrum concordance and context-dependent substitution rates, the metrics of the 2026 ONT benchmark.
- **Exit:**
  - the self-consistency loop passes at v=13 for both heads, and the documented failure at small v is reproduced as a guarded error or warning;
  - the baseline harness runs on at least one Illumina and one long-read dataset, and the report states where the native generator beats, matches or trails each baseline.

### Phase 4: exporters, wave 1 (ART/art_modern, InSilicoSeq, Badread, PBSIM3)
- Each exporter:
  - implements `--q-policy`, where the simulator couples errors to Q;
  - supports a base-profile fallback;
  - emits a fidelity report.
- Badread gets both the error model and the Q model. PBSIM3 gets QSHMM (Q-bearing) and ERRHMM.
- **Round-trip tests.** Export, run the real simulator, run skiver and FASTQ profiling, refit, and compare both the error rate and the Q marginals to the spec within tolerance. Simulators come from bioconda in a separate `export-roundtrip` workflow (weekly + `workflow_dispatch` + PRs touching `export/`), so the core Testing workflow stays fast.
- Unit tests validate generated files against each simulator's own loader where one is importable (InSilicoSeq `KDErrorModel`, Badread model parser). Confirm how InSilicoSeq couples substitutions to Q.
- **Exit:** all four simulators run on exported models; under the chosen `--q-policy`, the preserved quantity is within 10% of the spec, and the violated one is quantified in the fidelity report.

### Phase 5: exporters, wave 2, and wrapper recipes
- NEAT v4 (quality Markov + error model; optional dependency), Mason2 (including correct/wrong-base Q parameters) and wgsim parameter sets, and NanoSim error and quality model parts on a base model.
- Candidates, pending a format check: **ReSeq** (its Q process and Q-conditioned errors map closely onto both heads, a better Illumina target than ART if its stats file can be written) and **CycSim** (k-mer errors and error-state transitions for long reads). Drop either if its model file is not documented or stable.
- Documented CAMISIM (`art`/`nanosim3`/`wgsim` config) and MeSS recipes pointing at exported profiles.
- **Exit:** the same round-trip criterion (NanoSim: best-effort, documented gaps).

### Phase 6: additional evidence sources
In the order of §12 decision 6: PE overlap, self-reference BAM, ONT duplex, then importers. PE overlap is the first source that identifies `QualityWindow` effects reference-free; consider pulling it forward to directly after phase 3.

- `sources/pe_overlap.py` (Illumina): observation tuples carrying both mates' quality windows at each overlapping template base, with symmetric EM for mate attribution.
- `sources/bam.py`, a self-reference and spike-in source:
  - site masking by minor-allele frequency;
  - optional polishing-aware contig filtering;
  - per-read error *and quality* trajectories for `Latent(S)`.

  Ship a documented recipe (assemble, polish, minimap2) rather than wrapping assemblers.
- `sources/ont_duplex.py` (dorado duplex pairs, simplex Q tracks).
- `sources/importers.py` for GATK BQSR recalibration tables, DADA2 `learnErrors` matrices and InterOp error metrics.
- A cross-validation report comparing skiver-derived and alignment-derived models, per head and per component, on the same high-coverage genomes (real data, not synthetic). This is the project's main redundancy test (§4.4, risk 1). A minimal version, using an off-the-shelf assemble-and-align recipe and existing profilers rather than `sources/bam.py`, should run as early as phase 3.
- Multi-source joint fit with per-source disagreement diagnostics.
- **Exit:**
  - each source recovers a known model, including `QualityWindow` effects, on synthetic data;
  - on at least one real Illumina and one real ONT metagenome, the report quantifies skiver vs other-source agreement per component;
  - `QualityWindow` beats the centre-Q-only head E on held-out likelihood, or the report shows it doesn't.

### Phase 7: enhanced mode
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
- **Exit:** enhanced models beat default models on held-out likelihood for both heads on real ONT and Illumina datasets; PBSIM3/Badread exports gain the latent-state and error-conditioned Q models.

### Phase 8: reports, packaging, release
- `sem report`: a single HTML page (plotly loaded once, no MathJax) with:
  - coverage/filter diagnostics;
  - both heads' fitted effects;
  - calibration (empirical vs reported Q, by context and position);
  - identifiability flags;
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
| Exporter round-trips with real simulators (errors + qualities, per `--q-policy`) | `export-roundtrip` workflow | weekly, dispatch, PRs touching `export/` |
| Full-scale recovery on real datasets | manual workflow or HPC | before releases |

---

## 11. CI/CD and PR policy (in place)

The standard matches `EBI-Metagenomics/mimicc-ena-submission-assistant`: uv + Python 3.12, ruff format/check, mypy, pytest, and pre-commit with the same pinned revs. That reference repo has no branch protection configured, so this repo goes further:

- **Linting** (`lint`): pre-commit hygiene hooks, ruff format --check, ruff check, mypy (strict).
- **Testing** (`test`): pytest.
- **Release**: a `v*` tag checks that the tag matches the version, then tests, builds and creates a GitHub release with sdist + wheel.
- **Ruleset on `main`**: PR required, required checks `lint` + `test` (strict, up to date), no force-push, no deletion. Approvals are not required while there is a single maintainer; raise to 1 when collaborators join.
- To be added in later phases: `skiver-compat` (phase 1), `export-roundtrip` (phase 4), PyPI trusted publishing (phase 8).

---

## 12. Open decisions

1. **Name.** `sequencing-error-model` is a working name. The package and import names (`sequencing_error_model`, CLI `sem`) should change together before v0.1.0.
2. **Licence.** `pyproject.toml` declares MIT to match upstream skiver (MIT); there is no LICENSE file yet. Confirm before ported code lands.
3. **Fitting stack.** numpy/scipy core with a torch extra (recommended), or torch throughout like the fork.
4. **Priority of wave-1 exporters.** Recommended order: Badread → InSilicoSeq → ART → PBSIM3 (context-rich first, then most-used Illumina tools).
5. **Fork strategy.** Make a new fork of GZHoffie/skiver for enhanced outputs, base some of the changes on the current `timrozday-mgnify/skiver` fork but try to minimize modifications and if possible just extend with distinct module(s).
6. **Evidence-source priority.** Recommended order: paired-end overlap, then self-reference BAM, ONT duplex, FASTQ quality, importers. Phase 6 could move ahead of phase 5 if long-read realism matters more than second-wave exporters.
7. **Is skiver "primary"?** If the self-reference source turns out to dominate in accuracy, reposition skiver as the fast, reference-free default and alignment sources as the high-fidelity path. The cross-validation report (minimal version from phase 3, full version in phase 6) settles this with data. It is also the redundancy question (§4.4): if alignment sources dominate *and* existing profilers (ReSeq, Badread, CycSim) match the native generator on the same BAM, the project's value reduces to the joint Q model and the exporter hub, and scope should shrink accordingly.
8. **Factorisation.** Q-then-errors (§5.2, recommended: head E is directly the reported error profile and allows two-sided quality windows) vs errors-then-Q (Badread-style; simpler for long-read exports). The alternative is derivable from the joint model either way, so this only decides which one is fitted.
9. **Default `--q-policy`** for Q-coupled simulators: `preserve-quality` (recommended; realistic FASTQs for QC and aligner benchmarking) or `preserve-errors` (true error rates for assembly and variant benchmarking). It may be worth defaulting per use case.
10. **Quality windows vs neural window.** Start with additive `QualityWindow(m)` + low-rank `QualityxContext` and add `NeuralWindow` only if phase 6 shows held-out gains large enough to justify an export-unfriendly model.
11. **Baselines.** Which external simulators the baseline harness (phase 3) must include. Recommended: ReSeq for Illumina; Badread, PBSIM3 and CycSim for long reads. Also whether beating them should become a blocking exit criterion, rather than a reported one, before v0.1.0.

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
