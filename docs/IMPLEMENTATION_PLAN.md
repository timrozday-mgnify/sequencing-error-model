# Implementation plan: sequencing-error-model

Status: **draft, written 2026-09-15**. Phase 0 (repo, CI and PR policy) is done; everything else is planned.

## 1. Goal

A Python package and CLI that:

1. **learns** a sequencing error model from skiver outputs, with no reference genome or alignment required. Skiver is the primary evidence source; optional additional sources (§4.3) fill the gaps skiver can't cover;
2. **generates** reads with that model, keeping the `FASTA in → FASTQ + CIGAR out` contract that genome-blender already uses; and
3. **exports** the model to the native formats of popular, maintained read simulators, so it plugs into existing pipelines (including the metagenome wrappers CAMISIM and MeSS).

It runs in two modes:

| Mode | Input | skiver build |
|---|---|---|
| **default** | `skiver analyze` CSVs | unmodified upstream release, pinned (currently **v0.3.2**) |
| **enhanced** | `skiver dump` per-observation outputs (plus the default CSVs) | modified skiver: today the `timrozday-mgnify/skiver` fork, later new outputs such as homopolymer indel lengths |

Out of scope: signal-level simulation (squigulator, seq2squiggle), variant/mutation models, abundance/community modelling (CAMISIM and MeSS already do this), and error *correction*.

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
| `summary_phred.csv` | correct/error counts per reported Q | Marginal P(error \| Q) calibration |
| `summary_read_position.csv` | correct/error counts per index from read start and from read end | Marginal error-rate curve along the read |
| `summary_gc_content.csv` | error rate per GC bin | GC dependence |
| **`kvmer.csv`** | per key: `key`, `consensus_value`, `passes_filter`, `homopolymer_length`, `consensus_count`, `neighbor_count`, `total_count`, **20 per-op counts**, `consensus_count_up_to_v1..v` | **Per-locus (k+v)-base sequence context with per-op error counts.** The edit *position* inside the value is not recorded |

**`kvmer.csv` is the key to default-mode feasibility.** It gives, for every sampled locus, the full k+v consensus sequence plus counts of each single-edit error type. For each key:

- the survival columns fix how many observations survive to each t;
- the op counts are a *mixture over the unknown position t* of the edit.

A context model P(op | preceding L bases) can be fit by marginalising that latent position. The likelihood of op *e* at a key is Σₜ S(t−1) · h_e(ctx_t) · S(v−t), fit by EM or direct gradient on the mixture. With the reference-free consensus as truth, the context at each t is known exactly. So a context-dependent substitution/indel model is learnable from unmodified skiver, without `dump`.

### 2.2 Modified skiver (fork): what enhanced mode adds

The fork adds `skiver dump` (`src/dump.rs`) with `--raw`, `--base`, `--survival` and `--windows`. `ValueInfo` also gains `read_id` and `fragment_id`, so R1/R2 share an id.

| Output | Adds |
|---|---|
| `base_observations.tsv` | Per-base rows with the known edit position, Phred, `read_pos`, `dist_to_end`, strand, `fragment_id`. This allows joint models such as P(op \| context, Q, position, strand) (~32 GB for hq-illumina) |
| `raw_observations.tsv` | Key string per observation, which seeds context when v < L |
| `survival_observations.tsv` | Per-observation first-error time |
| `windows.bin` | Read-ordered 13-byte records (~430 MB, 75× smaller than the TSV). Enables read-level latent-state HMM fitting |

Candidate *new* fork outputs for later phases (§8, phase 6): homopolymer run length at the edit and indel length (> 1 bp), multi-edit (2+) values aligned to consensus, an explicit R1/R2 flag, and per-position Q histograms.

---

## 3. Lessons from developing the skiver fork

These come from `skiver/docs/hmm_error_model.md`, `docs/performance.md`, `examples/*synthetic_recovery*`, commit history and the report-integration work.

1. **The value length v controls recoverability.** Amplicon recovery test (source rate ≈ 3.4e-4, BAM truth 3.43e-4):
   - `v=13`: retrained 3.48e-4 ✓
   - `v=6`: 2.76e-4 (≈ −20%)
   - `v=1`: **0** (no errors observable)

   Keep v ≈ 13 and fit context from key bases when v < L. Reject inputs with tiny v.
2. **Mask impossible categories.** Destination-base substitutions (`sub_to_A` when the true base is A) must be masked in the likelihood, or fitted mass leaks into categories the generator discards and the marginal rate drops. Bake this into the core likelihood.
3. **Additive context scales and full tables don't.** `BaseContext(L)` is 4ᴸ rows, which is memory-bound past L≈6–8. AIC kept selecting longer context up to the L=8 screen cap. Default to `AdditiveContext(L)` and materialise dense k-mer tables only at export time, where Badread uses 7-mers.
4. **Not every fitted component is generative.** `PhredContext` and `Position` condition on quality and position, which are *outputs* at simulation time; `FragmentOverdispersion` only affects the training likelihood. Keep training-only and generative components explicitly separated in the model spec.
5. **Synthetic recovery needs `--use-all -l 0`** (the outlier filter hides generated errors) and **reference-guided dump** to isolate model error from consensus error. Low-error Illumina needs large synthetic sets before metrics mean anything.
6. **Data size dominates engineering.** The 32 GB TSV needed an HDF5 cache and a single-aggregation model search (all candidates are slices of one count tensor). Plan for streaming aggregation into count tensors from day one.
7. **Reports.** Quarto `-P` parameters only reach a `#| tags: [parameters]` cell. Embedding plotly.js per figure plus MathJax OOM'd pandoc at ~6 GB. Emit plotly once and skip MathJax.
8. **Determinism.** Consensus ties used to follow HashMap order; now broken by smallest value. Golden tests on exact outputs were what caught this.
9. **`is_forward` is strand, not R1/R2.** Paired-end R1/R2 differences need the fork's `fragment_id`/`read_id` or a new flag.
10. **Code structure.** The fork's Python grew as flat `scripts/` with a 2.5 k-line library, pickled torch artifacts and hard-coded platform reports. This repo starts as a proper package with a safe, versioned artifact format.

---

## 4. Landscape

### 4.1 Error *estimation* tools (context)

| Tool | Needs reference/alignment | Output |
|---|---|---|
| **skiver** | No (optionally reference-guided) | rate, hazard, spectrum, Q calibration, position, GC |
| samtools stats / Qualimap / alignment-based (e.g. BEST) | Yes | mismatch/indel rates by cycle and quality. Biased by incomplete references and strain divergence (skiver paper) |
| k-mer spectrum tools (GenomeScope-style) | No | a single error rate, no spectrum |
| Simulator "profilers" (`iss model`, `badread error_model`, NanoSim `read_analysis.py`, `art_profile_builder`, `neat model-seq-err`) | Yes, BAM/PAF against a reference (ART and NEAT can use FASTQ for qualities) | simulator-specific models |

**Our niche:** simulator-ready models for samples where no good reference exists (metagenomes, novel isolates), learned from skiver. We produce the same model files the profilers would, without the alignment step.

### 4.2 Read simulators: export targets

Maintenance data is from the GitHub API on 2026-09-15.

| Simulator | Reads | Last push / ★ | Error model format | Where it's used | Export fidelity from our model |
|---|---|---|---|---|---|
| **ART** / **art_modern** | Illumina (+ MGI, AVITI, Onso profiles) | art_modern 2026-07 / 39 (original ART unmaintained, still ubiquitous) | Text quality profiles (per-position Q distributions, optionally per base); indel rates as CLI flags (`-ir/-dr/-ir2/-dr2`); substitutions derived from Q | MeSS, CAMISIM, countless pipelines | Medium. Rates and position curve map well; context is lost (ART has none); Q distributions need a base profile (§5.3) |
| **InSilicoSeq** | Illumina (MiSeq/HiSeq/NextSeq/NovaSeq) | 2026-09 / 227 | `.npz` with `read_length`, `insert_size`, `mean_count_{forward,reverse}`, `quality_hist_{forward,reverse}`, `subst_choices_*` (per position × base), `ins_*`, `del_*` | Metagenome benchmarking | Medium–high. Per-position substitution/indel tables are populated from the position curve × spectrum; forward/reverse arrays come from the strand spectrum |
| **NEAT v4** | Illumina | 2026-08 / 72 | gzip-pickle dict `{error_model1, error_model2, qual_score_model1, qual_score_model2}` of `SequencingErrorModel(read_length, variant_probs, indel_len_model, insertion_model, transition_matrix, avg_seq_error)` | Variant-calling benchmarks | Medium. Needs NEAT's classes to pickle, so export runs through NEAT as an optional dependency |
| **Mason2** (SeqAn) | Illumina/454/Sanger | seqan 2026-08 / 502 | CLI parameters only (mismatch/indel probabilities, begin/end ramps, quality mean/sd) | Aligner benchmarks | Low (parametric), but trivial to emit |
| **wgsim** | Illumina | 2021 / 288 (unmaintained) | CLI base error rate | CAMISIM | Low. Emit a rate only, for CAMISIM compatibility |
| **Badread** | ONT/PacBio | 2026-07 / 301 | Error model: text, one line per 7-mer `KMER,p;ALT,p;...` (≤ 25 alternatives); Q model: `CIGAR;count;q:p,...` keyed by local CIGAR window | Long-read tool benchmarks | **High for context.** 7-mer → alternatives is a direct materialisation of our context model. Q model needs error-conditioned Q (enhanced mode, or a base model) |
| **PBSIM3** | PacBio/ONT | 2025-04 / 117 | ERRHMM text: per accuracy level `IP`/`EP`(match,sub,ins,del)/`TP` lines; QSHMM for quality | MeSS; top of the 2026 ONT simulator benchmark for length/Q realism | Medium. Sequence-independent, so context is lost; our `HMM(S)` layer (enhanced) or β-derived clustering maps onto states |
| **NanoSim** | ONT | 2026-03 / 311 | Directory of pickled KDEs + text Markov/histogram files (`_error_markov_model`, `_match_markov_model`, `*_hist`, `_error_rate.tsv`, ...) from alignments | CAMISIM (`nanosim3`) | Low–medium. Only the error-rate/Markov parts are estimable; length KDEs must come from a base model. Lowest priority |
| **genome-blender** | all | internal | `skiver-generate` subprocess contract | This project's current consumer | Full: it uses our native generator |

Metagenome wrappers: **MeSS** wraps ART + PBSIM3, and **CAMISIM** wraps ART, NanoSim3 and wgsim. Covering **ART, InSilicoSeq, Badread and PBSIM3** therefore reaches most short-read, long-read and metagenome-simulation users. NEAT, Mason2, wgsim and NanoSim are a second wave.

The 2026 ONT simulator benchmark (Badread, LongISLND, lrsim, NanoSim, PBSIM3, SimLoRD) found that no simulator reproduces all metrics. PBSIM3 was best for length and Q; only Badread and LongISLND captured context-dependent substitutions spanning ~2 orders of magnitude; homopolymer errors were poorly reproduced by most. This supports two choices:

- exporting context to Badread;
- treating homopolymer modelling as the headline enhanced-mode feature.

### 4.3 Evidence sources beyond skiver

Skiver's gaps (§5.2) are mostly about *which axes are observed together*. It sees context and op, but only marginal Q, only marginal position, no R1/R2, no indel length, and no read-level structure. Other tools and data observe different axes with different truth assumptions. Many of them need no external reference.

The design consequence: **every source reduces to the same count tensors** (op × context × Q × read position × strand × mate × read id). Each source declares which axes it observes and what it uses as "truth". Sources can be fitted alone or jointly, and disagreement between them is a first-class diagnostic.

| Source | Truth used | Reference needed? | Adds what skiver lacks | Limitations | Maintained |
|---|---|---|---|---|---|
| **Raw FASTQ quality profiling** (built in; same idea as `art_profile_builder`) | none, Q only | No | **Per-position, per-mate, per-base Q distributions**, read-length distribution. Removes the need for a "base profile" in ART/InSilicoSeq export | No error information by itself; must be combined with an error calibration (skiver `summary_phred.csv` or the sources below) | n/a |
| **Paired-end overlap discordance**: ErrorProfiler (yaotianran/ErrorProfiler), or our own implementation over NGmerge/fastp merges | Mate agreement in the overlap | No | Exact **cycle position × Q × substitution × dinucleotide context × R1/R2**, jointly. PCR errors cancel because both mates share them, so it isolates sequencer error | Illumina only; needs short inserts (fully overlapping for ErrorProfiler); substitutions only; which mate erred is ambiguous (resolve by Q or symmetric EM); covers the read ends of both mates | ErrorProfiler: 2025 preprint, code on GitHub; NGmerge 2025-11; fastp 2026-09 |
| **Self-reference ("assemble, then align")**: metaFlye / myloasm / hifiasm-meta (long), metaSPAdes (short) → minimap2 → pileup | Polished consensus of high-coverage contigs | No external one (the sample assembles its own) | **Everything alignment gives**: indel lengths, homopolymer run-length errors, joint Q × position × context, R1/R2, read-level error regimes (per-read error sequences for `HMM(S)` / PBSIM3 ERRHMM), insert size. The simulators' own profilers (`iss model`, `badread error_model`/`qscore_model`, NanoSim `read_analysis.py`) can then run directly on the BAM | Assembly errors become false "truth" (ONT homopolymers especially), so polish (medaka) and mask sites with minor-allele frequency above a threshold, like BQSR known sites. Strain mixtures and repeats need masking. Aligner scoring bias (minimap2 favours mismatches over indels; skiver paper). Compute-heavy. Only high-coverage genomes, as with skiver | Flye 2026-04, myloasm 2026-09, hifiasm-meta 2025-11, SPAdes 2026-09, minimap2 2026-05, medaka 2026-05 |
| **ONT duplex vs simplex** (dorado duplex) | Duplex read, as near-truth for its two simplex parents | No | **Long-read truth per molecule**: homopolymer indels, Q calibration, per-read error trajectories, without assembly | Needs duplex-capable chemistry and runs; only the duplex-paired fraction of reads; residual duplex errors are correlated with simplex errors in homopolymers | dorado 2026-09 |
| **UMI / duplex consensus** (fgbio) | Molecular consensus | No | Per-molecule truth for amplicon and targeted libraries; separates PCR from sequencer error | Only UMI libraries | fgbio 2026-08 |
| **Spike-in and control references**: Illumina PhiX, ONT lambda/DCS, mock communities (e.g. ZymoBIOMICS) | Known reference | Yes, but a free, exact one | Truth-grade alignment profiles for the *same run*, including Illumina InterOp per-cycle error metrics | The control's library prep differs from the sample's; PhiX reads often sit in `Undetermined` or are removed; mock-community strain drift | Illumina/interop 2026-04 |
| **Existing error tables (importers)**: GATK BQSR recalibration tables; DADA2 `learnErrors` matrices | Reference + known sites (BQSR); denoised ASVs (DADA2) | BQSR yes; DADA2 no | BQSR: **joint read group × Q × cycle × context** mismatch counts in a stable text format. DADA2: reference-free substitution × Q error matrices for amplicons | BQSR needs known sites, which rarely exist for metagenomes (use self-reference masking instead). DADA2 is substitutions only, amplicon only, and its matrices are smoothed fits rather than raw counts | GATK 2026-09, dada2 2026-07 |
| **k-mer spectrum** (GenomeScope2) | k-mer histogram model | No | An independent global rate check | Single rate, isolate-like genomes only (poor for metagenomes) | 2025-10 |
| **Error-corrector diffs** (Lighter, Rcorrector; BFC unmaintained since 2016) | Corrected read | No | Per-read substitution calls with Q, position and mate | Miscorrection and under-correction are coverage-dependent and unmeasured; mostly substitutions. **Not recommended** except as a sanity check | Lighter/Rcorrector 2026-01 |

**Recommendation**, ordered by gap closed per unit effort:

1. **Raw FASTQ quality profiling.** Trivial and streaming. It removes the biggest default-mode export limitation (the quality profiles). Fold it into phase 1.
2. **Self-reference BAM source.** The most general gap-filler, especially for long reads (homopolymers, indel lengths, read-level states). It also gives a real-data **cross-check of skiver-derived models** on the same high-coverage genomes, which the fork only ever validated on synthetic data.
3. **Paired-end overlap source.** A reference-free joint Q × cycle × context × R1/R2 source for Illumina, and cheap.
4. **ONT duplex/simplex source.** The best reference-free long-read truth when duplex data exist.
5. **Importers**: BQSR tables, DADA2 matrices, InterOp metrics, spike-in BAMs (a special case of item 2 with a known reference).

---

## 5. Feasibility and limitations

### 5.1 Learnable in default mode (unmodified skiver)

| Model element | Source | Confidence |
|---|---|---|
| Per-base error rate, clustering (λ, β) | `summary_error_rate.csv` | High (the published method) |
| Error-type composition (12 subs + 4 ins + 4 del), strand asymmetry | `summary_error_spectrum.csv` | High |
| Context dependence P(op \| L preceding bases), L ≤ k | `kvmer.csv` with latent edit position (§2.1) | Medium. Identifiable because position varies across keys; needs a recovery test |
| Marginal error rate vs read position (start/end) | `summary_read_position.csv` | High for rate; composition along the read assumed constant |
| Marginal Q calibration P(error \| Q) | `summary_phred.csv` | High, but marginal only |
| GC dependence | `summary_gc_content.csv` | Medium |

### 5.2 Not learnable, or only weakly, in default mode

- **Joint effects** (context × Q × position × strand). Each CSV is a separate marginal; we model them as independent multiplicative (log-additive) factors. That assumption is untestable without `dump`.
- **Quality score generation.** skiver reports error counts per Q, not the Q distribution per read position, so full ART/ISS quality profiles can't be built from skiver alone. The fix is to profile Q directly from the raw FASTQ (§4.3 source 1), which needs no truth.
- **Indels longer than 1 bp and multi-error windows.** 2+-edit values are dropped, so the indel length distribution is unobserved. Assume 1 bp (fine for Illumina, wrong for ONT homopolymers).
- **Homopolymer effects.** `kvmer.csv` has only the longest run in the value, and the position of an indel inside a run isn't identifiable. Only a coarse run-length covariate is possible.
- **Read-level heterogeneity** (the latent quality regimes of PBSIM3 ERRHMM). Only the aggregate β is available.
- **R1 vs R2** differences are unobservable.
- **Read length and insert-size distributions** are not skiver's job. Take them from the user or the simulator's defaults.

### 5.3 Structural limitations (both modes)

- **Reference-free truth is the consensus.** Coverage ≥ ~20× on at least some genomes is needed. Strain heterogeneity, repeats and low-coverage metagenomes inflate or filter estimates, so the outlier filter matters. Document minimum-coverage guidance and surface `passes_filter` statistics in reports.
- **FracMinHash subsampling** (`-c`) trades precision for memory. Low-error platforms need a low `c` or large inputs for stable context tables.
- **Simulator formats are mostly position-based, not context-based.** Exporting to ART, InSilicoSeq, NEAT or PBSIM3 marginalises context away. Only Badread and our native generator keep it. Every exporter reports what was dropped.
- **Quality profiles need Q distributions.** The preferred path is to profile per-position, per-mate Q from the raw FASTQ (§4.3) and pair it with skiver's P(error | Q) calibration and position curve. When the FASTQ isn't available, ART and InSilicoSeq exporters fall back to a *base profile* (a built-in one such as `HiSeq`/`NovaSeq`, or a user file), **recalibrated** to the skiver-estimated curves. The fallback is recorded in the model provenance and the fidelity report.
- **Simulator quirks.** ART derives substitutions from Q, so its substitution rate and quality are coupled. PBSIM3 ERRHMM outputs `!` qualities. NEAT's model files are Python pickles of NEAT classes, so there is a version-coupling risk. NanoSim's model directory mixes KDE pickles tied to its scikit-learn version.
- **skiver upstream drift.** The CSV schema changed between releases (v0.3.x added GC and survival files). Parsers validate headers and pin supported version ranges; a CI matrix runs released skiver binaries.

### 5.4 Verdict

Feasible. Default mode can produce a context-aware substitution/indel model, a position curve, strand asymmetry and Q calibration from unmodified skiver. That is enough for credible Illumina exports (ART, InSilicoSeq) via base-profile recalibration, and for a context-preserving Badread export. Long-read realism (homopolymers, indel lengths, read-level regimes) genuinely needs enhanced mode. Plan for enhanced mode to be where ONT/PacBio quality comes from.

---

## 6. Architecture

```
src/sequencing_error_model/
  sources/skiver_analyze.py  # parse + validate v0.3.x analyze CSVs → count tables (default mode)
  sources/skiver_dump.py     # streaming reader for fork dump TSV / windows.bin (enhanced mode)
  sources/fastq_quality.py   # per-position/mate/base Q + read-length profile from raw FASTQ
  sources/bam.py             # self-reference / spike-in alignments → pileup counts (site masking)
  sources/pe_overlap.py      # paired-end overlap discordance
  sources/ont_duplex.py      # simplex-vs-duplex alignments
  sources/importers.py       # GATK BQSR tables, DADA2 error matrices, InterOp metrics
  counts.py                  # count tensors with an observed-axes mask + truth label per source
  spec.py                # ErrorModelSpec: versioned JSON manifest + .npz arrays (no pickle)
  fit/                   # likelihoods + fitters: marginal, latent-position context, covariates, HMM
  select.py              # greedy criterion-based component search (port of skiver model_selection)
  generate.py            # native generator: FASTA → FASTQ + CIGAR (genome-blender contract)
  export/{art,iss,badread,pbsim3,neat,mason,wgsim,nanosim}.py
  cli.py                 # sem fit | select | generate | export <sim> | report | inspect
```

Principles:

- **One simulator-neutral model spec** (`spec.py`), the only thing exporters read. It holds:
  - provenance: skiver version, mode, k, v, c, input hashes;
  - global rate plus λ and β;
  - op composition (strand-aware);
  - the context model (additive parameters, with a materialisation helper for a dense L-mer table);
  - the position curve, Q calibration and optional GC curve;
  - optional enhanced components (Phred, Position, Strand, Homopolymer, indel length, `HMM(S)`);
  - a per-component `generative: bool` flag.

  Stored as JSON + `.npz`: safe to load, diffable, language-agnostic. It is not a torch pickle.
- **Fitting stack.** numpy/scipy (L-BFGS on count tensors) for the default mode: small installs, easy bioconda packaging. torch (and optionally pyro for VI uncertainty) as an `enhanced` extra for large dump-driven models and the HMM. Revisit if one stack proves sufficient.
- **Counts first.** Both skiver modes and every other evidence source reduce to count tensors, so model search is cheap (lesson 6). Joint fitting sums per-source log-likelihoods of the model marginalised onto each source's observed axes. Per-source nuisance terms (skiver consensus misses, assembly errors) are optional. The report shows per-source disagreement before anything is pooled.
- **Exporters are thin, tested adapters.** Each declares:
  - required spec fields;
  - an optional base profile;
  - a machine-readable *fidelity report* of what was dropped or approximated.
- **Optional simulator dependencies** (NEAT classes, InSilicoSeq's loader) are imported lazily and only needed for export or round-trip tests.

---

## 7. Mode capability matrix

| Component | Default | Enhanced | Generative | Exporters that can use it |
|---|---|---|---|---|
| Global rate, op composition | ✓ | ✓ | ✓ | all |
| Strand asymmetry | ✓ | ✓ | ✓ | InSilicoSeq (fwd/rev), native |
| Context `AdditiveContext(L)` | ✓ (latent position) | ✓ (exact position) | ✓ | Badread, native |
| Read-position curve | ✓ (marginal) | ✓ (joint) | ✓ | ART, InSilicoSeq, NEAT, Mason (ramps), native |
| Q calibration | ✓ (marginal) | ✓ (joint with op) | recalibration / Q sampling | ART, InSilicoSeq, Badread Q model (enhanced), native |
| GC | ✓ | ✓ | ✓ | native |
| Clustering β | ✓ | ✓ | approx. | PBSIM3 (2-state approximation), native |
| `HMM(S)` read-level states | – | ✓ (`windows.bin`) | ✓ | PBSIM3 ERRHMM, native |
| Homopolymer / indel length | coarse | ✓ (needs new fork outputs) | ✓ | Badread (via k-mers), native |
| Fragment overdispersion, R1/R2 | – | ✓ | training only / R1–R2 profiles | ART/ISS/NEAT R2 profiles |

---

## 8. Phases

Each phase lands as one or more PRs with green CI. Exit criteria are testable.

### Phase 0: repository, CI and PR policy ✓
uv/ruff/mypy/pytest, pre-commit (revs standardised with the MIMICC/ENA repos), Linting, Testing and Release workflows, PR template, and a `main` ruleset requiring PRs with passing `lint` + `test`.

### Phase 1: default-mode inputs
- `sources/skiver_analyze.py` with typed, header-validated parsers for all v0.3.x analyze CSVs.
- `sources/fastq_quality.py`: a streaming per-position, per-mate, per-base Q histogram and read-length distribution from raw FASTQ(.gz).
- **Fixtures.** A tiny synthetic genome plus reads with known injected errors, analysed by the *released* skiver binary. Fixture files are committed (< 1 MB).
- **CI job `skiver-compat`.** Downloads skiver release binaries (matrix: v0.3.1, v0.3.2), regenerates the fixtures and diffs them against the committed parsers. It runs on PR only when `sources/skiver_*` changes, plus a weekly schedule to catch new releases.
- **Exit:** every analyze file round-trips into count tables; unsupported versions fail with a clear error.

### Phase 2: model spec and default-mode fitting
- `spec.py` (JSON + npz, schema-versioned, with provenance).
- Fitters:
  - marginal composition and strand;
  - latent-position `AdditiveContext(L)` from `kvmer.csv`, with the true-base mask;
  - position curve (smoothed, monotone-optional);
  - Q calibration;
  - GC;
  - Weibull passthrough.
- Criterion-based L selection (AIC/BIC on held-out keys).
- **Exit:** on simulated data with known context effects, the recovered context log-odds correlate with truth (target r > 0.9 for dominant effects) and the marginal rate is within 5%.

### Phase 3: native generator and recovery harness
- Port `error_application.py` ideas to `generate.py`: vectorised per read, CIGAR in the header, and the same CLI contract as `skiver-generate` so genome-blender can switch without code changes.
- A recovery harness: generate, run skiver (released binary), fit, compare. Metrics: TV/KL on op probabilities, marginal rate, position curve.
- A small version runs in CI (a few seconds of reads); the full-scale version is a manual `workflow_dispatch` workflow.
- **Exit:** the self-consistency loop passes at v=13, and the documented failure at small v is reproduced as a guarded error or warning.

### Phase 4: exporters, wave 1 (ART/art_modern, InSilicoSeq, Badread, PBSIM3)
- Each exporter supports a base profile and emits a fidelity report.
- **Round-trip tests.** Export, run the real simulator, run skiver, refit and compare rates to the spec within a tolerance. Simulators come from bioconda in a separate `export-roundtrip` workflow (weekly + `workflow_dispatch` + PRs touching `export/`), so the core Testing workflow stays fast.
- Unit tests validate generated files against each simulator's own loader where one is importable (InSilicoSeq `KDErrorModel`, Badread model parser).
- **Exit:** all four simulators run on exported models, and the round-trip error rate is within 10% of the spec (position-curve correlation reported).

### Phase 5: exporters, wave 2, and wrapper recipes
- NEAT v4 (optional dependency), Mason2 and wgsim parameter sets, and NanoSim error-model parts on a base model.
- Documented CAMISIM (`art`/`nanosim3`/`wgsim` config) and MeSS recipes pointing at exported profiles.
- **Exit:** the same round-trip criterion (NanoSim: best-effort, documented gaps).

### Phase 6: additional evidence sources
- `sources/bam.py`, a self-reference and spike-in pileup source:
  - site masking by minor-allele frequency;
  - optional polishing-aware contig filtering;
  - per-read error sequences for `HMM(S)`.

  Ship a documented recipe (assemble, polish, minimap2) rather than wrapping assemblers.
- `sources/pe_overlap.py` (Illumina) and `sources/ont_duplex.py` (dorado duplex pairs).
- `sources/importers.py` for GATK BQSR recalibration tables, DADA2 `learnErrors` matrices and InterOp error metrics.
- A cross-validation report comparing skiver-derived and alignment-derived models on the same high-coverage genomes (real data, not synthetic).
- Multi-source joint fit with per-source disagreement diagnostics.
- **Exit:** each source recovers a known model on synthetic data. On at least one real Illumina and one real ONT metagenome, the report quantifies the skiver vs self-reference agreement per model component.

### Phase 7: enhanced mode
- `sources/skiver_dump.py`: a streaming aggregator from TSV and `windows.bin` into the same count tensors (no 32 GB in memory). Consider asking for a Parquet or binary `--base` in the fork.
- Port the composable components: `PhredContext`, `Position`, `Strand`, `FragmentOverdispersion`, `Homopolymer`, `HMM(S)`, plus greedy selection, keeping the generative/training-only split.
- **Proposed fork additions** (separate skiver PRs, each gated by a dump format version header):
  - homopolymer run length and indel length at edits;
  - 2+-edit values aligned to consensus;
  - R1/R2 flag;
  - per-position Q histogram;
  - error-conditioned Q, for the Badread Q model and ART/ISS profiles without a base.
- Offer the non-invasive ones upstream to GZHoffie/skiver so more of the default mode improves over time.
- **Exit:** enhanced models beat default models on held-out likelihood on real ONT and Illumina datasets, and PBSIM3/Badread exports gain the HMM and Q model.

### Phase 8: reports, packaging, release
- `sem report`: a single HTML page (plotly loaded once, no MathJax) with coverage/filter diagnostics, fitted curves, the fidelity report per export, and selection traces.
- PyPI and bioconda recipes, plus a container image for Nextflow (nf-core-style module for the synthetic-metagenomic-benchmark pipeline).
- **Exit:** v0.1.0 tagged via the Release workflow; the bioconda PR is open.

---

## 9. Testing strategy

| Layer | Where | Runs |
|---|---|---|
| Unit + property tests (parsers, likelihood masks, spec I/O) | `tests/` | every PR (Testing) |
| Golden fixtures from pinned skiver releases | `tests/fixtures/skiver-vX.Y.Z/` | every PR; regenerated by `skiver-compat` |
| Small statistical recovery (fixed seeds, tolerance bands) | `tests/recovery/` | every PR, under ~2 min |
| Exporter round-trips with real simulators | `export-roundtrip` workflow | weekly, dispatch, PRs touching `export/` |
| Full-scale recovery on real datasets | manual workflow or HPC | before releases |

---

## 10. CI/CD and PR policy (in place)

The standard matches `EBI-Metagenomics/mimicc-ena-submission-assistant`: uv + Python 3.12, ruff format/check, mypy, pytest, and pre-commit with the same pinned revs. That reference repo has no branch protection configured, so this repo goes further:

- **Linting** (`lint`): pre-commit hygiene hooks, ruff format --check, ruff check, mypy (strict).
- **Testing** (`test`): pytest.
- **Release**: a `v*` tag checks that the tag matches the version, then tests, builds and creates a GitHub release with sdist + wheel.
- **Ruleset on `main`**: PR required, required checks `lint` + `test` (strict, up to date), no force-push, no deletion. Approvals are not required while there is a single maintainer; raise to 1 when collaborators join.
- To be added in later phases: `skiver-compat` (phase 1), `export-roundtrip` (phase 4), PyPI trusted publishing (phase 8).

---

## 11. Open decisions

1. **Name.** `sequencing-error-model` is a working name. The package and import names (`sequencing_error_model`, CLI `sem`) should change together before v0.1.0.
2. **Licence.** `pyproject.toml` declares MIT to match upstream skiver (MIT); there is no LICENSE file yet. Confirm before ported code lands.
3. **Fitting stack.** numpy/scipy core with a torch extra (recommended), or torch throughout like the fork.
4. **Priority of wave-1 exporters.** Recommended order: Badread → InSilicoSeq → ART → PBSIM3 (context-rich first, then most-used Illumina tools).
5. **Fork strategy.** Keep enhanced outputs in the `timrozday-mgnify/skiver` fork, or upstream the dump subcommand to GZHoffie/skiver.
6. **Evidence-source priority.** Recommended order: FASTQ quality (phase 1), then self-reference BAM, paired-end overlap, ONT duplex, importers. Phase 6 could move ahead of phase 5 if long-read realism matters more than second-wave exporters.
7. **Is skiver "primary"?** If the self-reference source turns out to dominate in accuracy, reposition skiver as the fast, reference-free default and alignment sources as the high-fidelity path. The phase 6 cross-validation report settles this with data.

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
