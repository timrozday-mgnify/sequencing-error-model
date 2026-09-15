# `pe-overlap` mode

Fits an error model from Illumina read pairs whose mates overlap. Where both mates read the same template
base, their agreement is the truth for error labels. No reference is needed. Code:
[`sources/pe_overlap.py`](../src/sequencing_error_model/sources/pe_overlap.py). Design and evidence: plan
§6.4 and phase 4 in [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md).

## Usage

```bash
uv run python -m sequencing_error_model.sources.pe_overlap R1.fastq.gz R2.fastq.gz --output spec/ [--max-pairs N]
```

The command writes an `ErrorModelSpec` directory (`spec.json` + `arrays.npz`) and prints the pair statistics,
which are also stored in the spec's provenance:

| stat | meaning |
|---|---|
| `pairs` | pairs read |
| `no_overlap` | no confident overlap: none ≥ `--min-overlap` bases, the best offset above 20% mismatches, or another offset also matching well |
| `indel` | overlap implies an indel (dropped; gapless evidence can't label indels) |
| `agree`, `disagree` | overlapping template bases where the mates agree or disagree |

The spec can be passed straight to `sem-generate --model spec/`.

| option | default | |
|---|---|---|
| `--error-tokens` | `QualityWindow(1) Context(1,1) Position(4) Mate` | head E components |
| `--quality-tokens` | `QualityMarkov(1) Position(48) Mate Context(1,1)` | head Q components |
| `--smooth` | 30 | curvature prior across reported Q on head E's `QualityWindow` |
| `--window-l2` | 0.01 | L2 precision on `QualityWindow` weights (the other weights use 1) |
| `--min-overlap` | 20 | minimum comparable overlap bases |
| `--iterations` | 30 | maximum EM iterations |
| `--max-pairs` | all | stop after this many pairs |

## What it identifies

- **Substitutions only.** The spec records `identified_ops: ["substitution"]`, and its indel logits are −∞, so
  reads generated from it have no indels.
- **Sequencer errors only.** Errors both mates share (PCR, library, cluster) cancel.
- **Head E, per mate row:** substitution by true base, two-sided context, centre and neighbouring Q of that
  mate, cycle position, mate and strand. Position effects are fitted conditional on the positions the overlap
  covers.
- **Which mate erred.** When mates disagree, soft EM with head E over both mates' rows attributes the error.
  It starts from a coin flip, never from the higher Q.
- **Head Q** is fitted from every base of every pair, placed or not, conditioned on observed bases. Fitting it
  only on kept overlaps would bias it, because the dropped pairs are the low-Q ones.

Reported quality is a feature, never a label: placement uses bases only, and no rate is derived from Q.

## Placement

1. Each pair's gapless offset is the one with the lowest mismatch fraction over at least `min_overlap`
   comparable bases, then the longest. Adapter read-through needs no adapter list.
2. The pair is **ambiguous**, and dropped, if any offset more than 3 away also has a mismatch fraction below
   0.35. This catches repeats and low-complexity sequence. The threshold is absolute: a margin over the best
   offset would reject error-rich true overlaps.
3. An **indel** is implied, and the pair dropped, if splitting the overlap and shifting one side by ≤ 3
   removes ≥ 3 mismatches, and the split alignment passes the mismatch cap.
4. Otherwise the pair is dropped if the best offset has more than 20% mismatches. The cap is loose on
   purpose: NGmerge's 10% dropped the error-rich pairs and biased rates down 8% in recovery.

## Fitting defaults and why

Overlap evidence has few errors at high Q. With a single L2 penalty on every weight, the rarest Q bins are
pulled toward the overall error rate, and high-Q rates come out several times too high. The defaults give
`QualityWindow` a weak L2 (`--window-l2 0.01`) plus a second-order random-walk prior across reported Q
(`--smooth 30`). The prior penalises curvature, so a sparse Q bin follows its neighbours' trend.

On a 20,000-pair 2×125 simulation with ~10% of pairs merging, mean |log2(fitted / true rate)|:

| window_l2 | smooth | all Q | Q ≥ 30 | marginal rate |
|---|---|---|---|---|
| 1 | 0 | 0.83 | 1.54 | 1.073× |
| 1 | 30 | 0.84 | 1.68 | 1.082× |
| 0.01 | 0 | 1.44 | 1.10 | 1.017× |
| **0.01** | **30** | **0.60** | **0.97** | **1.022×** |

Either change alone doesn't help. A weak L2 alone collapses sparse bins, and smoothing alone can't move a
level that L2 holds. A first-difference prior was also tried; it flattened the curve.

Head Q's `Position(48)` default was chosen on held-out reads from a real MiSeq run. It cuts the median
per-position Q histogram TV from 0.177 (`Position(4)`) to 0.060. A second Markov lag and linear knot spacing
did not help.

## Choosing data

| library | use? | observed |
|---|---|---|
| Amplicons, short-insert libraries | yes | MiSeq 2×150 mock amplicons (SRR5240881): 97.7% of pairs placed; mismatch rate 0.01% where both mates have Q ≥ 30 |
| ~10% of pairs overlapping | yes, but high-Q rates are weakly identified | simulated 2×125, insert 294 ± 50: 10.5% placed; Q38–39 still ~5.5× too high |
| Long-insert genomic libraries | **no** | HiSeq 2×125 gut metagenome (ERR10889147): 3% placed, and ~0.3% of all pairs overlap divergent repeat copies longer than the read (1.6% false mismatches at Q ≥ 30) |

For long-insert data, use the `reference` mode (planned).

Aim for ~10% or more of pairs overlapping. For normal insert sizes with mean μ and sd σ and reads of length
L, the overlapping fraction is P(`min_overlap` ≤ insert ≤ 2L − `min_overlap`). For 2×125 reads and σ = 50,
μ ≈ 294 gives 10%.

## Cost

Reading pairs and building rows takes about 1 ms per pair. Fitting dominates: 20,000 MiSeq pairs took 12 min,
and 20,000 2×125 pairs 9.6 min, mostly head Q over every base. `--max-pairs` bounds a run.

## Known limits

- **Error head position.** With fixed-length reads, `Position` start and end splines are collinear and not
  identified. On SRR5240881 the fitted curve swings 0–6% while the raw disagreement rate is flat at ~0.5%.
- **Q38–39** rates stay ~5.5× too high when only ~10% of pairs merge. A real Q34 dip below its neighbours
  can't be followed by a smoothness prior.
- **Head Q** misfits narrow cycle-specific dips (cycles 40–41 on MiSeq, TV 0.32).
- **Nuisances not modelled.** Both mates erring to the same base is ignored. Divergent repeat copies longer
  than the read can't be told from true overlaps by bases alone.
