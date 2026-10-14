# Where `task test` spends its time

Measured 2026-10-14 on an unloaded M-series laptop, `pytest -p no:randomly` (sequential, 176 tests,
**370 s**). The plan's testing strategy budgets the per-PR recovery layer at "under ~2 min", so the suite is
about **3x over its own budget**. Re-measure with `--durations=20` before and after any change here; wall
times move by ~30 % with background load (the same suite took 492 s while two fits were running).

## It is not the tests, it is two functions

Profiled per test with `cProfile`, restricted to the package:

| test | wall | where it goes |
|---|---|---|
| `test_pe_overlap_recovers_substitution_head` | 66 s | **59 s in `fit.error.objective`** (1,247 calls, 24 ms each, over 14 EM refits) |
| `test_aligner_bias[minimap2-1.0-False]` | 56 s | **35 s in `fit.error.objective`** (544 calls) + **25 s in `generate.realign`** (200 calls). minimap2 itself: one `fork_exec`, negligible |

Read those two numbers as shares, not as wall time: `cProfile` inflates a Python-heavy row loop far more than
it inflates numpy work, so it overstated `realign` (25 s profiled against 22 s measured by `perf_counter`, in
a test whose honest total is 70 s). **Always A/B a change on the same machine in the same session** - the
absolute times in this document drift by ~25 % between sessions.

So the suite's cost is almost entirely two hot spots, both of them library code on the critical path of real
fits as well:

1. **`fit.error.objective`** - head E's L-BFGS objective and gradient. It dominates every slow test
   (`test_pe_overlap*`, `test_fit_error`, `test_sites`, `test_compare`, `test_latent`, `test_kmer_source`),
   because all of them fit head E. Already vectorised and chunked; the cost is the *number of evaluations*,
   not the per-call work.
2. **`generate.realign`** - the phase 5 soft-clip realignment, pure Python at 127 ms per read. It only shows
   up where reads are long and clipped (the `minimap2` aligner row), but it is paid by every real
   `reference`-mode fit run with `--unclip`, so it is not a test-only cost.

Scaffolding, fixtures, FASTQ I/O, subprocesses and the external aligners are all noise by comparison.

## What to do, by payoff over risk

1. ~~**Cap the inner optimiser on intermediate EM refits** (generalised EM).~~ **Tried, measured, rejected**
   (2026-10-14). `sources.pe_overlap.fit` and `fit.kmer.fit` warm-start each refit (`error.fit(init=...)`)
   but run L-BFGS to convergence every time - ~89 objective evaluations per refit, 15 refits, 1,323
   evaluations for one test - so capping all but the last looked like a several-fold win. It is not: the
   outer loop simply needs proportionally more iterations.

   | `maxiter` on intermediate refits | wall | `error.fit` calls | objective evaluations |
   |---|---|---|---|
   | none (today) | 85 s | 15 | 1,323 |
   | 20 | 92 s | 31 | 1,110 |
   | 50 | 90 s | 23 | 1,192 |
   | 100 | 80 s | 14 | 1,214 |

   16 % fewer evaluations for twice the E-steps, i.e. the total optimisation work is conserved and the extra
   E-steps make it slower. **Do not try this again.** What is left in this direction is smaller and real:
   `error.fit` rebuilds the design matrix from scratch on every refit (`_labelled` 5.4 s, `CountTable`
   validation 2 s of 82 s) although only the *counts* change between EM steps, so caching the rows and the
   matrix across refits is worth ~10 % with no change to the result.
2. ✓ **`generate.realign` done** (2026-10-14, 2.8x). It is a banded affine-gap DP whose row loop ran ~12
   numpy ops over the *whole* band, but with free template ends the band is `2w+1` wide while only the
   `m+1` cells with `0 <= col <= m` are ever valid - and `col` is contiguous along the band, so the valid
   part of a row is a slice. Computing only that slice, with a sentinel column so the `d+1` reads need no
   bounds check and the gap ramp hoisted out of the loop, gives **38.8 ms against 110.1 ms** per call on the
   aligner-bias workload (200 calls, 22.0 s to 7.8 s) and takes the minimap2 row from **69.7 s to 56.3 s**
   (A/B, same session). Verified exactly: 0 mismatches on 120 captured real calls and on 3,000 randomised
   small cases against the committed implementation.

   What is left here needs a decision, not a refactor: the DP is full because `free_template_ends` lets the
   read start anywhere in the window, but the caller always knows roughly where the aligner placed it. A band
   around that offset would be another large factor, and would change which alignments tie - which is exactly
   what `realign_observable`'s docstring warns about, so it is a semantics change and belongs in phase 5.
3. **Cut sample sizes that exceed their own stated requirement.**
   `test_pe_overlap_recovers_substitution_head` runs 10,000 pairs while its own comment says "Q37, the
   lowest-rate bin, needs 5,000 pairs here". **Measured: 5,000 pairs passes in 31 s against 66 s.** Before
   taking it, run it over several seeds - a test sitting exactly at its discrimination threshold trades
   runtime for flakiness, which is a worse deal than it looks.
4. **Split the genuinely end-to-end rows out of the per-PR suite.** `tests/test_aligners.py` (56 s for the
   minimap2 row) and `tests/test_sites.py`'s variation scenarios (28 + 14 + 13 s) are integration runs, and
   the plan's testing table already separates "every PR" from "dispatch". A `slow` marker plus a second task
   would need a CI change, since CI runs exactly `task lint` and `task test` - so this is a policy decision,
   not a refactor, and it should come *after* 1 and 2, which may make it unnecessary.
5. **Cheap and local**: `observations.CountTable.__post_init__` validates every key on construction (2 s over
   46 tables in one test), so intermediate tables pay O(rows) validation repeatedly. Validate once at the
   source boundary, or sample the keys, if it still shows up after 1 and 2.

## The slowest 20, for reference

```
66.31s  test_pe_overlap.py::test_pe_overlap_recovers_substitution_head
55.71s  test_aligners.py::test_aligner_bias[minimap2-1.0-False]
34.76s  test_compare.py::test_pe_overlap_vs_reference_estimates_library_errors
29.23s  test_fit_error.py::test_recovers_known_head
28.11s  test_sites.py::test_the_joint_model_removes_the_reference_modes_variant_inflation
21.70s  test_bias.py::test_pe_overlap_is_variant_immune_where_the_reference_mode_is_not
14.36s  test_sites.py::test_the_conservative_mask_removes_most_of_the_reference_modes_variant_inflation
13.20s  test_fit_error.py::test_select_keeps_true_effects_only
13.09s  test_latent.py::test_em_recovers_latent_classes
12.60s  test_sites.py::test_q_diagnostic_and_fit_follow_miscalibrated_q_through_variation
11.89s  test_kmer_source.py::test_level_from_rescales_the_fitted_level_and_is_recorded
 6.54s  test_kmer_source.py::test_per_mate_runs_identify_a_mate_term
 6.06s  test_kmer_source.py::test_cli_composes_a_spec_the_generator_draws_from
 5.64s  test_latent.py::test_compare_scores_a_latent_truth
 5.55s  test_recovery.py::test_cigar_mode_recovers_example_spec
 4.45s  test_indel.py::test_recovers_indel_lengths
 4.28s  test_mask_cost.py::test_mask_cost_on_clonal_reads
 4.18s  test_compare.py::test_cli_reports_both_levels
 3.64s  test_bias.py::test_grid_records_a_missing_skiver_and_still_tabulates
 2.77s  test_pe_overlap.py::test_cli_fits_spec_from_fastq
```

The 20 above are 340 s of the 370 s; the other 156 tests are 30 s between them. Nothing in the long tail is
worth touching.
