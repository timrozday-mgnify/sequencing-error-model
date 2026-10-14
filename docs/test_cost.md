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
| `test_aligner_bias[minimap2-1.0-False]` | 56 s | **35 s in `fit.error.objective`** (544 calls) + **25 s in `generate.realign`** (200 calls, 127 ms per read). minimap2 itself: one `fork_exec`, negligible |

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

1. **Cap the inner optimiser on intermediate EM refits** (generalised EM). `sources.pe_overlap.fit` and
   `fit.kmer.fit` warm-start each refit (`error.fit(init=...)`) but then run L-BFGS to full convergence every
   time - ~89 objective evaluations per refit, 14 refits, 1,247 evaluations for one test. Only the *last*
   refit needs to be converged; the outer loop re-optimises whatever the intermediate ones leave. An
   `error.fit(maxiter=...)` passed on all but the final iteration should cut evaluations several-fold.
   **Projected, not measured**: capping intermediates at ~15 evaluations would be ~5x fewer, and EM may need
   a few more outer iterations in exchange. This speeds up real fits too, and weakens no test.
2. **Speed up `generate.realign`** - 127 ms per read is a Python inner loop over the clip (`generate.py:343`
   is a genexpr called 1.2 M times in one test). Worth one pass of vectorisation; same argument as above,
   it is runtime cost, not test cost.
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
