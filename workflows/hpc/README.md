# workflows/hpc: long runs on Slurm

Runs that take more than ~10 minutes locally go here. One Nextflow workflow, two arms:

- **real**: SRR24523812 (the phase 5/6 isolate), phase 6's subset (pairs 1,020,001-1,100,000). Aligns it to the
  shovill assembly with minibwa, then fits `reference` with `Latent(0/2/3)` and `pe-overlap`, runs skiver and fits
  `kmer` at each `--k` on the subset and on a 1,000,000-pair depth slice, and compares the three modes
  (`compare --skiver --kmer-spec`) per k. This is plan phase 6b item 1's rerun at k >= 17. Outputs:
  `results/specs/`, `results/kmer/<slice>_k<k>/`, `results/compare/<slice>_k<k>.json`.
- **sweep**: one Slurm task per row of `sweeps/phase6b.csv` (`id,module,args,truth`), running
  `python -m sequencing_error_model.<module> --skiver ... <args> [--spec <truth>] --output <id>.json` into
  `results/sweep/`. `truth` names a spec from the real arm (`latent2`, `latent3`); those rows wait for it.
  The 105 rows are phase 6b item 2: seeds 0-9 at k 11-31 (`pool_*`), a 200 kb genome (`big_*`) and skiver's
  `-l` floor at 10/30/100x (`lb*`). The `Latent` truth rows (`latent{2,3}_k{11..31}_s{0..4}`) are not in yet:
  `recovery --skiver --spec` refuses a `Latent` truth (its compare has no `latent` field) until that is
  implemented. A row whose tolerances fail
  still publishes its report; a row that crashes is skipped and logged.

The image (`Dockerfile`) holds skiver v0.3.2 (the release binary), minibwa, minimap2 and the locked Python
dependencies. The repo's own code is bound in from this checkout (`<repo>/src` over `/app/src`), so a `git pull` on HPC
updates it; rebuild only when `uv.lock` or a tool version changes.

## Locally

```bash
nextflow run workflows/hpc -profile test --r1 SRR24523812_1.fastq.gz --r2 SRR24523812_2.fastq.gz --assembly contigs.fa
```

`test` uses the local `.venv` and `skiver-0.3.2`, no container, on the real subset with 1,000 training pairs. `-profile docker` runs the
full workflow in the image instead.

## HPC

Layout, as in kmer-functional-profiler's `hpc/kfp-ablations`:

```
<base>/sequencing-error-model/              this repo (git clone; no venv needed)
<base>/sequencing-error-model_runs/sem-hpc/
    site.config                            copied from hpc/site.config, edited for the site
    containers/sem-tools.sif
    data/SRR24523812_contigs.fa            the shovill assembly
    phase6b/                               a run directory: results/, work/, the head job's log
```

1. On the laptop: `bash workflows/hpc/build.sh --sif`, then copy `sem-tools.tar` and the assembly over and,
   on HPC, `singularity build containers/sem-tools.sif docker-archive://sem-tools.tar`.
2. Submit the head job from the run directory:

   ```bash
   cd <base>/sequencing-error-model_runs/sem-hpc/phase6b && sbatch ../../../sequencing-error-model/hpc/run.sh
   ```

   Extra options pass through to Nextflow: `--real false` for the sweep alone, `--sweep ''` for the real arm
   alone, `--k 21`. Re-submitting resumes.

The head job stages the two FASTQs from ENA (~1.5 GB) once into `work/stage-*`, so compute nodes need no
internet. If the head job's node has none either, download them on a login node and pass `--r1`/`--r2` paths.
