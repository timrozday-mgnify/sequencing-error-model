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
  The 155 rows are phase 6b item 2: seeds 0-9 at k 11-31 (`pool_*`), a 200 kb genome (`big_*`), the
  `Latent(2)`/`Latent(3)` truths at k 11-31 (`latent*`) and skiver's `-l` floor at 10/30/100x (`lb*`). A row whose tolerances fail
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

Set these once per shell, on both sides:

```bash
HPC=codon-login                                   # your ssh host for the cluster's login node
BASE=/hps/nobackup/<group>/<user>                 # shared storage the compute nodes can read
RUNS=$BASE/sequencing-error-model_runs/sem-hpc
```

1. **Build the image on the laptop** (Docker, amd64; `--sif` also writes `sem-tools.tar`, since there is no
   Singularity here to convert it):

   ```bash
   cd /Users/timrozday/Documents/sequencing-error-model && bash workflows/hpc/build.sh --sif
   ```

   Rebuild only when `uv.lock` or a tool version in the `Dockerfile` changes; the tag is their hash, so a stale
   `.sif` is easy to spot. The repo's own code is *not* in the image: it is bound in from the checkout at run time.

2. **Make the layout and copy over** the tar, the assembly and the site config:

   ```bash
   ssh $HPC "mkdir -p $RUNS/containers $RUNS/data $RUNS/phase6b"
   scp sem-tools.tar $HPC:$RUNS/containers/
   scp ~/Downloads/SRR24523812_shovill/run/contigs.fa $HPC:$RUNS/data/SRR24523812_contigs.fa
   scp hpc/site.config $HPC:$RUNS/site.config
   ```

3. **Build the .sif on HPC** and drop the tar:

   ```bash
   ssh $HPC
   cd $RUNS/containers
   module load singularity            # or apptainer, per the site
   singularity build sem-tools.sif docker-archive://sem-tools.tar && rm sem-tools.tar
   ```

   `site.config` already points at `../containers/sem-tools.sif` and `../data/SRR24523812_contigs.fa`, so the
   names above matter. Edit it for the site's `process.queue` / `process.clusterOptions` before submitting.

4. **Clone the repo** next to the runs tree (no venv needed; `git pull` is how code updates reach a run):

   ```bash
   git clone git@github.com:timrozday-mgnify/sequencing-error-model.git $BASE/sequencing-error-model
   ```

5. **Submit the head job** from the run directory:

   ```bash
   cd $RUNS/phase6b
   module load nextflow               # run.sh refuses to start without it on PATH
   sbatch ../../../sequencing-error-model/hpc/run.sh
   ```

   Extra options pass through to Nextflow: `--real false` for the sweep alone, `--sweep ''` for the real arm
   alone, `--k 21`. Re-submitting resumes.

The head job stages the two FASTQs from ENA (~1.5 GB) once into `work/stage-*`, so compute nodes need no
internet. If the head job's node has none either, download them on a login node and pass `--r1`/`--r2` paths.
