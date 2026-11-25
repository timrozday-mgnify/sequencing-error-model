# workflows/pfam: P2, the pfam mode against truth (Slurm + Singularity)

Runs phase P2 of [docs/pfam_mode_plan.md](../../docs/pfam_mode_plan.md): seven libraries whose genomes are known, and the
go / no-go on the `pfam` mode. You run it; nothing here is submitted by an agent.

| Run | What | Pairs kept |
|---|---|---|
| ERR2984773, SRR12324253, SRR12829160 | ZymoBIOMICS D6300 / D6305 even mix, 151, 150, 151 bp | 4.4, 4.7, 4.5 M |
| ERR2935805 | D6310 log mix, 101 bp | 4.8 M |
| SRR23926920, SRR32629731 | D6331 gut standard, 150 and 99 bp | 4.7, 3.6 M |
| SRR24523812 | the isolate of phases 5 and 6, against its shovill assembly | 4.0 M |

Every library gives undosed reads and two standard additions (+0.002 and +0.01 errors per base, added before ORF calling). Each is searched
with kfpfam's Mode A and HMMER (E <= 1000), then fitted: both anchors at each dose, and each anchor alone undosed. The undosed reads are aligned to the
genomes with minimap2, and `pfam_truth score` classifies every counted stop and compares every fit with the genomes. `VERDICT` applies the rule
(plan §7, P2).

## What it needs on the HPC

Most of it exists from earlier work:

| Input | From |
|---|---|
| `sem-tools.sif` | `workflows/hpc` (phase 6b). Rebuild only if `uv.lock` changed: `bash workflows/hpc/build.sh --sif` |
| `kfpfam-filter.sif`, `prep/`, `index/k4T6.5.tidx` (Pfam 39.0) | B11's `kfpfam_setup.sh`. The image must hold pyhmmer (kfpfam `3ccf01c` or later): see step 2 |
| `Pfam-A.hmm` 39.0 | B11's `pfam_setup.sh`: the file `prep` was built from |
| `recruit/pairs.py` | a kfpfam checkout (bound as a file; it is not in the image) |
| `truth/refs/{std,d6331}/refs.fna` | B9's `truth.sh refs` |
| `SRR24523812_contigs.fa` | `~/Downloads/SRR24523812_shovill/run/contigs.fa`, already copied for phase 6b |

The local prototype used Pfam 38.2; this uses 39.0. `pfam.read_hits` refuses hits whose accessions are not in `prep`.

## Steps

```bash
HPC=codon-login
BASE=/hps/nobackup/<group>/<user>
RUNS=$BASE/sequencing-error-model_runs/sem-pfam
```

1. Update the checkout on HPC (the image holds no repo code; it is bound from here):

   ```bash
   ssh $HPC "git -C $BASE/sequencing-error-model pull"
   ```

2. Make the run directory and the site config, then edit the config's `CHANGE` paths, queue and account:

   ```bash
   ssh $HPC "mkdir -p $RUNS/p2 && cp $BASE/sequencing-error-model/workflows/pfam/site.config $RUNS/site.config"
   ```

   Two checks on what it names, from a node with Singularity:

   ```bash
   singularity exec kfpfam-filter.sif python3 -c "import pyhmmer; print(pyhmmer.__version__)"
   ```

   ```bash
   singularity exec sem-tools.sif minimap2 --version
   ```

3. **A trial first**, one library and 20,000 pairs. It should finish in minutes and end with a `verdict.md` that says no power:

   ```bash
   cd $RUNS/p2 && sbatch $BASE/sequencing-error-model/workflows/pfam/run.sh --only ERR2984773 --max_pairs 20000 --outdir trial
   ```

   Check in `trial/score/ERR2984773.json`: `stats.reads_aligned` close to `stats.reads_exposed`, and `stops.unaligned` small.

4. **The run:**

   ```bash
   cd $RUNS/p2 && sbatch $BASE/sequencing-error-model/workflows/pfam/run.sh
   ```

   Re-submitting resumes. `--only A,B` runs some libraries; the verdict then counts those only.

Nextflow may warn about an `undefined parameter` (`container`, `filter_container`); it is harmless, the images named in `site.config` are used.

If the head job cannot reach ENA (`Can't stage file`; see [workflows/hpc/README.md](../hpc/README.md)), fetch the FASTQs where that works and set
`params.fastq_dir` to their directory. The names are ENA's, and only the files named in `runs.csv` are read.

## Cost

From the laptop: 1.2 core-ms per read for the search, ORFs included. Seven libraries x three doses x about 9 M reads = 190 M reads, so 65 core-hours there. **kfpfam's
README measured 13 times its laptop figure on the HPC's x86 nodes** before its kernels were rewritten, so expect 100-900 core-hours of `SEARCH`. A fit took 1.8 GB and 3 minutes at 0.3 M
reads of a simulated mock (half of its reads place, against 15 % in ERR10889147), so about 35 GB and 90 minutes for 9 M: 64 GB is asked for. Scoring holds every exposed codon: 128 GB is asked
for, unmeasured at that size. A task killed for memory retries with more.
Work directories hold the SAM (about 5 GB a library) and the hit tables (1-2 GB): budget about 150 GB, removable after the run.

## Bring back

Megabytes: `results/verdict.*`, `results/score/*.json` and `results/fits/` (specs and their reports). Nothing else.

```bash
rsync -av $HPC:$RUNS/p2/results/ results-pfam-p2/
```

## On a laptop

```bash
nextflow run workflows/pfam -profile test --runs my_runs.csv -c my.config \
    --filter_bin ~/Documents/kfpfam/filter/target/release/filter --pyhmmer_python <venv with pyhmmer>/bin/python
```

`my_runs.csv` has the columns of `runs.csv` with local FASTQ paths, and `my.config` sets `prep`, `index`, `hmm`, `pairs_py` and `genomes` as
`site.config` does. `test` keeps 2,000 pairs and one dose.
