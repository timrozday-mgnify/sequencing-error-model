#!/usr/bin/env bash
#SBATCH --job-name=sem-pfam
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Nextflow head job for workflows/pfam. Submit from a run directory beside site.config (README.md):
#   cd <runs>/sem-pfam/p2 && sbatch <repo>/workflows/pfam/run.sh [nextflow options]   # e.g. --only ERR2984773
# Re-submitting resumes. Reads come from FIRE (the head job is at EBI); --fire false fetches the URLs instead. REPO defaults to the layout in the README.
set -euo pipefail
REPO=${REPO:-../../../sequencing-error-model}
[[ -f ../site.config ]] || { echo "submit from <runs>/sem-pfam/<run>/, next to ../site.config" >&2; exit 1; }
[[ -f $REPO/workflows/pfam/main.nf ]] || { echo "no checkout at $REPO: set REPO" >&2; exit 1; }
command -v nextflow >/dev/null || { echo "nextflow not on PATH" >&2; exit 1; }
export NXF_OPTS='-Xms1g -Xmx4g'
nextflow run "$REPO/workflows/pfam" -profile slurm,singularity -c ../site.config --fire --outdir results -resume "$@"
