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
# Re-submitting resumes. REPO defaults to the layout in the README.
set -euo pipefail
REPO=${REPO:-../../../sequencing-error-model}
[[ -f ../site.config ]] || { echo "submit from <runs>/sem-pfam/<run>/, next to ../site.config" >&2; exit 1; }
[[ -f $REPO/workflows/pfam/main.nf ]] || { echo "no checkout at $REPO: set REPO" >&2; exit 1; }
command -v nextflow >/dev/null || { echo "nextflow not on PATH" >&2; exit 1; }
# With CONFIGS (a pipelines-nextflow-configs checkout, as for B11), reads come from FIRE with the FIRE_* Nextflow secrets.
fire=()
if [[ -n ${CONFIGS:-} ]]; then
    [[ -f $CONFIGS/profile-bases/codon/base.config ]] || { echo "no $CONFIGS/profile-bases/codon/base.config" >&2; exit 1; }
    fire=(-c "$CONFIGS/profile-bases/codon/base.config" -c "$REPO/workflows/pfam/fire.config")
fi
export NXF_OPTS='-Xms1g -Xmx4g'
nextflow run "$REPO/workflows/pfam" -profile slurm,singularity "${fire[@]}" -c ../site.config --outdir results -resume "$@"
