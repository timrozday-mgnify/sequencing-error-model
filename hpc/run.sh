#!/usr/bin/env bash
#SBATCH --job-name=sem-phase6b
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Nextflow head job for workflows/hpc. Submit from a run directory <runs>/sem-hpc/<run>/ (README.md):
#   cd <runs>/sem-hpc/phase6b && sbatch ../../../sequencing-error-model/hpc/run.sh [nextflow options]
# e.g. `--real false` for the sweep alone, `--sweep ''` for the real run alone. Re-submitting resumes.
set -euo pipefail
REPO=../../../sequencing-error-model  # <runs>/../sequencing-error-model
[[ -f ../site.config ]] || { echo "submit from <runs>/sem-hpc/<run>/, next to ../site.config" >&2; exit 1; }
command -v nextflow >/dev/null || { echo "nextflow not on PATH" >&2; exit 1; }
export NXF_OPTS='-Xms1g -Xmx4g'
nextflow run "$REPO/workflows/hpc" -profile slurm,singularity -c ../site.config --outdir results -resume "$@"
