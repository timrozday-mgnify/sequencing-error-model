#!/usr/bin/env bash
#SBATCH --job-name=sem-pfam-diagnose
#SBATCH --output=%x.%A_%a.log
#SBATCH --array=1-7
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=160G
#SBATCH --export=ALL
# P2's diagnosis (docs/pfam_mode_plan.md §7 P2, "Next"), one library of runs.csv per array task, on results/keep/:
#   docs/pfam_overlap_diag.py  what the mate-overlap anchor can see: errors in one mate or both, by Q and insert size
#   docs/pfam_p2_cells.py      the locus label and the within-locus contrast per op, Q, mate and position, against truth
# With the exports of README step 2 (BASE, B9_DIR, KF_DIR; REPO, SEM_HPC as the preflight), from <runs>/sem-pfam/p2:
#   cd $RUNS/p2 && sbatch $REPO/workflows/pfam/diagnose.sh        # or --array=2 for one library
# Writes diagnose/<run>.overlap.json, <run>.cells.json and <run>.cells.npz (the per-base table, reused on a rerun;
# a table from before its truth was masked is rebuilt). CELLS_ONLY=1 skips the overlap diagnostic.
set -euo pipefail
: "${BASE:?set BASE}" "${B9_DIR:?set B9_DIR}" "${KF_DIR:?set KF_DIR}"
REPO=${REPO:-$BASE/sequencing-error-model}
SEM_HPC=${SEM_HPC:-$BASE/sequencing-error-model_runs/sem-hpc}
SIF=$SEM_HPC/containers/sem-tools.sif
[[ -d results/keep ]] || { echo "submit from <runs>/sem-pfam/p2, beside results/keep" >&2; exit 1; }
row=$(awk -F, -v i="${SLURM_ARRAY_TASK_ID:?run with sbatch}" 'NR == i + 1' "$REPO/workflows/pfam/runs.csv")
run=${row%%,*} set=${row##*,}
case $set in
    std) genomes=$B9_DIR/truth/refs/std/refs.fna ;;
    d6331) genomes=$B9_DIR/truth/refs/d6331/refs.fna ;;
    isolate) genomes=$SEM_HPC/data/SRR24523812_contigs.fa ;;
    *) echo "$run: no genomes for '$set'" >&2; exit 1 ;;
esac
keep=results/keep/$run out=diagnose
mkdir -p "$out"
py() { singularity exec -B "$REPO/src:/app/src" -B "$BASE" -B "$B9_DIR" -B "$KF_DIR" "$SIF" python "$@"; }
echo "$run against $genomes"
[[ ${CELLS_ONLY:-0} == 1 ]] || py "$REPO/docs/pfam_overlap_diag.py" "$keep/aligned.sam" "$genomes" --output "$out/$run.overlap.json"
py "$REPO/docs/pfam_p2_cells.py" "$keep/$run.kept_1.fastq.gz" "$keep/$run.kept_2.fastq.gz" --pairs "$keep/0/pairs.tsv" \
    --prep "$KF_DIR/prep" --sam "$keep/aligned.sam" --reference "$genomes" --table "$out/$run.cells.npz" \
    --output "$out/$run.cells.json"
