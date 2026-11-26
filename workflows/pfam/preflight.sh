#!/usr/bin/env bash
#SBATCH --job-name=sem-pfam-preflight
#SBATCH --output=%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Check everything P2 reads, then write $RUNS/site.config with the paths found. Changes nothing else, apart from
# a .fai beside each genome FASTA (pysam would otherwise write it from several tasks at once).
#
#   BASE=/hps/nobackup/<group>/<user> B9_DIR=... B11_DIR=... KF_DIR=... bash preflight.sh
#
# BASE     holds sequencing-error-model/ (this checkout) and sequencing-error-model_runs/
# B9_DIR   the B9 nobackup directory (truth/refs/{std,d6331}/refs.fna from B9's `truth.sh refs`)
# B11_DIR  the B11 nobackup directory (pfam/39.0/Pfam-A.hmm from B11's pfam_setup.sh)
# KF_DIR   the kfpfam/ directory B11's kfpfam_setup.sh wrote (beside that script): the image, prep/ and index/
# KFPFAM   a kfpfam checkout, for recruit/pairs.py (default ~/src/kfpfam, as B11 used)
# QUEUE, ACCOUNT  optional: Slurm partition and account, written into site.config
#
# Run it where `singularity` works (a compute node: `srun --pty bash`, or sbatch this file). It prints every check
# and ends with READY or with what to fix.
set -uo pipefail
: "${BASE:?set BASE}" "${B9_DIR:?set B9_DIR}" "${B11_DIR:?set B11_DIR}" "${KF_DIR:?set KF_DIR}"
REPO=${REPO:-$BASE/sequencing-error-model}
RUNS=${RUNS:-$BASE/sequencing-error-model_runs/sem-pfam}
SEM_HPC=${SEM_HPC:-$BASE/sequencing-error-model_runs/sem-hpc}
KFPFAM=${KFPFAM:-$HOME/src/kfpfam}
SIF=$SEM_HPC/containers/sem-tools.sif
FILTER_SIF=$KF_DIR/kfpfam-filter.sif
PREP=$KF_DIR/prep
INDEX=$KF_DIR/index/k4T6.5.tidx
HMM=$B11_DIR/pfam/39.0/Pfam-A.hmm
PAIRS_PY=$KFPFAM/recruit/pairs.py
STD=$B9_DIR/truth/refs/std/refs.fna
D6331=$B9_DIR/truth/refs/d6331/refs.fna
ISOLATE=$SEM_HPC/data/SRR24523812_contigs.fa

bad=0
ok() { echo "  ok    $*"; }
no() { echo "  FIX   $*"; bad=$((bad + 1)); }
run() { singularity exec -B "$REPO/src:/app/src" "$@"; }  # as the workflow runs it (nextflow.config)

echo "== tools"
command -v singularity >/dev/null && ok "singularity $(singularity --version 2>&1 | head -1)" || no "no singularity on PATH: run on a compute node"
if command -v nextflow >/dev/null; then
    v=$(nextflow -version 2>&1 | grep -oE 'version [0-9]+\.[0-9]+' | head -1 | cut -d' ' -f2)
    [[ $(printf '%s\n25.04\n' "$v" | sort -V | head -1) == 25.04 ]] && ok "nextflow $v" || no "nextflow $v: need 25.04 or later"
else
    no "no nextflow on PATH (the head job needs it too)"
fi

echo "== code"
if git -C "$REPO" rev-parse >/dev/null 2>&1; then
    git -C "$REPO" fetch -q origin main 2>/dev/null || echo "  note  could not fetch origin: checking the local main only"
    here=$(git -C "$REPO" rev-parse HEAD) main=$(git -C "$REPO" rev-parse origin/main)
    [[ $(git -C "$REPO" branch --show-current) == main ]] || no "$REPO is not on main: git -C $REPO switch main"
    [[ $here == "$main" ]] && ok "$REPO at origin/main ${here:0:7}" || no "$REPO at ${here:0:7}, origin/main is ${main:0:7}: git -C $REPO pull"
    [[ -z $(git -C "$REPO" status --porcelain -- src workflows) ]] || no "$REPO has local changes under src/ or workflows/"
else
    no "no checkout at $REPO: set REPO"
fi

echo "== inputs"
for f in "$SIF" "$FILTER_SIF" "$PREP/models.tsv" "$PREP/offs.u32" "$PREP/profile.f32" "$INDEX" "$HMM" "$PAIRS_PY" "$STD" "$D6331" "$ISOLATE"; do
    [[ -s $f ]] && ok "$f" || no "missing $f"
done

echo "== the HMM is the one prep was built from"
if [[ -s $HMM && -s $PREP/models.tsv ]]; then
    if [[ -s $PREP/hmm.sha256 ]]; then
        [[ $(zcat -f "$HMM" | sha256sum | cut -d' ' -f1) == $(cat "$PREP/hmm.sha256") ]] \
            && ok "sha256 matches prep/hmm.sha256" || no "$HMM is not the file prep was built from (prep/hmm.sha256)"
    fi
    # what HMMER and prep must agree on: the same accessions and names, in the same order
    if cmp -s <(zcat -f "$HMM" | awk '$1 == "NAME" { n = $2 } $1 == "ACC" { split($2, a, "."); print a[1] "\t" n }') \
        <(awk -F'\t' 'NR > 1 { print $1 "\t" $2 }' "$PREP/models.tsv"); then
        ok "accessions and names match prep/models.tsv ($(($(wc -l < "$PREP/models.tsv") - 1)) models)"
    else
        no "$HMM and prep/models.tsv hold different models: another Pfam release (see B11's diagnose_183.sh)"
    fi
fi

echo "== images"
if command -v singularity >/dev/null; then
    [[ -s $SIF ]] && {
        run "$SIF" python -m sequencing_error_model.pfam_locus --help >/dev/null 2>&1 \
            && run "$SIF" python -m sequencing_error_model.pfam_truth --help >/dev/null 2>&1 \
            && ok "sem-tools runs this checkout's pfam_locus and pfam_truth" \
            || no "sem-tools cannot run this checkout's code: rebuild it (bash workflows/hpc/build.sh --sif)"
        v=$(run "$SIF" minimap2 --version 2>/dev/null) && ok "sem-tools minimap2 $v" || no "no minimap2 in sem-tools"
    }
    [[ -s $FILTER_SIF ]] && {
        v=$(singularity exec "$FILTER_SIF" python3 -c "import pyhmmer; print(pyhmmer.__version__)" 2>/dev/null) \
            && ok "kfpfam image pyhmmer $v" || no "no pyhmmer in $FILTER_SIF: kfpfam 3ccf01c or later (README step 2)"
        singularity exec "$FILTER_SIF" sh -c 'command -v filter' >/dev/null 2>&1 && ok "kfpfam image has filter" \
            || no "no filter binary in $FILTER_SIF"
    }
fi
[[ -s $PAIRS_PY ]] && echo "  note  pairs.py at kfpfam $(git -C "$KFPFAM" log -1 --format=%h 2>/dev/null || echo '?'), sha256 $(sha256sum "$PAIRS_PY" | cut -c1-16) (the laptop rehearsal: 9833cf6, 9ff396bc4d4fd99b)"

echo "== genome indexes"
for pair in "std:$STD" "d6331:$D6331" "isolate:$ISOLATE"; do
    name=${pair%%:*} f=${pair#*:}
    [[ -s $f ]] || continue
    if [[ -s $f.fai && $f.fai -nt $f ]]; then
        ok "$name: $f.fai"
    elif command -v singularity >/dev/null &&
        singularity exec -B "$(dirname "$f")" "$SIF" python -c "import pysam, sys; pysam.faidx(sys.argv[1])" "$f" 2>/dev/null; then
        ok "$name: built $f.fai"
    else
        no "$name: cannot write $f.fai"
    fi
done

echo "== reads and disk"
url=$(awk -F, 'NR == 2 { print $3 }' "$REPO/workflows/pfam/runs.csv")
curl -sfI --max-time 30 "$url" >/dev/null && ok "ENA reachable from $(hostname)" \
    || echo "  note  ENA not reachable from $(hostname): if the head job's node is the same, fetch the FASTQs and set fastq_dir (README)"
free=$(df -BG --output=avail "$BASE" 2>/dev/null | tail -1 | tr -dc 0-9)
[[ ${free:-0} -ge 250 ]] && ok "${free} GB free under $BASE" || no "${free:-?} GB free under $BASE: the run needs about 250 GB (150 work, 80 keep)"

echo "== site.config"
mkdir -p "$RUNS/p2"
new=$(mktemp)
cat > "$new" <<EOF
// Written by workflows/pfam/preflight.sh on $(date -I) from BASE=$BASE.
params {
    container        = '$SIF'
    filter_container = '$FILTER_SIF'
    prep     = '$PREP'
    index    = '$INDEX'
    hmm      = '$HMM'
    pairs_py = '$PAIRS_PY'
    genomes  = [
        std    : '$STD',
        d6331  : '$D6331',
        isolate: '$ISOLATE',
    ]
    // fastq_dir = '/path/to/fastq'   // if the head job cannot reach ENA (README)
}
process.cache = 'lenient'
${QUEUE:+process.queue = '$QUEUE'}
${ACCOUNT:+process.clusterOptions = '--account=$ACCOUNT'}
singularity.cacheDir = '../singularity_cache'
// Also retry node failures and preemption.
process.errorStrategy = {
    if (task.previousException?.toString()?.contains('terminated for an unknown reason')) return 'retry'
    if (task.exitStatus in ((130..145) + 104) && task.attempt <= 2) return 'retry'
    return 'finish'
}
EOF
if [[ -s $RUNS/site.config ]] && ! cmp -s <(tail -n +2 "$RUNS/site.config") <(tail -n +2 "$new"); then
    echo "  note  $RUNS/site.config exists and differs; the new one is $RUNS/site.config.new:"
    diff <(tail -n +2 "$RUNS/site.config") <(tail -n +2 "$new") | sed 's/^/        /'
    mv "$new" "$RUNS/site.config.new"
else
    mv "$new" "$RUNS/site.config" && ok "$RUNS/site.config"
fi

echo
if ((bad)); then
    echo "NOT READY: $bad to fix (FIX lines above)"
    exit 1
fi
echo "READY. Trial:  cd $RUNS/p2 && sbatch $REPO/workflows/pfam/run.sh --only ERR2984773 --max_pairs 20000 --outdir trial"
