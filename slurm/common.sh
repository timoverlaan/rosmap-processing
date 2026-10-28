# Shared setup for the job scripts in slurm/. Source it from a job script
# submitted from the repo root:
#
#   source slurm/common.sh
#
# Provides: CONTAINER, V2, SEAAD_H5AD, in_container, rp, skip_if_exists.
# Set EXTRA_BINDS (space-separated --bind options) before calling in_container
# to mount more paths.

set -euo pipefail
/usr/bin/scontrol show job -d "${SLURM_JOB_ID:-}" 2>/dev/null || true

CONTAINER="${CONTAINER:-./container_pixi_0-2-1.sif}"
V2="data/v2"

# SeaAD counts. run_all_seaad.sh writes a column-mapped copy into data/v2;
# without it, the original download (already column-mapped in place in 2025)
# is used.
SEAAD_ORIGINAL="data/seaAD/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"
SEAAD_V2="$V2/counts/SeaAD/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"
if [ -z "${SEAAD_H5AD:-}" ]; then
    if [ -f "$SEAAD_V2" ]; then SEAAD_H5AD="$SEAAD_V2"; else SEAAD_H5AD="$SEAAD_ORIGINAL"; fi
fi

# Network storage needs a Kerberos ticket, which the job only gets if there
# was a valid one at submit time. Stop here rather than mid-run.
if ! klist -s 2>/dev/null; then
    echo "no valid Kerberos ticket in this job -- before sbatch: klist; kinit -R"
    exit 1
fi
[ -f "$CONTAINER" ] || { echo "missing container $CONTAINER"; exit 1; }
mkdir -p slurm/out "$V2"

# The container has no .git, so record the commit here and pass it in.
GIT_COMMIT="$(git rev-parse HEAD 2>/dev/null || echo unknown)"
GIT_BRANCH="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
GIT_DIRTY="$( (git status --porcelain 2>/dev/null || true) | wc -l | tr -d ' ')"
echo "code: commit $GIT_COMMIT on $GIT_BRANCH, $GIT_DIRTY uncommitted change(s)"

# A node-local /tmp for the container (--containall's own is small).
JOBTMP="/tmp/$USER/rosmap_${SLURM_JOB_ID:-$$}"
mkdir -p "$JOBTMP"
trap 'rm -rf "$JOBTMP"' EXIT

EXTRA_BINDS="${EXTRA_BINDS:-}"

# Run a command inside the container, from /opt/app (the repo's src/ and data/
# are mounted there).
in_container() {
    # shellcheck disable=SC2086
    apptainer exec --writable-tmpfs --pwd /opt/app --containall \
        --bind "$JOBTMP":/tmp \
        --bind src/:/opt/app/src/ \
        --bind data/:/opt/app/data/ \
        $EXTRA_BINDS \
        --env PYTHONPATH=/opt/app/src \
        --env GIT_COMMIT="$GIT_COMMIT" \
        --env GIT_BRANCH="$GIT_BRANCH" \
        --env GIT_DIRTY="$GIT_DIRTY" \
        --env SLURM_JOB_ID="${SLURM_JOB_ID:-none}" \
        "$CONTAINER" "$@"
}

# python -m rosmap_processing ... inside the container
rp() {
    in_container pixi run python -m rosmap_processing "$@"
}

# skip_if_exists OUTPUT cmd args...: run the command unless OUTPUT exists, so a
# failed job can be resubmitted without redoing finished steps.
skip_if_exists() {
    local out="$1"; shift
    if [ -e "$out" ]; then
        echo "exists, skipping: $out"
    else
        "$@"
    fi
}

# Inputs and gene maps (gene maps from slurm/gene_maps.sh)
MIT_COUNTS="$V2/counts/ROSMAP_MIT/combined.h5ad"
GREEN_COUNTS="$V2/counts/ROSMAP/combined.h5ad"
REF="$V2/reference"
OUT="$V2/processed"
# Green 2024 has no projid of its own; its individualID is complete.
GREEN_DONOR=individualID

# _mit NAME SHARED_GENES [args...]: ROSMAP-MIT, genes chosen within SHARED_GENES
_mit() {
    local name="$1" shared="$2"; shift 2
    echo "======================================== ROSMAP-MIT: $name"
    skip_if_exists "$OUT/rosmap_mit_${name}.h5ad" \
        rp pipeline scanpy "$MIT_COUNTS" "$OUT/rosmap_mit_${name}.h5ad" \
            --gene-map "$REF/rosmap_mit_gene_map.tsv" --restrict-genes "$shared" \
            --k-neighbors 30 "$@"
}

# _aligned seaad|green NAME: exactly the genes of rosmap_mit_NAME, same order
# (fails if any is missing)
_aligned() {
    local dataset="$1" name="$2"
    local ref_h5ad="$OUT/rosmap_mit_${name}.h5ad"
    echo "======================================== $dataset: $name"
    if [ "$dataset" = seaad ]; then
        skip_if_exists "$OUT/seaad_${name}.h5ad" \
            rp pipeline scanpy "$SEAAD_H5AD" "$OUT/seaad_${name}.h5ad" --layer UMIs \
                --gene-map "$REF/seaad_gene_map.tsv" --import-genes "$ref_h5ad" --k-neighbors 30
    else
        skip_if_exists "$OUT/rosmap_green_${name}.h5ad" \
            rp pipeline scanpy "$GREEN_COUNTS" "$OUT/rosmap_green_${name}.h5ad" \
                --gene-map "$REF/rosmap_green_gene_map.tsv" --donor-column "$GREEN_DONOR" \
                --import-genes "$ref_h5ad" --k-neighbors 30
    fi
}

# process_pair NAME [gene selection args...]
# ROSMAP-MIT -> rosmap_mit_NAME.h5ad (genes chosen among those SeaAD also has),
# then SeaAD -> seaad_NAME.h5ad with exactly the same genes.
process_pair() {
    local name="$1"; shift
    _mit "$name" "$REF/shared_gene_ids.txt" "$@"
    _aligned seaad "$name"
}

# process_triplet NAME [gene selection args...]
# Like process_pair, but genes are chosen among those all three datasets have,
# and Green 2024 -> rosmap_green_NAME.h5ad gets the same genes too.
process_triplet() {
    local name="$1"; shift
    _mit "$name" "$REF/shared_gene_ids_3way.txt" "$@"
    _aligned seaad "$name"
    _aligned green "$name"
}

# process_green NAME [gene selection args...]
# Green 2024 on its own -> rosmap_green_NAME.h5ad, with its own HVGs.
process_green() {
    local name="$1"; shift
    echo "======================================== Green 2024: $name"
    skip_if_exists "$OUT/rosmap_green_${name}.h5ad" \
        rp pipeline scanpy "$GREEN_COUNTS" "$OUT/rosmap_green_${name}.h5ad" \
            --gene-map "$REF/rosmap_green_gene_map.tsv" --donor-column "$GREEN_DONOR" \
            --k-neighbors 30 "$@"
}
