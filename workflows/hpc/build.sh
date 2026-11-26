#!/usr/bin/env bash
# Build the workflow's image (tools + locked Python deps; the code comes from the repo checkout at run time).
#   bash workflows/hpc/build.sh                 # docker image, for local runs
#   bash workflows/hpc/build.sh --sif           # also sem-tools.tar: copy it to HPC and run there
#                                               #   singularity build sem-tools.sif docker-archive://sem-tools.tar
#   bash workflows/hpc/build.sh --push REGISTRY # also push REGISTRY/sem-tools:<tag>
# Rebuild only when uv.lock or a tool version in the Dockerfile changes. The tag is uv.lock's hash plus the
# Dockerfile's, so a stale .sif is easy to spot.
set -euo pipefail
repo="$(cd "$(dirname "$0")/../.." && pwd)"
tag=$(cat "$repo/uv.lock" "$repo/workflows/hpc/Dockerfile" | shasum | cut -c1-10)
image="sem-tools:$tag"
docker build --platform linux/amd64 -t "$image" -f "$repo/workflows/hpc/Dockerfile" "$repo"
echo "built $image"
case "${1:-}" in
    --sif)
        docker save "$image" -o sem-tools.tar
        sing=$(command -v singularity || command -v apptainer || true)
        if [ -n "$sing" ]; then "$sing" build sem-tools.sif docker-archive://sem-tools.tar && rm sem-tools.tar; fi
        ;;
    --push)
        docker tag "$image" "$2/$image"
        docker push "$2/$image"
        ;;
esac
