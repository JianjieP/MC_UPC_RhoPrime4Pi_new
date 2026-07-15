#!/usr/bin/env bash

set -euo pipefail

VERSION="v0p1"
ARCH="${SCRAM_ARCH:-unknown_arch}"
CMSSW="${CMSSW_VERSION:-unknown_cmssw}"

NAME="rhoprime_PbPb5360_NoTag_${VERSION}_${ARCH}_${CMSSW}"
OUTPUT="${PWD}/${NAME}.tgz"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAGE="$(mktemp -d)"

cleanup() {
    rm -rf "$STAGE"
}
trap cleanup EXIT

mkdir -p \
    "$STAGE/bin" \
    "$STAGE/data" \
    "$STAGE/cards" \
    "$STAGE/metadata"

install -m 755 \
    "$PROJECT_ROOT/packaging/runcmsgrid.sh" \
    "$STAGE/runcmsgrid.sh"

install -m 755 \
    "$PROJECT_ROOT/build/rhoprime_lhe" \
    "$STAGE/bin/rhoprime_lhe"

install -m 644 \
    "$PROJECT_ROOT/grids/grid_PbPb5360_NoTag_prod.root" \
    "$STAGE/data/grid_PbPb5360_NoTag_prod.root"

install -m 644 \
    "$PROJECT_ROOT/cards/production/rhoprime_PbPb5360_NoTag.card" \
    "$STAGE/cards/rhoprime_PbPb5360_NoTag.card"

git -C "$PROJECT_ROOT" rev-parse HEAD \
    > "$STAGE/metadata/git_commit.txt"

{
    echo "SCRAM_ARCH=${SCRAM_ARCH:-unset}"
    echo "CMSSW_VERSION=${CMSSW_VERSION:-unset}"
    echo "Compiler=$(c++ --version | head -1)"
    echo "ROOT=$(root-config --version)"
} > "$STAGE/metadata/build_info.txt"

cat > "$STAGE/README_RUNTIME.md" <<'EOF'
Run:
  ./runcmsgrid.sh NEVENTS RANDOM_SEED NCPU

Output:
  cmsgrid_final.lhe
EOF

(
    cd "$STAGE"

    sha256sum \
        bin/rhoprime_lhe \
        data/grid_PbPb5360_NoTag_prod.root \
        cards/rhoprime_PbPb5360_NoTag.card \
        metadata/git_commit.txt \
        metadata/build_info.txt \
        > SHA256SUMS
)

# 注意：打包 STAGE 里面的内容，而不是 STAGE 目录本身
tar -czf "$OUTPUT" -C "$STAGE" .

echo "Created:"
echo "$OUTPUT"

echo
echo "Tarball content:"
tar -tzf "$OUTPUT"