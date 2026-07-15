#!/usr/bin/env bash

set -euo pipefail

fail() {
    echo "[make_tarball] ERROR: $*" >&2
    exit 1
}

VERSION="v0p1"
ARCH="${SCRAM_ARCH:-unknown_arch}"
CMSSW="${CMSSW_VERSION:-unknown_cmssw}"

NAME="rhoprime_PbPb5360_NoTag_${VERSION}"
OUTPUT="${PWD}/${NAME}.tgz"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-${PROJECT_ROOT}/build}"
BINARY="${BUILD_DIR}/rhoprime_lhe"
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

[[ -x "$BINARY" ]] || fail "Generator executable not found: $BINARY"

if readelf -d "$BINARY" | grep -q '/usr/lib64/root'; then
    fail "Refusing to package $BINARY because it is linked to system ROOT (/usr/lib64/root). Rebuild in cmsenv and/or set BUILD_DIR to a CMSSW ROOT build."
fi

install -m 755 \
    "$PROJECT_ROOT/packaging/runcmsgrid.sh" \
    "$STAGE/runcmsgrid.sh"

install -m 755 \
    "$BINARY" \
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
    echo "BUILD_DIR=${BUILD_DIR}"
    echo "Compiler=$(c++ --version | head -1)"
    echo "ROOT=$(root-config --version)"
    echo "BinaryRpath=$(readelf -d "$BINARY" | grep -E 'RPATH|RUNPATH' || true)"
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