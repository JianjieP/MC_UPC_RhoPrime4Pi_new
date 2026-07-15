#!/usr/bin/env bash

set -euo pipefail

fail() {
    echo "[rhoprime] ERROR: $*" >&2
    exit 1
}

if [[ $# -lt 3 ]]; then
    fail "Usage: $0 NEVENTS RANDOM_SEED NCPU"
fi

NEVENTS="$1"
RANDOM_SEED="$2"
NCPU="$3"

[[ "$NEVENTS" =~ ^[0-9]+$ ]] || fail "Invalid NEVENTS: $NEVENTS"
[[ "$RANDOM_SEED" =~ ^[0-9]+$ ]] || fail "Invalid RANDOM_SEED: $RANDOM_SEED"
[[ "$NCPU" =~ ^[0-9]+$ ]] || fail "Invalid NCPU: $NCPU"

# tarball 解压后的根目录
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BINARY="${HERE}/bin/rhoprime_lhe"
GRID="${HERE}/data/grid_PbPb5360_NoTag_prod.root"
CARD="${HERE}/cards/rhoprime_PbPb5360_NoTag.card"

OUTPUT="${PWD}/cmsgrid_final.lhe"
TMP_OUTPUT="${PWD}/cmsgrid_final.lhe.tmp"

[[ -x "$BINARY" ]] || fail "Generator executable not found: $BINARY"
[[ -r "$GRID" ]] || fail "Grid file not found: $GRID"
[[ -r "$CARD" ]] || fail "Card file not found: $CARD"

rm -f "$OUTPUT" "$TMP_OUTPUT"

echo "========================================"
echo " rho-prime LHE generation"
echo "========================================"
echo "Events : $NEVENTS"
echo "Seed   : $RANDOM_SEED"
echo "CPUs   : $NCPU"
echo "Grid   : $GRID"
echo "Card   : $CARD"
echo "Host   : $(hostname)"
echo "Date   : $(date -u)"
echo "========================================"

# 私有动态库存在时使用
if [[ -d "${HERE}/lib" ]]; then
    export LD_LIBRARY_PATH="${HERE}/lib:${LD_LIBRARY_PATH:-}"
fi

"$BINARY" \
    --grid "$GRID" \
    --config "$CARD" \
    --events "$NEVENTS" \
    --seed "$RANDOM_SEED" \
    --threads "$NCPU" \
    --output "$TMP_OUTPUT"

[[ -s "$TMP_OUTPUT" ]] || fail "Generator produced an empty LHE file"

# 检查实际 event 数
ACTUAL_EVENTS="$(
    grep -c '^[[:space:]]*<event>' "$TMP_OUTPUT" || true
)"

if [[ "$ACTUAL_EVENTS" -ne "$NEVENTS" ]]; then
    fail "Requested $NEVENTS events, but found $ACTUAL_EVENTS"
fi

# XML 完整性检查
if command -v xmllint >/dev/null 2>&1; then
    xmllint --stream --noout "$TMP_OUTPUT" \
        || fail "xmllint validation failed"
else
    echo "[rhoprime] WARNING: xmllint is unavailable"
fi

# 所有检查通过后再正式命名
mv "$TMP_OUTPUT" "$OUTPUT"

echo "Generated $ACTUAL_EVENTS events"
echo "Output: $OUTPUT"
echo "SHA256: $(sha256sum "$OUTPUT" | awk '{print $1}')"
echo "rho-prime generation completed successfully"