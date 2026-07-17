#!/usr/bin/env bash
set -euo pipefail

fail() {
    echo "[generate_grid_opt_crab] ERROR: $*" >&2
    exit 1
}

TOTAL_JOBS="${TOTAL_JOBS:-200}"
BASE_OUTPUT="grid_PbPb5360_Notag_smoke.root"
JOBINFO="grid_PbPb5360_Notag_smoke_jobinfo.txt"
PACKAGE="grid_PbPb5360_Notag_smoke_crab_job.tar"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

[[ "$TOTAL_JOBS" =~ ^[0-9]+$ ]] || fail "Invalid TOTAL_JOBS: $TOTAL_JOBS"
[[ "$TOTAL_JOBS" -gt 0 ]] || fail "TOTAL_JOBS must be positive"

RAW_CRAB_JOB_ID="${1:-${CRAB_Id:-}}"
if [[ -z "$RAW_CRAB_JOB_ID" && -n "${_CONDOR_JOB_AD:-}" && -r "${_CONDOR_JOB_AD:-}" ]]; then
    RAW_CRAB_JOB_ID="$(grep -i '^CRAB_Id =' "$_CONDOR_JOB_AD" | tr -d '"' | tail -1 | awk '{print $NF}')"
fi
[[ "$RAW_CRAB_JOB_ID" =~ ^[0-9]+$ ]] || fail "Cannot determine numeric CRAB job id from CRAB_Id=${CRAB_Id:-unset}"
[[ "$RAW_CRAB_JOB_ID" -ge 1 ]] || fail "CRAB job id must be 1-based and positive: $RAW_CRAB_JOB_ID"
[[ "$RAW_CRAB_JOB_ID" -le "$TOTAL_JOBS" ]] || fail "CRAB job id $RAW_CRAB_JOB_ID is outside [1, $TOTAL_JOBS]"
JOB_ID=$((RAW_CRAB_JOB_ID - 1))

if [[ -f "./generate_grid_opt_job" ]]; then
    chmod +x "./generate_grid_opt_job"
    BINARY="./generate_grid_opt_job"
elif [[ -f "./build_cmssw_1511/generate_grid_opt_job" ]]; then
    chmod +x "./build_cmssw_1511/generate_grid_opt_job"
    BINARY="./build_cmssw_1511/generate_grid_opt_job"
else
    fail "Cannot find generate_grid_opt_job in the CRAB sandbox"
fi

printf -v JOB_PAD "%03d" "$JOB_ID"
SPLIT_ROOT="grid_PbPb5360_Notag_smoke_split${TOTAL_JOBS}_job${JOB_PAD}of${TOTAL_JOBS}.root"

rm -f "$BASE_OUTPUT" "$SPLIT_ROOT" "$JOBINFO" "$PACKAGE"

{
    echo "========================================"
    echo " generate_grid_opt CRAB split job"
    echo "========================================"
    echo "Host       : $(hostname)"
    echo "Date UTC   : $(date -u)"
    echo "CRAB_Id    : ${CRAB_Id:-unset}"
    echo "Raw job id : $RAW_CRAB_JOB_ID"
    echo "JOB_ID     : $JOB_ID"
    echo "TOTAL_JOBS : $TOTAL_JOBS"
    echo "Binary     : $BINARY"
    echo "Base output: $BASE_OUTPUT"
    echo "Split root : $SPLIT_ROOT"
    echo "========================================"

    "$BINARY" \
        PbPb 5360 NoTag "$BASE_OUTPUT" \
        400 100 1 4 1 100 400 1 \
        "$JOB_ID" "$TOTAL_JOBS"

    [[ -s "$SPLIT_ROOT" ]] || fail "Expected split output was not produced: $SPLIT_ROOT"

    echo "========================================"
    echo "Produced: $SPLIT_ROOT"
    echo "Size    : $(stat -c '%s' "$SPLIT_ROOT") bytes"
    echo "SHA256  : $(sha256sum "$SPLIT_ROOT" | awk '{print $1}')"
    echo "========================================"
} 2>&1 | tee "$JOBINFO"

tar -cf "$PACKAGE" "$SPLIT_ROOT" "$JOBINFO"
[[ -s "$PACKAGE" ]] || fail "Failed to create CRAB output package: $PACKAGE"

cat > FrameworkJobReport.xml <<'EOF'
<FrameworkJobReport>
  <ReadBranches>
  </ReadBranches>
  <PerformanceReport>
  </PerformanceReport>
  <GeneratorInfo>
  </GeneratorInfo>
  <RunReport>
  </RunReport>
  <ExitCode Value="0"/>
</FrameworkJobReport>
EOF
