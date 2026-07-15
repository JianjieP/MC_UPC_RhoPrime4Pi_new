#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR=${BUILD_DIR:-build_cmssw}

./${BUILD_DIR}/generate_bose build/grid_PbPb5360_Notag_smoke_split20_merged.root build/events_nonbose_5000000.root 5000000 12345 500 0
./${BUILD_DIR}/analyze_nonbose build/events_nonbose_5000000.root build/analyze_nonbose_5000000.root 5000000

./${BUILD_DIR}/compare_bose_modulations build/events_bose_5000000.root build/events_nonbose_5000000.root --no-bell
