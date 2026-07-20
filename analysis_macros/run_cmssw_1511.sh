#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cmssw_src="${repo_dir}/../CMSSW_15_1_1/src"

if [[ ! -d "${cmssw_src}" ]]; then
  echo "Cannot find CMSSW_15_1_1/src next to ${repo_dir}" >&2
  exit 1
fi

cd "${cmssw_src}"
eval "$(scram runtime -sh)"
cd "${repo_dir}"

exec "$@"
