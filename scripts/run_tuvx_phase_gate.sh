#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_tuvx_phase_gate.sh <phase> <history_nc> [--compare <history_nc>] [--fallback <history_nc>] [--j-vars <vars>] [--allow-missing-j]

Phases:
  0  LNOx-O3 baseline       -> nonnegative + verify_ox_conservation.py
  1  Solar geometry         -> phase 0 + night-jzero
  2  TUV-x coupling         -> phase 1 + transition-smooth + decomp-compare + fallback-compare

Options:
  --compare <history_nc>    Comparison file for phase 2 decomposition check.
  --fallback <history_nc>   Reference file for phase 2 fallback-compare check.
  --j-vars <vars>           Comma-separated photolysis variable names (default: j_no2).
  --allow-missing-j         Allow missing j variables for exploratory/manual runs.

Notes:
  - Phase 1/2 night and transition checks require `coszr` and the requested
    j-variable(s) unless `--allow-missing-j` is used.
  - Override checker path with CHECK_TUVX_PHASE_SCRIPT (or CHECK_TUVX_SCRIPT).
  - Override Ox checker path with VERIFY_OX_SCRIPT.
EOF
}

if [[ $# -lt 2 ]]; then
  usage
  exit 2
fi

phase="$1"
history_nc="$2"
shift 2

compare_nc=""
fallback_nc=""
j_vars="${J_VARS:-j_no2}"
allow_missing_j="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --compare)
      compare_nc="${2:-}"
      shift 2
      ;;
    --fallback)
      fallback_nc="${2:-}"
      shift 2
      ;;
    --j-vars)
      j_vars="${2:-}"
      shift 2
      ;;
    --allow-missing-j)
      allow_missing_j="true"
      shift
      ;;
    *)
      echo "[ERROR] unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
checker="${CHECK_TUVX_PHASE_SCRIPT:-${CHECK_TUVX_SCRIPT:-${script_dir}/check_tuvx_phase.py}}"
ox_checker="${VERIFY_OX_SCRIPT:-${script_dir}/verify_ox_conservation.py}"

if [[ ! -x "${checker}" ]]; then
  echo "[ERROR] checker script not executable: ${checker}" >&2
  exit 2
fi

if [[ ! -f "${ox_checker}" ]]; then
  echo "[ERROR] Ox checker script not found: ${ox_checker}" >&2
  exit 2
fi

j_missing_flags=()
if [[ "${allow_missing_j}" == "true" ]]; then
  j_missing_flags=(--allow-missing)
fi

common_nonnegative() {
  "${checker}" nonnegative -i "${history_nc}" --vars qNO,qNO2,qO3 --tol 1.0e-20
}

common_ox() {
  python3 "${ox_checker}" -i "${history_nc}" --max-drift-pct 0.01 --warn-drift-pct 0.1
}

common_night() {
  "${checker}" night-jzero -i "${history_nc}" \
    --coszr-var coszr \
    --j-vars "${j_vars}" \
    --night-threshold 0.0 \
    --abs-tol 1.0e-20 \
    "${j_missing_flags[@]}"
}

common_transition() {
  "${checker}" transition-smooth -i "${history_nc}" \
    --coszr-var coszr \
    --j-vars "${j_vars}" \
    --transition-coszr 0.08 \
    --day-coszr 0.20 \
    --dt-seconds 3.0 \
    --j-floor 1.0e-12 \
    --max-p99-jump 0.35 \
    --max-max-jump 0.75 \
    --max-p99-curvature 0.50 \
    "${j_missing_flags[@]}"
}

common_decomp() {
  if [[ -z "${compare_nc}" ]]; then
    echo "[ERROR] phase ${phase} requires --compare <history_nc> for decomp-compare" >&2
    exit 2
  fi
  "${checker}" decomp-compare \
    -a "${history_nc}" -b "${compare_nc}" \
    --vars qNO,qNO2,qO3 \
    --max-rel-l2 1.0e-3 \
    --max-final-mean-rel 1.0e-3
}

common_fallback() {
  if [[ -z "${fallback_nc}" ]]; then
    echo "[ERROR] phase ${phase} requires --fallback <history_nc> for fallback-compare" >&2
    exit 2
  fi
  "${checker}" fallback-compare \
    -a "${history_nc}" -b "${fallback_nc}" \
    --vars qNO,qNO2,qO3 \
    --max-rel-l2 1.0e-12 \
    --max-final-mean-rel 1.0e-12 \
    --max-abs-diff 1.0e-18
}

case "${phase}" in
  0)
    common_nonnegative
    common_ox
    ;;
  1)
    common_nonnegative
    common_ox
    common_night
    ;;
  2)
    common_nonnegative
    common_ox
    common_night
    common_transition
    common_decomp
    common_fallback
    ;;
  *)
    echo "[ERROR] unsupported phase '${phase}'" >&2
    usage
    exit 2
    ;;
esac

echo "[PASS] phase ${phase} checks completed."
