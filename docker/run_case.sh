#!/usr/bin/env bash

set -euo pipefail

MPAS_ROOT="${CHEMPAS_ROOT:-/mpas}"
DATA_ROOT="${CHEMPAS_DATA_ROOT:-/data/CheMPAS}"
MPI_RANKS=8
DOWNLOAD=1
FORCE_INIT="${CHEMPAS_FORCE_INIT:-0}"
PYTHON="${PYTHON:-python3}"

usage() {
    cat <<'EOF'
Usage:
  chempas-run-case CASE [options]

Cases:
  supercell-abba       Run the Supercell ABBA chemistry example
  supercell-lnox       Run the Supercell LNOx + O3 chemistry example
  chapman-nox-global   Run the Global Chapman + NOx chemistry example

Options:
  --setup-only         Prepare data/configuration but do not run MPAS
  --no-download        Do not download missing NCAR test-case tarballs
  --force-init         Regenerate init NetCDF files
  --help               Show this help text

Environment:
  CHEMPAS_ROOT         Source/runtime root in the image (default: /mpas)
  CHEMPAS_DATA_ROOT    Mounted data root (default: /data/CheMPAS)
  CHEMPAS_RUN_DURATION Override config_run_duration in copied namelists
  CHEMPAS_FORCE_INIT   Set to 1 to regenerate init files
EOF
}

log() {
    printf '[chempas] %s\n' "$*" >&2
}

die() {
    printf '[chempas] ERROR: %s\n' "$*" >&2
    exit 1
}

require_file() {
    local path="$1"
    [[ -f "${path}" ]] || die "required file not found: ${path}"
}

require_dir() {
    local path="$1"
    [[ -d "${path}" ]] || die "required directory not found: ${path}"
}

run_mpi() {
    local -a cmd=(mpiexec)
    if [[ "$(id -u)" == "0" ]]; then
        cmd+=(--allow-run-as-root)
    fi
    cmd+=(-n "${MPI_RANKS}" "$@")
    log "running: ${cmd[*]}"
    "${cmd[@]}" 2>&1 | tee run.out
}

copy_tree_contents() {
    local src="$1"
    local dst="$2"
    require_dir "${src}"
    mkdir -p "${dst}"
    cp -a "${src}/." "${dst}/"
}

download_case_data() {
    local case_name="$1"
    local required_file="$2"
    local case_dir="${DATA_ROOT}/${case_name}"
    local archive="${DATA_ROOT}/${case_name}.tar.gz"
    local url="https://www2.mmm.ucar.edu/projects/mpas/test_cases/v7.0/${case_name}.tar.gz"

    if [[ -f "${case_dir}/${required_file}" ]]; then
        return
    fi

    [[ "${DOWNLOAD}" == "1" ]] || die "missing ${case_dir}/${required_file} and downloads are disabled"

    mkdir -p "${DATA_ROOT}"
    log "downloading ${url}"
    curl -fL "${url}" -o "${archive}"
    log "extracting ${archive}"
    tar xzf "${archive}" -C "${DATA_ROOT}"
    rm -f "${archive}"

    require_file "${case_dir}/${required_file}"
}

install_executables() {
    local run_dir="$1"
    require_file "${MPAS_ROOT}/init_atmosphere_model"
    require_file "${MPAS_ROOT}/atmosphere_model"
    ln -sf "${MPAS_ROOT}/init_atmosphere_model" "${run_dir}/init_atmosphere_model"
    ln -sf "${MPAS_ROOT}/atmosphere_model" "${run_dir}/atmosphere_model"
}

install_physics_tables() {
    local run_dir="$1"
    local path
    require_dir "${MPAS_ROOT}/physics_tables"
    shopt -s nullglob
    for path in "${MPAS_ROOT}/physics_tables"/*; do
        [[ -f "${path}" ]] && ln -sf "${path}" "${run_dir}/$(basename "${path}")"
    done
    shopt -u nullglob
}

install_tuvx_data() {
    local run_dir="$1"
    require_dir "/usr/local/share/musica/tuvx_data"
    ln -sfnT "/usr/local/share/musica/tuvx_data" "${run_dir}/data"
}

ensure_stream_vars() {
    local run_dir="$1"
    shift
    local stream="${run_dir}/stream_list.atmosphere.output"
    local name
    require_file "${stream}"
    for name in "$@"; do
        grep -Fxq "${name}" "${stream}" || printf '%s\n' "${name}" >> "${stream}"
    done
}

copy_chemistry_files() {
    local run_dir="$1"
    shift
    local name
    for name in "$@"; do
        require_file "${MPAS_ROOT}/micm_configs/${name}"
        cp "${MPAS_ROOT}/micm_configs/${name}" "${run_dir}/"
    done
}

archive_outputs() {
    local run_dir="$1"
    local ts
    ts="$(date -u +%Y%m%dT%H%M%SZ)"
    shopt -s nullglob
    local files=()
    local path
    for path in "${run_dir}/output.nc" "${run_dir}/run.out"; do
        [[ -e "${path}" ]] && files+=("${path}")
    done
    for path in "${run_dir}"/log.atmosphere.*.out "${run_dir}"/log.atmosphere.*.err; do
        [[ -e "${path}" ]] && files+=("${path}")
    done
    if (( ${#files[@]} == 0 )); then
        shopt -u nullglob
        return
    fi
    mkdir -p "${run_dir}/archive/${ts}"
    mv "${files[@]}" "${run_dir}/archive/${ts}/"
    shopt -u nullglob
    log "archived stale run outputs under ${run_dir}/archive/${ts}"
}

write_active_chemistry_blocks() {
    local namelist="$1"
    local block="$2"
    require_file "${namelist}"
    NAMELIST_PATH="${namelist}" CHEMISTRY_BLOCKS="${block}" "${PYTHON}" - <<'PY'
from pathlib import Path
import os

path = Path(os.environ["NAMELIST_PATH"])
block = os.environ["CHEMISTRY_BLOCKS"].rstrip() + "\n"
lines = path.read_text().splitlines(keepends=True)
records = {"musica", "chemistry", "photolysis", "lnox"}

out = []
in_active_record = False
for line in lines:
    stripped = line.lstrip()
    if not in_active_record and stripped.startswith("&"):
        record = stripped[1:].split(None, 1)[0].lower()
        if record in records:
            in_active_record = True
            continue
    if in_active_record:
        if stripped.startswith("/"):
            in_active_record = False
        continue
    out.append(line)

while out and out[-1].strip() == "":
    out.pop()
out.append("\n\n")
out.append(block)
path.write_text("".join(out))
PY
}

override_run_duration_if_requested() {
    local namelist="$1"
    [[ -n "${CHEMPAS_RUN_DURATION:-}" ]] || return 0
    require_file "${namelist}"
    NAMELIST_PATH="${namelist}" RUN_DURATION="${CHEMPAS_RUN_DURATION}" "${PYTHON}" - <<'PY'
from pathlib import Path
import os
import re

path = Path(os.environ["NAMELIST_PATH"])
duration = os.environ["RUN_DURATION"]
pattern = re.compile(r"^(\s*)config_run_duration\s*=.*$")
lines = path.read_text().splitlines()
updated = False
for i, line in enumerate(lines):
    match = pattern.match(line)
    if match:
        lines[i] = f"{match.group(1)}config_run_duration = '{duration}'"
        updated = True
        break
if not updated:
    raise SystemExit(f"config_run_duration not found in {path}")
path.write_text("\n".join(lines) + "\n")
PY
    log "overrode config_run_duration in ${namelist}: ${CHEMPAS_RUN_DURATION}"
}

ensure_uniform_tracer() {
    local init_file="$1"
    shift
    require_file "${init_file}"
    INIT_FILE="${init_file}" TRACERS="$*" "${PYTHON}" - <<'PY'
from netCDF4 import Dataset
import numpy as np
import os

path = os.environ["INIT_FILE"]
pairs = [item.split("=", 1) for item in os.environ["TRACERS"].split()]
with Dataset(path, "r+") as ds:
    n_times = ds.dimensions["Time"].size
    n_cells = ds.dimensions["nCells"].size
    n_levels = ds.dimensions["nVertLevels"].size
    dtype = ds.variables["qv"].dtype if "qv" in ds.variables else "f8"
    for name, raw_value in pairs:
        value = float(raw_value)
        if name in ds.variables:
            var = ds.variables[name]
        else:
            var = ds.createVariable(name, dtype, ("Time", "nCells", "nVertLevels"))
            var.units = "kg kg^{-1}"
            var.long_name = name
        var[:] = np.full((n_times, n_cells, n_levels), value, dtype=dtype)
    ds.sync()
PY
}

verify_output() {
    local run_dir="$1"
    shift

    require_file "${run_dir}/output.nc"
    require_file "${run_dir}/log.atmosphere.0000.out"
    grep -Eq "Critical error messages[[:space:]]*=[[:space:]]*0" "${run_dir}/log.atmosphere.0000.out" \
        || die "log does not report zero critical errors: ${run_dir}/log.atmosphere.0000.out"

    OUTPUT_FILE="${run_dir}/output.nc" EXPECTED_VARS="$*" "${PYTHON}" - <<'PY'
from netCDF4 import Dataset
import numpy as np
import os

path = os.environ["OUTPUT_FILE"]
expected = os.environ["EXPECTED_VARS"].split()
with Dataset(path) as ds:
    if "Time" not in ds.dimensions:
        raise SystemExit("missing Time dimension")
    n_time = len(ds.dimensions["Time"])
    if n_time < 2:
        raise SystemExit(f"expected at least two Time records, found {n_time}")
    missing = [name for name in expected if name not in ds.variables]
    if missing:
        raise SystemExit(f"missing expected variables: {', '.join(missing)}")
    for name in expected:
        var = ds.variables[name]
        sample = np.ma.asarray(var[0, ...])
        if sample.size == 0:
            raise SystemExit(f"{name} has an empty first record")
        if np.ma.is_masked(sample) and bool(np.all(sample.mask)):
            raise SystemExit(f"{name} first record is fully masked")
        filled = np.asarray(sample.filled(np.nan) if np.ma.isMaskedArray(sample) else sample)
        if not np.isfinite(filled).any():
            raise SystemExit(f"{name} first record has no finite values")
print(f"verified {path}: Time={n_time}, variables={','.join(expected)}")
PY
}

ensure_supercell_base() {
    download_case_data supercell supercell_grid.nc
    require_file "${DATA_ROOT}/supercell/supercell.graph.info.part.${MPI_RANKS}"
}

prepare_supercell_run_dir() {
    local suffix="$1"
    local run_dir="${DATA_ROOT}/supercell_${suffix}"
    local base_dir="${DATA_ROOT}/supercell"

    ensure_supercell_base
    mkdir -p "${run_dir}"
    cp "${base_dir}/supercell_grid.nc" "${run_dir}/"
    cp "${base_dir}/supercell.graph.info.part.${MPI_RANKS}" "${run_dir}/"
    [[ -f "${base_dir}/supercell.graph.info" ]] && cp "${base_dir}/supercell.graph.info" "${run_dir}/"
    copy_tree_contents "${MPAS_ROOT}/test_cases/supercell" "${run_dir}"
    install_executables "${run_dir}"
    install_physics_tables "${run_dir}"
    override_run_duration_if_requested "${run_dir}/namelist.atmosphere"

    if [[ "${FORCE_INIT}" == "1" || ! -f "${run_dir}/supercell_init.nc" ]]; then
        log "initializing supercell in ${run_dir}"
        (cd "${run_dir}" && run_mpi ./init_atmosphere_model)
        require_file "${run_dir}/supercell_init.nc"
    fi

    PREPARED_RUN_DIR="${run_dir}"
}

run_supercell_abba() {
    local setup_only="$1"
    local run_dir
    prepare_supercell_run_dir abba
    run_dir="${PREPARED_RUN_DIR}"

    copy_chemistry_files "${run_dir}" abba.yaml
    write_active_chemistry_blocks "${run_dir}/namelist.atmosphere" "\
&chemistry
    config_micm_file = 'abba.yaml'
/"
    ensure_uniform_tracer "${run_dir}/supercell_init.nc" qA=0.0 qB=0.0 qAB=1.0
    (cd "${run_dir}" && "${PYTHON}" "${MPAS_ROOT}/scripts/init_tracer_sine.py" \
        -i supercell_init.nc -t qAB --create --waves-x 2 --amplitude 0.4 --offset 0.6)

    log "prepared Supercell ABBA run directory: ${run_dir}"
    [[ "${setup_only}" == "1" ]] && return

    archive_outputs "${run_dir}"
    (cd "${run_dir}" && run_mpi ./atmosphere_model)
    verify_output "${run_dir}" qA qB qAB
}

run_supercell_lnox() {
    local setup_only="$1"
    local run_dir
    prepare_supercell_run_dir lnox
    run_dir="${PREPARED_RUN_DIR}"

    copy_chemistry_files "${run_dir}" lnox_o3.yaml tuvx_no2.json tuvx_upper_atm.csv
    install_tuvx_data "${run_dir}"
    ensure_stream_vars "${run_dir}" j_jNO2
    write_active_chemistry_blocks "${run_dir}/namelist.atmosphere" "\
&chemistry
    config_micm_file = 'lnox_o3.yaml'
/

&photolysis
    config_tuvx_config_file = 'tuvx_no2.json'
    config_tuvx_top_extension = .true.
    config_tuvx_upper_column_mode = 'legacy_static'
    config_tuvx_extension_file = 'tuvx_upper_atm.csv'
    config_tuvx_update_interval = 600.0
    config_j_no2_max = 0.01
    config_chemistry_latitude = 35.86
    config_chemistry_longitude = -97.93
/

&lnox
    config_lnox_gating_mode = 'altitude'
    config_lnox_source_rate = 0.5
    config_lnox_w_threshold = 5.0
    config_lnox_w_ref = 10.0
    config_lnox_z_min = 5000.0
    config_lnox_z_max = 12000.0
    config_lnox_nox_tau = 0.0
/"
    (cd "${run_dir}" && "${PYTHON}" "${MPAS_ROOT}/scripts/init_lnox_o3.py" -i supercell_init.nc)

    log "prepared Supercell LNOx run directory: ${run_dir}"
    [[ "${setup_only}" == "1" ]] && return

    archive_outputs "${run_dir}"
    (cd "${run_dir}" && run_mpi ./atmosphere_model)
    verify_output "${run_dir}" qNO qNO2 qO3 j_jNO2
}

ensure_jw_init() {
    local jw_dir="${DATA_ROOT}/jw_baroclinic_wave"
    download_case_data jw_baroclinic_wave x1.40962.grid.nc
    require_file "${jw_dir}/x1.40962.graph.info.part.${MPI_RANKS}"

    copy_tree_contents "${MPAS_ROOT}/test_cases/jw_baroclinic_wave" "${jw_dir}"
    install_executables "${jw_dir}"
    install_physics_tables "${jw_dir}"

    if [[ "${FORCE_INIT}" == "1" || ! -f "${jw_dir}/x1.40962.init.nc" ]]; then
        log "initializing JW baroclinic-wave base data in ${jw_dir}"
        (cd "${jw_dir}" && run_mpi ./init_atmosphere_model)
        require_file "${jw_dir}/x1.40962.init.nc"
    fi
}

run_chapman_nox_global() {
    local setup_only="$1"
    local jw_dir="${DATA_ROOT}/jw_baroclinic_wave"
    local run_dir="${DATA_ROOT}/chapman_nox_global"

    ensure_jw_init
    mkdir -p "${run_dir}"
    copy_tree_contents "${MPAS_ROOT}/test_cases/chapman_nox_global" "${run_dir}"
    install_executables "${run_dir}"
    install_physics_tables "${run_dir}"
    copy_chemistry_files "${run_dir}" chapman_nox.yaml tuvx_chapman_nox.json
    require_file "${run_dir}/tuvx_upper_atm.csv"
    install_tuvx_data "${run_dir}"
    ensure_stream_vars "${run_dir}" j_jO2 j_jO3_O j_jO3_O1D j_jNO2
    ln -sf "${jw_dir}/x1.40962.init.nc" "${run_dir}/x1.40962.init.nc"
    ln -sf "${jw_dir}/x1.40962.graph.info.part.${MPI_RANKS}" "${run_dir}/x1.40962.graph.info.part.${MPI_RANKS}"
    override_run_duration_if_requested "${run_dir}/namelist.atmosphere"

    if [[ "${FORCE_INIT}" == "1" || ! -f "${run_dir}/x1.40962.chapman_nox_init.nc" ]]; then
        (cd "${run_dir}" && "${PYTHON}" "${MPAS_ROOT}/scripts/init_chapman_nox.py" \
            -i x1.40962.init.nc \
            -o x1.40962.chapman_nox_init.nc)
        require_file "${run_dir}/x1.40962.chapman_nox_init.nc"
    fi

    log "prepared Global Chapman + NOx run directory: ${run_dir}"
    [[ "${setup_only}" == "1" ]] && return

    archive_outputs "${run_dir}"
    (cd "${run_dir}" && run_mpi ./atmosphere_model)
    verify_output "${run_dir}" qO2 qO3 qO qNO qNO2 j_jO2 j_jO3_O j_jO3_O1D j_jNO2
}

main() {
    local case_name="${1:-}"
    local setup_only=0

    if [[ -z "${case_name}" || "${case_name}" == "--help" || "${case_name}" == "-h" ]]; then
        usage
        exit 0
    fi
    shift || true

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --setup-only)
                setup_only=1
                ;;
            --no-download)
                DOWNLOAD=0
                ;;
            --force-init)
                FORCE_INIT=1
                ;;
            --help|-h)
                usage
                exit 0
                ;;
            *)
                die "unknown option: $1"
                ;;
        esac
        shift
    done

    mkdir -p "${DATA_ROOT}"
    require_dir "${MPAS_ROOT}/test_cases"
    require_dir "${MPAS_ROOT}/micm_configs"
    require_dir "${MPAS_ROOT}/scripts"

    case "${case_name}" in
        supercell-abba)
            run_supercell_abba "${setup_only}"
            ;;
        supercell-lnox)
            run_supercell_lnox "${setup_only}"
            ;;
        chapman-nox-global)
            run_chapman_nox_global "${setup_only}"
            ;;
        *)
            usage
            die "unknown case: ${case_name}"
            ;;
    esac
}

main "$@"
