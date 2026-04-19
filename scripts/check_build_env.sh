#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  scripts/check_build_env.sh
  eval "$(scripts/check_build_env.sh --export)"

Checks the local build prerequisites used by CheMPAS-A.
Auto-detects the Fortran compiler toolchain (LLVM/flang or GCC/gfortran).

Options:
  --export   Print shell exports for the detected working environment.
  --help     Show this help text.
EOF
}

mode="report"
if [[ $# -gt 1 ]]; then
    usage >&2
    exit 2
fi
if [[ $# -eq 1 ]]; then
    case "$1" in
        --export) mode="export" ;;
        --help|-h) usage; exit 0 ;;
        *)
            usage >&2
            exit 2
            ;;
    esac
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

status=0
override_pkgconfig=""
pkg_config_path_value="${PKG_CONFIG_PATH:-}"
netcdf_root=""
netcdff_root=""
pnetcdf_root=""
pio_root=""
musica_status=""
preferred_pkgconfig_dir="${HOME}/software/lib/pkgconfig"

# Auto-detect compiler toolchain
fc_compiler=""
cc_compiler=""
cxx_compiler=""
make_target=""

if command -v flang >/dev/null 2>&1; then
    fc_compiler="flang"
    cc_compiler="clang"
    cxx_compiler="clang++"
    make_target="llvm"
elif command -v gfortran >/dev/null 2>&1; then
    fc_compiler="gfortran"
    cc_compiler="gcc"
    cxx_compiler="g++"
    make_target="gfortran"
fi

report() {
    if [[ "${mode}" == "report" ]]; then
        printf '%s\n' "$*"
    fi
}

fail() {
    status=1
    report "[fail] $*"
}

pass() {
    report "[ok] $*"
}

resolve_root() {
    local required_path_a="$1"
    local required_path_b="${2:-}"
    shift 2
    local candidate

    for candidate in "$@"; do
        [[ -n "${candidate}" ]] || continue
        if [[ -e "${candidate}/${required_path_a}" ]]; then
            if [[ -z "${required_path_b}" || -e "${candidate}/${required_path_b}" ]]; then
                printf '%s\n' "${candidate}"
                return 0
            fi
        fi
    done

    return 1
}

pkg_config_env() {
    if [[ -n "${pkg_config_path_value}" ]]; then
        printf '%s\n' "${pkg_config_path_value}"
    fi
}

pkg_config_has_valid_musica() {
    local env_path="$1"
    local includedir=""
    local libdir=""

    if ! includedir="$(PKG_CONFIG_PATH="${env_path}" pkg-config --variable=includedir musica-fortran 2>/dev/null)"; then
        return 1
    fi
    if ! libdir="$(PKG_CONFIG_PATH="${env_path}" pkg-config --variable=libdir musica-fortran 2>/dev/null)"; then
        return 1
    fi

    [[ -f "${includedir}/musica_micm.mod" ]] || return 1

    if [[ -e "${libdir}/libmusica-fortran.a" || -e "${libdir}/libmusica-fortran.dylib" || -e "${libdir}/libmusica-fortran.so" ]]; then
        if [[ -e "${libdir}/libmusica.a" || -e "${libdir}/libmusica.dylib" || -e "${libdir}/libmusica.so" ]]; then
            return 0
        fi
    fi

    return 1
}

write_musica_override() {
    local build_root="$1"
    local yaml_libs=""

    if [[ -d /opt/homebrew/opt/yaml-cpp/lib ]]; then
        yaml_libs="-L/opt/homebrew/opt/yaml-cpp/lib "
    fi

    override_pkgconfig="/tmp/chemps-musica-fortran.pc"
    cat > "${override_pkgconfig}" <<EOF
prefix=${build_root}
exec_prefix=\${prefix}
libdir=\${prefix}/lib
includedir=\${prefix}/mod_fortran

yaml_lib=yaml-cpp

Name: musica-fortran
Description: Fortran wrapper for the MUSICA library for modeling atmospheric chemistry
Version: 0.14.5
Cflags: -I\${includedir}
Libs: -L\${libdir} ${yaml_libs}-lmusica-fortran -lmusica -lmechanism_configuration -l\${yaml_lib}
EOF

    if [[ -n "${pkg_config_path_value}" ]]; then
        pkg_config_path_value="/tmp:${pkg_config_path_value}"
    else
        pkg_config_path_value="/tmp"
    fi
}

probe_musica_link() {
    local cflags
    local libs
    local tmp_src
    local tmp_exe
    local -a cflags_array
    local -a libs_array

    cflags="$(PKG_CONFIG_PATH="${pkg_config_path_value}" pkg-config --cflags musica-fortran)"
    libs="$(PKG_CONFIG_PATH="${pkg_config_path_value}" pkg-config --libs musica-fortran)"
    read -r -a cflags_array <<<"${cflags}"
    read -r -a libs_array <<<"${libs}"

    tmp_src="$(mktemp /tmp/check_musica_fortran.XXXXXX)"
    mv "${tmp_src}" "${tmp_src}.f90"
    tmp_src="${tmp_src}.f90"
    tmp_exe="${tmp_src%.f90}.x"
    cat > "${tmp_src}" <<'EOF'
program test_musica_fortran
  use musica_util, only : string_t
  use musica_micm, only : get_micm_version
  type(string_t) :: version_string
  version_string = get_micm_version()
  print *, trim(version_string%value_)
end program test_musica_fortran
EOF

    if OMPI_FC="${fc_compiler}" mpifort "${cflags_array[@]}" "${tmp_src}" -o "${tmp_exe}" "${libs_array[@]}" -lstdc++ >/dev/null 2>&1; then
        rm -f "${tmp_src}" "${tmp_exe}"
        return 0
    fi

    rm -f "${tmp_src}" "${tmp_exe}"
    return 1
}

report "CheMPAS-A build environment preflight"

if [[ -n "${fc_compiler}" ]]; then
    pass "Fortran compiler: ${fc_compiler} ($(command -v ${fc_compiler}))"
    pass "Make target: ${make_target}"
else
    fail "No Fortran compiler found (expected flang or gfortran in PATH)"
fi

if command -v mpifort >/dev/null 2>&1; then
    pass "mpifort: $(command -v mpifort)"
else
    fail "mpifort not found in PATH"
fi

if command -v pkg-config >/dev/null 2>&1; then
    pass "pkg-config: $(command -v pkg-config)"
else
    fail "pkg-config not found in PATH"
fi

# NetCDF: check for .dylib (macOS) or .so (Linux) or .a
netcdf_lib_pattern=""
if netcdf_root="$(resolve_root "include/netcdf.mod" "lib/libnetcdff.dylib" "${NETCDF:-}" "/opt/homebrew" "${CONDA_PREFIX:-}")"; then
    netcdf_lib_pattern="dylib"
elif netcdf_root="$(resolve_root "include/netcdf.mod" "lib/libnetcdff.so" "${NETCDF:-}" "${CONDA_PREFIX:-}")"; then
    netcdf_lib_pattern="so"
elif netcdf_root="$(resolve_root "include/netcdf.mod" "lib/libnetcdff.a" "${NETCDF:-}" "${CONDA_PREFIX:-}")"; then
    netcdf_lib_pattern="a"
fi
if [[ -n "${netcdf_root}" ]]; then
    pass "NETCDF=${netcdf_root}"
else
    fail "NETCDF not found; expected include/netcdf.mod and lib/libnetcdff.*"
fi

# NETCDFF override for flang: Homebrew's libnetcdff is gfortran-mangled and
# cannot resolve flang-mangled __QMnetcdfP* symbols from libmusica.a (pulled in
# via TUV-x). MUSICA-LLVM ships a flang-built netcdf-fortran at flang-deps/
# that matches. Detect and export NETCDFF so the Makefile links against it.
if [[ "${fc_compiler}" == "flang" ]]; then
    if netcdff_root="$(resolve_root "include/netcdf.mod" "lib/libnetcdff.a" "${NETCDFF:-}" "${HOME}/EarthSystem/MUSICA-LLVM/flang-deps/netcdf-fortran-install")"; then
        pass "NETCDFF=${netcdff_root} (flang-built netcdf-fortran for MUSICA+TUV-x)"
    else
        report "[note] No flang-built netcdf-fortran found; MUSICA+TUV-x builds will fail to link."
        report "       Build one at ~/EarthSystem/MUSICA-LLVM/flang-deps/netcdf-fortran-install,"
        report "       or set NETCDFF to an existing flang-built tree."
    fi
fi

if pnetcdf_root="$(resolve_root "include/pnetcdf.h" "lib/libpnetcdf.a" "${PNETCDF:-}" "${HOME}/software" "${CONDA_PREFIX:-}")" || \
   pnetcdf_root="$(resolve_root "include/pnetcdf.h" "lib/libpnetcdf.so" "${PNETCDF:-}" "${HOME}/software" "${CONDA_PREFIX:-}")"; then
    pass "PNETCDF=${pnetcdf_root}"
else
    fail "PNETCDF not found; expected include/pnetcdf.h and lib/libpnetcdf.*"
fi

if pio_root="$(resolve_root "include/pio.mod" "lib/libpiof.a" "${PIO:-}" "${HOME}/software")"; then
    pass "PIO=${pio_root}"
else
    fail "PIO not found; expected include/pio.mod and lib/libpiof.a"
fi

# MUSICA detection: check installed pkg-config first, then fall back to build trees
musica_build_root=""
if [[ "${fc_compiler}" == "flang" ]]; then
    musica_build_root="$(resolve_root "mod_fortran/musica_micm.mod" "lib/libmusica-fortran.a" "${MUSICA_BUILD_DIR:-}" "${repo_root}/../MUSICA-LLVM/build" "${HOME}/EarthSystem/MUSICA-LLVM/build" || true)"
elif [[ "${fc_compiler}" == "gfortran" ]]; then
    musica_build_root="$(resolve_root "mod_fortran/musica_micm.mod" "lib/libmusica-fortran.a" "${MUSICA_BUILD_DIR:-}" "${repo_root}/../MUSICA/build" "${HOME}/EarthSystem/MUSICA/build" || true)"
fi

if command -v pkg-config >/dev/null 2>&1; then
    if [[ -f "${preferred_pkgconfig_dir}/musica-fortran.pc" ]]; then
        if [[ -n "${pkg_config_path_value}" ]]; then
            pkg_config_path_value="${preferred_pkgconfig_dir}:${pkg_config_path_value}"
        else
            pkg_config_path_value="${preferred_pkgconfig_dir}"
        fi
    fi

    if PKG_CONFIG_PATH="${pkg_config_path_value}" pkg-config --exists musica-fortran 2>/dev/null && \
       pkg_config_has_valid_musica "${pkg_config_path_value}"; then
        musica_status="pkg-config"
        pass "musica-fortran resolved from current pkg-config search path"
    elif [[ -n "${musica_build_root}" ]]; then
        write_musica_override "${musica_build_root}"
        if pkg_config_has_valid_musica "${pkg_config_path_value}"; then
            musica_status="override"
            pass "musica-fortran override written to ${override_pkgconfig}"
        else
            fail "wrote ${override_pkgconfig}, but it still does not resolve a valid MUSICA module/lib tree"
        fi
    else
        fail "musica-fortran not found, and no local MUSICA build tree was detected"
    fi
fi

if [[ ${status} -eq 0 && -n "${musica_status}" ]]; then
    if probe_musica_link; then
        pass "MUSICA link probe succeeded with OMPI_FC=${fc_compiler} mpifort"
    else
        fail "MUSICA link probe failed; check yaml-cpp and the musica-fortran.pc library paths"
    fi
fi

if [[ "${mode}" == "export" ]]; then
    if [[ ${status} -ne 0 ]]; then
        exit ${status}
    fi

    printf 'export NETCDF=%q\n' "${netcdf_root}"
    if [[ -n "${netcdff_root}" ]]; then
        printf 'export NETCDFF=%q\n' "${netcdff_root}"
    fi
    printf 'export PNETCDF=%q\n' "${pnetcdf_root}"
    printf 'export PIO=%q\n' "${pio_root}"
    printf 'export PKG_CONFIG_PATH=%q\n' "${pkg_config_path_value}"
    printf 'export OMPI_FC=%q\n' "${fc_compiler}"
    printf 'export OMPI_CC=%q\n' "${cc_compiler}"
    printf 'export OMPI_CXX=%q\n' "${cxx_compiler}"
    exit 0
fi

if [[ ${status} -eq 0 ]]; then
    report ""
    report "Suggested next command:"
    report "  eval \"\$(scripts/check_build_env.sh --export)\""
    if [[ -n "${netcdff_root}" ]]; then
        report "  make -j8 ${make_target} CORE=atmosphere PIO=${pio_root} NETCDF=${netcdf_root} NETCDFF=${netcdff_root} PNETCDF=${pnetcdf_root} PRECISION=double MUSICA=true"
    else
        report "  make -j8 ${make_target} CORE=atmosphere PIO=${pio_root} NETCDF=${netcdf_root} PNETCDF=${pnetcdf_root} PRECISION=double MUSICA=true"
    fi
fi

exit ${status}
