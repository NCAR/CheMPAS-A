#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"

eval "$("${repo_root}/scripts/check_build_env.sh" --export)"

src="${script_dir}/chapman_nox_level2_reproducer.F90"
config="${repo_root}/micm_configs/chapman_nox.yaml"
build_dir="${TMPDIR:-/tmp}/chempas_chapman_nox_reproducer"
exe="${build_dir}/chapman_nox_level2_reproducer"

mkdir -p "${build_dir}"

cflags="$(pkg-config --cflags musica-fortran)"
libs="$(pkg-config --libs musica-fortran)"

# shellcheck disable=SC2206
cflags_array=(${cflags})
# shellcheck disable=SC2206
libs_array=(${libs})

OMPI_FC="${FC:-flang}" mpifort \
  "${cflags_array[@]}" \
  "${src}" \
  -o "${exe}" \
  "${libs_array[@]}"

args=("${config}")
if [[ $# -gt 0 ]]; then
  args+=("$@")
fi

"${exe}" "${args[@]}"
