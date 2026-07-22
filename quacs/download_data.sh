#!/usr/bin/env bash
# Download MPAS-A meshes and idealized test case archives.
# Usage: ./quacs/download_data.sh [DATA_DIR [--all]]
#   DATA_DIR:   directory in which to store the data (default: data)
#   Without --all: downloads only the 480-km mesh + static + JW baroclinic wave
#   With --all:    also downloads the 240-km mesh + static
set -euo pipefail

DATA_DIR="${1:-data}"
mkdir -p "${DATA_DIR}"
cd "${DATA_DIR}"

MESH_BASE="https://www2.mmm.ucar.edu/projects/mpas/atmosphere_meshes"
IDEAL_BASE="https://www2.mmm.ucar.edu/projects/mpas/test_cases/v7.0"

download() {
    local url="$1"
    local basename
    basename="$(basename "${url}")"
    if [ -f "${basename}" ]; then
        echo "  Already exists: ${basename}"
    else
        echo "  Downloading: ${basename}"
        wget -q --show-progress "${url}"
    fi
}

echo "=== 480-km quasi-uniform mesh (x1.2562) ==="
download "${MESH_BASE}/x1.2562.tar.gz"
download "${MESH_BASE}/x1.2562_static.tar.gz"

echo "=== Jablonowski-Williamson baroclinic wave test case ==="
download "${IDEAL_BASE}/jw_baroclinic_wave.tar.gz"

if [[ "${2:-}" == "--all" ]]; then
    echo "=== 240-km quasi-uniform mesh (x1.10242) ==="
    download "${MESH_BASE}/x1.10242.tar.gz"
    download "${MESH_BASE}/x1.10242_static.tar.gz"
fi

echo ""
echo "=== Extracting archives ==="
for f in *.tar.gz; do
    echo "  Extracting: ${f}"
    tar xzf "${f}"
done

echo ""
echo "Done. Data in: $(pwd)"
ls -lh
