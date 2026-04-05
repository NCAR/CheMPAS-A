#!/usr/bin/env bash
#
# Download WRF physics lookup tables from MPAS-Data repository.
# Called by the physics Makefile during build.
#

set -euo pipefail

MPAS_DATA_TAG="v8.2"
FILES_DIR="physics_wrf/files"
COMPAT_FILE="${FILES_DIR}/COMPATIBILITY"

# Skip if already present and compatible
if [[ -f "${COMPAT_FILE}" ]] && grep -q "8.2" "${COMPAT_FILE}" 2>/dev/null; then
    echo "Physics data files already present (${MPAS_DATA_TAG})."
    exit 0
fi

echo "Downloading WRF physics data files (${MPAS_DATA_TAG})..."

mkdir -p "${FILES_DIR}"

TARBALL_URL="https://codeload.github.com/MPAS-Dev/MPAS-Data/tar.gz/${MPAS_DATA_TAG}"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "${TMPDIR}"' EXIT

if command -v curl >/dev/null 2>&1; then
    curl -sL "${TARBALL_URL}" | tar xz -C "${TMPDIR}"
elif command -v wget >/dev/null 2>&1; then
    wget -qO- "${TARBALL_URL}" | tar xz -C "${TMPDIR}"
else
    echo "ERROR: Neither curl nor wget found. Cannot download physics data files." >&2
    exit 1
fi

EXTRACTED="$(ls -d "${TMPDIR}"/MPAS-Data-*)/atmosphere/physics_wrf/files"
if [[ ! -d "${EXTRACTED}" ]]; then
    echo "ERROR: Expected directory not found in tarball." >&2
    exit 1
fi

cp -v "${EXTRACTED}"/* "${FILES_DIR}/"
echo "Physics data files installed to ${FILES_DIR}/"
