# QUACS Development Instructions

This document describes how to get up-and-running as a QUACS developer. It covers:
* Setting up your VSCode to develop in the CheMPAS-A container environment
* Building and running CheMPAS-A in your development container
* Alternative command-line option for building and running CHeMPAS-A in a container

---

## Prerequisites

* [Docker](https://docs.docker.com/get-docker/) or [Podman](https://podman.io/getting-started/installation) (with the `docker` CLI alias)
* [VS Code](https://code.visualstudio.com/) with the [Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension
* Git

---

## Setup
### Option 1: Developing inside the VS Code Dev Container (recommended)

The `.devcontainer/` configuration mounts the repository source at `/mpas` inside a pre-built container that has all build dependencies (MPI, NetCDF, PnetCDF, MUSICA, etc.) already installed.

1. **Open the repo in VS Code** and, when prompted, click **Reopen in Container**.  
   Alternatively, open the Command Palette (`Ctrl+Shift+P`) and run **Dev Containers: Reopen in Container**.

2. VS Code will pull the `dev` image, start the container, and open a terminal at `/mpas`.

3. **Build the init_atmosphere core**:
   ```bash
   make -j$(nproc) gnu CORE=init_atmosphere AUTOCLEAN=true MUSICA=true
   ```

4. **Build the atmosphere core**:
   ```bash
   make -j$(nproc) gnu CORE=atmosphere AUTOCLEAN=true MUSICA=true
   ```
   Both executables (`init_atmosphere_model`, `atmosphere_model`) are placed in `/mpas`.

---

### Option 2: Building the container image locally (alternative approach)

0. This is only needed if you want to build the `dev` Docker image from scratch (e.g., to modify the `Containerfile`).

   ```bash
   docker build --target dev -t chempas-deps -f docker/Containerfile .
   ```

1. To build the full runtime image with MPAS executables baked in (the `build` stage clones physics externals and downloads WRF lookup tables automatically):
   ```bash
   docker build -t chempas --target build -f docker/Containerfile .
   ```

2. You can then run the container and proceed to the next section:
   ```bash
   docker run -it chempas bash
   ```
---

## Downloading test data and running the JW baroclinic wave test

These steps assume the executables have already been built (see above).

1. **Download the 480-km mesh and test case** (into the `data/` directory):
   ```bash
   quacs/download_data.sh data
   ```
   Pass `--all` as a second argument to also fetch the 240-km mesh:
   ```bash
   quacs/download_data.sh data --all
   ```

2. **Run the Jablonowski-Williamson baroclinic wave test**:
   ```bash
   quacs/run_jw_test.sh [NPROCS]
   ```
   `NPROCS` is the number of MPI ranks (default: 1). For example, with 2 ranks:
   ```bash
   quacs/run_jw_test.sh 2
   ```
   Output is written to `data/jw_480km/output.nc`.


---

## CI/CD

The GitHub Actions workflow (`.github/workflows/chempas.yml`) runs on every push and pull request to the `proto-quacs` branch. It:

1. **Builds the `dev` dependency image** and pushes it to the GitHub Container Registry (`ghcr.io`), using a SHA-tagged image for reproducibility and a `latest` tag on the default branch.
2. **Builds both MPAS cores** (`init_atmosphere` and `atmosphere`) inside that image.
3. **Runs the JW baroclinic wave test** on the 480-km mesh with 2 MPI ranks and verifies that `output.nc` is produced.

