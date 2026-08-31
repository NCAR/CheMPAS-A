# CheMPAS-A Docker Workflow

This directory provides a Docker build and run harness for the three CheMPAS-A
chemistry demonstration cases:

- `supercell-abba`
- `supercell-lnox`
- `chapman-nox-global`

The image builds `init_atmosphere_model` and `atmosphere_model` with
MUSICA/MICM support. Its Python environment contains only NumPy and netCDF4,
which are required by the three case initializers and the built-in output
checks. Large NCAR input data and generated outputs live in a mounted data
directory, not in the image.

The build pins MUSICA to commit `1403e3d22717bc87f3bf9d0aa591caf039c92bbc`
(MUSICA-Fortran 0.16.5 with MICM, TUV-x, and MIEM), the known-good revision
qualified for these three container cases. The runtime image creates the small
`chempas-runtime` conda environment from `docker/environment.mpas.yml` via
micromamba.

The supported container target is Linux/AMD64. Docker Desktop can run it on
Apple Silicon through its x86-64 emulation support.

## Build

From the repository root:

```bash
docker build --platform linux/amd64 \
  -f docker/Containerfile --target run -t chempas-a:docker .
```

## Data Volume

Use a host directory for input data and outputs:

```bash
mkdir -p "$HOME/Data/CheMPAS"
```

The container sees that directory as `/data/CheMPAS`.

## Run Cases

Prepare and run Supercell ABBA:

```bash
docker run --rm \
  --platform linux/amd64 \
  --shm-size=1g \
  -v "$HOME/Data/CheMPAS:/data/CheMPAS" \
  chempas-a:docker supercell-abba
```

Prepare and run Supercell LNOx + O3:

```bash
docker run --rm \
  --platform linux/amd64 \
  --shm-size=1g \
  -v "$HOME/Data/CheMPAS:/data/CheMPAS" \
  chempas-a:docker supercell-lnox
```

Prepare and run Global Chapman + NOx:

```bash
docker run --rm \
  --platform linux/amd64 \
  --shm-size=1g \
  -v "$HOME/Data/CheMPAS:/data/CheMPAS" \
  chempas-a:docker chapman-nox-global
```

The runner downloads missing NCAR MPAS v7.0 tarballs when needed. To require
preexisting data, add `--no-download`. The `--shm-size=1g` allocation is
required by the eight-rank OpenMPI runs; Docker's default shared-memory volume
is too small.

## Setup Only

To prepare a run directory without launching MPAS:

```bash
docker run --rm \
  --platform linux/amd64 \
  --shm-size=1g \
  -v "$HOME/Data/CheMPAS:/data/CheMPAS" \
  chempas-a:docker supercell-abba --setup-only
```

## Output Locations

The runner uses isolated run directories:

```text
$HOME/Data/CheMPAS/supercell_abba/
$HOME/Data/CheMPAS/supercell_lnox/
$HOME/Data/CheMPAS/chapman_nox_global/
```

Stale `output.nc`, `run.out`, and atmosphere logs are moved under
`archive/<timestamp>/` before each model run.

## Useful Environment Overrides

```bash
CHEMPAS_FORCE_INIT=1         # regenerate init NetCDF files
CHEMPAS_RUN_DURATION=00:10:00  # override copied namelist run duration
```

Example smoke run with a shorter duration:

```bash
docker run --rm \
  --platform linux/amd64 \
  --shm-size=1g \
  -e CHEMPAS_RUN_DURATION=00:10:00 \
  -v "$HOME/Data/CheMPAS:/data/CheMPAS" \
  chempas-a:docker supercell-abba
```

For bounded qualification runs, use `00:10:00` for each supercell case and
`01:00:00` for `chapman-nox-global`. The global case writes hourly output, so
the one-hour duration is required for the verifier to see both initial and
final records.

## Verification

After a run, the harness checks:

- `output.nc` exists.
- `log.atmosphere.0000.out` exists.
- the log reports `Critical error messages = 0`.
- expected chemistry variables are present and finite in `output.nc`.

Expected variables:

```text
supercell-abba:     qA qB qAB
supercell-lnox:     qNO qNO2 qO3 j_jNO2
chapman-nox-global: qO2 qO3 qO qNO qNO2 j_jO2 j_jO3_O j_jO3_O1D j_jNO2
```

## Compose

The optional compose file builds the `run` target and mounts
`$HOME/Data/CheMPAS`:

```bash
docker compose -f docker/compose.yaml run --rm chempas supercell-abba
```
