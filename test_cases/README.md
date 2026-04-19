# Idealized Test Cases

Reference namelists and streams configurations for MPAS idealized test cases.
These files are tracked in the repo; the actual run data lives in `~/Data/MPAS/`.

## Data Download

Download and extract the three test cases from NCAR:

```bash
mkdir -p ~/Data/MPAS && cd ~/Data/MPAS
curl -LO https://www2.mmm.ucar.edu/projects/mpas/test_cases/v7.0/supercell.tar.gz
curl -LO https://www2.mmm.ucar.edu/projects/mpas/test_cases/v7.0/mountain_wave.tar.gz
curl -LO https://www2.mmm.ucar.edu/projects/mpas/test_cases/v7.0/jw_baroclinic_wave.tar.gz
for f in *.tar.gz; do tar xzf "$f"; done
rm *.tar.gz
```

## Test Cases

| Case | Init Case | Mesh | Cells | Levels | dt (s) | Duration | Chemistry Use |
|------|-----------|------|-------|--------|--------|----------|---------------|
| `supercell/` | 5 | ~500 m | ~40k | 40 | 3 | 2 hours | LNOx-O3 dev, TUV-x |
| `mountain_wave/` | 6 | ~577 m | ~2k | 70 | 6 | 5 hours | — |
| `jw_baroclinic_wave/` | 2 | 120 km | 40,962 | 26 | 450 | 16 days | Chapman + NOx (planned) |

### Notes

- **supercell**: tropospheric domain (0–20 km), doubly-periodic. Primary dev case for LNOx-O3 mechanism and TUV-x photolysis integration.
- **jw_baroclinic_wave**: global domain (no lateral BCs). The standard v7.0 config has 26 levels with a damping layer at 22 km, which provides only shallow stratospheric coverage. The vertical grid will be extended to higher altitude with increased resolution in the stratosphere for Chapman + NOx chemistry development.

## Setup and Initialization

After downloading, copy the tracked configs and initialize each case:

```bash
# Copy reference configs (with io_type="netcdf") to run directories
for case in supercell mountain_wave jw_baroclinic_wave; do
  cp test_cases/$case/* ~/Data/MPAS/$case/
  ln -sf ~/EarthSystem/CheMPAS/init_atmosphere_model ~/Data/MPAS/$case/
  ln -sf ~/EarthSystem/CheMPAS/atmosphere_model ~/Data/MPAS/$case/
done

# Initialize each case
for case in supercell mountain_wave jw_baroclinic_wave; do
  cd ~/Data/MPAS/$case
  mpiexec -n 4 ./init_atmosphere_model
  cd -
done
```

## Available Partition Files (MPI ranks)

| Case | Partitions |
|------|-----------|
| supercell | 2, 4, 8, 12, 16, 24, 32 |
| mountain_wave | 2, 4, 6, 8 |
| jw_baroclinic_wave | 2, 4, 6, 8, 12, 16, 24 |

## Running

See [RUN.md](../RUN.md) for run instructions and verification steps.

## Files Tracked Here

For each test case:
- `namelist.init_atmosphere` — initialization namelist
- `streams.init_atmosphere` — initialization I/O streams
- `namelist.atmosphere` — model run namelist
- `streams.atmosphere` — model run I/O streams (with `io_type="netcdf"`)
- `stream_list.atmosphere.output` — output variable list

## Data Source

All test case data from [NCAR MPAS idealized test cases (v7.0)](https://www2.mmm.ucar.edu/projects/mpas/site/access_code/idealized.html).
