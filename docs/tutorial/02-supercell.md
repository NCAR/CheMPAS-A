# Chapter 2: Supercell with ABBA and LNOx

## 2.1 What you'll learn

```{admonition} Work in progress
:class: warning

Section content coming.
```

By the end of this chapter you will:

- Run the supercell idealized convection case in CheMPAS-A.
- Switch between two MUSICA/MICM chemistry mechanisms — the toy ABBA
  reactant set and the LNOx + O3 lightning-NOx tropospheric setup —
  by editing the `&musica` namelist block.
- Compare what tracer transport looks like in isolation (ABBA) against
  what it looks like coupled to active gas-phase chemistry and a
  lightning-NOx source (LNOx + O3).

## 2.2 The supercell case

```{admonition} Work in progress
:class: warning

Section content coming.
```

The supercell test case is an idealized deep-convection setup adapted
from Weisman and Klemp (1982): a horizontally homogeneous, strongly
sheared environment perturbed by a warm bubble that triggers a
long-lived rotating storm. CheMPAS-A's tracked configuration
uses 60 stretched vertical levels spanning 0–50 km (~300 m at the
surface, ~1 km near the lid), a 3-second dynamics timestep, Kessler
warm-rain microphysics, and a 2-hour integration.

It is a useful chemistry testbed for two reasons. First, the rotating
updraft produces strong, well-organized vertical transport that lifts
boundary-layer tracers into the upper troposphere; the cold-pool
outflow then spreads species horizontally near the surface. Second,
the dynamics are deterministic and the chemistry adds no feedback to
the dynamics — which means a chemistry change shows up cleanly as a
chemistry-only signal, rather than confounded with a different storm
evolution.

**[Figure 2.1: Supercell initial state — potential temperature and
moisture cross-section. To be added.]**

## 2.3 Setup checklist

```{admonition} Work in progress
:class: warning

Section content coming.
```

Before you run anything, confirm:

```bash
# 1. The atmosphere_model executable exists and is from this branch.
ls -la ~/EarthSystem/CheMPAS-A/atmosphere_model

# 2. The supercell run directory is set up.
cd ~/Data/CheMPAS/supercell
ls namelist.atmosphere streams.atmosphere supercell.graph.info.part.8

# 3. The vertical-level file is in place (read by init_atmosphere_model).
ls supercell_zeta_levels.txt

# 4. The conda env is active for plotting.
conda activate mpas
python -c "import netCDF4, numpy, matplotlib; print('ok')"
```

If any of these fail, see `BUILD.md` in the repository root for the
model build, `RUN.md` in the repository root for the run-directory
layout, and `test_cases/README.md` for the data download.

## 2.4 Initialization

```{admonition} Work in progress
:class: warning

Section content coming.
```

If `supercell_init.nc` is already in the run directory, you can skip
ahead to running the model. To regenerate the initial condition from
scratch:

```bash
cd ~/Data/CheMPAS/supercell
mpiexec -n 8 ~/EarthSystem/CheMPAS-A/init_atmosphere_model
```

This produces `supercell_init.nc`, which contains the prognostic state
on the unstructured Voronoi mesh, the stretched vertical levels read
from `supercell_zeta_levels.txt`, and the Kessler microphysics
variables.

Always run with 8 MPI ranks for the supercell case — the partition
file in the run directory (`supercell.graph.info.part.8`) is keyed to
that rank count, and a mismatched partition file causes a segfault in
the dynamics solver. See `RUN.md` in the repository root for the full
rank-vs-partition table.

## 2.5 Run with the ABBA mechanism

```{admonition} Work in progress
:class: warning

Section content coming.
```

The ABBA mechanism (`micm_configs/abba.yaml`) is a toy reactant set —
qA, qB, and qAB — used as a transport-only sandbox: chemistry advances
the species, but the rate constants are chosen so that on the run
duration the field is dominated by advection rather than reaction.
This makes it a clean way to verify that tracers are being carried by
the dynamics correctly before adding real chemistry on top.

**Edit the namelist.** Open `namelist.atmosphere` in
`~/Data/CheMPAS/supercell/` and ensure the `&musica` block reads:

```fortran
&musica
    config_micm_file = 'abba.yaml'
/
```

ABBA does not need TUV-x photolysis or the lightning-NOx options, so
the block is minimal.

**Archive any prior output and run.** MPAS defaults to
`clobber_mode = never_modify`; if `output.nc` already exists the run
will silently skip all output writes. Always move the previous output
aside before re-running:

```bash
cd ~/Data/CheMPAS/supercell

timestamp=$(date +%Y%m%d_%H%M%S)
[ -f output.nc ] && mv output.nc output.${timestamp}.nc
[ -f log.atmosphere.0000.out ] && \
    mv log.atmosphere.0000.out log.atmosphere.0000.${timestamp}.out

mpiexec -n 8 ~/EarthSystem/CheMPAS-A/atmosphere_model
```

**Verify the run completed cleanly.** The tail of
`log.atmosphere.0000.out` should report zero critical errors:

```
Critical error messages = 0
```

**Plot.** With the conda env active:

```bash
~/miniconda3/envs/mpas/bin/python \
    ~/EarthSystem/CheMPAS-A/scripts/plot_chemistry.py -o supercell_abba.png
```

**[Figure 2.2: qA, qB, qAB at t = 2 h, ABBA mechanism. To be added.]**

What to look for: qA and qB are conserved tracers transported by the
flow; qAB is produced where qA and qB co-exist. In the supercell, the
strongest gradients sit along the updraft and in the cold-pool
outflow.

## 2.6 Run with the LNOx + O3 mechanism

```{admonition} Work in progress
:class: warning

Section content coming.
```

The LNOx + O3 setup (`micm_configs/lnox_o3.yaml`) is a tropospheric
gas-phase configuration with three prognostic species — NO, NO₂, and
O₃ — and a parameterized lightning-NOx source term tied to the model's
vertical velocity. It is the smallest realistic chemistry case in
CheMPAS-A: enough species and reactions to exercise the MICM solver,
TUV-x photolysis, and the LNOx source coupling, without the cost of a
full tropospheric mechanism.

**Initialize the LNOx tracers.** The supercell init file does not
contain NO / NO₂ / O₃; populate them with a one-time script:

```bash
cd ~/Data/CheMPAS/supercell
~/miniconda3/envs/mpas/bin/python \
    ~/EarthSystem/CheMPAS-A/scripts/init_lnox_o3.py -i supercell_init.nc
```

This sets NO = 0, NO₂ = 0, and O₃ = 50 ppbv (background) throughout
the domain.

**Edit the namelist.** Replace the `&musica` block in
`namelist.atmosphere` with the full LNOx tropospheric setup:

```fortran
&musica
    config_micm_file = 'lnox_o3.yaml'
    config_tuvx_config_file = 'tuvx_no2.json'
    config_tuvx_top_extension = .true.
    config_tuvx_extension_file = 'tuvx_upper_atm.csv'
    config_lnox_source_rate = 0.5
    config_lnox_w_threshold = 5.0
    config_lnox_w_ref = 10.0
    config_lnox_z_min = 5000.0
    config_lnox_z_max = 12000.0
    config_lnox_j_no2 = 0.01
    config_lnox_nox_tau = 0.0
    config_chemistry_latitude = 35.86
    config_chemistry_longitude = -97.93
/
```

The lightning-NOx source injects NO into grid cells where the vertical
velocity exceeds `config_lnox_w_threshold` and the height falls
between `config_lnox_z_min` and `config_lnox_z_max`. See
[docs/guides/TUVX_INTEGRATION.md](../chempas/guides/TUVX_INTEGRATION.md)
for the TUV-x configuration files referenced in the block.

**Archive prior output and run.** Same pattern as the ABBA run:

```bash
timestamp=$(date +%Y%m%d_%H%M%S)
[ -f output.nc ] && mv output.nc output.${timestamp}.nc
[ -f log.atmosphere.0000.out ] && \
    mv log.atmosphere.0000.out log.atmosphere.0000.${timestamp}.out

mpiexec -n 8 ~/EarthSystem/CheMPAS-A/atmosphere_model
```

**Plot.** The dedicated LNOx plotting script produces the standard
diagnostic set (vertical cross-sections, time series, NO₂ partitioning
ratio):

```bash
~/miniconda3/envs/mpas/bin/python \
    ~/EarthSystem/CheMPAS-A/scripts/plot_lnox_o3.py
```

**[Figure 2.3: NO, NO₂, O₃ at t = 2 h, LNOx + O3 mechanism. To be added.]**

What to look for: a localized NO source in the updraft column where
the vertical-velocity threshold is exceeded, downwind transport of
NO + NO₂ along the anvil, and an O₃ depletion signature in the freshly
emitted plume (titration by NO).

## 2.7 Comparing the two runs

```{admonition} Work in progress
:class: warning

Section content coming.
```

**[Figure 2.4: Side-by-side comparison of ABBA tracer transport and
LNOx + O3 chemistry at t = 2 h. To be added.]**

## 2.8 Verifying numerically

```{admonition} Work in progress
:class: warning

Section content coming.
```

## 2.9 Next steps

```{admonition} Work in progress
:class: warning

Section content coming.
```
