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

**[Figure 2.2: qA, qB, qAB at t = 2 h, ABBA mechanism. To be added.]**

## 2.6 Run with the LNOx + O3 mechanism

```{admonition} Work in progress
:class: warning

Section content coming.
```

**[Figure 2.3: NO, NO₂, O₃ at t = 2 h, LNOx + O3 mechanism. To be added.]**

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
