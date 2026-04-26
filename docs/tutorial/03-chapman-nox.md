# Chapter 3: Chapman + NOx Photostationary State

```{admonition} Work in progress
:class: warning

This chapter is being actively written. Commands and expected output
are provisional; figures are placeholders.
```

The Chapman + NOx photostationary-state (PSS) tutorial walks through a
small-domain integration of the Chapman ozone cycle plus NOx, with
TUV-x photolysis driven by an extended atmosphere column that reaches
above the model lid. The analytical Leighton expression for [NO]/[NO₂]
under steady-state photolysis is a clean numerical sanity check on the
coupled MICM + TUV-x configuration.

## 3.1 What you'll learn

```{admonition} Work in progress
:class: warning

Section content coming.
```

By the end of this chapter you will:

- Run the Chapman + NOx idealized stratospheric-chemistry case in
  CheMPAS-A on the supercell mesh.
- Generate the TUV-x upper-atmosphere extension CSV and understand
  why TUV-x needs photons from above the model lid.
- Verify the chemistry against the analytical Leighton photostationary
  state and the regression suite.

## 3.2 The Chapman + NOx case

```{admonition} Work in progress
:class: warning

Section content coming.
```

The Chapman cycle is the canonical four-reaction pure-oxygen
photochemistry that maintains a stratospheric ozone column:

$$
\begin{aligned}
\mathrm{O_2} + h\nu &\rightarrow 2\,\mathrm{O} \\
\mathrm{O} + \mathrm{O_2} + \mathrm{M} &\rightarrow \mathrm{O_3} + \mathrm{M} \\
\mathrm{O_3} + h\nu &\rightarrow \mathrm{O} + \mathrm{O_2} \\
\mathrm{O} + \mathrm{O_3} &\rightarrow 2\,\mathrm{O_2}
\end{aligned}
$$

Adding NOx introduces the catalytic
NO–NO₂–O₃ cycle (NO + O₃ → NO₂ + O₂; NO₂ + hν → NO + O), which
modulates ozone titration by tying its evolution to NOx photolysis.
On any timescale longer than a few seconds, [NO] / [NO₂] in sunlight
relaxes to the **Leighton photostationary state**, the analytical
target of section 3.7.

The Chapman cycle is global-stratospheric physics, but
`scripts/init_chapman.py` seeds a 1-D AFGL mid-latitude-summer ozone
profile uniformly across the supercell mesh, and the chemistry has no
feedback on dynamics. This chapter therefore uses the small
(~85 km × 85 km × 50 km top) supercell grid as a column-like sandbox
— what matters is the vertical structure of the photolysis driver and
the chemistry's ability to settle into the PSS, both of which TUV-x
sees through the column extension introduced in section 3.3.
Horizontal dynamics are present but largely irrelevant to the PSS
demonstration.

**[Figure 3.1: AFGL mid-latitude-summer O₃ profile interpolated to the
supercell vertical grid (the initial state qO3 produces). To be added.]**

## 3.3 The TUV-x column extension

```{admonition} Work in progress
:class: warning

Section content coming.
```

The MPAS atmosphere lid for the supercell case sits at 50 km — that's
roughly stratopause level, but Chapman-cycle photolysis depends on UV
radiation that has already been attenuated by the entire ozone column
*above* the photolysis cell, including the ~50–100 km region MPAS
itself does not simulate. Without an extension, TUV-x sees vacuum
above 50 km, jO₃ and jNO₂ are off by a non-trivial factor at high
altitudes, and the Chapman steady state never establishes properly.

The fix is `micm_configs/tuvx_upper_atm.csv`: a tracked CSV carrying
temperature, air number density, and ozone number density on a
uniform 5-km grid from 50 to 100 km. The temperature and air values
come from the US Standard Atmosphere 1976 tables; the ozone values
come from the AFGL mid-latitude-summer constituent profile. At
runtime, `mpas_tuvx.F::load_extension_csv` stitches MPAS midpoint
values (lower slice) and CSV midpoint values (upper slice) into a
single radiator column for TUV-x, blending across the boundary so
the profile is continuous.

The stitch lives in `src/core_atmosphere/chemistry/mpas_tuvx.F`; for
the broader integration story, see
[docs/chempas/guides/TUVX_INTEGRATION.md](../chempas/guides/TUVX_INTEGRATION.md).

## 3.4 Generating and verifying the extension CSV

```{admonition} Work in progress
:class: warning

Section content coming.
```

**Generate the CSV.** The generator is parameterized but defaults to
the configuration the runtime expects (50–100 km, 10 layers, 5-km
spacing). Run:

```bash
cd ~/EarthSystem/CheMPAS-A
~/miniconda3/envs/mpas/bin/python \
    scripts/gen_tuvx_upper_atm.py --out micm_configs/tuvx_upper_atm.csv
```

The script emits a header line followed by one row per edge with
columns `z_km, T_K, n_air_molec_cm3, n_O3_molec_cm3`. The output path
must match the `config_tuvx_extension_file` value in the namelist
(set in section 3.6 below).

**Verify the stitched column.** The companion plotter overlays the
MPAS region with the extension-CSV region as TUV-x actually sees
them, including the edge-blending the runtime applies at the 50-km
boundary:

```bash
cd ~/Data/CheMPAS/supercell
~/miniconda3/envs/mpas/bin/python \
    ~/EarthSystem/CheMPAS-A/scripts/plot_extension_profiles.py \
    -i output.nc
```

Note: this verification requires an `output.nc` from a CheMPAS-A run.
You can run it now against any prior supercell `output.nc` (for
example, the LNOx run from Chapter 2) just to inspect the column
shape, and then run it again after the Chapman + NOx run lands in
section 3.6 below to see the stitched column the actual run used.

**[Figure 3.2: Stitched T, n_air, and n_O₃ vertical profiles from
mpas_tuvx.F. MPAS region (below 50 km) and extension-CSV region
(above 50 km) overplotted. To be added.]**

## 3.5 Initializing the Chapman tracers

```{admonition} Work in progress
:class: warning

Section content coming.
```

## 3.6 Run with the Chapman + NOx mechanism

```{admonition} Work in progress
:class: warning

Section content coming.
```

**[Figure 3.3: Vertical profiles of qO3, qNO, qNO2, and the NO/NO₂
ratio at t = 2 h, mid-domain column, Chapman + NOx mechanism. To be
added.]**

## 3.7 The photostationary-state diagnostic

```{admonition} Work in progress
:class: warning

Section content coming.
```

## 3.8 Verifying numerically

```{admonition} Work in progress
:class: warning

Section content coming.
```

**[Figure 3.4: Simulated vs. analytical Leighton [NO]/[NO₂] ratio vs.
height at the final timestep. To be added.]**

## 3.9 Next steps

```{admonition} Work in progress
:class: warning

Section content coming.
```
