# TUV-x Integration Summary

This note summarizes the CheMPAS TUV-x integration work completed to date and
records the development test case used to validate it. It is intended as the
stable overview document for TUV-x work; detailed implementation planning lives
in `docs/plans/2026-03-06-tuvx-photolysis-integration.md`, and recorded run
results live in `docs/results/TEST_RUNS.md`.

## Scope

CheMPAS currently uses TUV-x to compute `j_no2` for the idealized supercell
LNOx-O3 chemistry mechanism. The integration was developed in phases:

1. Phase 0: fixed-rate LNOx-O3 chemistry and lightning-NOx source
2. Phase 1: solar-zenith-angle fallback photolysis
3. Phase 2: clear-sky TUV-x photolysis from host atmospheric profiles
4. Phase 3: cloud attenuation in TUV-x using host `qc` and `qr`

The current implementation target beyond this summary is later Phase 3
hardening and follow-on photolysis realism, not the original Phase 1-2 bring-up.

## What Was Added

### Core chemistry changes

- `src/core_atmosphere/chemistry/mpas_lightning_nox.F`
  - operator-split lightning NO source before the MICM solve
- `src/core_atmosphere/chemistry/mpas_solar_geometry.F`
  - fallback solar geometry for idealized chemistry runs
- `src/core_atmosphere/chemistry/mpas_tuvx.F`
  - TUV-x initialization, per-column profile updates, and `j_no2` extraction
- `src/core_atmosphere/chemistry/mpas_atm_chemistry.F`
  - provider selection between fallback photolysis and TUV-x
  - MPAS host-state extraction for `rho_air`, temperature, `qv`, `qO3`, `qc`,
    and `qr`
  - `j_no2` diagnostic writing
- `src/core_atmosphere/chemistry/musica/mpas_musica.F`
  - scalar photolysis setter for Phase 1
  - field photolysis setter for Phase 2 and Phase 3

### Configuration and metadata

- `src/core_atmosphere/Registry.xml`
  - `config_tuvx_config_file`
  - `j_no2` diagnostic field
  - updated `config_lnox_j_no2` semantics as daytime `j_max`
- `micm_configs/tuvx_no2.json`
  - minimal TUV-x NO2 photolysis configuration using host profiles
- `micm_configs/lnox_o3.yaml`
  - tropospheric NO-NO2-O3 chemistry used for TUV-x development

### Validation tooling

- `scripts/check_tuvx_phase.py`
  - Phase-gate checks adapted to `qNO`, `qNO2`, `qO3`, and `j_no2`
- `scripts/run_tuvx_phase_gate.sh`
  - wrapper for Phase 0-2 gate checks
- `scripts/verify_ox_conservation.py`
  - Ox conservation verification for the fixed-rate prerequisite case

## Implementation Summary By Phase

### Phase 0: LNOx-O3 prerequisite

Purpose:
- replace the earlier placeholder chemistry with a scientifically meaningful
  tropospheric NO-NO2-O3 cycle
- establish the operator-split lightning source path before adding TUV-x

Delivered:
- `NO + O3 -> NO2 + O2`
- `NO2 + hv -> NO + O3`
- lightning NO source controlled by updraft threshold, reference velocity, and
  altitude window
- optional sink path via separate MICM config

Outcome:
- stable supercell runs
- correct O3 titration in the storm core
- exact Ox conservation in the controlled no-source/no-transport diagnostic

### Phase 1: fallback solar geometry

Purpose:
- validate photolysis-time dependence before introducing radiative transfer

Delivered:
- Spencer-style solar geometry helper in `mpas_solar_geometry.F`
- `j_no2 = j_max * max(0, cos_sza)` fallback path
- `j_no2` diagnostic output in MPAS history files

Outcome:
- daytime `j_no2` scaled exactly with `cos_sza`
- nighttime `j_no2` went to zero
- fallback path remained available when `config_tuvx_config_file = ''`

### Phase 2: clear-sky TUV-x

Purpose:
- replace scalar `j_no2` with column-dependent clear-sky photolysis computed
  from host atmospheric profiles

Delivered:
- `config_tuvx_config_file` namelist switch
- TUV-x module initialization and per-column execution
- extraction of height, density, temperature, water vapor, and ozone from MPAS
- `j_no2(level, cell)` write-back into MICM rate parameters

Outcome:
- physically plausible clear-sky vertical structure
- surface `j_no2` in the expected order of magnitude
- stable 15-minute idealized supercell chemistry runs

### Phase 3: cloud attenuation

Purpose:
- let TUV-x respond to storm cloud structure rather than treating all columns
  as clear-sky

Delivered:
- host-driven cloud radiator in TUV-x
- cloud optical depth derived from Kessler `qc` and `qr`
- cloud-aware `j_no2` attenuation and above-cloud enhancement

Outcome:
- strong `j_no2` suppression inside cloudy columns
- clear-sky columns preserved Phase 2 behavior
- above-cloud enhancement consistent with cloud albedo feedback

## Development Test Case

The TUV-x development case is the idealized supercell already used for the
LNOx-O3 chemistry work.

### Location and runtime

| Item | Value |
|------|-------|
| Run directory | `~/Data/CheMPAS/supercell` |
| Grid / case | Idealized supercell, 60 stretched levels (0–50 km) |
| Start time | `0000-01-01_18:00:00` |
| Nominal run duration | `00:30:00` for Phase 1, `00:15:00` for Phase 2/3 comparisons |
| Dynamics timestep | `3.0 s` |
| MPI layout | `mpiexec -n 8` |
| Decomposition prefix | `supercell.graph.info.part.` |
| Coordinates used for fallback SZA | `35.86 N`, `97.93 W` (Kingfisher, Oklahoma) |
| Physics suite | `none` |
| Microphysics | `mp_kessler` |

### Chemistry configuration

The tracked `test_cases/supercell/namelist.atmosphere` development settings are:

| Namelist option | Value | Meaning |
|-----------------|-------|---------|
| `config_micm_file` | `lnox_o3.yaml` | Tropospheric development mechanism |
| `config_lnox_source_rate` | `0.5` | NO source in `ppbv/s` when `w - w_threshold = w_ref` |
| `config_lnox_w_threshold` | `5.0` | Updraft threshold in `m/s` |
| `config_lnox_w_ref` | `10.0` | Excess updraft used for source scaling |
| `config_lnox_z_min` | `5000.0` | Minimum source altitude in `m` |
| `config_lnox_z_max` | `12000.0` | Maximum source altitude in `m` |
| `config_lnox_j_no2` | `0.01` | Phase 1 daytime `j_max` in `s^-1` |
| `config_lnox_nox_tau` | `0.0` | Sink disabled |
| `config_chemistry_latitude` | `35.86` | Fallback SZA latitude |
| `config_chemistry_longitude` | `-97.93` | Fallback SZA longitude |
| `config_tuvx_config_file` | `tuvx_no2.json` | Enable TUV-x for Phase 2/3 |

Additional setup used during development:

- initialize `qO3 = 50 ppbv`, `qNO = 0`, `qNO2 = 0` with
  `scripts/init_lnox_o3.py`
- copy `lnox_o3.yaml` and, for Phase 2/3, `tuvx_no2.json` into the run
  directory
- use `io_type="netcdf"` in `streams.atmosphere` on the macOS/LLVM build path

Why the synthetic start time:
- Phase 1 used a fixed daytime UTC that gives daylight at the Kingfisher test
  coordinates without requiring grid coordinates in the chemistry path
- the DC3-specific timestamp remains later validation work once grid-aware
  coordinates are available in the workflow

## Recorded Development Results

### Phase 1

- `cos_sza = 0.508` at `18:00 UTC`
- `cos_sza = 0.516` by `18:30 UTC`
- `j_no2 = 0.00508 -> 0.00516 s^-1`
- midnight test produced `j_no2 = 0`

### Phase 2

- surface `j_no2 = 7.2e-3 s^-1`
- top-of-domain `j_no2 = 1.2e-2 s^-1`
- `NO` peak `29.9 ppbv`
- `NO2` peak `6.5 ppbv`
- `O3` minimum `43.5 ppbv`
- empty `config_tuvx_config_file` correctly fell back to Phase 1 behavior

### Phase 3

- cloud develops after about 5 minutes in the supercell
- `346` cloudy columns and `27,734` clear columns at `t = 15 min`
- cloudy minimum `j_no2 = 6.27e-5 s^-1`
- clear-sky surface `j_no2 = 7.16e-3 s^-1`
- above-cloud `j_no2` reaches about `1.5e-2 s^-1`
- `NO2` remains higher in cloudy columns where photolysis recycling is reduced

For the complete tables, plots, and pass/fail notes, see
`docs/results/TEST_RUNS.md`.

## Multi-Photolysis Plumbing

CheMPAS couples TUV-x to MICM through a single-name convention: for every
photolysis reaction, the same string is used in four places —

```
<rxn_name>  == TUV-x JSON "name:" field
            == MICM yaml reaction "name:" field
            == Registry diagnostic variable "j_<rxn_name>"
            == MICM rate parameter key "PHOTO.<rxn_name>"
```

At chemistry init, `tuvx_init` enumerates every photolysis reaction the
TUV-x config registers and caches names + indices. The chemistry driver
then calls `musica_cache_photo_indices(names)` so MICM looks up
`PHOTO.<name>` for each one. At each step, `tuvx_compute_photolysis`
fills a `(n_rates, nVertLevels)` slab per cell, and a single call to
`musica_set_photolysis_rates` writes the full `(n_rates, nVertLevels,
nCells)` array into both the coupled and reference MICM states.

### Supported mechanisms (in `micm_configs/`)

| MICM yaml | TUV-x JSON | Rates | Purpose |
|-----------|------------|-------|---------|
| `lnox_o3.yaml` | `tuvx_no2.json` | jNO2 | Tropospheric LNOx-O3 development |
| `chapman_full.yaml` | `tuvx_chapman.json` | jO2, jO3_O, jO3_O1D | Stratospheric Chapman cycle |
| `chapman_nox.yaml` | `tuvx_chapman_nox.json` | jO2, jO3_O, jO3_O1D, jNO2 | Chapman + NOx catalytic O3 destruction |

Pick a paired MICM/TUV-x config via `&musica` in `namelist.atmosphere`:

```
&musica
    config_micm_file = 'chapman_full.yaml'
    config_tuvx_config_file = 'tuvx_chapman.json'
    config_tuvx_top_extension = .true.
    config_tuvx_extension_file = 'tuvx_upper_atm.csv'
    config_chemistry_latitude = 35.86
    config_chemistry_longitude = -97.93
/
```

Chapman and Chapman+NOx runs require realistic initial O3 (AFGL
mid-latitude-summer) — seed with `scripts/init_chapman.py`, optionally
`--seed-nox` for the NOx setup.

The Phase-1 `cos(SZA)` fallback is still available: when
`config_tuvx_config_file` is empty, chemistry runs with a single
synthesised `jNO2` rate equal to `config_lnox_j_no2 * max(0, cos_sza)`.
The fallback is single-rate only; Chapman and Chapman+NOx mechanisms
require TUV-x.

## Column Extension Above MPAS Top

TUV-x can see a climatology column above the MPAS domain so UV attenuation
by above-domain air and O3 is captured. Enable via the `&musica` namelist:

```
config_tuvx_top_extension  = .true.
config_tuvx_extension_file = 'tuvx_upper_atm.csv'
```

The CSV file lives in the run directory. Format (one edge per row, bottom-up):

```
z_km,T_K,n_air_molec_cm3,n_O3_molec_cm3
50.00,270.65,2.13e16,6.45e11
55.00,260.77,1.18e16,2.02e11
...
100.00,195.08,1.19e13,2.85e07
```

The shipped default (`micm_configs/tuvx_upper_atm.csv`) covers 50–100 km at
5-km spacing (10 layers) using US Standard Atmosphere 1976 for T and air
density, and AFGL mid-latitude-summer O3. Regenerate with
`scripts/gen_tuvx_upper_atm.py`.

At init, TUV-x is built on `nVertLevels + n_ext_layers` grid sections
(70 for the default supercell setup). During each photolysis call the
per-column profile is stitched:

- MPAS layers 1..nVertLevels: from host state (rho, T, qO3).
- Extension layers nVertLevels+1..nVertLevels+n_ext: midpoint values
  averaged from CSV edges.
- Edge values are the usual neighbour-averages; the MPAS/extension join
  edge averages the top MPAS midpoint with the first extension midpoint.
- Cloud OD is zero in the extension.

Extension-layer `j` values are computed by TUV-x but discarded — only the
first `nVertLevels` slice flows back to MICM.

The extension CSV must start at the MPAS top (zero gap). For the 50-km
supercell this means `z_km = 50.0` is the first CSV row.

## Current Limits And Follow-On Work

- **No cloud shadows.** TUV-x runs as independent 1D columns (plane-parallel).
  At SZA ~59° the `j_no2` cross-section shows purely vertical structure — a
  cloud attenuates photolysis directly below it but has zero effect on
  neighboring columns. In reality the solar beam enters at an angle, so clouds
  should cast shadows to one side, reducing photolysis in adjacent clear-sky
  columns. A slant-column approximation (trace the beam through neighboring
  columns at the geometric SZA, accumulate their cloud OD as above-column
  attenuation) would capture the dominant effect without full 3D radiative
  transfer.
- fallback SZA still uses namelist coordinates; grid-aware chemistry geometry is
  deferred
- current TUV-x work targets `j_no2` in the tropospheric LNOx-O3 mechanism, not
  the full Chapman photolysis set
- aerosols, earth-sun-distance refinement, and ice hydrometeors are later work
- clear-sky/cloudy column split optimization is deferred
- some extended gate checks remain future work, including longer transition
  tests and MPI decomposition comparison

## Pointers

- Detailed plan:
  `docs/plans/2026-03-06-tuvx-photolysis-integration.md`
- Recorded test results:
  `docs/results/TEST_RUNS.md`
- Run instructions:
  `RUN.md`
