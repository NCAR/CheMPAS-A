# Chapter 8: Chemistry Coupling

This chapter documents the runtime chemistry coupling between MPAS and the
MUSICA stack — MICM (chemistry solver) and TUV-x (photolysis solver) — as
implemented in `src/core_atmosphere/chemistry/`. The chapter assumes the
runtime tracer infrastructure described in [Chapter 7](07-runtime-tracers.md);
namelist details are in Appendix B. For the upstream chemistry-package
documentation see <https://musica.readthedocs.io/>,
<https://micm.readthedocs.io/>, and <https://tuv-x.readthedocs.io/>.

## 8.1 Overview

The chemistry pipeline lives in `src/core_atmosphere/chemistry/` and consists
of five Fortran modules:

| Module | File | Role |
|--------|------|------|
| `mpas_atm_chemistry`  | `mpas_atm_chemistry.F`        | Generic init/step/finalize driver, gather/scatter, photolysis-update gating |
| `mpas_musica`         | `musica/mpas_musica.F`        | MICM solver instance, coupled and reference states, species table, unit conversion, photolysis rate parameters |
| `mpas_tuvx`           | `mpas_tuvx.F`                 | TUV-x setup with from-host height grid, air/T/O3/O2 profiles, and cloud radiator |
| `mpas_solar_geometry` | `mpas_solar_geometry.F`       | Spencer (1971) per-cell `cos(SZA)` |
| `mpas_lightning_nox`  | `mpas_lightning_nox.F`        | Vertical-velocity-gated NO source |

The driver is `chemistry_step` (`mpas_atm_chemistry.F:330`), called once per
MPAS dynamics step. It is operator-split into five phases:

1. Optional photolysis-rate update (TUV-x or `cos(SZA)` fallback), gated by
   `config_tuvx_update_interval`.
2. Lightning NOx injection into `qNO` (no-op when the mechanism does not
   include `qNO`).
3. MPAS → MICM state gather, with conversion to mol m⁻³.
4. MICM solve, optionally split into `config_chem_substeps` sub-steps, with
   an optional uncoupled reference solve in parallel.
5. MICM → MPAS state scatter, with conversion back to mass mixing ratio.

Every routine that touches MICM is wrapped in `#ifdef MPAS_USE_MUSICA`, so
the entire chemistry pipeline compiles out without the `MUSICA=true` build
flag (see [Section 3.6](03-building.md#36-building-with-chemistry-musica-support)).

## 8.2 Initialization

`chemistry_init` (`mpas_atm_chemistry.F:73`) drives chemistry startup, in
order:

1. Read the `&musica` namelist record (paths, lightning NOx, TUV-x, solver
   controls, solar geometry).
2. `musica_init` (`musica/mpas_musica.F:74`) — instantiate the persistent
   `micm_t` solver from `config_micm_file` using the Rosenbrock standard
   ordering, allocate two `state_t` instances (coupled `state` and reference
   `state_ref`), populate the `chem_species` table from MICM's species
   ordering, read each species' `__molar mass` and `__initial concentration`
   properties, and seed both states.
3. `resolve_mpas_indices` (`musica/mpas_musica.F:181`) — for each MICM
   species name `X`, look up the MPAS scalar index from the pool dimension
   `index_qX`. The MPAS tracer pool was already extended at startup by
   `atm_extend_scalars_for_chemistry` (Chapter 7).
4. `chemistry_seed_chem` (`mpas_atm_chemistry.F:895`) — copy the MICM
   initial state back to the MPAS scalars on both time levels via
   `micm_to_mpas_chem` (`musica/mpas_musica.F:447`). Skipped per-tracer when
   the MPAS array already shows spatial gradients (i.e., the tracer was
   initialized from the input file).
5. `chemistry_lightning_nox_init` (`mpas_atm_chemistry.F:1027`) — read
   `config_lnox_*`; remains a no-op when `qNO` is not in the mechanism.
6. `tuvx_init` (`mpas_tuvx.F:118`) — register from-host grids and profiles,
   optionally extend the column with an upper-atmosphere CSV climatology,
   construct the TUV-x solver against a 102-bin CAM wavelength grid, and
   cache every photolysis reaction TUV-x reports. Skipped when
   `config_tuvx_config_file` is empty.
7. `musica_cache_photo_indices` (`musica/mpas_musica.F:722`) — for every
   reaction name TUV-x cached (or the single fallback name `jNO2`), look up
   `PHOTO.<name>` in MICM's `rate_parameters_ordering` and store the stride
   index for later writes.

`assign_rate_parameters` (`musica/mpas_musica.F:654`) populates MICM's
`rate_parameters` array at `musica_init` time. By convention, `USER.*`
parameters default to `1.0` so the YAML mechanism's `scaling_factor`
defines the effective rate; `PHOTO.*` parameters default to `0`. The
pseudo-first-order NOx loss `LOSS.no_loss` and `LOSS.no2_loss` parameters
are wired to `1 / config_lnox_nox_tau` when that namelist value is
positive.

Module-level state cached at init, used by `chemistry_step`:

- Solar geometry: `chem_lat`, `chem_lon`, `chem_use_grid_coords`,
  `chem_j_no2_max`.
- TUV-x update interval: `tuvx_update_interval`, `tuvx_time_since_last`.
- TUV-x activation and tracer indices: `use_tuvx`, `idx_qO3`, `idx_qc`,
  `idx_qr`.
- MICM controls: `chem_substeps_val`, `micm_relative_tolerance_val`,
  `use_ref_solve`.

## 8.3 The Chemistry Time Step

`chemistry_step(dt, currTime, mesh, state, diag, dimensions, time_lev)`
(`mpas_atm_chemistry.F:330`) is the per-step entry point.

**Phase 0 — Photolysis update.** A module accumulator
`tuvx_time_since_last` tracks simulated seconds since the photolysis block
last fired. When the accumulator reaches `config_tuvx_update_interval`,
the block fires:

- Compute per-cell `cos(SZA)` (broadcast from `config_chemistry_latitude`/
  `_longitude` when `config_chemistry_use_grid_coords = .false.`, otherwise
  per-cell from `latCell`/`lonCell`).
- Allocate `photo_rates(n_rates, nVertLevels, nCells)` and either
  call `tuvx_compute_photolysis` per cell or fill the single-rate fallback
  `jNO2 = j_max * max(0, cos(SZA))`.
- Write the rates into MICM via `musica_set_photolysis_rates`
  (`musica/mpas_musica.F:782`) and into the diag pool via
  `chemistry_set_photolysis_diag` (`mpas_atm_chemistry.F:631`).
- Reset the accumulator.

On skipped steps the accumulator just ticks forward and MICM keeps using
the rates last written. Setting `config_tuvx_update_interval = 0` (the
default) updates every step.

**Phase 1 — Lightning NOx injection.** `lightning_nox_inject`
(`mpas_lightning_nox.F:149`) modifies `scalars(idx_qNO, :, :)` in place.
Operator-split: the increment is added before the MICM call, not as a
tendency.

**Phase 2 — MPAS → MICM gather.** `chemistry_from_MPAS`
(`mpas_atm_chemistry.F:703`) reconstructs per-cell ρ, T, p (Section 8.4),
then calls `MICM_from_chemistry` (`musica/mpas_musica.F:544`) which writes
MICM's `state%conditions` and converts each species' mixing ratio to
mol m⁻³.

**Phase 3 — MICM solve.** The MPAS dt is divided into
`config_chem_substeps` calls to `musica_step` (`musica/mpas_musica.F:290`),
each advancing the coupled state by `dt / N`. TUV-x rates set in Phase 0
remain frozen across sub-steps. When `config_chemistry_ref_solve = .true.`,
`musica_step_ref` (`musica/mpas_musica.F:374`) runs in lockstep on the
reference state.

**Phase 4 — MICM → MPAS scatter.** `chemistry_to_MPAS`
(`mpas_atm_chemistry.F:816`) reconstructs ρ and calls `MICM_to_chemistry`
(`musica/mpas_musica.F:607`), which converts the integrated mol m⁻³
state back to mass mixing ratio and writes it into the scalars pool at
`time_lev`.

The phases are wrapped in MPAS timer regions: `chemistry`, `chem MPAS->MICM`,
`chem MICM solve`, `chem MICM ref solve` (when enabled), and
`chem MICM->MPAS`. They appear in the standard MPAS performance log.

## 8.4 State Transfer

MPAS stores chemistry tracers in the `scalars` pool as mass mixing ratio
(kg of species per kg of moist air). MICM works in volumetric concentration
(mol m⁻³). Conversion happens once per direction per chemistry step.

**Air state reconstruction.** `chemistry_from_MPAS`
(`mpas_atm_chemistry.F:703`) computes:

- ρ = `zz` × `rho_zz` × (1 + qv) [kg m⁻³]
- T = (`theta_m` / (1 + Rv/Rd × qv)) × `exner` [K]
- p = `pressure_p` + `pressure_base` (with an ideal-gas-law fallback) [Pa]

following the MPAS physics-interface idiom. These per-cell, per-level
quantities feed `state%conditions(:)%temperature` /`pressure` /`air_density`
on line 583 of `mpas_musica.F`. `air_density` is reported in mol m⁻³
(ρ / M_AIR with M_AIR = 0.0289644 kg mol⁻¹), as MICM expects.

**MICM cell layout.** MICM's working state is a flat array indexed by a
single grid-cell index that runs over the full (column, level) product:
`micm_cell = (iCell - 1) * nVertLevels + k` (line 579 of `mpas_musica.F`).
Per-species concentration storage is
`1 + (micm_cell - 1) * cell_stride + (micm_index - 1) * var_stride`,
where the strides are read from `state%species_strides` and the species
indices come from `state%species_ordering`.

**Forward conversion** (`MICM_from_chemistry`, line 544):

```
conc[mol m⁻³] = q[kg/kg] · ρ[kg/m³] / M[kg/mol]
```

`M` (the species molar mass) is read once at init from the MICM property
`__molar mass` and cached in `chem_species(:)%molar_mass`.

**Reverse conversion** (`MICM_to_chemistry`, line 607):

```
q[kg/kg] = conc[mol m⁻³] · M[kg/mol] / ρ[kg/m³]
```

The same ρ cache is used in both directions; `MICM_to_chemistry` does not
re-compute T or p — only concentrations are scattered back.

`micm_to_mpas_chem` (`musica/mpas_musica.F:447`) is a seed-only variant
used by `chemistry_seed_chem` to broadcast each species' MICM initial
concentration across all MPAS columns at startup. It does not write
environmental conditions.

## 8.5 Photolysis

Photolysis rates are external rate parameters from MICM's perspective:
MICM does not solve photolysis itself but reads
`state%rate_parameters(PHOTO.<name>)` for every photolysis reaction.
CheMPAS-A supplies these from one of two sources, set per build/run by
the namelist.

**TUV-x with from-host atmosphere and cloud radiator** (preferred). When
`config_tuvx_config_file` is non-empty, `tuvx_init` (`mpas_tuvx.F:118`)
constructs the solver using:

- A from-host height grid in km, sized for the composite column (MPAS
  layers + optional extension layers).
- From-host air, temperature, O3, and O2 profiles in molecule cm⁻³.
- A cloud radiator on a 102-bin CAM wavelength grid.

`tuvx_init` enumerates every photolysis reaction the JSON config declares,
caches the names and stride indices in `photo_names`/`tuvx_indices`, and
hands the same name list to `musica_cache_photo_indices` so MICM and TUV-x
agree on the reaction set.

`tuvx_compute_photolysis` (`mpas_tuvx.F:325`) runs once per cell per update:

1. Build the composite column = MPAS layers + extension layers.
2. Convert MPAS mid-layer values to TUV-x units (number density in
   molecule cm⁻³, layer densities in molecule cm⁻²).
3. Update the from-host profiles (`air_profile%set_*`, `temp_profile`,
   `o3_profile`, `o2_profile`) and call
   `air_profile%calculate_exo_layer_density(7.0_dk, ...)` to populate the
   above-top layer for spherical-geometry slant-path computations.
4. Compute per-layer cloud optical depth and write it to the cloud
   radiator's `set_optical_depths`. Single-scattering albedo (0.999999)
   and asymmetry parameter (0.85) are constants.
5. Call `tuvx_solver%run(sza_rad, esd, photo_rates, heating_rates, ...)`.
6. Average TUV-x's edge-defined rates to mid-layer values for the MPAS
   slice; extension-layer rates are computed but discarded.
7. Skip the TUV-x call entirely at night (cos(SZA) ≤ 0) and return zeros.

**Solar zenith angle.** `solar_cos_sza` (`mpas_solar_geometry.F:24`)
implements the Spencer (1971) declination + equation-of-time formula.
With `config_chemistry_use_grid_coords = .true.`, every cell uses its
own `latCell`/`lonCell`; otherwise all cells share the namelist
`config_chemistry_latitude` and `config_chemistry_longitude`, intended for
idealized cases where solar geometry should be uniform across the domain.

**Cloud optical depth.** `compute_cloud_optical_depth`
(`mpas_tuvx.F:732`) uses

```
τ = 3 · LWC · dz / (2 · r_eff · ρ_water)
```

with `r_eff = 10 µm` for cloud water (`qc`) and `r_eff = 500 µm` for rain
water (`qr`). LWC is reconstructed from `q · ρ_air` in each layer. The
total τ is the sum of cloud and rain contributions; by the choice of
`r_eff`, rain contributes ~50× less per unit mass than cloud water.

**Upper-atmosphere column extension.** When
`config_tuvx_top_extension = .true.`, `load_extension_csv`
(`mpas_tuvx.F:636`) reads `(z_km, T_K, n_air_cm⁻³, n_O3_cm⁻³)` edge values
from `config_tuvx_extension_file` and stitches them onto the top of the
MPAS column. The extension's temperature profile is anchored to the MPAS
top mid-layer T so the join is continuous (the MPAS sounding may lack a
stratopause, while the climatology does not). Heights must be strictly
increasing; the file format is one header row plus N data rows.

**Phase 1 fallback.** When `config_tuvx_config_file` is empty, the driver
falls back to a single rate `jNO2 = config_lnox_j_no2 · max(0, cos(SZA))`,
filling slot 1 of `photo_rates`. The MICM mechanism must then declare
`PHOTO.jNO2` for `musica_cache_photo_indices` to wire up — otherwise
`chemistry_init` aborts with an error.

**Update interval.** `config_tuvx_update_interval` gates the entire
photolysis block in seconds. The default (0) updates every chemistry
step; positive values hold rates between updates and MICM reuses the
last-set values. Useful when TUV-x dominates the chemistry-step cost
and rates change slowly relative to the MPAS dt.

## 8.6 Lightning NOx

`mpas_lightning_nox.F` is a stand-alone source-only module: it injects
`qNO` mass mixing ratio in cells where the layer-midpoint vertical
velocity exceeds a threshold and the layer height lies inside a
configured range. The chemistry response (NO + O3 → NO2, NO2 + hν →
NO + O3) is handled by MICM.

`lightning_nox_init` (`mpas_lightning_nox.F:61`) looks up `index_qNO` in
the scalars pool. If `qNO` is not in the active mechanism (e.g., the
placeholder `abba.yaml`), the module sets `lnox_active = .false.` and
the inject hook becomes a no-op. The same disabled state results when
`config_lnox_source_rate ≤ 0` or `config_lnox_w_ref ≤ 0`.

`lightning_nox_inject` (`mpas_lightning_nox.F:149`) is called from
`chemistry_step` before the MPAS → MICM gather (operator-split, not a
tendency). For each (cell, level) pair:

- Layer-mid w = ½ · (`w(k)` + `w(k+1)`).
- Layer-mid z = ½ · (`zgrid(k)` + `zgrid(k+1)`).
- Source increment in mass mixing ratio:

```
Δq = config_lnox_source_rate · max(0, w_mid - config_lnox_w_threshold) / config_lnox_w_ref
     · dt · 1e-9 · (M_NO / M_AIR)
```

where `1e-9` converts ppbv to mole fraction and `M_NO / M_AIR` ≈
0.030 / 0.029 (kg mol⁻¹) converts mole fraction to mass mixing ratio.
Active only when `w_mid > config_lnox_w_threshold` and
`config_lnox_z_min ≤ z_mid ≤ config_lnox_z_max`. With the default
`config_lnox_source_rate = 0`, the entire injection is disabled.

## 8.7 Solver Controls

Three namelist options shape MICM's per-step behavior. They exist because
stiff mechanisms with fast null cycles (e.g., the Chapman O / O₃ pair)
can drive Backward Euler onto non-physical implicit fixed points at the
full MPAS dt; sub-stepping or tighter tolerances avoid the bad root.

**`config_chem_substeps`** (default 1). The MPAS dt is divided into N
sub-steps and `musica_step` is called N times with `dt / N`
(`mpas_atm_chemistry.F:548`). Photolysis rates are set once per outer
step and frozen across sub-steps.

**`config_micm_relative_tolerance`** (default 1e-6, MICM's own default).
Tightening (e.g., to 1e-9) forces MICM's adaptive controller to subdivide
its internal steps more aggressively before accepting a step. Passed to
`musica_init` and applied at solver-state construction.

**Re-entrant solve loop.** `musica_step` (`musica/mpas_musica.F:290`)
wraps `micm%solve` in a sub-call loop with `MAX_SUB_CALLS = 100`. MICM's
adaptive controller has its own internal step budget
(`max_number_of_steps`, default 11) and returns early when it exhausts
that budget without reaching the requested interval. The wrapper inspects
`solver_stats%final_time()`, computes the remaining duration, and
resubmits until the full interval is covered. It aborts if MICM advances
0 s on any sub-call (genuine stall) or if the 100-sub-call cap is hit.
This pattern is borrowed from MUSICA-Fortran's `column_model.F90` example
and is required for stiff transients (e.g., Chapman chemistry at sunrise).

**Reference solve.** When `config_chemistry_ref_solve = .true.`,
`musica_step_ref` (`musica/mpas_musica.F:374`) runs in lockstep on
`state_ref`, an independent MICM state seeded once from the coupled state
on the first chemistry step (`copy_state_to_ref`,
`musica/mpas_musica.F:966`). The reference state is not coupled to
MPAS — it sees the same photolysis rates and the same time step, but no
advection or per-step environmental updates. Subsequent divergence
between coupled and reference states is therefore attributable to
transport. The reference solve doubles chemistry cost and is off by
default.

## 8.8 Diagnostics and Logging

**Photolysis diag fields.** `chemistry_set_photolysis_diag`
(`mpas_atm_chemistry.F:631`) writes the photolysis rate array into
`j_<reaction>` diag-pool fields after every photolysis update. Reactions
whose `j_<reaction>` field is not declared in `Registry.xml` are silently
skipped, so the chemistry driver remains agnostic about which mechanism
is loaded. `zero_photolysis_diag` (`mpas_atm_chemistry.F:672`) clears
candidate fields (`j_jNO2`, `j_jO2`, `j_jO3_O`, `j_jO3_O1D`) at init,
before the first history write.

**MICM solver statistics.** `log_solver_stats`
(`musica/mpas_musica.F:891`) emits MICM's per-step counters: function
calls, Jacobian updates, total internal steps, accepted and rejected
steps, LU decompositions, linear solves, and the final time the solver
reached. Logged after every `musica_step` outer call.

**Re-entrant sub-call summary.** `musica_step` logs a single line when
MICM completed the requested interval in more than one sub-call,
indicating the adaptive controller's internal budget was exhausted but
the wrapper recovered. Recurring sub-call activity is a signal to
tighten `config_micm_relative_tolerance` or raise
`config_chem_substeps`.

**Coupled vs reference column comparison.** When
`config_chemistry_ref_solve = .true.`, `log_column_comparison`
(`musica/mpas_musica.F:997`) is called from `chemistry_step` after each
MICM solve. For a probe cell (`nCells / 2`), it logs coupled and
reference concentrations every `nVertLevels / 4` levels for every
species. Differences quantify advection effects.

**Per-species seed log.** `micm_to_mpas_chem` (`musica/mpas_musica.F:447`)
logs `min`/`max` of each tracer after the initial-state seed, and the
species-table init in `musica_init` logs each species' resolved molar
mass.

For full namelist documentation including units and default values, see
Appendix B.
