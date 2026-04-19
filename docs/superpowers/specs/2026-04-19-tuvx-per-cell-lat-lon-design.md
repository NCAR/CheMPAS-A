# TUV-x Per-Cell Solar Geometry — Design

Date: 2026-04-19
Status: Design (awaiting user review)
Target files: `src/core_atmosphere/Registry.xml`,
              `src/core_atmosphere/chemistry/mpas_atm_chemistry.F`

## Goal

Enable per-cell solar geometry in the chemistry path by reading
`latCell` / `lonCell` from the mesh pool, gated by a new
`config_chemistry_use_grid_coords` namelist option (default `.false.`).
When the switch is enabled, `cos_sza` is computed per cell from grid
coordinates; when disabled, the existing namelist scalars
(`config_chemistry_latitude` / `config_chemistry_longitude`) are
broadcast to every cell, preserving exact current behavior for all
existing test cases (supercell, mountain wave, baroclinic wave).

## Scope

**In scope:**
- New namelist option in `&musica`: `config_chemistry_use_grid_coords`
- Per-cell `cos_sza_cell(:)` array in `chemistry_step`, replacing the
  scalar `cos_sza`
- Downstream substitution of `cos_sza` with `cos_sza_cell(iCell)` in the
  TUV-x call and in the fallback `j_no2` calculation

**Out of scope:**
- TUV-x sub-stepping / call-frequency throttling (separate spec, lands
  next as a sequential PR)
- Init-side seeding (no current need for per-cell lat/lon at init)
- Lightning_nox path (does not use lat/lon)
- Per-cell time zone or earth-sun-distance corrections
- Vectorizing the `solar_cos_sza` API itself (deferred; would be a wider
  refactor with no measurable benefit at current TUV-x cost)

## Pre-existing State

- `latCell` and `lonCell` are required mesh fields in
  `src/core_atmosphere/Registry.xml` (lines 1380–1384), dimension
  `nCells`, units `rad`.
- `mpas_solar_geometry.F::solar_cos_sza` takes scalar `lat_deg` /
  `lon_deg` (degrees) and returns scalar `cos_sza`.
- `chemistry_step` in `mpas_atm_chemistry.F` currently computes a single
  scalar `cos_sza` at line 379 from `chem_lat` / `chem_lon` (namelist
  values stored at module scope) and uses it everywhere downstream:
  - Line 434: TUV-x call (`tuvx_compute_photolysis(cos_sza, ...)`).
  - Lines 453–454: fallback `j_no2` calc
    (`j_no2_value = chem_j_no2_max * max(0.0, cos_sza)`).
- `chem_lat` / `chem_lon` are populated in `chemistry_init` from
  `config_chemistry_latitude` / `config_chemistry_longitude` (lines
  159–164).
- TUV-x integration documentation flags this gap explicitly:
  `docs/guides/TUVX_INTEGRATION.md`, "Current Limits And Follow-On Work"
  section: *"fallback SZA still uses namelist coordinates; grid-aware
  chemistry geometry is deferred"*.

## Design

### Namelist & Registry

Add one new namelist option to the existing `&musica` record in
`src/core_atmosphere/Registry.xml`, immediately after
`config_chemistry_longitude`:

```xml
<nml_option name="config_chemistry_use_grid_coords" type="logical" default_value=".false."
     description="If .true., compute per-cell solar geometry from latCell/lonCell. If .false., broadcast namelist config_chemistry_latitude/longitude to every cell (preserves idealized-case behavior)."
     possible_values=".true. or .false."/>
```

Default `.false.` preserves exact current behavior for every test case
and existing run. The two existing scalars
(`config_chemistry_latitude`, `config_chemistry_longitude`) are kept,
unchanged in semantics; they are simply ignored on a per-cell basis when
the new switch is `.true.`. Their description strings stay as-is — they
remain meaningful for the idealized-case path.

No edits to MICM YAML configs or TUV-x JSON configs. No new fields,
dimensions, or pool entries.

### `chemistry_step` data flow

Inside `chemistry_step` in
`src/core_atmosphere/chemistry/mpas_atm_chemistry.F`:

1. **Read the new switch and grid coords.** Add a new module-scope
   `logical, save :: chem_use_grid_coords = .false.`, populated in
   `chemistry_init` from `config_chemistry_use_grid_coords` (mirroring
   the existing `chem_lat` / `chem_lon` pattern at lines 159–164).
   In `chemistry_step`, retrieve `latCell` / `lonCell` from the `mesh`
   pool when the switch is `.true.`:

   ```fortran
   real(kind=RKIND), dimension(:), pointer :: latCell, lonCell
   if (chem_use_grid_coords) then
       call mpas_pool_get_array(mesh, 'latCell', latCell)
       call mpas_pool_get_array(mesh, 'lonCell', lonCell)
   end if
   ```

2. **Replace scalar `cos_sza` with a per-cell array.** Allocate once per
   call (deallocated at end), declared near the existing `cos_sza`
   variable:

   ```fortran
   real(kind=RKIND), allocatable :: cos_sza_cell(:)
   real(kind=RKIND), parameter   :: RAD2DEG = 57.29577951308232_RKIND
   ...
   allocate(cos_sza_cell(nCells))
   if (chem_use_grid_coords) then
       do iCell = 1, nCells
           cos_sza_cell(iCell) = solar_cos_sza(DoY, hour_utc, &
                                               latCell(iCell) * RAD2DEG, &
                                               lonCell(iCell) * RAD2DEG)
       end do
   else
       cos_sza_cell(:) = solar_cos_sza(DoY, hour_utc, chem_lat, chem_lon)
   end if
   ```

   The single existing scalar `cos_sza` declaration is removed; the
   existing `cos_sza = solar_cos_sza(...)` line at 379 is removed.

3. **Substitute downstream uses:**
   - Line 434, TUV-x call: replace `cos_sza` with `cos_sza_cell(iCell)`.
   - Lines 453–454, fallback path: replace
     ```fortran
     j_no2_value = chem_j_no2_max * max(0.0_RKIND, cos_sza)
     photo_rates(1, :, :) = j_no2_value
     ```
     with
     ```fortran
     do iCell = 1, nCells
         photo_rates(1, :, iCell) = chem_j_no2_max * max(0.0_RKIND, cos_sza_cell(iCell))
     end do
     ```

4. **Deallocate** `cos_sza_cell` at the end of the chemistry_step body
   alongside the other locals.

Memory cost: 8 bytes × `nCells` per call (≈800 KB at 100K cells,
negligible compared to `photo_rates` at `n_rates × nVertLevels ×
nCells`).

### Edge cases

- **Switch `.false.` on idealized cases:** Broadcasting the namelist
  scalar across all cells produces a `cos_sza_cell` array where every
  entry is identical and equal to the current scalar `cos_sza` value.
  All downstream arithmetic is bit-for-bit identical to the current
  code. The supercell, mountain wave, and baroclinic wave runs must
  produce identical output to the pre-change baseline.
- **Switch `.true.` on a Cartesian/idealized mesh:** Allowed — the user
  opted in. `latCell` / `lonCell` will reflect whatever the init step
  produced, likely small or zero values. Result is meaningful only if
  the user knows what they're doing. No defensive check; this is a
  user-configuration concern.
- **`latCell` / `lonCell` absent from mesh pool:** Cannot happen — they
  are required mesh fields. No defensive check needed.
- **Polar / extreme latitudes:** `solar_cos_sza` already handles
  arbitrary lat/lon; per-cell call adds nothing new.

### Validation

1. **Reproducibility check (switch `.false.`):** Run the existing
   supercell Phase 2/3 case with the new code and
   `config_chemistry_use_grid_coords = .false.`. Diff output `j_no2`,
   `qNO`, `qNO2`, `qO3` against the pre-change reference run. Expected:
   bit-for-bit identical (no operation reordering is introduced).
2. **Per-cell sanity check (switch `.true.`):** Run the
   Jablonowski–Williamson baroclinic wave case
   (`test_cases/jw_baroclinic_wave/`, real spherical 120-km mesh) with
   `config_chemistry_use_grid_coords = .true.`. Plot a horizontal slice
   of `j_no2` at a single time step and confirm it shows a clear
   longitudinal gradient (terminator) and latitudinal SZA dependence
   consistent with the simulation time of day.
3. **Phase-gate:** No change to `scripts/check_tuvx_phase.py` or
   `run_tuvx_phase_gate.sh` — they continue to gate the existing
   supercell case.

## Files Touched

**Modified (2):**
- `src/core_atmosphere/Registry.xml`
- `src/core_atmosphere/chemistry/mpas_atm_chemistry.F`

**Not touched:**
- `src/core_atmosphere/chemistry/mpas_solar_geometry.F` (signature
  unchanged)
- `src/core_atmosphere/chemistry/mpas_tuvx.F` (signature unchanged)
- `src/core_atmosphere/chemistry/musica/mpas_musica.F`,
  `src/core_atmosphere/chemistry/mpas_lightning_nox.F` (don't use
  lat/lon)
- `src/core_atmosphere/mpas_atm_core.F` (doesn't touch lat/lon)
- MICM YAML configs, TUV-x JSON configs, namelist defaults in
  `test_cases/`

## Commit Shape

Single commit on `develop`:

```
feat(chemistry): add per-cell solar geometry switch for TUV-x
```

Body explains the new `config_chemistry_use_grid_coords` switch,
default `.false.`, and the bit-reproducibility guarantee for existing
test cases.

## Doc Updates

- `docs/guides/TUVX_INTEGRATION.md` "Current Limits And Follow-On Work":
  remove the bullet *"fallback SZA still uses namelist coordinates;
  grid-aware chemistry geometry is deferred"* (now resolved).
- Add a brief note to the same doc describing the new namelist option
  and when to enable it.

Doc updates can land in the same commit as the code change or as a
small follow-on commit, at the implementer's discretion.

## Sequencing

This is the first of two sequential TUV-x PRs (per the user's choice in
brainstorming, option B). The next PR will add TUV-x sub-stepping /
call-frequency throttling on top of this baseline.
