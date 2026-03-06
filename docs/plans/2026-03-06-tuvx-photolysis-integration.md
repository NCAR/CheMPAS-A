# TUV-x Photolysis Integration Plan

## Document Status

- `Historical Context:` Adapted from ancestor project plan
  (`MPAS-Model-ACOM-dev/PLAN_TUVx.md`). The original plan assumed static
  Registry.xml tracers; this version accounts for CheMPAS Phase 2 runtime
  tracer allocation.
- `Current State:` Planning stage. Not yet started.
- `Use This As:` Primary reference for TUV-x integration work.

## Goal

Integrate TUV-x into MPAS-MUSICA so MICM photolysis rates are computed from
radiative transfer using MPAS atmospheric state, replacing hardcoded constants.

## Key Adaptation Notes (vs Ancestor Plan)

The ancestor plan was written before CheMPAS Phase 2 (runtime tracer
allocation). Key differences:

1. **Chemistry tracers are runtime-allocated, not in Registry.xml.** Switching
   from ABBA to Chapman requires only changing `config_micm_file` to point to
   `chapman_simplified.json`. The tracers `qO3`, `qO2`, `qO` will be
   discovered automatically from the MICM config — no Registry.xml edits for
   tracers.

2. **MICM config is the source of truth for chemistry species**, not
   Registry.xml. The ancestor plan's "make Registry.xml authoritative" applies
   only to namelist options, not tracers.

3. **ABBA code removal is simpler.** The generic species coupling (Phase 1) and
   runtime allocation (Phase 2) already eliminated most ABBA-specific code.
   The remaining ABBA artifacts are the `abba.yaml` config file and the
   `state_ref` reference state path.

4. **The `tend_` prefix convention** for `scalars_tend` constituent names is
   already in place.

## Scope Decisions (From Ancestor, Still Valid)

- Target real-world meshes first (correctness over early optimization).
- Keep `state_ref` chemistry path (useful for advection diagnostics).
- Phase 1 uses per-column photolysis (no SZA-only binning).
- Primary test domain top is ~20 km. Full stratospheric Chapman validation
  deferred to Phase 6 high-top/offline tests.
- Optional supercell lightning-NOx experiment deferred to Phase 8.

## Phase 1 Mechanism: Simplified Chapman

5 reactions, no O1D:

| Reaction | Type | Notes |
|----------|------|-------|
| O + O3 -> 2 O2 | ARRHENIUS | |
| O + O2 + M -> O3 + M | ARRHENIUS | termolecular |
| O2 + hv -> 2 O | PHOTOLYSIS (`jo2_b`) | |
| O3 + hv -> O + O2 | PHOTOLYSIS (`jo3_a`) | |
| O3 + hv -> O + O2 | PHOTOLYSIS (`jo3_b`) | redirected O1D channel |

Runtime-discovered MPAS tracers: `qO3`, `qO2`, `qO`

`M` is air number density from `state%conditions(:)%air_density`, not a
transported tracer.

## Architecture (Phase 1)

```
MPAS Atmosphere Core
====================

Dynamics -> Physics -> Chemistry (mpas_atm_chemistry.F)
                           |
                           v
                     mpas_musica.F
                     =============
                     1) Extract MPAS state (T, P, rho, tracers, zgrid, lat/lon, time)
                     2) Build/update per-column TUV-x profiles
                     3) Run TUV-x per column -> j-values
                     4) Map j-values -> MICM rate_parameters
                     5) MPAS tracers -> MICM concentrations
                     6) MICM solve (adaptive internal sub-steps)
                     7) MICM concentrations -> MPAS tracers
```

## Phases

### Phase 0: Baseline and Code Hardening

**Goal:** Clean baseline for Chapman chemistry without TUV-x.

**Implementation:**
1. Create `chapman_simplified.json` MICM config with fixed (constant)
   photolysis rates — no TUV-x yet. Switch `config_micm_file` to point to
   it. Runtime tracer discovery handles `qO3/qO2/qO` automatically.
2. Verify the coupled run produces physically plausible Chapman chemistry
   with constant photolysis (O3 photolysis, O recombination, steady state).
3. Add remaining MUSICA namelist options to Registry.xml:
   - `config_tuvx_config_file`
   - `config_chemistry_latitude` / `config_chemistry_longitude`
   - `config_musica_internal_max_step_s`
   - `config_musica_solver_fallback`
4. Keep `state_ref` for now (useful for diagnosing advection effects).
5. Fix chemistry loops to use `nCellsSolve` instead of `nCells` (exclude halo
   cells in MPI runs).
6. Code hardening:
   - Eliminate per-timestep allocate/deallocate churn for work arrays.
   - Add warning log on pressure fallback.
   - Add explicit `musica_finalize` cleanup.
   - Cache species/rate-parameter indices at init.

**Exit criteria:**
- Build passes, initializes with Chapman tracers via runtime discovery.
- Chapman chemistry with fixed photolysis rates produces physically plausible
  O3/O2/O steady-state behavior.
- Tracer fields remain physically bounded (no negative qO3, qO2, qO).

### Phase 1: Solar Geometry and Day/Night Physics

**Goal:** Physically correct solar zenith angle for chemistry.

**Implementation:**
1. Per-cell SZA from model time + `latCell/lonCell`.
2. Idealized fallback using `config_chemistry_latitude/longitude`.
3. Nighttime mask (SZA >= 90 deg).

**Exit criteria:**
- Plausible diurnal SZA cycle across latitudes/seasons.
- Nighttime photolysis rates numerically near zero.

### Phase 2: TUV-x Radiative Transfer (Diagnostic Mode)

**Goal:** Validate TUV-x photolysis physics before coupling to chemistry.

**Implementation:**
1. Add `mpas_tuvx_setup.F`, compile/link in MUSICA makefiles.
2. Initialize TUV-x with MPAS vertical grid per column.
3. Update profiles every chemistry step from MPAS state.
4. Run TUV-x per column, log/store j-values only (no MICM update yet).

**Exit criteria:**
- j-value relative error <= 10% vs standalone MUSICA column results.
- No unphysical j discontinuities.

### Phase 3: Photolysis-to-MICM Coupling

**Goal:** Mapped photolysis rates drive correct oxygen chemistry.

**Implementation:**
1. Cache `PHOTO.jo2_b`, `PHOTO.jo3_a`, `PHOTO.jo3_b` indices at init.
2. Inject per-column j-values into `state%rate_parameters`.
3. Remove default `rate_parameters = 1.0` for photolysis entries.

**Exit criteria:**
- Oxygen budget drift <= 0.5% in closed-window tests.
- Correct sign of diurnal tendencies (dO3/dt < 0 in photolysis regime).
- No top-level O/O3 artifact growth.

### Phase 4: Solver Robustness

**Goal:** Stable integration of stiff oxygen chemistry.

**Implementation:**
1. Adaptive sub-stepping with configurable max internal step.
2. Non-convergence and zero-progress guards.
3. Fallback solver path (BackwardEulerStandardOrder).

**Exit criteria:**
- Chemistry convergence >= 99.9% of calls.
- No infinite-loop sub-step events.

### Phase 5: Real-World Robustness

**Goal:** Reproducibility across decomposition and runtime settings.

**Exit criteria:**
- Cross-decomposition differences within tolerance.
- Phase-gate diagnostics pass in automated regression.

### Phase 6: Full Chapman (O1D)

**Goal:** Reintroduce O1D pathway with high-top/offline validation.

### Phase 7: Performance Optimization

**Goal:** Reduce wall-clock cost after physics is validated.

### Phase 8 (Optional): Lightning-NOx Experiment

**Goal:** Localized NOx production and O3 titration in supercell test.

## Key Constraints

1. **Vertical grid from MPAS** — TUV-x height edges from `zgrid(:, iCell)`,
   no hardcoded grids.
2. **Profiles from MPAS state** — No static atmosphere data files.
3. **Photolysis mapping must match MICM names exactly** — see mapping table
   in ancestor plan.
4. **Domain-top limitation** — 20 km top only samples UTLS; full Chapman
   validation requires high-top or offline.
5. **Source terms through MICM** — Lightning-NOx via MICM rate parameters,
   not direct tracer tendency hacks.

## Phase Gate Runbook

See ancestor plan (`MPAS-Model-ACOM-dev/PLAN_TUVx.md`) for:
- Required namelist block
- Runtime settings for gate runs
- Required output variables
- Log tag conventions
- Pass/fail script specifications
- Phase-by-phase gate matrix
- Artifact storage convention

These will be adapted to CheMPAS conventions as each phase begins.

## Dependencies

- MUSICA-Fortran with TUV-x support
- TUV-x v5.4 data files (cross sections, quantum yields, solar data)
- TUV-x config file (e.g., `tuv_5_4.json`)
- MICM simplified Chapman config (Phase 1) and full Chapman config (Phase 6)
- Python `netCDF4` for phase-gate scripts

## Planned New Files

| File | Purpose |
|------|---------|
| `src/core_atmosphere/chemistry/musica/mpas_tuvx_setup.F` | TUV-x init and profile update |
| `scripts/check_tuvx_phase.py` | Phase-gate pass/fail checks |
| `scripts/run_tuvx_phase_gate.sh` | Phase wrapper for gate checks |

## Planned Modified Files

| File | Changes |
|------|---------|
| `mpas_musica.F` | TUV-x lifecycle, rate-parameter mapping, adaptive sub-stepping |
| `mpas_atm_chemistry.F` | Extract runtime fields for TUV-x, fix nCellsSolve loops |
| `Registry.xml` | Add MUSICA namelist options (not tracers — those are runtime) |
| `chemistry/musica/Makefile` | Compile/link mpas_tuvx_setup.o |
