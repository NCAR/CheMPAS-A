# Photolysis and Tropospheric Chemistry Integration Plan

## Document Status

- `Historical Context:` Adapted from ancestor project plans — the TUV-x
  photolysis plan (`MPAS-Model-ACOM-dev/PLAN_TUVx.md`) and the DAVINCI
  lightning-NOx/O3 mechanism (`DAVINCI-MPAS/PLAN.md` Phase 6, `SCIENCE.md`).
- `Current State:` Phase 0 complete (2026-03-06). LNOx-O3 mechanism runs
  end-to-end with correct unit conversion, prescribed photolysis (j_NO2),
  configurable NOx sink (tau), and w_ref safety guard. Case B (storm chemistry)
  passes with 0.5 ppbv/s source: O3 titration, NO2 production, photolysis
  recycling all verified. Case A (equilibrium) runs but rigorous Ox conservation
  requires transport-disabled test. Docs updated.
- `Use This As:` Primary reference for post-ABBA chemistry development.

## Locked Decisions (2026-03-06)

1. **NOx sink path:** Carry over NOx loss via MICM-configured sink terms (not
   hard-coded Fortran sink tendencies).
2. **Lightning source path:** Operator-split pre-MICM injection via standalone
   `mpas_lightning_nox.F` module. MICM doesn't support spatially varying source
   terms, so injection happens before solve.
3. **Verification:** Include a dedicated controlled-NOx pulse/equilibration
   case for Leighton-ratio validation.
4. **MICM units:** Concentrations in mol/m³. Arrhenius `A` must be in m³/mol/s
   (convert from cm³/molecule/s via `A × Nₐ × 10⁻⁶`). Photolysis rate
   parameters in s⁻¹, set externally.
5. **Initial conditions:** Set via MPAS init file (kg/kg), not MICM config.
   MICM `__initial concentration` is required by the parser but overwritten
   by state transfer; set to 0.

## Strategic Direction

The ABBA mechanism validated the MPAS-MICM coupling infrastructure and the
Phase 2 runtime tracer allocation. It served its purpose brilliantly. Now we
need a mechanism that is *scientifically meaningful* in our supercell test
domain.

**Key insight from the DAVINCI sister project:** The supercell domain is
tropospheric (0–20 km). Chapman oxygen photolysis (`j_O2`) is negligible below
the stratosphere — it would be a plumbing test, not real science. But
lightning-NOx-driven O3 photochemistry is *exactly* what happens in supercell
thunderstorms. This is the natural next mechanism for CheMPAS.

The DAVINCI project has already implemented and validated this system through
Phase 6 (LNOx-O3 tropospheric photochemistry, code complete on the same
supercell domain). That work is gold — it gives us a proven 3-species
mechanism, analytical verification targets, and DC3 observational validation
data to aim for.

## Why LNOx-O3 Instead of Chapman

We love this direction for five reasons:

1. **Physically meaningful in the supercell domain.** O3 titration where
   lightning injects NO, O3 recovery downwind as NO2 photolyzes — these are
   real, observable phenomena in thunderstorms.

2. **Only 3 prognostic species with compact chemistry.** Simple enough to
   verify analytically (Leighton photostationary ratio), but scientifically
   interesting.

3. **Runtime tracer discovery handles the species registration.** Create an
   `lnox_o3.yaml` MICM config, point `config_micm_file` at it, and `qNO`,
   `qNO2`, `qO3` appear automatically.

4. **DC3 validation targets exist.** The Deep Convective Clouds and Chemistry
   campaign (Barth et al., 2015) provides real observational data for exactly
   this scenario — see `DAVINCI-MPAS/DC3.md` for the full reference.

5. **Proven in the sister project.** DAVINCI Phase 6 implemented this system
   with an internal ODE solver on the same supercell mesh. The chemistry,
   Jacobian, conservation properties, and verification criteria are all worked
   out.

## Phase 0: LNOx-O3 with Fixed Rates (No TUV-x Yet)

**Goal:** Replace ABBA with a 3-species NO/NO2/O3 tropospheric photochemistry
mechanism using constant (prescribed) photolysis and reaction rates. Validate
that the runtime tracer machinery works end-to-end with the new species.

### Chemistry

The tropospheric photostationary system (from DAVINCI Phase 6):

```
NO + O3 → NO2 + O2     (k: Arrhenius, temperature-dependent)
NO2 + hv → NO + O3     (j_NO2: prescribed constant, ~0.01 s⁻¹ daytime)
```

The net O3 reaction is `NO2 + hv → NO + O + (O + O2 + M → O3)`, written as a
single step because atomic O reaches steady state in microseconds in the
troposphere.

### ODE System

```
d[NO]/dt  =  j*[NO2] - k*[NO]*[O3] + S_ltg - [NO]/tau
d[NO2]/dt = -j*[NO2] + k*[NO]*[O3] - [NO2]/tau
d[O3]/dt  =  j*[NO2] - k*[NO]*[O3]
```

where `S_ltg` is the lightning NO source term and `tau` is a configurable NOx
sink timescale represented in the MICM mechanism.

**Conservation:** Ox = [O3] + [NO2] is conserved when both lightning source and
NOx sink are disabled.

### Source/Sink Coupling Strategy (Decided)

1. **Operator-split lightning source** — A standalone Fortran module
   (`mpas_lightning_nox.F`) injects NO into MPAS scalars before MICM runs each
   timestep. Source scales linearly with updraft strength:
   `S = rate * max(0, w - w_thr) / w_ref` in cells where altitude is in range.
   The configured `rate` is reached when `w - w_thr = w_ref`.
   This is the same approach used in DAVINCI-MPAS, adapted to work with MPAS
   pools directly (no DAVINCI state types). The module is a complete no-op if
   `qNO` is not present in the mechanism or `config_lnox_source_rate = 0`.

2. **Chemistry in MICM** — The NO/NO2/O3 reactions (Arrhenius titration +
   photolysis) are defined in the MICM config. MICM handles the ODE solve.

3. **Sink in MICM config** — NOx loss is represented by mechanism-defined
   first-order sink terms with configurable timescale `tau`. The current
   `lnox_o3.yaml` includes `FIRST_ORDER_LOSS` reactions for NO and NO2, and
   `mpas_musica.F` sets their rates to `1/tau` when `config_lnox_nox_tau > 0`.

4. **Namelist control** — Seven parameters in `&musica`:
   `config_lnox_source_rate`, `config_lnox_w_threshold`, `config_lnox_w_ref`,
   `config_lnox_z_min`, `config_lnox_z_max`, `config_lnox_j_no2`, and
   `config_lnox_nox_tau`.

### Implementation Progress

- [x] `mpas_lightning_nox.F` — lightning source module (init + inject)
- [x] Registry.xml — namelist parameters for lightning source
- [x] `mpas_atm_chemistry.F` — hook lightning init and inject into chemistry pipeline
- [x] `chemistry/Makefile` — build integration with dependency ordering
- [x] Build passes with MUSICA=true
- [x] `micm_configs/lnox_o3.yaml` — MICM config (YAML v1 format) with:
   - Species: NO, NO2, O3 (molar masses, `__initial concentration: 0`)
   - Arrhenius reaction: NO + O3 → NO2 (A=1.084e6 m³/mol/s, C=-1370)
   - Photolysis reaction: NO2 + hv → NO + O3 (rate set externally)
- [x] Runtime tracer discovery: `qNO`, `qNO2`, `qO3` created automatically
- [x] `scripts/init_lnox_o3.py` — initialize tracers in supercell_init.nc
- [x] `scripts/plot_lnox_o3.py` + `scripts/style.py` — visualization suite
- [x] `mpas_musica.F` — default photolysis rate parameters to 0 (was 1.0)
- [x] Arrhenius A parameter corrected from cm³/molecule/s to m³/mol/s
- [x] `micm_configs/lnox_o3.yaml` — added FIRST_ORDER_LOSS reactions for NO
  and NO2 with rate = 1/tau (set externally via `config_lnox_nox_tau`)
- [x] **Case B (storm chemistry):** 15-min supercell run with lightning source
   produces visible O3 titration in updraft core. NO peaks ~3500 ppbv (at
   artificially high 50 ppbv/s source rate), NO2 produced via Arrhenius, O3
   depleted to near zero where NO is injected.
- [x] **Case A (equilibrium diagnostic):** Ran with j_NO2=0.01, source off,
  uniform NOx+O3 init. Chemistry directionally correct; rigorous Ox
  conservation requires transport-disabled test.
- [x] Tuned source rate to 0.5 ppbv/s — physically realistic, ~28 ppbv NO peak
- [x] ppbv conversions verified: q [kg/kg] × (M_air/M_species) × 1e9
- [x] Corrected lightning source: delta_q = rate × scale × dt × 1e-9 × (M_NO/M_AIR)
- [x] Added `w_ref <= 0` guard to prevent divide-by-zero
- [x] Updated `RUN.md` and `TEST_RUNS.md` with LNOx-O3 workflow and results

### Verification Criteria (from DAVINCI Phase 6)

| Check | Criterion |
|-------|-----------|
| Non-negativity | qNO >= 0, qNO2 >= 0, qO3 >= 0 everywhere |
| O3 background | O3 ≈ 50 ppbv away from storm |
| O3 titration | O3 depressed in updraft core where NO is injected |
| O3 recovery | O3 rebounds away from fresh NO source (NO2 photolysis) |
| Leighton ratio | In Case A, φ = j[NO2]/(k[NO][O3]) approaches 1 in photostationary regions |
| Ox conservation | In Case A with source off and very large tau, [O3]+[NO2] conserved to solver tolerance |
| Sink behavior | With finite tau and source off, NOx decays consistently with configured sink timescale |
| Unit consistency | Verification metrics reported in ppbv use documented q↔concentration conversion |

### Exit Criteria

- Build passes, initializes with LNOx-O3 tracers via runtime discovery.
- Case A passes Leighton + Ox checks with controlled NOx initialization.
- 30-minute supercell run (Case B) produces physically plausible O3 titration.
- Tracer fields remain non-negative.
- Lightning source activates only in updraft cores within altitude range.

## Review Findings Incorporated (2026-03-06)

1. **Sink-path gap:** Plan and ODE include NOx sink, but current mechanism
   config does not yet implement sink reactions. Keep sink architecture in MICM
   and add missing mechanism terms.
2. **Photolysis-zero blocker:** `j` defaults to 0 globally; Case A Leighton and
   Ox conservation checks cannot pass until nonzero controlled `j_NO2` is
   available for the diagnostic case.
3. **Architecture consistency fix:** Source path is operator-split pre-MICM;
   sink path is MICM-configured. Plan wording now reflects this split
   consistently.
4. **Source-unit correctness:** Current source path needs explicit conversion
   review/fix (`ppbv/s` to `kg/kg/s`) before quantitative interpretation.
5. **`state_ref` caveat:** With operator-split source active each chemistry
   step, advection-vs-chemistry attribution using `state_ref` must account for
   source contamination (or disable source for pure transport diagnostics).
6. **Numerical safety:** Add `w_ref > 0` validation to prevent singular source
   scaling.
7. **Documentation drift:** User-facing run/test docs must be updated to match
   LNOx-O3 workflow and verification outputs.

## Phase 1: Solar Geometry and Day/Night Physics

(Unchanged from previous plan — per-cell SZA from model time + lat/lon,
nighttime mask for j_NO2 = 0.)

## Phase 2: TUV-x Radiative Transfer (Diagnostic Mode)

(Unchanged — validate TUV-x j-values before coupling.)

## Phase 3: TUV-x Coupled Photolysis

Replace fixed j_NO2 with TUV-x-computed values. Add j_O3 channels when
ready for full Chapman extension.

## Later Phases

- Phase 4: Solver robustness under real-world forcing
- Phase 5: Real-world robustness and reproducibility
- Phase 6: Full Chapman (O1D) reintroduction
- Phase 7: Performance optimization
- Phase 8 (Optional): Extended NOx chemistry (PAN, HNO3, organic nitrates)

## Key Constraints

1. **Vertical grid from MPAS** — TUV-x height edges from `zgrid(:, iCell)`.
2. **Profiles from MPAS state** — No static atmosphere data files.
3. **Domain-top limitation** — 20 km top only samples troposphere/UTLS.
4. **Source/sink representation split** — Lightning-NOx source is operator-split
   pre-MICM; NOx sink remains mechanism-defined within MICM.
5. **Keep `state_ref`** — Useful for diagnosing advection effects on chemistry.

## Reference Material

### DAVINCI Sister Project

The DAVINCI project (`~/EarthSystem/DAVINCI-MPAS/`) contains:

- `SCIENCE.md` — Lightning NOx physics, Leighton framework, DC3 findings
- `PLAN.md` Phase 6 — LNOx-O3 mechanism details, ODE system, Jacobian,
  verification criteria
- `DC3.md` — Deep Convective Clouds and Chemistry campaign reference
  (Barth et al., 2015), observational validation targets
- `TUV.md` — TUV-x algorithm and data file reference

### Ancestor TUV-x Plan

`MPAS-Model-ACOM-dev/PLAN_TUVx.md` contains:
- Full 9-phase TUV-x integration plan with physical verification gates
- Phase gate runbook (namelist, runtime settings, pass/fail scripts)
- Fortran implementation sketches
- Photolysis-to-MICM rate parameter mapping

## Dependencies

- MUSICA-Fortran with MICM support (already linked)
- MICM LNOx-O3 mechanism config (`micm_configs/lnox_o3.yaml`, base config done;
  sink extension pending)
- TUV-x support in MUSICA-Fortran (Phase 2+)
- Python `netCDF4` for verification scripts
