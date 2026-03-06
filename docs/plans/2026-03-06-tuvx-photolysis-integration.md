# Photolysis and Tropospheric Chemistry Integration Plan

## Document Status

- `Historical Context:` Adapted from ancestor project plans — the TUV-x
  photolysis plan (`MPAS-Model-ACOM-dev/PLAN_TUVx.md`) and the DAVINCI
  lightning-NOx/O3 mechanism (`DAVINCI-MPAS/PLAN.md` Phase 6, `SCIENCE.md`).
- `Current State:` Planning stage. Phase 0 next.
- `Use This As:` Primary reference for post-ABBA chemistry development.

## Locked Decisions (2026-03-06)

1. **NOx sink path:** Carry over NOx loss via MICM-configured sink terms (not
   hard-coded Fortran sink tendencies).
2. **Lightning source path:** Keep source integrated in the chemistry solver
   path (no pre-solver operator-split source update).
3. **Verification:** Include a dedicated controlled-NOx pulse/equilibration
   case for Leighton-ratio validation.

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
   `lnox_o3.json` MICM config, point `config_micm_file` at it, and `qNO`,
   `qNO2`, `qO3` appear automatically. (Source/sink coupling details still
   require focused chemistry-coupler implementation work.)

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

1. **Integrated source in solver path** — Implement lightning source as a MICM
   source reaction with per-cell rate parameter updates each chemistry step.
2. **Sink in MICM config** — Represent NOx loss through mechanism-defined
   first-order sink terms (e.g., NO and NO2 sink pathways) with configurable
   timescale `tau`.
3. **No operator-split source path in Phase 0** — Do not inject source as a
   pre-MICM tracer tendency update.

### Implementation Steps

1. Create `lnox_o3.json` MICM mechanism config with:
   - Species: NO, NO2, O3
   - Arrhenius reaction: NO + O3 → NO2 + O2
   - Photolysis reaction: NO2 + hv → NO + O3 (fixed rate)
   - Lightning source reaction(s) for NO with per-cell rate parameter input
   - NOx sink reaction(s) with configurable `tau`
   - Molar masses for each species

2. Point `config_micm_file` at the new config. Runtime tracer discovery
   gives us `qNO`, `qNO2`, `qO3` automatically.

3. **Case A (equilibrium diagnostic):** initialize O3 to uniform background
   (~50 ppbv) and apply a small controlled NOx pulse/seed (non-zero NO and NO2),
   with lightning source off, to verify Leighton behavior and Ox conservation.

4. **Case B (storm chemistry):** initialize O3 background (~50 ppbv), NO=NO2=0,
   lightning source on, to verify O3 titration/recovery patterns.

5. Add/verify diagnostic conversions for ppbv reporting from MPAS tracer mixing
   ratios so verification thresholds are unit-consistent.

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
- Source/sink behavior matches configured mechanism terms (no operator-split path).

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
4. **Source/sink terms through MICM solver path** — Lightning-NOx source and
   NOx sink both represented through MICM mechanism/rate-parameter pathways.
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
- MICM LNOx-O3 mechanism config (to be created)
- TUV-x support in MUSICA-Fortran (Phase 2+)
- Python `netCDF4` for verification scripts
