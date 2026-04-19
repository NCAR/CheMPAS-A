# Chapman + NOx explodes qO3 on chem_box — open issue

## Problem

On the 64-cell `chem_box` mesh, running `chapman_nox.yaml` +
`tuvx_chapman_nox.json` (all four photolysis rates) produces an
explosion in qO3 within the first chemistry step:

```
t=0    qO3 max=23,981 ppb   qO max=0     (AFGL init)
t=5m   qO3 max=93,929,448 ppb   qO max=1.9e7 ppb
t=10m  qO3 max=105,028,599 ppb  qO max=1.9e7 ppb
```

qO3 saturates near 0.15 kg/kg (≈ 10 % mass fraction of air) — about
four orders of magnitude too large. qO pegs at ≈ 0.01 kg/kg.

Per-cell mass balance is preserved (total oxygen atoms conserved),
and the saturation value stays the same across several yaml tweaks,
which looks like a fixed point of whatever the mechanism is
actually solving — not a blow-up.

## Test case (64-cell chem_box)

All files are tracked in the repo so the run is reproducible:

- `test_cases/chem_box/namelist.atmosphere` — `config_dt = 3.0 s`,
  10–60 min run, 8-MPI-rank partition, `config_micm_file =
  chapman_nox.yaml`, `config_tuvx_config_file = tuvx_chapman_nox.json`,
  `config_tuvx_top_extension = .true.`,
  `config_tuvx_extension_file = tuvx_upper_atm.csv`.
- `test_cases/chem_box/streams.atmosphere`,
  `streams.init_atmosphere`, `stream_list.atmosphere.output`,
  `supercell_zeta_levels.txt`, `namelist.init_atmosphere`.
- `micm_configs/chapman_nox.yaml` — full Chapman + NOx mechanism (7
  species incl. M, O1D; 9 reactions; 4 photolysis channels).
- `micm_configs/chapman_nox_noO1D.yaml` — reduced variant with O1D
  folded into jO3_O1D → O (see hypothesis #2 below). Same 4-photolysis
  wiring so `tuvx_chapman_nox.json` works unchanged.
- `micm_configs/tuvx_chapman_nox.json` — TUV-x config: 4 photolysis
  reactions (jO2, jO3_O1D, jO3_O, jNO2), air/O2/O3 radiators all
  driven from MPAS host state.
- `micm_configs/tuvx_upper_atm.csv` — 10-layer 50→100 km US-Std-Atm
  extension (z, T, n_air, n_O3) appended above the MPAS column by
  `mpas_tuvx.F`.

Run directory: `~/Data/CheMPAS/chem_box/`. Init file
(`chem_box_init.nc`) has AFGL qO3 (max ≈ 24 ppb at surface, ≈ 8 ppm
at the stratopause), qO2 = 0.23141 kg/kg (dry air VMR 0.2095), qNO
and qNO2 with small smooth profiles, qO = qO1D = 0.

### Mesh and grid build (how the 64-cell case was constructed)

The chem_box mesh is a doubly-periodic planar hexagonal mesh, 8×8
cells at 500 m centre-to-centre spacing (domain 4000 m × 3464 m,
`on_a_sphere = "NO"`, `is_periodic = "YES"`). Vertical grid is 60
stretched levels to 50 km, reusing the supercell zeta profile. All
tools come from the `mpas` conda env (`mpas-tools` + `metis`
packages on conda-forge).

Exact reproduction from scratch, run in the chem_box run directory:

```bash
# 1. Build the 64-cell periodic hex mesh (8x8 at 500 m).
#    planar_hex writes an MPAS-format grid directly — no conversion
#    needed for atmosphere_model input.
~/miniconda3/envs/mpas/bin/planar_hex --nx 8 --ny 8 --dc 500 \
    -o chem_box_grid.nc

# 2. Derive the METIS graph description from the grid. MpasMeshConverter
#    writes a graph.info file as a side effect of converting; we use
#    the graph.info and discard the "converted" grid file since the
#    planar_hex output already matches MPAS format exactly.
~/miniconda3/envs/mpas/bin/MpasMeshConverter.x \
    chem_box_grid.nc /tmp/_converted.nc
mv graph.info chem_box.graph.info
rm -f /tmp/_converted.nc

# 3. Partition the graph for the MPI rank counts you'll run.
#    gpmetis is from conda-forge's `metis` package. The partition is
#    randomized, so the .part.N files won't byte-match another run,
#    but any valid 8-way partition works identically at runtime.
~/miniconda3/envs/mpas/bin/gpmetis chem_box.graph.info 8
~/miniconda3/envs/mpas/bin/gpmetis chem_box.graph.info 4

# 4. Copy the reference configs from the repo.
cp ~/EarthSystem/CheMPAS/test_cases/chem_box/* .

# 5. Initialize. init_case = 5 is the supercell (Weisman-Klemp)
#    thermodynamic sounding; chem_box reuses it for its moisture
#    column and as a background state for chemistry. Runs with
#    any of the partitions you generated (use -n matching a
#    .part.N file).
ln -sf ~/EarthSystem/CheMPAS/init_atmosphere_model .
mpiexec -n 4 ./init_atmosphere_model

# 6. Stage chemistry config + atmosphere executable.
ln -sf ~/EarthSystem/CheMPAS/atmosphere_model .
ln -sf ~/EarthSystem/CheMPAS/LANDUSE.TBL .
cp ~/EarthSystem/CheMPAS/micm_configs/chapman_nox.yaml .
cp ~/EarthSystem/CheMPAS/micm_configs/tuvx_chapman_nox.json .
cp ~/EarthSystem/CheMPAS/micm_configs/tuvx_upper_atm.csv .
```

Run the model with `mpiexec -n 8 ./atmosphere_model`. Outputs land
in `output.nc` (one snapshot per `output_interval` in
`streams.atmosphere`; set to `00:00:03` when debugging step-by-step
photolysis behaviour, `00:05:00` for normal runs).

Global attributes on `chem_box_grid.nc` captured by `ncdump -h` for
verification: `dc = 500.0 m`, `nx = 8`, `ny = 8`, `x_period = 4000 m`,
`y_period ≈ 3464.1 m`, 64 cells, 192 edges, 128 vertices.

## What was verified (facts, not hypotheses)

**Rate constants are correct.** Every ARRHENIUS rate, at T = 220 K,
matches JPL 2020 within 1 % after the cgs → SI conversion
(`×N_A × 10⁻⁶` for bimolecular, `×N_A² × 10⁻¹²` for termolecular).
Compared against the known-good `lnox_o3.yaml` `NO + O3` rate as a
ground-truth calibration point.

**MICM Arrhenius form** (from
`MUSICA/build/_deps/micm-src/include/micm/process/rate_constant/arrhenius_rate_constant.hpp`):
`k = A × exp(C/T) × (T/D)^B × (1 + E·P)` with defaults `D = 300`,
`E = 0`. Omitting D/E from the YAML is fine.

**Third-body M is handled correctly by MICM.** Earlier notes here
claimed `is_third_body` was parsed but never consumed. That was
wrong. Trace through the installed MUSICA library:

- `MUSICA/src/micm/v1_parse.cpp:38-41` → calls
  `s.SetThirdBody()` when the YAML sets `is third body: true`.
- `micm-src/include/micm/system/species.hpp:162` →
  `SetThirdBody` installs `parameterize_ = [](c){ return c.air_density_; }`.
- `micm-src/include/micm/system/phase.hpp:66-90` → parameterized
  species are filtered out of `StateSize`, `UniqueNames`, and
  `SpeciesNames`.
- `micm-src/include/micm/process/chemical_reaction.hpp:127-132` →
  during rate-constant evaluation, parameterized reactants are
  multiplied into `fixed_reactants` via `parameterize_(conditions)`,
  so R2 = A·(T/300)^B · air_density · [O] · [O2] (correct termolecular).
- Runtime log (`~/Data/CheMPAS/chem_box/log.atmosphere.0000.out:583-590`)
  confirms `state%species_ordering` lists **6** species (NO2, NO, O,
  O3, O2, O1D) — M is correctly excluded.
- `air_density_` units are mol m⁻³ (`conditions.hpp:14`), matching
  what the coupler sets (`rho / M_AIR`).

**Mesh & dynamics are fine.** `chem_box` with `lnox_o3.yaml` +
`tuvx_no2.json` runs cleanly for 10 min (qO3 stable at 50 ppb,
NO at 0). The same mesh with `chapman_nox.yaml` explodes in 5 min.

**TUV-x photolysis rates look reasonable by altitude** (extracted
from the `j_*` diagnostic fields in `output.nc` at t=5 min):

```
    z (km)          jO2        jO3_O      jO3_O1D         jNO2
     0.1    0.000e+00    2.758e-11    0.000e+00    2.786e-06
     9.2    0.000e+00    3.877e-10   3.723e-225    3.722e-05
    21.5    0.000e+00    2.393e-08    1.494e-53    1.390e-03
    35.4    0.000e+00    4.564e-06    2.150e-09    5.973e-03
    49.5    1.736e-10    4.407e-04    4.242e-04    8.438e-03
```

`jO2` is effectively zero everywhere except the MPAS top, as expected
(Schumann-Runge bands attenuated out). This is the **first** CheMPAS
run that has O2 photolysis in the mechanism; prior runs had `jNO2`
only via `lnox_o3.yaml`.

**Variants tried (all still explode):**
1. `A: 8.018e-17` (original, cgs units) — R2 too slow, qO3 DECAYS
   exponentially to ~0.4 ppb in 1 h.
2. `A: 6.0e-34` + M as reactant — same behaviour as #1 (wrong scale).
3. `A: 2.18e2` + M as reactant — qO3 EXPLODES as above.
4. Copied from canonical `MUSICA/configs/v1/chapman/config.yaml`
   (`A: 217.6, B: -2.4, D: 300, E: 0`) — same explosion as #3.
5. Removed `__tracer type: CONSTANT` from O2 — qO2 now depletes to
   25 % of initial, confirming mass balance, but qO3 still at same
   saturation value.
6. `chapman_nox_noO1D.yaml` — O1D dropped, jO3_O1D channel wired to
   produce O(3P) directly (R6 quenching collapsed into photolysis).
   **Output bit-for-bit identical to the full mechanism**; rules out
   O1D stiffness / zero-initial-value as the trigger.

## What the standalone MICM column model shows

`MUSICA/fortran/examples/column_model.F90` ships as a coupled TUV-x
+ MICM 1-D Chapman driver (120 layers, 0–120 km, Boulder CO, June
21, 15-min steps over 24 h). Rebuilt with `MUSICA_ENABLE_TESTS=ON`
and run against the canonical `MUSICA/configs/v1/chapman/config.yaml`:

```
t=0h : O3_max = 6.902e-06 mol/m³ at z = 21.5 km (AFGL-like init)
t=6h : O3_max = 6.902e-06 mol/m³              O_max = 5.5e-9
t=12h: O3_max = 6.929e-06 mol/m³              O_max = 5.7e-8
t=18h: O3_max = 6.995e-06 mol/m³              O_max = 1.1e-7
t=24h: O3_max = 6.995e-06 mol/m³              O_max = 6.3e-8
```

Stratopause O3 stays pinned at a physical ~15 ppm over 24 h. **No
explosion.** So the bug is NOT in MICM's Chapman mechanism itself —
it is specific to the CheMPAS + MUSICA coupling.

## The smoking gun: step-1 photolysis rates

Setting `output_interval = 00:00:03` (one snapshot per chemistry
step) and printing `[MUSICA] jO3_O min/max` at every step reveals a
dramatic anomaly at **step 1 only**:

```
step  1: jO3_O min= 3.634e-04  max= 5.627e-04   (nearly UNIFORM, TOA value everywhere)
step  2: jO3_O min= 1.661e-09  max= 4.518e-04   (normal column attenuation)
step  5: jO3_O min= 7.097e-10  max= 4.404e-04
step 10: jO3_O min= 3.433e-10  max= 4.404e-04
step 50: jO3_O min= 5.767e-11  max= 4.405e-04
step 200: jO3_O min= 1.375e-11  max= 4.409e-04   (correct surface attenuation)
```

At step 1, **every cell and level** is receiving the top-of-
atmosphere rate (no O3 column attenuation). From step 2 onward, the
profile is attenuated correctly. The column in `j_jO3_O` diagnostics
in `output.nc` at t=3 s matches the `[MUSICA]` log output for step 1
(3.6e-4 uniform), confirming this is not a logging artifact.

Given the rate constant R2 ≈ 45 000 × [O2] [M] s⁻¹ at surface, a
bootstrap pulse of O atoms from the step-1 UV blast is enough to
produce ~0.02 kg/kg of qO3 in one 3-s step (observed). On step 2
the now-huge O3 column heavily attenuates UV (proper operation of
TUV-x), but the damage has been done — the system has already
landed in the wrong attractor and cannot recover.

## Debug instrumentation confirmed the input to TUV-x is correct

Temporary instrumentation in `tuvx_compute_photolysis` (reverted
after use) logged the profile values handed to TUV-x at the first
few calls:

```
[DEBUG tuvx] call=0  air_nd(1)=2.36e19   temp_mid(1)=299 K
                     o3_nd(1)=6.34e11 molec/cm³  (≈ 20 ppb surface)
                     o3_ld col_sum=9.18e18 molec/cm²  (≈ 340 DU — realistic!)
```

Readback via `o3_profile%get_layer_densities` after
`set_layer_densities` returned **exactly** the same numbers we
passed in, so TUV-x's internal profile storage is populated
correctly on the first call.

**Warmup test** — calling `tuvx_solver%run()` twice at the first
chemistry step produced identical photolysis rates both times, so
this is NOT a lazy-initialization / first-call-is-stale issue
inside TUV-x.

Something downstream of the correctly-set profile is producing the
wrong radiative transfer result at step 1. Expected surface OD for
a 340 DU column at ~300 nm is ~3.7 → transmission ~2.5 %. TUV-x is
giving transmission ~26 % at the surface (OD ~1.3), about 12× less
attenuation than physics would predict.

## Hypotheses still on the table

1. **TUV-x radiator update ordering.** `radiator_from_host_t::update_state`
   in `MUSICA/build/_deps/tuvx-src/src/radiative_transfer/radiators/from_host.F90`
   is a no-op. The "base" radiator (used by air/O2/O3 in our TUV-x
   config) has an `update_state` in `radiator.F90:169-230` that reads
   `radiator_profile%layer_dens_` to compute `layer_OD_`, and is
   called from `radiative_transfer%calculate:170-177`. In principle
   this should run on every `tuvx_solver%run()`. Worth confirming
   that on step 1 the base radiator's `layer_OD_` actually reflects
   our provided `layer_dens_`.
2. **MICM photolysis scaling_factor.** `UserDefinedRateConstant::Calculate`
   returns `custom_parameter × scaling_factor_`. Default is 1.0, but
   if our YAML's photolysis reactions inherit a different factor, the
   effective rate would be amplified. Worth grepping `scaling_factor`
   in the running build and logging the actual `k` MICM computes
   per photolysis reaction at step 1.
3. **Profile edge values never set for O3/O2.** `tuvx_compute_photolysis`
   calls `set_midpoint_values` + `set_layer_densities` for O3 and O2,
   but not `set_edge_values`. `edge_val_` stays at the constructor
   default of zero. The canonical `profile_from_host_t::update()`
   (not exposed via C bindings) REQUIRES edge_values. Possibly the
   radiative transfer solver reads edge_val_ somewhere.
4. **`exo_layer_dens_(N+1)` never updated.** `internal_set_layer_densities`
   populates `layer_dens_(1:N)` and `exo_layer_dens_(1:N)` but leaves
   `exo_layer_dens_(N+1)` at the constructor zero. The canonical
   `update()` sets it from `edge_val_(N+1) * hscale_`. Used by
   spherical geometry air-mass in `radiative_transfer%calculate`
   (air profile only). Probably not the bug — air column has zero
   above the extension anyway — but worth ruling out.
5. **Rosenbrock stiffness failure.** With `__absolute tolerance: 1e-12`
   on O/O1D and an initial state of exactly zero, the first Jacobian
   could be degenerate. Looser tolerances (`1e-10`) and seeding a
   small non-zero qO haven't been tried yet. But removing O1D
   entirely (variant #6 above) didn't help, which argues against
   pure solver stiffness as the culprit.

## Current state shipped

- `chapman_nox.yaml` in canonical MICM-v1 form with SI rate
  constants (correct per JPL); the explosion is reproducible.
- `chapman_nox_noO1D.yaml` reduced mechanism; also explodes.
- `chem_box` test case can be used for multi-day LNOx diurnal-cycle
  work with `lnox_o3.yaml` while this is being debugged.
- All multi-photolysis plumbing (mpas_tuvx.F multi-rate,
  `musica_set_photolysis_rates`, `j_<name>` diagnostic fields) is
  working correctly — confirmed by TUV-x rates appearing in the
  output with sensible altitude dependence from step 2 onward.
- Standalone MICM column model (`MUSICA/fortran/examples/`) is
  built and runs the canonical Chapman config successfully — use
  this as a reference oracle when debugging.

## Path forward

Two parallel lines of investigation:

**Line A — confirm MICM rate constants at step 1.** Briefly
instrument MICM (or dump `state%rate_constants_`) to confirm the
effective `k` for each reaction at step 1 is what we expect from
the j-values we fed in. Rules in/out hypothesis #2.

**Line B — accept step-1 photo rates as the trigger; ask why
TUV-x produces TOA rates for a 340-DU column.** Candidates:
hypothesis #1 (radiator refresh), #3 (missing edge values), #4
(exo layer density). Could also be reproduced in a pure-TUV-x
standalone driver fed the same profile.

Lower-risk workaround if neither A nor B resolves quickly: skip
chemistry on the very first step (let dynamics run once, then
allow chem to start). That sidesteps the bootstrap without
addressing the root cause. Not a real fix — flag it as mitigation
only.
