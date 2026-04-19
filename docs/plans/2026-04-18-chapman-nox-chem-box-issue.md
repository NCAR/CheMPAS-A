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

## What was verified

**Rate constants are correct.** Every ARRHENIUS rate, at T = 220 K,
matches JPL 2020 within 1 % after the cgs → SI conversion
(`×N_A × 10⁻⁶` for bimolecular, `×N_A² × 10⁻¹²` for termolecular).
Compared against the known-good `lnox_o3.yaml` `NO + O3` rate as a
ground-truth calibration point.

**MICM Arrhenius form** (from
`MUSICA/build/_deps/micm-src/include/micm/process/rate_constant/arrhenius_rate_constant.hpp`):
`k = A × exp(C/T) × (T/D)^B × (1 + E·P)` with defaults `D = 300`,
`E = 0`. Omitting D/E from the YAML is fine.

**Third-body handling.** Searching the micm-src tree for
`is_third_body` shows the flag is parsed by
`mechanism_configuration` but never consumed inside MICM's rate
calculation. `M` is treated as a plain reactant species; putting it
on both sides of R2 does NOT double-count.

**Mesh & dynamics are fine.** `chem_box` with `lnox_o3.yaml` +
`tuvx_no2.json` runs cleanly for 10 min (qO3 stable at 50 ppb,
NO at 0). The same mesh with `chapman_nox.yaml` explodes in 5 min.
So the issue is specifically in the Chapman+NOx mechanism
definition.

**TUV-x photolysis rates look reasonable by altitude:**
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
run that has O2 photolysis in the mechanism; prior runs had
`jNO2` only via `lnox_o3.yaml`.

**Variants tried (all still explode):**
1. `A: 8.018e-17` (original, cgs units) — R2 too slow, qO3 DECAYS
   exponentially to ~0.4 ppb in 1 h
2. `A: 6.0e-34` + M as reactant — same behaviour as #1 (wrong scale)
3. `A: 2.18e2` + M as reactant — qO3 EXPLODES as above
4. Copied from canonical `MUSICA/configs/v1/chapman/config.yaml`
   (`A: 217.6, B: -2.4, D: 300, E: 0`) — same explosion as #3
5. Removed `__tracer type: CONSTANT` from O2 — qO2 now depletes to
   25 % of initial, confirming mass balance, but qO3 still at same
   saturation value

## Hypotheses (not yet tested)

1. **MICM Rosenbrock integrator stiffness failure.** With
   `__absolute tolerance: 1e-12` on O/O1D and an initial state of
   exactly zero, the first Jacobian could be degenerate, pushing the
   step estimator into a wrong attractor. Worth testing looser
   tolerances (`1e-10`) or seeding a small non-zero qO initial value.
2. **jO2 underflow interactions.** The `jO2 = 0` at all but the top
   level is a physical result (Schumann-Runge) but may interact
   weirdly with MICM's PHOTOLYSIS reaction when many cells have
   exactly `rate_parameter = 0`. Worth setting a small floor
   (e.g. `max(jO2, 1e-30)`) in `mpas_atm_chemistry.F`.
3. **Species ordering mismatch.** When the yaml declares 7 species
   but MICM reports them in a different order than the name list
   iteration, there could be a subtle species-index vs state-index
   confusion specific to chapman (not lnox). Requires walking
   through `MICM_from_chemistry` / `MICM_to_chemistry` stride
   calculations with 7 species vs 3.
4. **`__atoms: 3` field on O3.** Canonical config has this; our
   `chapman_nox.yaml` now has it too. If MICM uses it as a
   stoichiometric multiplier somewhere it could affect rate scaling.

## Current state shipped

- `chapman_nox.yaml` is in canonical MICM-v1 form with SI rate
  constants (correct per JPL); the explosion is reproducible
- `chem_box` test case works cleanly with `lnox_o3.yaml` and can
  be used for multi-day LNOx diurnal-cycle work while the Chapman
  issue is being debugged
- All other multi-photolysis plumbing (mpas_tuvx.F multi-rate,
  musica_set_photolysis_rates, j_<name> diagnostic fields) is
  working correctly — confirmed by TUV-x rates appearing in the
  output with sensible altitude dependence

## Path forward

Debug the chapman explosion on a **single-column standalone MICM
test** (no MPAS) using `MUSICA/configs/v1/chapman/config.yaml`.
If that canonical config also explodes on this machine, the bug is
in our MICM build / environment, not the yaml. If it runs correctly,
diff our `chapman_nox.yaml` against the canonical until the
explosion is isolated to one line.

The `MUSICA/fortran/examples/` tree has standalone MICM driver
programs that can be pointed at any config file.
