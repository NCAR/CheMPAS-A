# MPAS-MUSICA Development TODO

## Chemistry Integration

- [ ] **Add namelist option for tracer seeding control**
  - Add `config_micm_seed_tracers` to `&musica` namelist (default `.true.`)
  - When `.false.`, skip MICM->MPAS tracer seeding to preserve custom initial conditions
  - Current workaround: seeding is skipped if tracers have spatial gradients

## MICM Issues

- [ ] **Report FIRST_ORDER_LOSS rate=0 bug to MUSICA team**
  - MICM applies nonzero loss even when FIRST_ORDER_LOSS rate parameters are set to 0
  - Caused 8.67% Ox drift and 48% NOx loss in 2-minute Case A (should be zero)
  - Workaround: omit LOSS reactions from config when sink is disabled
  - Reproduce: run `lnox_o3.yaml` with LOSS reactions, `config_lnox_nox_tau = 0`
  - Compare with config without LOSS reactions — conservation is perfect (0.0000%)
  - File issue at https://github.com/NCAR/musica or contact MUSICA team

## Future Enhancements

- [ ] **Use uxarray for proper unstructured mesh visualization** (low priority)
  - Current approach uses Delaunay triangulation of cell centers
  - uxarray respects actual MPAS Voronoi mesh topology
  - Would eliminate edge artifacts and improve accuracy
  - https://uxarray.readthedocs.io/

- [x] Solar geometry / day-night physics for j_NO2 (Phase 1)
- [x] TUV-x photolysis rate calculations (Phase 2)
- [ ] **TUV-x single-column optimization** — Currently solves TUV-x for all
  28,080 cells every timestep. With fixed SZA and nearly identical profiles
  away from the storm, a single representative column would give ~28,000x
  speedup with negligible accuracy loss. Significant runtime impact observed.
- [x] Domain-integrated Ox conservation test for Case A verification
- [x] Cloud opacity in TUV-x (Phase 3) — from-host radiator with LWC-based OD
- [ ] **Slant-column cloud shadows in TUV-x** — Currently each column runs
  independent 1D radiative transfer (plane-parallel). At SZA ~59° the j_NO2
  cross-section shows purely vertical structure with no cloud shadow cast
  sideways. In reality, the solar beam passes through neighboring columns at
  an angle, so a cloud in one column should reduce actinic flux in adjacent
  columns on the shadow side. Options:
  - *Slant-column approximation:* trace the solar beam path at the geometric
    SZA through neighboring columns, accumulate their cloud OD, and apply as
    an above-cloud attenuation factor before the 1D TUV-x solve. Moderate
    complexity, physically motivated.
  - *3D radiative transfer:* full treatment (e.g., SHDOM). Accurate but very
    expensive and likely impractical for online chemistry.
  - *Effective cloud fraction:* parameterize shadow effects statistically
    using cloud fraction and SZA without explicit ray tracing. Simplest but
    least physical.
- [ ] Aerosol opacity in TUV-x — extend from-host radiator with aerosol OD
- [ ] Aerosol chemistry integration
- [ ] Parallel processing optimization for MICM solver
