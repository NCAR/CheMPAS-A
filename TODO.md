# MPAS-MUSICA Development TODO

## Chemistry Integration

- [ ] **Add namelist option for tracer seeding control**
  - Add `config_micm_seed_tracers` to `&musica` namelist (default `.true.`)
  - When `.false.`, skip MICM->MPAS tracer seeding to preserve custom initial conditions
  - Current workaround: seeding is skipped if tracers have spatial gradients

- [ ] **Output reference state to NetCDF**
  - Add runtime chemistry reference tracers (e.g., `q*_ref`) without hardcoding species in `Registry.xml`
  - Modify `mpas_musica.F` to copy reference-state concentrations to MPAS pools each timestep
  - Add reference variables to output stream in `streams.atmosphere`
  - Update `plot_chemistry.py` to compare coupled vs reference states
  - This enables visualization of advection effects (coupled - reference = advection contribution)

## Future Enhancements

- [ ] **Use uxarray for proper unstructured mesh visualization** (low priority)
  - Current approach uses Delaunay triangulation of cell centers
  - uxarray respects actual MPAS Voronoi mesh topology
  - Would eliminate edge artifacts and improve accuracy
  - https://uxarray.readthedocs.io/

- [ ] Solar geometry / day-night physics for j_NO2 (Phase 1)
- [ ] TUV-x photolysis rate calculations (Phase 2-3)
- [x] Domain-integrated Ox conservation test for Case A verification
- [ ] Aerosol chemistry integration
- [ ] Parallel processing optimization for MICM solver
