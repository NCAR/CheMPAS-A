# MPAS-MUSICA Development TODO

## Chemistry Integration

- [ ] **Add namelist option for tracer seeding control**
  - Add `config_micm_seed_tracers` to `&musica` namelist (default `.true.`)
  - When `.false.`, skip MICM->MPAS tracer seeding to preserve custom initial conditions
  - Current workaround: seeding is skipped if tracers have spatial gradients

- [ ] **Output reference state to NetCDF**
  - Add `qAB_ref`, `qA_ref`, `qB_ref` variables to `Registry.xml`
  - Modify `mpas_musica.F` to copy reference state concentrations to MPAS pools each timestep
  - Add reference variables to output stream in `streams.atmosphere`
  - Update `plot_chemistry.py` to compare coupled vs reference states
  - This enables visualization of advection effects (coupled - reference = advection contribution)

## Future Enhancements

- [ ] **Use uxarray for proper unstructured mesh visualization** (low priority)
  - Current approach uses Delaunay triangulation of cell centers
  - uxarray respects actual MPAS Voronoi mesh topology
  - Would eliminate edge artifacts and improve accuracy
  - https://uxarray.readthedocs.io/

- [ ] Full atmospheric chemistry mechanisms (tropospheric ozone, etc.)
- [ ] TUV-x photolysis rate calculations
- [ ] Aerosol chemistry integration
- [ ] Species-specific molar mass lookup from MICM config
- [ ] Parallel processing optimization for MICM solver
