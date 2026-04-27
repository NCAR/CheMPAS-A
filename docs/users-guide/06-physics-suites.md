# Chapter 6: Physics Suites

Beginning with version 4.0, MPAS-Atmosphere introduces a new way of selecting the physics schemes to be used in a simulation. Rather than selecting individual parameterization schemes for different processes (e.g., convection, microphysics, etc.), the preferred method is for the user to select a *suite* of parameterization schemes that have been tested together. The selection of a physics suite is made via the namelist option `config_physics_suite` in the `&physics` namelist record. Each of the available suites are described in the sections that follow.

Although the preferred method for selecting the schemes in a simulation is via the choice of a suite, the need to enable or disable individual schemes, or to substitute alternative schemes for the suite default, is recognized. Accordingly, it is possible to override the choice of any individual parameterization scheme through the namelist options described in [Appendix B](0B-model-namelist.md). This is useful, e.g., to disable all parameterizations except for microphysics when running some idealized simulations. The details of selecting individual physics parameterizations are explained in Section 6.4.

## 6.1 Suite: mesoscale_reference

The default physics suite in CheMPAS-A is the `mesoscale_reference` suite, which contains the schemes listed in Table 6.1. This suite has been tested for mesoscale resolutions (> 10 km cell spacing), and is not appropriate for convective-scale simulations because the Tiedtke scheme will remove convective instability before resolved-scale motions (convective cells) can respond to it.

*Table 6.1: The set of parameterization schemes used by the `mesoscale_reference` physics suite.*

| Parameterization | Scheme |
|------------------|--------|
| Convection | New Tiedtke |
| Microphysics | WSM6 |
| Land surface | Noah |
| Boundary layer | YSU |
| Surface layer | Revised Monin-Obukhov |
| Radiation, LW | RRTMG |
| Radiation, SW | RRTMG |
| Cloud fraction for radiation | Xu-Randall |
| Gravity wave drag by orography | YSU |

## 6.2 Suite: convection_permitting

The `convection_permitting` physics suite is appropriate at spatial resolutions allowing for both explicitly resolved hydrostatic and nonhydrostatic motions. It has been tested for mesh spacings from several hundred kilometers down to 3 km in MPAS. The Grell-Freitas convection scheme transitions from a conventional parameterization of deep convection at hydrostatic scales (cell spacings of several tens of kilometers) to a parameterization of precipitating shallow convection at cell spacings less than 10 km. This is the recommended suite for any MPAS applications where convection-permitting meshes (dx < 10 km) are employed, including variable-resolution meshes spanning hydrostatic to nonhydrostatic resolutions.

*Table 6.2: The set of parameterization schemes used by the `convection_permitting` physics suite.*

| Parameterization | Scheme |
|------------------|--------|
| Convection | Grell-Freitas |
| Microphysics | Thompson (non-aerosol aware) |
| Land surface | Noah |
| Boundary layer | MYNN |
| Surface layer | MYNN |
| Radiation, LW | RRTMG |
| Radiation, SW | RRTMG |
| Cloud fraction for radiation | Xu-Randall |
| Gravity wave drag by orography | YSU |

## 6.3 Suite: none

The only other recognized physics suite in CheMPAS-A is the `none` suite, which sets all physics parameterizations to `off`. This suite is primarily intended for use with idealized simulations. For example, the idealized supercell test case makes use of the `none` suite, but with the microphysics scheme explicitly overridden:

```
config_physics_suite = 'none'
config_microp_scheme = 'mp_kessler'
```

## 6.4 Selecting Individual Physics Parameterizations

Selecting or disabling an individual physics parameterization may be accomplished by setting the appropriate namelist variable to one of its possible options; possible options for individual parameterizations, along with details of those options, are given in Table 6.3. Note that all parameterization options may be set to `'off'` to disable the parameterization of the associated process.

*Table 6.3: Possible options for individual physics parameterizations. Namelist variables should be added to the `&physics` namelist record.*

| Parameterization | Namelist variable | Possible options | Details |
|------------------|-------------------|-----------------|---------|
| Convection | `config_convection_scheme` | `cu_tiedtke` | Tiedtke (WRF 3.8.1) |
| | | `cu_ntiedtke` | New Tiedtke (WRF 4.5) |
| | | `cu_grell_freitas` | Modified version of scale-aware Grell-Freitas (WRF 3.6.1) |
| | | `cu_kain_fritsch` | Kain-Fritsch (WRF 3.2.1) |
| | | `off` | Disable convection |
| Microphysics | `config_microp_scheme` | `mp_wsm6` | WSM 6-class (WRF 4.5) |
| | | `mp_thompson` | Thompson non-aerosol aware (WRF 3.8.1) |
| | | `mp_thompson_aerosols` | Thompson aerosol-aware (WRF 4.1.4) |
| | | `mp_kessler` | Kessler |
| | | `off` | Disable cloud microphysics |
| Land surface | `config_lsm_scheme` | `sf_noah` | Noah (WRF 4.5) |
| | | `sf_noahmp` | Noah-MP 5.0.1 |
| | | `off` | Disable the land-surface scheme |
| Boundary layer | `config_pbl_scheme` | `bl_ysu` | YSU (WRF 4.5) |
| | | `bl_mynn` | MYNN (WRF 3.6.1) |
| | | `off` | Disable boundary-layer parameterization |
| Surface layer | `config_sfclayer_scheme` | `sf_monin_obukhov` | Monin-Obukhov (WRF 4.5) |
| | | `sf_monin_obukhov_rev` | Revised Monin-Obukhov (WRF 4.5) |
| | | `sf_mynn` | MYNN (WRF 3.6.1) |
| | | `off` | Disable the surface-layer scheme |
| Radiation, LW | `config_radt_lw_scheme` | `rrtmg_lw` | RRTMG (WRF 3.8.1) |
| | | `cam_lw` | CAM (WRF 3.3.1) |
| | | `off` | Disable long-wave radiation |
| Radiation, SW | `config_radt_sw_scheme` | `rrtmg_sw` | RRTMG (WRF 3.8.1) |
| | | `cam_sw` | CAM (WRF 3.3.1) |
| | | `off` | Disable short-wave radiation |
| Cloud fraction for radiation | `config_radt_cld_scheme` | `cld_fraction` | Xu and Randall (1996) |
| | | `cld_incidence` | 0/1 cloud fraction depending on q_c + q_i |
| | | `cld_fraction_thompson` | Thompson cloud fraction scheme |
| | | `off` | Disable cloud-fraction calculation; if radiation remains enabled, the model switches this to `cld_incidence` at runtime |
| Gravity wave drag by orography | `config_gwdo_scheme` | `bl_ysu_gwdo` | YSU (WRF 4.5) |
| | | `bl_ugwp_gwdo` | NOAA/GSL orographic gravity wave drag (see also `config_ngw_scheme`) |
| | | `off` | Disable orographic gravity wave drag |
