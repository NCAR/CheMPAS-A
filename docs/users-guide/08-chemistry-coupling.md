# Chapter 8: Chemistry Coupling

```{admonition} Under construction
:class: warning

This chapter describes the CheMPAS-A ↔ MUSICA/MICM chemistry coupling —
state transfer, unit conventions, external rate-parameter wiring,
photolysis input from TUV-x, and the chemistry sub-stepping controls.
A full draft is in preparation. Source notes for this chapter live in
``docs/chempas/musica/MUSICA_INTEGRATION.md`` and
``docs/chempas/guides/TUVX_INTEGRATION.md``.
```

## 8.1 Overview

CheMPAS-A advances chemistry via the MUSICA/MICM solver framework. On each
chemistry step, the model gathers temperature, pressure, and the chemistry
species from the MPAS `scalars` pool, hands MICM a per-cell state in
mol m⁻³, and writes the integrated state back into the scalar tracers.
Photolysis rate parameters are supplied either by TUV-x (when configured)
or by a cos(SZA) fallback for development cases. Lightning NOx and the
chemistry sub-stepping / tolerance controls are exposed through the
`&musica` namelist record.

For details of the MICM solver itself — supported mechanism families,
rate-constant forms, stiffness handling, and configuration-file syntax —
refer to the MUSICA documentation at
<https://musica.readthedocs.io/>.

## 8.2 State Transfer

*To be drafted.* Will cover the gather/scatter between MPAS scalars (kg
kg⁻¹ for moisture, mol mol⁻¹ for chemistry) and MICM's mol m⁻³ working
state, including the temperature/pressure inputs and the unit conversion
in `mpas_musica.F`.

## 8.3 External Rate Parameters and Photolysis

*To be drafted.* Will cover `musica_set_photolysis`, the TUV-x clear-sky
and cloud-radiator coupling, and the cos(SZA) fallback path.

## 8.4 Chemistry Sub-stepping and Tolerances

*To be drafted.* Will cover `config_chem_substeps` and
`config_micm_relative_tolerance`, including when sub-stepping is required
to avoid non-physical implicit roots in stiff mechanisms (e.g. Chapman).
