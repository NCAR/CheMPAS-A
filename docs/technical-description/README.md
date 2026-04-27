# Technical Description — port notes

This subtree is a lightly edited port of the draft NCAR Technical Note
*A Description of the Model for Prediction Across Scales, Atmosphere,
Version 8* by Skamarock, Duda, Klemp, and Fowler (29 May 2025).

## Port status

- **Text and equations:** ported from the upstream draft, with light
  editorial fixes and source paths normalized to this repository layout.
- **Figures:** deferred to a follow-on session. Placeholders of the form
  `**[Figure N.M: caption. To be added next session.]**` mark the
  intended figure location in each chapter.
- **Bibliography:** inline citations (e.g. *Ringler et al. 2008*) are
  kept as written; the entries are collected in
  [`0F-bibliography.md`](0F-bibliography.md).

## Figures plan (Option C, next session)

Three approaches are in play for the ~20 figures:

1. **PDF/raster extraction** for the mesh schematics, vertical-grid
   diagram, terrain-following coordinate visualization, tangent-plane
   schematics, and the variable-resolution globe.
2. **Mermaid diagrams** for the flowcharts (Figure 1.2 top-level flow
   chart, Figures 3.1 and 3.3 RK3 timestep pseudocode, Figure 7.3 LBC
   pseudocode). These render natively in Sphinx with
   ``sphinxcontrib-mermaid``.
3. **matplotlib regeneration** for Figure 3.2 (Runge–Kutta response
   curves).

## Build

From ``docs/``:

```
make html
```

Dollar-math and AMS-math are already enabled in ``conf.py``; no
additional Sphinx configuration is required for this subtree.
