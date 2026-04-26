# Chapter 3: Chapman + NOx Photostationary State

```{admonition} Work in progress
:class: warning

This chapter is being actively written. Commands and expected output
are provisional; figures are placeholders.
```

The Chapman + NOx photostationary-state (PSS) tutorial walks through a
small-domain integration of the Chapman ozone cycle plus NOx, with
TUV-x photolysis driven by an extended atmosphere column that reaches
above the model lid. The analytical Leighton expression for [NO]/[NO₂]
under steady-state photolysis is a clean numerical sanity check on the
coupled MICM + TUV-x configuration.

## 3.1 What you'll learn

```{admonition} Work in progress
:class: warning

Section content coming.
```

By the end of this chapter you will:

- Run the Chapman + NOx idealized stratospheric-chemistry case in
  CheMPAS-A on the supercell mesh.
- Generate the TUV-x upper-atmosphere extension CSV and understand
  why TUV-x needs photons from above the model lid.
- Verify the chemistry against the analytical Leighton photostationary
  state and the regression suite.

## 3.2 The Chapman + NOx case

```{admonition} Work in progress
:class: warning

Section content coming.
```

The Chapman cycle is the canonical four-reaction pure-oxygen
photochemistry that maintains a stratospheric ozone column:

$$
\begin{aligned}
\mathrm{O_2} + h\nu &\rightarrow 2\,\mathrm{O} \\
\mathrm{O} + \mathrm{O_2} + \mathrm{M} &\rightarrow \mathrm{O_3} + \mathrm{M} \\
\mathrm{O_3} + h\nu &\rightarrow \mathrm{O} + \mathrm{O_2} \\
\mathrm{O} + \mathrm{O_3} &\rightarrow 2\,\mathrm{O_2}
\end{aligned}
$$

Adding NOx introduces the catalytic
NO–NO₂–O₃ cycle (NO + O₃ → NO₂ + O₂; NO₂ + hν → NO + O), which
modulates ozone titration by tying its evolution to NOx photolysis.
On any timescale longer than a few seconds, [NO] / [NO₂] in sunlight
relaxes to the **Leighton photostationary state**, the analytical
target of section 3.7.

The Chapman cycle is global-stratospheric physics, but
`scripts/init_chapman.py` seeds a 1-D AFGL mid-latitude-summer ozone
profile uniformly across the supercell mesh, and the chemistry has no
feedback on dynamics. This chapter therefore uses the small
(~85 km × 85 km × 50 km top) supercell grid as a column-like sandbox
— what matters is the vertical structure of the photolysis driver and
the chemistry's ability to settle into the PSS, both of which TUV-x
sees through the column extension introduced in section 3.3.
Horizontal dynamics are present but largely irrelevant to the PSS
demonstration.

**[Figure 3.1: AFGL mid-latitude-summer O₃ profile interpolated to the
supercell vertical grid (the initial state qO3 produces). To be added.]**

## 3.3 The TUV-x column extension

```{admonition} Work in progress
:class: warning

Section content coming.
```

## 3.4 Generating and verifying the extension CSV

```{admonition} Work in progress
:class: warning

Section content coming.
```

**[Figure 3.2: Stitched T, n_air, and n_O₃ vertical profiles from
mpas_tuvx.F. MPAS region (below 50 km) and extension-CSV region
(above 50 km) overplotted. To be added.]**

## 3.5 Initializing the Chapman tracers

```{admonition} Work in progress
:class: warning

Section content coming.
```

## 3.6 Run with the Chapman + NOx mechanism

```{admonition} Work in progress
:class: warning

Section content coming.
```

**[Figure 3.3: Vertical profiles of qO3, qNO, qNO2, and the NO/NO₂
ratio at t = 2 h, mid-domain column, Chapman + NOx mechanism. To be
added.]**

## 3.7 The photostationary-state diagnostic

```{admonition} Work in progress
:class: warning

Section content coming.
```

## 3.8 Verifying numerically

```{admonition} Work in progress
:class: warning

Section content coming.
```

**[Figure 3.4: Simulated vs. analytical Leighton [NO]/[NO₂] ratio vs.
height at the final timestep. To be added.]**

## 3.9 Next steps

```{admonition} Work in progress
:class: warning

Section content coming.
```
