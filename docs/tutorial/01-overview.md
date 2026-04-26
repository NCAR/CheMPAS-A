# Chapter 1: Overview

```{admonition} Work in progress
:class: warning

This chapter is being actively written. Commands and cross-references
are provisional and may change.
```

The CheMPAS Tutorial walks through CheMPAS-A's idealized chemistry test
cases, the same cases driven by the numerical regression suite
(`scripts/regression.py`). Where the [User's Guide](../users-guide/index.rst)
is reference-style — a verbatim port of the upstream MPAS-Atmosphere
documentation — this tutorial is narrative: run the case, look at the
output, understand what the chemistry is doing.

## What this tutorial assumes

- `atmosphere_model` is built. See `BUILD.md` in the repository root
  and [Chapter 3 of the User's Guide](../users-guide/03-building.md).
- `~/Data/CheMPAS/<case>/` is set up with the namelist, streams, graph
  partition, and (where needed) initial-condition files. See `RUN.md`
  in the repository root and `test_cases/README.md`.
- The conda environment `mpas` is available for plotting:
  `conda activate mpas`.

## Chapters

- [Chapter 2: Supercell with ABBA and LNOx](02-supercell.md) — idealized
  deep convection, run with two MUSICA/MICM mechanisms, with a
  side-by-side comparison.
- [Chapter 3: Chapman + NOx Photostationary State](03-chapman-nox.md) —
  small-domain Chapman cycle plus NOx, where the analytical PSS solution
  is a clean numerical sanity check. *(Placeholder; content coming.)*

## Verifying numerically

The tutorial focuses on what each run looks like. For numerical match
against tracked reference values, the source of truth is the regression
suite:

```bash
python scripts/regression.py run --case supercell
```

See `docs/superpowers/specs/2026-04-19-regression-suite-design.md` for
the regression suite's design.
