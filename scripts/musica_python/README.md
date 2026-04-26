# Standalone MUSICA-Python Examples

Pedagogical scripts that exercise the same MICM mechanism configs the
CheMPAS-A tutorial chapters use, *without* MPAS in the loop. Each
script mirrors a section of the tutorial:

| Script                  | Tutorial section            | Mechanism             |
| ----------------------- | --------------------------- | --------------------- |
| `abba_box.py`           | §2.10 (Chapter 2 supercell) | `abba.yaml`           |
| `lnox_box.py`           | §2.11 (Chapter 2 supercell) | `lnox_o3.yaml`        |
| `chapman_nox_column.py` | §3.10 (Chapter 3 Chapman)   | `chapman_nox.yaml`    |

Each script writes `<case>.nc` (xarray Dataset) and `<case>.png`
into this directory.

## Dependencies

`numpy`, `xarray`, `matplotlib`, `netCDF4` are already in the `mpas`
conda environment. Two extra installs:

```bash
~/miniconda3/envs/mpas/bin/pip install musica ussa1976 pvlib
```

(Equivalent: `pip install 'musica[tutorial]'` per MUSICA's README.)

## Running

```bash
~/miniconda3/envs/mpas/bin/python scripts/musica_python/abba_box.py
~/miniconda3/envs/mpas/bin/python scripts/musica_python/lnox_box.py
~/miniconda3/envs/mpas/bin/python scripts/musica_python/chapman_nox_column.py
```

## Plot style

All three scripts import `scripts/style.py` (via a `sys.path` shim)
and call `style.apply_ncar_style()`. Axis labels and titles use
`style.species_label()` and `style.format_title()` for consistent
NCAR formatting.
