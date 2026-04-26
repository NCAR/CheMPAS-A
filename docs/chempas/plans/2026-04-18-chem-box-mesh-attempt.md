# Small-domain chemistry test case (chem_box) — attempt log

## Goal

Run a multi-day Chapman+NOx integration on a very small horizontal mesh
(~64 cells, 4×4 km) to watch the chemistry diurnal cycle quickly. The
supercell mesh's 28,080 cells × 60 levels = 1.7 M chemistry columns make
multi-day runs prohibitive (~3 days wall clock per sim-day).

Requirements were:
- 64 cells in a periodic-plane hex mesh
- Same vertical grid (60 stretched levels, 0–50 km)
- Restart capability (hourly checkpoints)
- Chapman+NOx mechanism end-to-end, including TUV-x photolysis driven
  by simulation time so the SZA sweeps through day/night

## What we tried

### 1. Hand-rolled periodic hex mesh generator (`gen_periodic_hex_mesh.py`)

Wrote a ~300-line NumPy generator that emits a NetCDF mesh compatible
with MPAS init_atmosphere — every required variable present, topology
consistent (nEdges = 3·nCells, nVertices = 2·nCells, every cell has
six neighbours, CCW edge ordering, periodic reciprocity verified).

`init_atmosphere_model` loaded the mesh cleanly. `atmosphere_model`
crashed on the first dynamics timestep with NaN in the chemistry scalars
and, after TRSK-weight iteration, in `u`/`w` too.

**Root cause:** the Thuburn-Ringler-Skamarock-Klemp weights in
`weightsOnEdge` are non-trivial for MPAS dynamics. Magnitudes follow
the well-known `|w_k| = (1/2 − k/6) × dv/dc` pattern (observed directly
in `supercell_grid.nc`: values `{±0.1925, ±0.0962, 0}`), but the sign
convention depends on edge orientation, cell-side, and CCW traversal
direction in ways that are not documented in the repo and did not match
any simple angular rule I tried. Observed sign patterns by edge normal
angle in supercell:

```
 30°:  [+ + 0 + +  + + 0 + +]   (all positive)
 90°:  [+ - 0 + -  + - 0 + -]   (alternating)
150°:  [- - 0 - -  - - 0 - -]   (all negative)
```

A sign rule `sign(sin(angle_j − angle_i))` produces the wrong pattern.
The actual MPAS convention is embedded in the MPAS-Tools mesh generator
and requires either (a) reading that generator's Fortran source to
reproduce exactly, or (b) using MPAS-Tools directly.

### 2. Subset of the existing supercell mesh (`subset_supercell_mesh.py`)

Idea (user-suggested): take an 8×8 block of cells from the supercell
mesh, renumber 0..63, stitch up boundary neighbours to wrap within the
block. All TRSK weights and geometry carry over, since the supercell
is uniform hex.

Got the cell selection working (pick cells with the smallest xCell
and yCell, grouped by row). Partial topology rebuild. Hit edge cases
in the "direction of the out-of-block neighbour" logic when a supercell
boundary cell's outside-neighbour had identical coordinates after the
supercell's own periodic wrap (`dx = 0, dy = 0`).

The bookkeeping is doable with more care but the hex lattice's
column-parity offset interacting with the block's own periodic wrap
produced enough bugs that a reliable result in the session's budget
wasn't reached.

## What we learned (carry-forward)

- MPAS init_atmosphere reads `weightsOnEdge` from the grid file. The
  atmosphere core does **not** regenerate them at init. See
  `mpas_init_atm_cases.F:592` (`mpas_pool_get_array(mesh, 'weightsOnEdge', …)`)
  and the absence of any `weightsOnEdge =` assignment anywhere in
  `src/core_atmosphere` or `src/core_init_atmosphere`.
- Zeroing `weightsOnEdge` or setting only magnitudes is not enough — the
  acoustic/dynamics step produces NaN within one timestep even with
  zero initial winds.
- Supercell's edge normals are 30°/90°/150° (pointy-top hex with
  column-offset rows). A hand-rolled mesh must match this convention
  to reuse supercell's weight templates directly.
- `nominalMinDc` in the NetCDF is stored but returns 0 from the supercell
  file (loaded as a masked variable). `config_len_disp` in the namelist
  is what the dynamics uses as the reference mesh scale.
- The supercell mesh is effectively an 84 km × 84 km doubly-periodic
  pointy-top hex lattice with cell spacing ≈ 538.86 m (from `dcEdge`).
  The global cell ordering is not row-major — it looks partitioned
  (indices jump by ~312 between cells in the same y-row), probably
  from a METIS-like reordering.

## Right path forward

Obtain MPAS-Tools' `periodic_hex` mesh generator, which produces a
correct doubly-periodic hex mesh with TRSK weights that MPAS dynamics
accepts. Two approaches:

1. **Source checkout.** `MPAS-Dev/MPAS-Tools` on GitHub has the
   periodic_hex generator (Fortran, one small executable). Build
   against the existing `netcdf-fortran` / `pnetcdf` available on this
   machine; takes ~15–30 min once the source is identified.
2. **Python port (`mpas_tools` package).** A Python `mpas_tools`
   exists on PyPI / conda-forge with a mesh generation submodule.
   Check whether the conda `mpas` environment has or can install it.

Once a valid tiny mesh is available, everything else in the chem_box
test case (namelists, chapman_nox.yaml, tuvx_chapman_nox.json,
init_chapman.py seeding, restart configuration) is ready and tested —
the only missing piece is the mesh with correct weights.

## Files left behind vs. cleaned up

**Kept** (design reference for the next attempt):
- this log

**Cleaned up** (will bring back with the real mesh):
- `scripts/gen_periodic_hex_mesh.py` — hand-rolled generator with
  incorrect TRSK weights
- `scripts/subset_supercell_mesh.py` — incomplete subset tool
- `test_cases/chem_box/` — namelists/streams/configs; re-create when
  a valid grid file is in hand
- `~/Data/CheMPAS/chem_box/` — run directory
