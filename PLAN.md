# Plan: Generalize Chemistry Coupling — Remove Hardcoded ABBA

## Status: Phase 1 Complete

Phase 1 was implemented and merged to `develop` in commit `2dee808`.

## Context

All chemistry tracer coupling in CheMPAS was hardcoded for the 3-species ABBA mechanism (AB, A, B). Species names appeared as string literals in 7+ places each, molar masses were Fortran `parameter` constants, and tracer indices were fetched individually by name. This made it impossible to switch chemistry mechanisms without editing Fortran source code.

Phase 1 generalized the coupling layer to loop over species dynamically using MICM's runtime API. Registry.xml tracer definitions remain static for now — runtime scalars allocation (Phase 2) will be a separate future effort.

## Files Modified (Phase 1)

| File | Changes |
|------|---------|
| `src/core_atmosphere/chemistry/musica/mpas_musica.F` | Dynamic `species_entry_t` table, generic coupling loops |
| `src/core_atmosphere/chemistry/mpas_atm_chemistry.F` | Removed per-species index args, renamed `chemistry_seed_chem` |
| `src/core_atmosphere/mpas_atm_core.F` | Updated call sites: `chemistry_seed_chem`, `chem_seed_after_stream` |
| `src/core_atmosphere/Registry.xml` | Wrapped qAB/qA/qB in `#ifdef MPAS_USE_MUSICA` |
| `~/Data/MPAS/supercell/abba.yaml` | Added `__molar mass` to each species |

## What Was Done (Phase 1)

### Step 1: Add `__molar mass` to `abba.yaml` — DONE

Added `__molar mass` property (kg/mol) to each species in the MICM config:
- A: 0.029, B: 0.029, AB: 0.058

### Step 2: Add species table to `mpas_musica.F` — DONE

Replaced hardcoded `parameter` block (M_AB, M_A, M_B) with module-level derived type and allocatable array:

```fortran
type :: species_entry_t
    character(len=64) :: micm_name    ! e.g. 'AB'
    character(len=64) :: mpas_name    ! e.g. 'qAB'
    integer           :: micm_index   ! from state%species_ordering%index(name)
    integer           :: mpas_index   ! from MPAS pool dimension
    real(kind=8)      :: molar_mass   ! kg/mol from MICM config
end type
```

### Step 3: Populate species table in `musica_init` — DONE

Generic loop over `state%species_ordering`:
- Builds `mpas_name` via `'q' // trim(micm_name)` convention
- Reads `molar_mass` via `micm%get_species_property_double` (scalar return)
- Uses `state%species_ordering%index(name)` for correct stride-based `micm_index`
- Seeds both `state` and `state_ref` with initial concentrations

### Step 4: Add `resolve_mpas_indices` subroutine — DONE

Looks up MPAS pool dimensions (`index_qAB`, etc.) and fills `chem_species(i)%mpas_index`. Called from `chemistry_init` after `musica_init`.

### Step 5–6: Rewrite `MICM_from/to_chemistry` — DONE

Generic species loops replace explicit per-species blocks. No more index arguments in signatures.

### Step 7: Rename `micm_to_mpas_abba` → `micm_to_mpas_chem` — DONE

Generic seeding loop over `chem_species(:)`. Gradient detection uses first species.

### Step 8: Generalize `log_column_comparison` — DONE

Dynamic loop prints each species name and coupled/reference concentrations.

### Step 9: Update `mpas_atm_chemistry.F` — DONE

Removed per-species index declarations and pool lookups. Simplified call signatures. Renamed public interface.

### Step 10: Update `mpas_atm_core.F` — DONE

All call sites updated: `chemistry_seed_abba` → `chemistry_seed_chem`, `abba_seed_after_stream` → `chem_seed_after_stream`.

### Step 11: Wrap Registry.xml tracers — DONE

Used `#ifdef MPAS_USE_MUSICA` / `#endif` around qAB/qA/qB entries in both `scalars` and `scalars_tend` var_arrays. Non-MUSICA builds have no chemistry tracers.

## Verification Results

1. **Build:** MUSICA=true build succeeds (MUSICA-Fortran 0.13.0)
2. **Run:** Supercell test passes with 8 MPI ranks, ~32s wall time, 0 errors/warnings
3. **Log check:** Dynamic species discovery confirmed:
   - `[MUSICA] Species A molar_mass=0.029 kg/mol`
   - `[MUSICA] Species AB molar_mass=0.058 kg/mol`
   - `[MUSICA] Species B molar_mass=0.029 kg/mol`
   - `[MUSICA] Resolved A -> MPAS index 5`
   - `[MUSICA] Resolved AB -> MPAS index 4`
   - `[MUSICA] Resolved B -> MPAS index 6`
4. **Chemistry correctness:** Mass conservation holds (qA + qB + qAB ≈ 1.0)
5. **Net change:** -89 lines (244 insertions, 333 deletions) across 5 files

## Lessons Learned

- `state%species_ordering%index(name)` returns the stride-based index needed for concentration array access; the iteration counter from `species_ordering%name(i)` is NOT the same thing
- `micm%get_species_property_double` returns a scalar `real(8)`, not an array
- MPAS logging `realArgs` must use `kind=RKIND` (not `kind=4`) when `-fdefault-real-8` is active
- Running with wrong MPI rank count (e.g., 1 rank on 8-partition mesh) causes segfaults in RK3 dynamics — always test with correct partition

## Phase 2: Future Work

- **Runtime scalars allocation** — Remove tracers from Registry.xml entirely; allocate at runtime based on MICM config
- **`__do advect` filtering** — Support non-advected species
- **Fallback molar mass table** — For community chemistry configs that lack `__molar mass`
- **Generic visualization** — Update `plot_chemistry.py` to discover species dynamically instead of hardcoding qA/qB/qAB
- **Non-MUSICA build verification** — Confirm clean compilation without MUSICA flag
