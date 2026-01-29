# MPAS–MICM tracer consistency notes

## Overview
Recent changes add a lightweight consistency check between MICM species (from the MICM configuration) and MPAS tracer names (from the `scalars` var_array metadata). The checker enforces a naming convention where an MPAS tracer corresponding to MICM species `X` is named `qX` (e.g., MICM `AB` → MPAS `qAB`). It logs every successful match, reports extra/missing entries, and aborts initialization if any MICM species cannot be found.

## Key code changes

### mpas_atm_chemistry
- `chemistry_init` now accepts the `structs` pool so it can fetch tracer names:

```fortran
use mpas_derived_types, only : field3dReal
use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_field
...
type (field3dReal), pointer     :: scalarsField
character(len=StrKIND), dimension(:), pointer :: tracer_names
...
call mpas_pool_get_subpool(structs, 'state', statePool)
if (associated(statePool)) then
    call mpas_pool_get_field(statePool, 'scalars', scalarsField, 1)
    if (associated(scalarsField)) tracer_names => scalarsField%constituentNames
end if
...
call musica_init(filepath_ptr, nVertLevels, error_code, error_message, tracer_names)
```

### mpas_musica
- `musica_init` accepts optional `registry_tracer_names` and invokes a new checker right after MICM state creation:

```fortran
subroutine musica_init(..., registry_tracer_names)
    ...
    call check_registry_tracer_consistency(state, registry_tracer_names, &
            check_error_code, check_error_message)
    if (check_error_code /= 0) then
        error_code = check_error_code
        error_message = check_error_message
        return
    end if
```

- The checker implements the `q` prefix rule, logs matches, and reports discrepancies:

```fortran
subroutine check_registry_tracer_consistency(state, registry_tracer_names, error_code, error_message)
    ...
    expected_tracer = 'q' // species_name
    if (expected_tracer == trim(registry_tracer_names(i_tracer))) then
        call mpas_log_write('[MUSICA] Matched MICM species ' // trim(species_name) // &
                            ' to MPAS tracer ' // trim(expected_tracer))
        registry_used(i_tracer) = .true.
    end if
    ...
    if (len_trim(missing_list) > 0) then
        error_code = 1
        error_message = '[MUSICA] MICM species missing in MPAS registry tracers:' // trim(adjustl(missing_list))
    end if
    if (len_trim(extra_list) > 0) then
        call mpas_log_write('[MUSICA] MPAS registry tracers not present in MICM config:' // trim(adjustl(extra_list)))
    end if
end subroutine
```

## Behavior
- **Match logging:** Each MICM species matched to an MPAS tracer produces a log line.
- **Missing MICM species:** Initialization aborts with an error; the missing list is logged.
- **Extra MPAS tracers:** Logged as informational; does not abort.
- **No tracer names available:** Checker is skipped (a log message notes this).

## Naming convention
- MICM species names are taken verbatim from the MICM configuration.
- MPAS expects the corresponding tracer names to be prefixed with `q` (mass mixing ratios), e.g., `A` → `qA`, `B` → `qB`, `AB` → `qAB`.

## Notes
- The AB/A/B initialization remains in place for testing.
- `mpas_atm_core.F` caller was updated to pass `blocklist%structs` into `chemistry_init`.
- No changes were made to `mpas_atm_core.F90` (it is a symlink). The changes are compiled via the `.F` source.
