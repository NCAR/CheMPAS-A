# MPAS-MICM tracer mapping notes

## Current behavior

Chemistry tracers are now allocated dynamically at runtime. They are no longer
hardcoded in `src/core_atmosphere/Registry.xml`.

Initialization flow:

1. `atm_extend_scalars_for_chemistry` (in `mpas_atm_core_interface.F`) queries
   MICM species names from the configured mechanism file.
2. The routine extends `scalars` and `scalars_tend` var_array metadata and adds
   `index_q*` dimensions in the state/tend pools.
3. `chemistry_init` calls `musica_init`, then `resolve_mpas_indices` to map
   MICM species to MPAS tracer indices.

## Naming convention

- MICM species names are used verbatim from the mechanism.
- MPAS tracer names are generated as `q` + MICM species name.
  - Example: `AB -> qAB`, `A -> qA`, `B -> qB`

## Error behavior

- If dynamic tracer discovery fails, initialization returns an error.
- If `resolve_mpas_indices` cannot find an `index_q*` dimension for any species,
  chemistry initialization fails.
- Runtime chemistry tracers are currently incompatible with
  `config_apply_lbcs=true`; this is guarded with a critical log message.

## Historical note

An earlier implementation used a registry-vs-MICM consistency checker
(`check_registry_tracer_consistency`). That checker was removed when runtime
tracer allocation became the source of truth.
