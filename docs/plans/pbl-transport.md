# Implementation plan: PBL transport of MICM tracers

Status: planned; transport implementation and qualification are pending.

- Tracking issue: [CheMPAS-A #40](https://github.com/NCAR/CheMPAS-A/issues/40).
- Branch: `chempas-a-pbl`, from `chempas-a` at
  `89008ed50ce8f9ffd10876deed7c9a48ec573a38`.
- Companion: [convective transport plan](https://github.com/NCAR/CheMPAS-A/blob/chempas-a-convect/docs/plans/convective-transport.md).
- Coordination: [MUSICA #1073](https://github.com/NCAR/musica/issues/1073).
- Source/issue audit: 2026-10-18. These are development plans, not supported
  capabilities. Standalone MICM contribution #39 is a separate effort.

Implement vertical gas transport for every prognostic MICM tracer discovered at
runtime, retaining the existing native MPAS moisture transport. Emissions,
lightning sources, deposition, scavenging, and aerosol transport are excluded
from this qualification. Keep those sources disabled in the transport fixtures.

## Findings at the branch base

| Location | Finding and consequence |
| --- | --- |
| `chemistry/musica/mpas_musica.F`: `species_entry_t`, `resolve_mpas_indices`, `host_binding` | `chem_species(:)%mpas_index` is the authoritative mapping. H2O is bound to native `qv`; prescribed O2/N2 can also use `qv` as an index anchor. An index alone does not identify an independently transported gas. |
| `physics/physics_wrf/module_bl_ysu.F`: `ysu` | The wrapper declares `nmix = 0`; its existing pack/unpack loops cannot transport MICM species. |
| MMM-physics `bl_ysu.F90`: `bl_ysu_run` | A generic `qmix`/`qmixtnp` loop already exists. It uses the moisture diffusion matrix, with its own `xkzq` and entrainment treatment. It does not copy the heat/water nonlocal source terms to passive gases. The solve uses `dt2 = 2*dt`, with tendencies divided by `dt2`. |
| `physics/mpas_atmphys_driver_pbl.F`: `pbl_to_MPAS` | The exported `kzq` is assigned `kzh_p`. It is not sufficient evidence of the scalar operator actually used inside YSU. |
| `physics/physics_wrf/module_bl_mynn.F` and MMM-physics `bl_mynn.F90` | Chemistry arguments are behind `WRF_CHEM`; the MPAS call supplies none. The dormant wrapper has inconsistent chemistry identifiers. The external chemistry writeback is also guarded by `WRF`, with a positive concentration floor. Enabling a macro is not a complete adapter. |
| MMM-physics `bl_mynn_subroutines.F90`: `mynn_mix_chem`, `DMP_mf` | Existing code combines diffusion and plume tracer fluxes. Its chemistry solve includes deposition arguments and a prescribed top value. Closed-column gas transport requires explicit zero-flux boundaries and no artificial tracer floor. |
| `physics/mpas_atmphys_interface.F`: `MPAS_to_physics` | `rho_p = zz*rho_zz*(1+qv)` includes water vapor; chemistry obtains dry density as `zz*rho_zz`. Moist pressure/density weights cannot be assumed to conserve dry-air tracer mass. |
| `physics/mpas_atmphys_manager.F` and `mpas_atmphys_driver.F` | `config_pbl_interval` creates an alarm and sets `dt_pbl`, but the driver is called every model step when PBL and surface layer are enabled. No PBL alarm check gates that call at this base. |
| `mpas_atm_core.F`: `atm_do_timestep` | Order is physics, chemistry on scalar time level 1, then dynamics. Tracer tendencies computed inside the earlier physics call can already be stale when dynamics uses them. |
| `physics/mpas_atmphys_todynamics.F`: `physics_get_tend_work` | Only fixed native scalar indices receive PBL tendencies, multiplied by `rho_zz` from time level 2. |
| `dynamics/mpas_atm_time_integration.F` | Ordinary and monotonic scalar advances handle source tendencies differently; the monotonic path applies them to old state and clears them. Verify both the mass time level and single application before enabling a transport mode. |

Paths in this table are relative to `src/core_atmosphere/`. The external source
audited was [NCAR/MMM-physics at a4baf7f3](https://github.com/NCAR/MMM-physics/tree/a4baf7f3243d1db0dbc5f63473f895bdbdc05c30),
selected by `Externals.cfg` tag `20250616-MPASv8.3`.

## Shared interface with the convection branch

Implement this foundation first in a separate commit on `chempas-a-pbl`; reuse
that commit on `chempas-a-convect`. Do not create two independent species maps,
budget implementations, or tendency application loops. The following names are
proposed interfaces, not existing routines or namelist options.

1. Build a compact `transport_to_scalar(:)` map after MICM resolves its indices.
   Include every non-host-bound, writable prognostic MICM species, regardless of
   mechanism name, ordering, or scalar contiguity. Validate range and uniqueness;
   fail on an unexpected collision with a native field. Log each species,
   scalar index, and exclusion reason. Do not filter species by initial value,
   chemical lifetime, or a fixed list. An empty map is a valid no-op.
2. Store the map in block-owned host data. Physics accepts indices and ordinary
   Fortran arrays, without importing MICM solver objects or `mpas_musica`.
   Populate the host data from chemistry initialization and supply an empty map
   without MUSICA. Preserve the current one-block-per-rank chemistry restriction.
   This avoids a physics/chemistry build cycle: these libraries currently build
   as independent prerequisites.
3. Add a physics-side helper, provisionally `mpas_atmphys_tracer_transport.F`,
   with compact buffers `q(nTransport,nVertLevels,nOwnedCells)` and separate PBL
   and convection tendency buffers. Register only needed persistent cache and
   diagnostic fields through existing pools/runtime registration. Allocate outside
   threaded loops; keep tile workspaces private and clear every inactive column.
   Transport all host vertical levels, including any outside a restricted MICM
   solve domain. Bound workspace memory with tiles/species chunks; do not rerun
   the meteorological scheme for each species.
4. Units are `q = kg species / kg dry air`, `R = dq/dt` in `s-1`, dry layer mass
   per area `m_d = rho_dry*dz`, and interface tracer flux in
   `kg species m-2 s-1`. For an upward-positive flux,
   `m_d(k)*R(k) = F(k-1/2) - F(k+1/2)`. The domain species budget is
   `sum_owned(areaCell*m_d*q)`. Establish the exact equivalence to the MPAS
   `rho_zz`, `zz`, and vertical metrics at the application time level.
5. Cache meteorological transport coefficients/plume properties during physics.
   Construct gas tendencies later, at the once-per-model-step tendency assembly
   after chemistry, using current authoritative `scalars` and actual model `dt`.
   Carry `dt` explicitly to the helper. Cache coefficients, not tracer-dependent
   increments that would be replayed after chemistry or another transport update.
6. Use a scratch column for the deterministic sequence PBL then convection:
   `q1 = PBL(q0,dt)`, `q2 = CU(q1,dt)`,
   `R_pbl = (q1-q0)/dt`, `R_cu = (q2-q1)/dt`.
   Add the sum once at each mapped scalar index in `physics_get_tend`, after
   tendency clearing, using the proven MPAS mass conversion. Do not also write
   these trial updates to the prognostic state. The normal scalar advance owns
   state mutation; native moisture tendencies retain their existing path.
7. Prove the resulting full-step increment through ordinary/monotonic transport,
   RK stages, and dynamics substeps. Include fresh post-chemistry state, chemistry
   skipped intervals, and the existing scalar halo exchanges. If the advection
   source treatment cannot preserve the required positivity/conservation, fix
   that treatment for mapped gas sources before declaring that mode supported;
   do not hide failures with post-update clipping.
8. Persist any held operator coefficients, active masks, and validity times needed
   across restart. Per-step gas tendencies are reconstructed and are not replayed.
   Fail explicitly on missing/incompatible enabled cache state unless a tested
   deterministic reconstruction is available. Free map/cache/workspace data at
   finalize and test repeated initialization.

Introduce `config_micm_pbl_transport = .false.` as an independent opt-in control.
Its disabled path must preserve the baseline. Empty `config_micm_file` and
`MUSICA=false` remain no-ops. Resolve `suite` before checking support; log the
resolved scheme and transport capability. An enabled, unsupported combination
must stop during initialization with an actionable message. `pbl_scheme='off'`
is a supported zero-transport case; an active PBL with surface layer off is not.

## Scheme implementation stages

| Stage/configuration | Planned treatment | Gate before claiming support |
| --- | --- | --- |
| YSU, default PBL interval | Replace the zero tracer interface with runtime-sized inputs/outputs. Expose the actual scalar exchange/entrainment operator so it can be evaluated on fresh gas state, independent of the earlier meteorology solve. Use conservative dry-air weights and closed boundaries. | Resolve the `2*dt` convention, top/interface indexing, and moist-pressure versus dry-mass weighting with column reference tests. |
| YSU stable/unstable and top-down mixing | Retain the selected scheme's scalar diffusivity and entrainment response. Document precisely which passive-gas nonlocal terms exist; zero surface tracer flux does not justify copying thermal/moisture source terms. | Compare normal and `config_ysu_pblmix` configurations against the exposed scheme operator; any omitted gas entrainment/nonlocal term must be addressed or the configuration excluded. |
| MYNN with `config_mynn_edmf=0`, `config_mynn_edmf_dd=0` | Add a clean MPAS passive-gas interface using the scheme's scalar diffusivity and an implicit dry-mass flux-divergence solve. | Uniform, zero, pulse, sharp-gradient, and nonuniform-density columns; zero flux at both external boundaries. |
| MYNN with updraft EDMF | Export/reuse plume geometry and entrainment coefficients; recompute plume concentrations from current gas profiles, include plume tracer flux and compensating environment transport together with diffusion. | Updraft-only and combined diffusion/EDMF budgets and nonnegativity. `edmf_a`, `edmf_w`, or `exch_h` alone are not the complete operator. |
| MYNN downdraft EDMF | Extend the adapter through the downdraft pathway and its environmental compensation. | Remains explicitly unsupported until downdraft tracer budgets pass; do not silently use the updraft-only operator. |
| Nondefault PBL intervals | Diagnose and repair the alarm/`dt_pbl` contract as a separately reviewable change before qualification. Specify coefficient refresh times and per-step gas application. | Initially reject nondefault intervals when gas PBL transport is enabled. Lift the restriction only after interval/restart tests and native-meteorology regressions establish the intended cadence. Disabled gas transport must not acquire an incidental cadence change. |

Do not enable WRF chemistry/emissions macros globally. Adapt passive-gas
algorithms with no deposition velocity, fire enhancement, or imposed positive
floor. Use a monotone implicit diffusion solve and conservative flux limiting
or substeps for nonlocal terms. Any roundoff repair must be conservative and
reported; a zero species must remain zero.

Required MMM-physics exports are an external dependency. Prefer a versioned
MMM-physics change with a recorded immutable revision in the external manifest.
If a temporary patch is needed, check it into CheMPAS-A with source hash checks,
deterministic application, and retained license notices. An edited ignored
`physics_mmm` checkout is not a reproducible implementation deliverable.

## Commit-sized work and tracking

- [x] Audit #40, branch base, native ownership, scheme hooks, and time ordering.
- [x] Record this branch plan and link the companion effort.
- [ ] P1: shared map, empty/no-MUSICA path, block lifecycle, capability checks,
  dry-mass budget helper, and diagnostic schema. Publish its commit ID here for
  reuse by the convection branch; include exact Makefile dependencies.
- [ ] P2: standalone conservative column kernels and independent reference
  fixtures; prove `rho_zz` conversion and tendency consumption for each supported
  advection mode before connecting active mixing.
- [ ] P3: YSU runtime interface/operator exports, generic tendency assembly,
  default-interval integration, and stable/unstable/top-down qualification.
- [ ] P4: MYNN local diffusion, then updraft EDMF, then downdraft EDMF in separate
  commits. Keep unqualified combinations rejected at initialization.
- [ ] P5: nondefault interval contract, cache restart behavior, MPI/halo checks,
  and coupled PBL/convection/chemistry qualification after #41 is available.
- [ ] P6: GCC/gfortran and Clang/Flang eight-rank results, supported-configuration
  table, reproducible fixtures, and user configuration documentation.

Likely files: `Registry.xml`; chemistry initialization/runtime registration in
`chemistry/musica/mpas_atm_chemistry.F`; the new physics helper; physics driver,
interface, manager, control, initialization/finalization and tendency assembly;
YSU/MYNN wrappers; necessary MMM exports; and affected Makefiles. Modify core
time integration or halo groups only where the verified application contract
requires it. Keep shared foundation commits separate from scheme commits.

## Verification and completion criteria

Create version-controlled generators, namelists, mechanism variants, checks,
and small column reference data under `testing/chemistry_transport/`. Record
external input URLs, checksums, exact revisions, compiler/dependency versions,
build commands, MPI command, and elapsed simulated time. Large model inputs and
NetCDF outputs belong in the run-data tree, with a manifest and summarized
results tracked alongside this plan.

| Test group | Required evidence |
| --- | --- |
| Mapping/ownership | Empty/inactive chemistry; 1, 3, and 64 independently transported gases; reordered mechanisms; synthetic noncontiguous indices interleaved with native/passive fields; bound H2O and prescribed O2/N2. Every eligible species participates exactly once; excluded fields and unrelated scalars do not receive gas tendencies. |
| Column physics | Nonuniform dry density and layer spacing; uniform profiles; isolated lower/middle/top pulses; sharp gradients; zero and near-zero species; diffusion/nonlocal components separately; zero coefficients and large stable timesteps. Check every species, not just a sum over gases. Compare diffusion against an independently assembled matrix or analytic constant-coefficient mode. |
| Active PBL | Stable and convective boundary-layer fixtures with YSU and each qualified MYNN mode; show nonzero exchange and redistribution of every gradient tracer. Keep convection and chemical reactions off initially, with an inert MICM mechanism still loaded. |
| Time/parallel/restart | Repeat many updates with chemistry intervals of 1, 3, and 7 model steps; test qualified physics intervals separately, including no-mixing transitions. Compare one and eight ranks by `indexToCellID`, include ranks with no active columns, and restart between refresh times. Check thread-private buffers if OpenMP is enabled. |
| Coupled/regression | Reactive ABBA and Chapman/NOx cases, then PBL+convection+MICM. Check transport-only budgets before chemistry and reaction invariants afterwards. Compare disabled transport to the unmodified base, including native meteorology. Exercise each supported scalar advection setting. |
| Toolchains/output | Separate double-precision GCC/gfortran and Clang/Flang builds of the full pinned CheMPAS-A/MUSICA stack; eight-rank NetCDF runs in both. A build of plain MPAS or a prior GCC-only chemistry run does not satisfy this gate. |

Initial double-precision acceptance thresholds (freeze them with fixtures before
implementation; justify any later changes):

- Isolated closed-column relative species-mass residual <= `1e-12` per call and
  <= `1e-10` over 1,000 calls; full integration accumulated transport-budget
  residual <= `1e-8`. Use both absolute and relative errors. Normalize using
  `max(initial species mass, sum(m_d)*1e-30)` for per-area columns, and include
  area for domain budgets. Exactly zero species must remain zero.
- Uniform-profile and permutation comparisons: absolute error <= `1e-20 kg/kg`
  plus `1e-12` times the reference magnitude. Nonnegativity allows roundoff no
  smaller than `-(1e-20 + 1e-12*max_column(q_initial)) kg/kg`; larger violations
  fail even if a later chemistry clamp would remove them.
- Serial/eight-rank, restart/continuous, and short compiler-comparison fixtures:
  `atol=1e-18 kg/kg`, `rtol=1e-10` after cell-ID alignment. Report observed
  differences separately; long nonlinear meteorology is not a bitwise test.
- Inactive/disabled tendency is exactly zero. Disabled transport reproduces the
  base fields on the same build/settings, preferably bitwise, excluding logs and
  timing metadata. All output is finite; expected duration and restart events
  are reached; no rank reports an error.

Write per-species PBL tendency, integrated increment, min/max, before/after dry
mass, boundary flux, conservation residual, active-column count, coefficient
refresh time, and applied `dt` to selectable diagnostics/NetCDF or a compact
budget stream. Check the realized state increment as well as the computed
tendency so omitted or duplicate application is visible.

The prepared local baseline is `/Users/dfillmor/Data/CheMPAS`, with results in
`TEST-RESULTS.txt` and `summary.json`. Its idealized cases have PBL disabled;
they are regression tests, not evidence of PBL transport. Reuse their input and
NetCDF validation infrastructure, add deterministic active-PBL fixtures, and
track portable recipes rather than making that absolute path a requirement.

Complete #40 only after the supported subset is explicitly listed, all required
tests for that subset pass on both toolchains, remaining configurations fail
clearly, and the result manifest and implementation commit IDs are recorded
here. No transport test results are claimed by this planning commit.
