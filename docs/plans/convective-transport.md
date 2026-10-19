# Implementation plan: convective transport of MICM tracers

Status: planned; transport implementation and qualification are pending.

- Tracking issue: [CheMPAS-A #41](https://github.com/NCAR/CheMPAS-A/issues/41).
- Branch: `chempas-a-convect`, from `chempas-a` at
  `89008ed50ce8f9ffd10876deed7c9a48ec573a38`.
- Companion/dependency: [PBL transport plan](https://github.com/NCAR/CheMPAS-A/blob/chempas-a-pbl/docs/plans/pbl-transport.md).
- Coordination: [MUSICA #1073](https://github.com/NCAR/musica/issues/1073).
- Source/issue audit: 2026-10-18. These are development plans, not supported
  capabilities. Standalone MICM contribution #39 is a separate effort.

Add parameterized convective redistribution for every runtime prognostic MICM
gas. Use the actual selected scheme's plume transport, with conservative dry-air
budgets. Emissions, lightning sources, scavenging, deposition, aqueous chemistry,
and aerosol transport remain separate work. Keep these processes disabled in
transport qualification; precipitation must not remove an inert gas.

## Findings at the branch base

| Location | Finding and consequence |
| --- | --- |
| `chemistry/musica/mpas_musica.F` | `chem_species(:)%mpas_index` supplies runtime indices. Host-bound H2O and prescribed O2/N2 must be excluded from the independent gas loop; the latter may use `qv` only as an index anchor. |
| `physics/mpas_atmphys_driver_convection.F` | Calls Grell-Freitas, Kain-Fritsch, original Tiedtke, and new Tiedtke. None of these MPAS calls carries a runtime gas array or generic gas tendency. |
| `physics/physics_wrf/module_cu_gf.mpas.F` | Deep `CUP_gf` and shallow `CUP_GF_SH` contain plume, entrainment/detrainment, environmental transport, and closure calculations. Host diagnostics expose `xmb_total`, `xmb_shallow`, and cloud levels, but not the full transport. `ishallow=1` is fixed in `mpas_atmphys_vars.F`. |
| `physics/physics_wrf/module_cu_kfeta.F` | Updraft/downdraft and exchange arrays (`UMF`, `DMF`, `UER`, `UDR`, `DER`, `DDR`) exist internally. They undergo closure adjustments and scaling. Exporting an early trial plume would give the wrong transport. |
| `physics/physics_wrf/module_cu_tiedtke.F` | Original Tiedtke has internal `cumastr`, up/down mass fluxes, and distinct convective types. The host wrapper does not export a complete tracer operator. |
| `physics/physics_wrf/module_cu_ntiedtke.F`, its pre/post wrappers, and MMM-physics `cu_ntiedtke.F90` | New Tiedtke is a separate implementation. `cumastrn`, `cuascn`, and downdraft routines carry mass/exchange information internally. It needs its own adapter and vertical-index audit. |
| `physics/mpas_atmphys_manager.F` | Convection is scheduled through `config_conv_interval`/`l_conv`, while the scheme calls receive `dt_dyn`; Kain-Fritsch also receives `cudt=dt_cu/60` and `n_cu`. These are distinct time conventions. |
| `physics/mpas_atmphys_driver_convection.F`: `update_convection_step1` | Kain-Fritsch's `nca` lifetime is decremented every step and native tendencies are cleared when it expires. Gas operator validity must follow this lifecycle too. |
| `mpas_atm_core.F` and `physics/mpas_atmphys_todynamics.F` | Physics precedes MICM, then dynamics assembles held native tendencies. Only fixed moisture indices receive convective contributions. A gas tendency built before chemistry must not be blindly reused afterwards. |
| `physics/mpas_atmphys_interface.F` | Physics density includes vapor (`zz*rho_zz*(1+qv)`), whereas authoritative gas mixing ratios use dry air. Flux units and weights need explicit conversion. |

Paths are relative to `src/core_atmosphere/`. The audited external revision is
[NCAR/MMM-physics a4baf7f3](https://github.com/NCAR/MMM-physics/tree/a4baf7f3243d1db0dbc5f63473f895bdbdc05c30),
selected by `Externals.cfg` tag `20250616-MPASv8.3`. Do not substitute a newer
WRF/GF implementation without recording and qualifying that scientific change.

## Shared integration contract

Reuse the foundation commit planned as P1 in `chempas-a-pbl`, then add the
convection adapter on this branch. Record that dependency's commit ID here when
available. Its species mapping, tendency assembly, diagnostics, and lifecycle
must remain shared; this branch does not depend on all PBL schemes being ready.
PBL-disabled tests can qualify convection independently.

The proposed interface has these requirements:

1. Map each non-host-bound, writable prognostic MICM species through its resolved
   `mpas_index`, validating uniqueness and bounds. Support empty maps, arbitrary
   species order/count, and noncontiguous indices. Reject unexplained native-field
   collisions. Record exclusions; never add a second moisture tendency.
2. Keep the map/cache in block-owned host data, populated after chemistry
   initialization. Physics consumes plain arrays and indices, without a dependency
   on MICM solver objects. `MUSICA=false` and empty chemistry produce no gas work.
   Retain the current one-block-per-rank chemistry limitation. Allocate reusable
   compact buffers outside threaded loops; use thread-private column scratch.
   Include all host vertical levels, even outside a restricted MICM solve domain.
   Bound workspace memory with tiles/species chunks and reuse meteorological
   coefficients across species.
3. Each scheme returns a species-independent transport description: finalized
   plume mass fluxes, entrainment/detrainment, source/termination levels,
   environmental compensation, and active/validity information. Build gas plume
   concentrations from current environmental gas state, rather than keeping old
   plume concentrations across chemistry or PBL updates.
4. Evaluate gas transport at the once-per-model-step tendency assembly after
   MICM. On scratch profiles compute `q1=PBL(q0,dt)` then `q2=CU(q1,dt)` and
   return `R_cu=(q2-q1)/dt`; if PBL is off, `q1=q0`. Add PBL and convection
   tendencies once through the shared mapped loop in `physics_get_tend`.
   The normal scalar advance alone changes authoritative state. Do not apply
   the scratch update separately or call it again at every RK stage.
5. Gas `q` is `kg species/kg dry air`; `R` is `dq/dt`; use dry layer mass per area
   `m_d=rho_dry*dz`. Verify the exact conversion into MPAS density-weighted
   `scalars_tend` against its `rho_zz` time level and scalar advance. The monotonic
   advance applies/clears sources differently from the ordinary advance; both
   require a realized-increment test before being listed as supported.
6. Clear all per-step tendencies, inactive columns, halo/padding slots, and expired
   cache entries. Compute only owned columns; use the established scalar halo
   exchanges at state consumption. Do not copy old gas buffers into halos or sum
   halo cells into budgets. Cache metadata must distinguish scheme-not-called,
   still-valid held convection, and a newly diagnosed inactive column.

Proposed opt-in control: `config_micm_convective_transport = .false.`. Resolve
the physics suite first and report the actual scheme. Scheme `off` is a valid
zero-transport case. Reject enabled but unqualified scheme/closure/cadence
combinations at initialization; never fall back to a generic diffusion surrogate.
Disabled gas transport must retain baseline native physics and chemistry.

## Conservative column operator

Use a plume/environment mass-exchange formulation consistent with the selected
scheme. Account for entrainment and detrainment at every participating layer,
including launch and final detrainment; include updrafts, available downdrafts,
and compensating environmental motion. Preserve the scheme's final closure,
ensemble weights, scale awareness, and active masks. Do not infer plume profiles
from precipitation, cloud-base mass flux alone, or a generic mixing depth.

For signed upward dry-air mass fluxes `M_u`, `M_d`, and `M_e`, require zero net
column air transport at an interface after environmental compensation:
`M_u+M_d+M_e=0`. The total gas flux is
`F_s=M_u*q_u,s+M_d*q_d,s+M_e*q_e,s`, including all plume contributions.
Apply `m_d(k)*R_s(k)=F_s(k-1/2)-F_s(k+1/2)` with zero external flux at the ground
and model top. Internal entrainment/detrainment must cancel between plume and
environment; do not count it again as an independent source after assembling the
total flux. A vertically uniform gas then has zero net transport tendency.

Document each scheme's sign, stagger, pressure unit, horizontal-area convention,
and moist/dry mass convention in the adapter. Audit vapor/condensate loading and
plume humidity before converting a moist plume flux into dry-air flux; a single
column-wide humidity factor is insufficient. Reconcile finalized dry-air
exchange continuity before transporting species, and compare against the scheme's
native exchange equations in its own coordinates.

Use conservative donor-mass substeps or a positive implicit exchange solve.
For explicit transfers, bound each substep by total outgoing dry mass from each
reservoir, using a documented safety factor (initial target <= 0.8 of available
mass), and recompute gas fluxes on the evolving profile. Share that transport
schedule across species. Never set negative concentrations to zero without
accounting for the lost mass, and never impose a positive floor on a zero gas.
When combined with resolved advection, verify actual state positivity; repair
the gas-source integration if necessary rather than relying on later MICM clamps.

## Scheme scope and implementation order

| Stage | Implementation | Qualification before enabling |
| --- | --- | --- |
| Grell-Freitas deep and shallow | Extend `CUP_gf` and `CUP_GF_SH` to export the final scaled exchange operator, including up/down mass entrainment and detrainment, closure/ensemble weights and environment terms. Combine deep and shallow contributions exactly once. | Separate deep-only, shallow-only, combined, downdraft, and inactive fixtures. Audit the current fixed `ishallow=1`; isolated component testing may use the column harness. Explicitly list validated `config_gfconv_closure_deep` and `config_gfconv_closure_shallow` values. |
| New Tiedtke (`cu_ntiedtke`) | Add exports through the MPAS pre/run/post wrappers and MMM routines; capture final updraft/downdraft exchange and supported convective types. | Deep/shallow cases and any reachable midlevel mode, vertical reversal/interface checks, and rain/no-rain inert budgets. This is the mesoscale suite's default and gets its own gate. |
| Kain-Fritsch (`cu_kain_fritsch`) | Expose the accepted, closure-adjusted up/down exchanges and integrate the operator over the actual cloud lifetime. Convert any cell-total flux with the correct cell area. | Triggered deep/shallow, downdraft, internal no-new-trigger/held-cloud, expiry, and retrigger cases; correct `nca`, seconds/minutes, and interval semantics. |
| Original Tiedtke (`cu_tiedtke`) | Add a separate adapter around its finalized `cumastr` transport and convective-type selection. | Deep/shallow/midlevel pathways that are reachable in this version, downdraft and inactive cases. Passing new Tiedtke does not qualify this implementation. |
| Unqualified variants | No gas transport support yet, including untested closures or combinations of MYNN EDMF and shallow convection. | Initialization must reject them while the gas option is on. Extend the support table only with recorded tests. |

Start with GF because the convection-permitting suite selects it; qualify its
deep and shallow paths together before advertising that scheme. Subsequent
stages extend scheme coverage without restricting the runtime species set.
For joint MYNN EDMF/shallow-cumulus runs, document which existing meteorological
parameterization owns each subgrid flux and ensure the gas adapters do not add
either flux twice. Do not change the meteorological closure to make gas budgets
pass; unsupported overlapping configurations remain explicit.

## Cadence, held operators, and restart

The meteorological refresh follows the existing convection alarm and each
scheme's internal activation/lifetime rules. Between refreshes, apply a valid
held *operator* to fresh gas profiles for the actual model step; do not reuse a
species tendency. Zero/replace a cache on a new inactive diagnosis. Kain-Fritsch
must retain its operator during a valid `nca` window and clear it when the cloud
expires, even if the next convection alarm has not rung.

Initially qualify `config_conv_interval='none'`, then integer intervals of 3 and
5 model steps, crossed with chemistry intervals of 1, 3, and 7 steps. Reject
shorter-than-step or nonintegral convection intervals until their scheduling is
defined and tested; do not silently round them through `n_cu`. Verify the first
step, skipped calls, last partial validity window, retrigger, and restarted clock.
Do not multiply a held tendency by `dt_cu` on every intervening dynamics step.

Restart the needed species-independent operator arrays, active flags, refresh
time, expiration/lifetime state, and configuration signature. Rebuild the species
map from the mechanism and validate compatibility. Reconstruct per-step gas
tendencies from current state. A restart taken between refreshes must match the
continuous run; silently discarding the held operator is not equivalent. Define
an explicit error or tested reconstruction for old restart files lacking cache
data. Coordinate this with the PBL cache without changing chemistry cadence.

## Commit-sized work and tracking

- [x] Audit #41, all four MPAS scheme selections, ownership, and time ordering.
- [x] Record this branch plan and its shared PBL dependency.
- [ ] C1: reuse the shared mapping/lifecycle/budget foundation commit from P1;
  add convection capability checks and capture/replay interfaces with no active
  gas transport. Record the shared commit ID here.
- [ ] C2: conservative plume/environment column kernel, independent exchange
  references, per-species diagnostics, and source-consumption tests.
- [ ] C3: GF deep/shallow final-operator exports and integration; active/inactive
  eight-rank fixtures with chemistry and PBL transport disabled.
- [ ] C4: new Tiedtke, Kain-Fritsch, and original Tiedtke in separately reviewable
  adapter/fixture commits; update support only after each gate passes.
- [ ] C5: held-operator lifetime, nondefault interval, restart, and halo/MPI
  qualification, including inactive transitions and no active columns on a rank.
- [ ] C6: coupled PBL/convection/MICM runs, GCC/gfortran and Clang/Flang eight-rank
  NetCDF results, supported configurations, and reproducible result manifest.

Likely changes are in `Registry.xml`, the shared physics helper, convection
driver/from/to/expiry routines, `mpas_atmphys_todynamics.F`, physics manager,
control and lifecycle, scheme wrappers and their internals, and Makefiles.
Core scalar integration changes are conditional on the source-consumption proof.
MMM changes for new Tiedtke require a versioned external revision or a tracked,
hash-checked, deterministically applied patch with retained license notices.
Do not rely on modifying an ignored external checkout. Preserve upstream
attribution and document any algorithmic differences from its passive tracers.

## Verification and completion criteria

Use the shared `testing/chemistry_transport/` fixture/checker layout from PBL.
Track small deterministic input generators, inert/reactive mechanisms, namelists,
streams, reference data, commands, and checksums. Keep large NetCDF results in a
run-data tree, with a version-controlled manifest and per-case/species summary.

| Test group | Required evidence |
| --- | --- |
| Species coverage | 0/inactive, 1, 3, and 64 transported gases; reordered mechanisms; noncontiguous maps with unrelated scalars between them; bound H2O and prescribed O2/N2. Distinct gradient/pulse profiles exercise every eligible scalar, including the last index. |
| Independent column physics | Hand-specified two/three-layer exchange cases and a separately assembled conservative transfer matrix; updraft only, downdraft, entrainment/detrainment, compensating environment, uniform and zero gases, nonuniform dry mass/layer spacing, and large-Courant substeps. Check units and reversed vertical coordinates independently of the implementation. |
| Real scheme activation | Active deep/shallow cases per qualified scheme, inert gases, no emissions/deposition/chemistry sources, and initial PBL gas transport off. Demonstrate nonzero finalized mass flux and species redistribution, not merely successful model completion or rainfall. Include neighboring inactive columns. |
| Cadence/lifecycle | Repeated refresh/hold cycles, chemistry refresh between convection calls, KF expiry/retrigger, disabled process, active-to-inactive transitions, and restart halfway through a held interval. Record expected/actual call counts and integrated increments. |
| Coupled physics | PBL+convection with inert tracers, then ABBA and Chapman/NOx chemistry. Separate transport mass residuals from chemical production/loss and check chemical invariants at the appropriate operator boundaries. |
| Parallel/toolchains | Align one/eight-rank outputs by cell ID; test eight-rank NetCDF runs using full pinned CheMPAS-A/MUSICA with both GCC/gfortran and Clang/Flang, separately built. Compare restart/continuous runs and short cross-compiler fixtures. |

Use the same initial double-precision thresholds as the PBL plan: closed-column
relative mass residual <= `1e-12` per call and <= `1e-10` over 1,000 calls;
full-run accumulated transport-budget residual <= `1e-8`. Normalize with
`max(initial species mass, sum(m_d)*1e-30)` per area, adding area for domain
budgets, and report absolute errors too. A zero species stays exactly zero.
Uniform/permutation profiles use `atol=1e-20 kg/kg`, `rtol=1e-12`;
nonnegativity permits only `-(1e-20 + 1e-12*max_column(q_initial)) kg/kg`
roundoff. Serial/eight-rank, restart, and short compiler comparisons use
`atol=1e-18 kg/kg`, `rtol=1e-10`. Freeze thresholds with fixtures before
implementation and document justified revisions; do not loosen them to hide a
systematic budget error. Inactive tendencies must be exactly zero and disabled
transport must reproduce the same-build baseline, preferably bitwise.

For every gas and component, report requested and realized increments,
before/after dry mass, boundary flux, conservation residual, min/max, substep
count, active columns, closure/type, operator refresh/expiry times, and applied
`dt`. Save selectable per-species convection tendencies and component budgets in
NetCDF or a compact diagnostic stream. Reject nonfinite output, missing fields,
incomplete duration, or errors on any rank.

The prepared local baseline at `/Users/dfillmor/Data/CheMPAS` supplies working
GCC binaries, data and NetCDF checks. Its idealized supercell and JW/Chapman
cases do not activate parameterized convection; resolved supercell overturning
does not validate this feature. Add deterministic active-column fixtures and a
documented active-convection model case. Stage only their meteorological inputs,
record URLs/checksums and trigger diagnostics, and keep recipes portable.

Complete #41 only when the declared scheme/closure/cadence subset passes both
toolchains, unsupported selections fail clearly, every transported species has
budget evidence, and the result manifest and implementation commits are recorded
here. No implementation or transport test results are claimed by this plan.
