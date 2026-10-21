# TMT versus MUSICA validation

Validated on 2026-10-20 (America/New_York), on branch `chempas-a-tmt` based on
`89008ed50ce8f9ffd10876deed7c9a48ec573a38`. All nine configured cases passed with
both backends, including the five full wiki/template cases. The largest
field-scaled difference across all compared outputs and restart files was
**3.663e-11**. All nonchemistry floating-point fields were bitwise identical.

The required TMT compatibility fixes are recorded in [TMT `0f9ce16`](https://github.com/NCAR/TMT/commit/0f9ce1655c6dd24bdb45aaf3266d3860638bc3b9).

## Comparison conditions

Both model executables use GCC/gfortran 16.2.0, double precision, `-O3`, eight
MPI ranks, and one OpenMP thread. They were built from the same CheMPAS-A source
patch with `MUSICA=true` or `TMT=true`. Initial conditions, namelists, streams,
and partitions are identical. All 230 input hash entries
(86 distinct files) were verified against the case manifests.
MUSICA 0.16.5 supplies MIEM in both configurations; `mpas_miem.F` is unchanged.

Every numerical field and time record in `output.nc` and every scheduled
restart file was compared. Checks include complete simulation duration,
finite/unmasked output, zero logged errors, expected chemistry and photolysis
fields, ABBA stoichiometry and domain element conservation, and restart counts.
Integer and character fields also match exactly.

The errors below are `max(abs(TMT - MUSICA)) / max(abs(TMT), abs(MUSICA))`,
where each maximum spans a whole field and all records. Each table entry is
the largest such error among the indicated fields and files. These scales
describe whole fields; relative errors at individual near-zero values can
be larger. The strict roundoff reporting threshold is 64 storage epsilons
(approximately 1.42e-14 for double precision). Longer integrations and the
default-library photolysis comparison can exceed that threshold while
remaining very close. Solver tolerances and model physics were unchanged.

## Model results

Durations use `days_HH:MM:SS` or `HH:MM:SS`. Record counts refer to `output.nc`;
the comparison also includes restart files. Photolysis errors are zero for
ABBA cases, which do not enable photolysis.
The one-minute LNOx case retains zero NO and NO2; the Chapman case supplies
the evolving NOx comparison.

| Case | Duration | Output records | Maximum chemistry error | Maximum photolysis error |
| --- | --- | ---: | ---: | ---: |
| `mountain_wave_abba` | 00:10:00 | 2 | 7.234e-15 | 0.000e+00 |
| `jw_abba` | 01:00:00 | 2 | 3.068e-15 | 0.000e+00 |
| `supercell_abba` | 00:10:00 | 2 | 6.375e-15 | 0.000e+00 |
| `supercell_lnox` | 00:01:00 | 3 | 0.000e+00 | 1.366e-14 |
| `chapman_nox_global` | 01:00:00 | 2 | 1.225e-11 | 3.628e-11 |
| `full_mountain_wave_abba` | 05:00:00 | 31 | 1.191e-12 | 0.000e+00 |
| `full_chapman_nox_global` | 1_00:00:00 | 25 | 1.573e-11 | 3.663e-11 |
| `full_jw_abba` | 16_00:00:00 | 17 | 3.731e-14 | 0.000e+00 |
| `full_supercell_abba` | 02:00:00 | 13 | 4.952e-14 | 0.000e+00 |

## Floating-point control

The installed MUSICA/TUV-x uses the compiler's default floating-point
contraction setting, while TurboTUV specifies `-ffp-contract=off`. A separate
diagnostic executable rebuilds only the 128 TUV-x Fortran objects with
contraction disabled and preserves the other 35 MUSICA archive objects byte
for byte. The original build's 128 TUV-x objects were verified to match the
installed archive exactly. Neither production installation was replaced.

Matched-setting comparisons are recorded in `fp-control/summary.json`:

| Case | Maximum chemistry error | Maximum photolysis error | Photolysis fields bitwise identical |
| --- | ---: | ---: | --- |
| `supercell_lnox` | 0.000e+00 | 0.000e+00 | yes |
| `chapman_nox_global` | 4.279e-14 | 0.000e+00 | yes |
| `full_chapman_nox_global` | 1.899e-11 | 5.379e-11 | no |

The short LNOx and one-hour Chapman photolysis comparisons become bitwise
identical with matching contraction settings. In the one-day coupled case,
photolysis first differs at the second hourly output. Chemistry has already
introduced small differences in ozone, which the host passes back into TUV-x
(`mpas_atm_chemistry.F` passes `scalars(idx_qO3, :, iCell)` to
`tuvx_compute_photolysis`). This feedback is consistent with the small
remaining differences; matching compiler settings does not make the complete
coupled trajectory bitwise identical.

## Compatibility fixes and library checks

TMT required two compatibility fixes discovered by the model tests:

- Preserve unquoted YAML `NO` as a species name and accept plain strings in
  phase species lists. Regression tests cover both cases and Boolean/numeric
  parsing.
- Accept solar zenith angles in radians at the public MUSICA-compatible API,
  then convert to degrees for TurboTUV. The previous degree interpretation
  caused large photolysis errors. C, Fortran, Python documentation, examples,
  and regression tests now use the public radian convention. The C regression
  compares native-degree and public-radian calls bitwise at daytime, twilight,
  and nighttime angles.

Library validation passed: all 30 MUSICA CTest entries, all 13 TMT CTest
entries, the TMT/MUSICA MIEM coexistence test, and relocated installed-package
checks through CMake, pkg-config, C, Fortran, Python, and coupled chemistry.
The coexistence test exercises both allocation/error systems and MIEM in the
same process. The backend selection changes conditional `use` imports; shared
library-call instrumentation was subsequently added to both configurations.

## Library timing instrumentation

Shared timers now measure the `micm%solve` and `tuvx_solver%run` library calls
in both configurations, excluding MPAS conversion and logging. They report
rank minimum/mean/maximum accumulated wall seconds and global call counts.
The optional reference-chemistry solve is reported separately. Both
instrumented executables compile successfully; the numerical results above
come from the preceding executables without this instrumentation.

Sequential measurements with 18 MPI ranks are prepared and will run after
this initial validation commit is pushed. No instrumented timing results are
claimed in this revision. Earlier validation runtimes shared the machine;
the final TMT supercell was also paused while MUSICA finished. Those elapsed
times are not controlled performance measurements. The case inputs and
partition provenance are staged in `timing18/` in the local artifact directory.

## Retained evidence

The local artifact directory is
`/Users/dfillmor/Data/CheMPAS/tmt-validation-20261020`.
It contains binary snapshots, library hashes, compiler settings, dependency
revisions, source patches and added files, all run logs, per-case verification,
per-field numerical comparisons, and the scripts used to run and compare the
cases. `before-angle-fix/` preserves the superseded diagnostic runs separately.
`manifest.json`, `summary.json`, and `input-audit.json` record provenance and
results. The original validation binary snapshots are retained separately from
the `_timed` executables and their hashes. The numerical suite shared the
machine; `timing18/` stages the subsequent sequential instrumented runs.
