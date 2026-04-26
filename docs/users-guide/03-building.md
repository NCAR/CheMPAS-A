# Chapter 3: Building MPAS

## 3.1 Prerequisites

To build MPAS, compatible C and Fortran compilers are required. Additionally, the MPAS software relies on the PIO parallel I/O library to read and write model fields, and the PIO library requires the standard netCDF library as well as the parallel-netCDF library from Argonne National Laboratory. All libraries must be compiled with the same compilers that will be used to build MPAS. Section 3.2 summarizes the basic procedure of installing the required I/O libraries for MPAS.

In order for the MPAS makefiles to find the PIO, parallel-netCDF, and netCDF include files and libraries, the environment variables `PIO`, `PNETCDF`, and `NETCDF` should be set to the root installation directories of the PIO, parallel-netCDF, and netCDF installations, respectively.

An MPI installation such as MPICH or OpenMPI is also required, and there is no option to build a serial version of the MPAS executables. MPAS-Atmosphere v5.0 introduces the capability to use hybrid parallelism using MPI and OpenMP; however, the use of OpenMP *should be considered experimental* and generally does not offer any performance advantage. The primary reason for releasing a shared-memory capability is to make this code available to collaborators for future development.

## 3.2 Compiling I/O Libraries

> **IMPORTANT NOTE:** *The instructions provided in this section for installing libraries have been successfully used by MPAS developers, but due to differences in library versions, compilers, and system configurations, it is recommended that users consult documentation provided by individual library vendors should problems arise during installation. The MPAS developers cannot assume responsibility for third-party libraries.*

Although most recent versions of the netCDF and parallel-netCDF libraries should work, the most tested versions of these libraries are netCDF 4.4.x and parallel-netCDF 1.8.x. Users are strongly encouraged to use either the latest PIO 2.x version from <https://github.com/NCAR/ParallelIO/>, or PIO versions 1.7.1 or 1.9.23, as other versions have not been tested or are known to not work with MPAS. The netCDF and parallel-netCDF libraries must be installed before building the PIO library.

### 3.2.1 NetCDF

Version 4.4.x of the netCDF library may be downloaded from <http://www.unidata.ucar.edu/downloads/netcdf/index.jsp>. The Unidata page provides detailed instructions for building the netCDF C and Fortran libraries; both the C and Fortran interfaces are needed by PIO. If netCDF-4 support is desired, the zlib and HDF5 libraries will need to be installed prior to building netCDF. *Before proceeding to compile PIO the `NETCDF` environment variable should be set to the netCDF root installation directory.*

### 3.2.2 Parallel-NetCDF

Version 1.8.x of the parallel-netCDF library may be downloaded from <https://trac.mcs.anl.gov/projects/parallel-netcdf/wiki/Download>. *Before proceeding to compile PIO the `PNETCDF` environment variable should be set to the parallel-netCDF root installation directory.*

### 3.2.3 PIO

Beginning with the MPAS v5.2 release, either of the PIO 1.x or 2.x library versions may be used. The two major versions have slightly different APIs; by default, the MPAS build system assumes the PIO 1.x API, but the PIO 2.x library versions may be used by adding the `USE_PIO2=true` option when compiling MPAS as described in Section 3.3.

If compiling with the PIO 1.x library versions, users are strongly encouraged to choose either PIO 1.7.1 or PIO 1.9.23, as other 1.x versions may not work; these two specific versions may be obtained from:
- <https://github.com/NCAR/ParallelIO/releases/tag/pio1_7_1>
- <https://github.com/NCAR/ParallelIO/releases/tag/pio1_9_23>

The PIO 2.x library versions support integrated performance timing with the GPTL library; however, the MPAS infrastructure does not currently provide calls to initialize this library when it is used in PIO 2.x. Therefore, it is recommended to add `-DPIO_ENABLE_TIMING=OFF` when running the cmake command to build PIO 2.x versions.

After PIO is built and installed the `PIO` environment variable should be set to the directory where PIO was installed. Recent versions of PIO support the specification of an installation prefix, while some older versions do not, in which case the `PIO` environment variable should be set to the directory where PIO was compiled.

## 3.3 Compiling MPAS

> **IMPORTANT NOTE:** *Before compiling MPAS, the `NETCDF`, `PNETCDF`, and `PIO` environment variables must be set to the library installation directories as described in the previous section.*

The MPAS code uses only the `make` utility for compilation. Rather than employing a separate configuration step before building the code, all information about compilers, compiler flags, etc., is contained in the top-level `Makefile`; each supported combination of compilers (i.e., a configuration) is included in the `Makefile` as a separate make target, and the user selects among these configurations by running make with the name of a build target specified on the command-line, e.g.,

```
make gfortran
```

to build the code using the GNU Fortran and C compilers. The available targets are listed in the table below, and additional targets can be added by editing the `Makefile` in the top-level directory.

| Target | Fortran compiler | C compiler | MPI wrappers |
|--------|-----------------|------------|--------------|
| `xlf` | xlf90 | xlc | mpxlf90 / mpcc |
| `pgi` | pgf90 | pgcc | mpif90 / mpicc |
| `ifort` | ifort | gcc | mpif90 / mpicc |
| `gfortran` | gfortran | gcc | mpif90 / mpicc |
| `llvm` | flang | clang | mpifort / mpicc |
| `bluegene` | bgxlf95_r | bgxlc_r | mpxlf95_r / mpxlc_r |

The MPAS framework supports multiple *cores* -- currently a shallow water model, an ocean model, a land-ice model, a non-hydrostatic atmosphere model, and a non-hydrostatic atmosphere initialization core -- so the build process must be told which core to build. This is done by either setting the environment variable `CORE` to the name of the model core to build, or by specifying the core to be built explicitly on the command-line when running make. For the atmosphere core, for example, one may run either:

```
setenv CORE atmosphere
make gfortran
```

or:

```
make gfortran CORE=atmosphere
```

If the `CORE` environment variable is set and a core is specified on the command-line, the command-line value takes precedence; if no core is specified, either on the command line or via the `CORE` environment variable, the build process will stop with an error message stating such. Assuming compilation is successful, the model executable, named `${CORE}_model` (e.g., `atmosphere_model`), should be created in the top-level MPAS directory.

In order to get a list of available cores, one can simply run the top-level `Makefile` without setting the `CORE` environment variable or passing the core via the command-line:

```
> make
( make error )
make[1]: Entering directory '/scratch/MPAS-Release'

Usage: make target CORE=[core] [options]

Example targets:
    ifort
    gfortran
    xlf
    pgi

Available Cores:
    atmosphere
    init_atmosphere
    landice
    ocean
    seaice
    sw
    test

Available Options:
    DEBUG=true       - builds debug version. Default is optimized version.
    USE_PAPI=true    - builds version using PAPI for timers. Default is off.
    TAU=true         - builds version using TAU hooks for profiling. Default is off.
    AUTOCLEAN=true   - forces a clean of infrastructure prior to build new core.
    GEN_F90=true     - Generates intermediate .f90 files through CPP, and builds with them.
    TIMER_LIB=opt    - Selects the timer library interface to be used for profiling the model.
                       TIMER_LIB=native - Uses native built-in timers in MPAS
                       TIMER_LIB=gptl   - Uses gptl for the timer interface
                       TIMER_LIB=tau    - Uses TAU for the timer interface
    OPENMP=true      - builds and links with OpenMP flags. Default is to not use OpenMP.
    USE_PIO2=true    - links with the PIO 2 library. Default is to use the PIO 1.x library.
    PRECISION=single - builds with default single-precision real kind. Default is double-precision.

************ ERROR ************
No CORE specified. Quitting.
************ ERROR ************
```

## 3.4 Selecting a Single-Precision Build

Beginning with version 2.0, MPAS-Atmosphere can be compiled and run in single-precision, offering faster model execution and smaller input and output files. Beginning with version 5.0, the selection of the model precision can be made on the command-line, with no need to edit the `Makefile`. To compile a single-precision CheMPAS-A executable, add `PRECISION=single` to the build command, e.g.,

```
make gfortran CORE=atmosphere PRECISION=single
```

Regardless of which precision the CheMPAS-A `init_atmosphere` and `atmosphere` cores were compiled with, either single- or double-precision input files may be used. In general, the MPAS infrastructure should correctly detect the precision of input files, but one may also explicitly specify the precision of files in an input stream by adding the `precision` attribute to the stream definition as described in [Section 5.2](05-configuring-io.md#52-optional-stream-attributes).

## 3.5 Building with Chemistry (MUSICA) Support

CheMPAS-A's chemistry features — runtime species discovery from MICM,
MUSICA/MICM coupling, and TUV-x photolysis — require an external
MUSICA-Fortran installation and an additional build flag. CheMPAS-A
supports two host platforms: **macOS** (LLVM/flang) and **Ubuntu**
(GCC/gfortran via conda).

A preflight script, `scripts/check_build_env.sh`, auto-detects the
toolchain and exports the required environment (`NETCDF`, `PNETCDF`,
`PIO`, `PKG_CONFIG_PATH`). Source it in the same shell as `make`:

```
eval "$(scripts/check_build_env.sh --export)"
```

`PKG_CONFIG_PATH` must be present in that shell, because the Makefile
invokes `pkg-config` at parse time; exporting it in a separate shell and
then running `make` later is not sufficient.

To compile with chemistry support, add `MUSICA=true` to the build
command. On **macOS** (LLVM/flang):

```
eval "$(scripts/check_build_env.sh --export)" && make -j8 llvm \
    CORE=atmosphere MUSICA=true \
    PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" PRECISION=double
```

On **Ubuntu** (GCC/gfortran via conda; requires `conda activate mpas`):

```
eval "$(scripts/check_build_env.sh --export)" && make -j8 gfortran \
    CORE=atmosphere MUSICA=true \
    PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" PRECISION=double
```

Without `MUSICA=true`, the chemistry hooks compile out and the `&musica`
namelist record is ignored at runtime. The MUSICA-Fortran library must
be built with the same Fortran compiler used for CheMPAS-A — flang and
gfortran `.mod` files are not interchangeable. Refer to the MUSICA
documentation at <https://musica.readthedocs.io/> for MUSICA-Fortran
build instructions.

For full preflight, dependency-build, and troubleshooting notes — including
the LLVM-on-macOS and conda-on-Ubuntu setup — see `BUILD.md` in the
repository root. For the runtime chemistry features this build flag
enables, see [Chapter 7](07-runtime-tracers.md) and
[Chapter 8](08-chemistry-coupling.md).

## 3.6 Cleaning

To remove all files that were created when the model was built, including the model executable itself, make may be run for the `clean` target:

```
make clean
```

As with compiling, the core to be cleaned is specified by the `CORE` environment variable, or by specifying a core explicitly on the command-line with `CORE=`.
