# Chapter 1: Introduction

The Model for Prediction Across Scales (MPAS) is a collaborative project for developing atmosphere, ocean and other earth-system simulation components for use in climate, regional climate and weather studies. The primary development partners are the climate modeling group at Los Alamos National Laboratory (COSIM) and the National Center for Atmospheric Research. Both primary partners are responsible for the MPAS framework, operators and tools common to the applications; LANL has primary responsibility for the ocean and land ice models, and NCAR has primary responsibility for the atmospheric model.

The defining features of MPAS are the unstructured Voronoi meshes and C-grid discretization used as the basis for many of the model components. The unstructured Voronoi meshes, formally Spherical Centroidal Voronoi Tesselations (SCVTs), allow for both quasi-uniform discretization of the sphere and local refinement. The C-grid discretization, where the normal component of velocity on cell edges is prognosed, is especially well-suited for higher-resolution, mesoscale atmosphere and ocean simulations.

## 1.1 MPAS Atmosphere

The atmospheric component of MPAS, as with all MPAS components, uses an unstructured centroidal Voronoi mesh (grid, or tessellation) and C-grid staggering of the state variables as the basis for the horizontal discretization in the fluid-flow solver. The unstructured variable resolution meshes can be generated having smoothly-varying mesh transitions as illustrated in Figure 1.1, and this capability ameliorates many issues associated with the traditional mesh refinement strategy of one-way and two-way grid nesting where the transitions are abrupt. Using the flexibility of the MPAS meshes, we are working towards applications in high-resolution numerical weather prediction (NWP) and regional climate, in addition to global uniform-resolution NWP and climate applications.

**[Figure 1.1: A variable resolution MPAS global mesh. To be added next session.]**

The MPAS atmosphere consists of an atmospheric fluid-flow solver (the dynamical core) and a subset of the Advanced Research WRF (ARW, Skamarock et al. 2021b) model atmospheric physics.

## 1.2 Technical Note Purpose and Structure

The purpose of this technical note is to provide a description of the MPAS-Atmosphere model. Importantly, this includes a complete description of the governing equations in the fluid-flow solver and their discretization. While aspects of the solver are described in the peer-reviewed literature, a comprehensive and detailed description appears only here in this technical note, along with updates to the techniques. It is our intent to keep this technical note up-to-date with the latest MPAS-Atmosphere release.

We begin in Chapter 2 with a description of the unstructured horizontal MPAS centroidal Voronoi mesh and with the structured vertical component that make up the full 3-dimensional MPAS mesh. An understanding of this mesh is a prerequisite to understanding the discretization of the dynamical core (the MPAS-Atmosphere fluid-flow solver). Chapter 3 introduces the continuous fluid flow equations used in MPAS-A and describes the time integration of these equations. MPAS-A uses a split-explicit time integration technique (similar to that used in the WRF model, Skamarock et al. 2021b) that requires defining mean and perturbation quantities for the prognostic variables and time-dependent equations for these variables that are used in specific sections of the solver. Chapter 4 describes the spatial discretization of the fluid-flow equations. Chapter 5 describes the discretization of the primary filters employed in MPAS-A, including 2nd and 4th order horizontal filters for cell-centered variables and for the horizontal velocity, and some aspects of filters applied to absorb vertically-propagating gravity waves. Chapter 6 describes the computation of the cell-centered zonal and meridional velocity using radial basis functions, and the horizontal frontogenesis function. Chapter 7 describes the regional configuration of MPAS focussing on the boundary conditions and how they are applied in the time integration. Chapter 8 describes MPAS-A initialization options, focussing on the real-data initialization but also providing an overview of available idealized cases. Initialization of an idealized or real atmospheric state that can be used in MPAS is accomplished by a separate MPAS model call *init_atmosphere_model*, in contrast to *atmosphere_model*. Chapter 9 provides a brief overview of physics parameterizations available in MPAS-Atmosphere. The physics schemes available in MPAS are evolving quickly both in their formulation and in how MPAS-A accesses them. The final chapter in the technical note, Chapter 10, outlines key aspects of the MPAS infrastructure that supports the atmosphere model.

In additional to describing the MPAS-A model and utilities, this technical note references the MPAS model code, and references to that code will appear in boxes with a pink background:

:::{admonition} MPAS code
:class: note

Here is an example of how a statement referencing the MPAS code appears in this technical note. File names in the MPAS source code directory are given in blue. The main text is italicized to further distinguish it from the main body of the technote. Within an **MPAS code** section, subroutine names, variable names and named Fortran loops are given in red.
:::

In future MPAS-Atmosphere releases we will include comments in the MPAS source code referencing this technical note, and we will be adding additional markers in the code, and references to them in this technical note, that are easily found using simple search techniques.

## 1.3 MPAS-Atmosphere Code Overview

MPAS-Atmosphere is comprised of two main executables that can be built from the system contained in the MPAS-A Github repository, and these executables are called *cores* in the MPAS system. The *init_atmosphere* core is responsible for producing atmospheric states that the *atmosphere* core can use as an initial state that it integrates forward in time.

### 1.3.1 The init_atmosphere core

The *init_atmosphere* core is responsible for

- Interpolating *static* fields to the MPAS horizontal mesh, including terrain heights, land-use categorizations, etc.
- Creating the 3-dimensional MPAS mesh. See Chapter 2.
- Interpolating analytic or an idealized initial state to the 3D MPAS mesh.
- Interpolating 3D analyses, or 3D states from other MPAS integrations, to the 3D MPAS mesh.
- Enforcing hydrostatic balance on the interpolated 3D state. See Chapter 8 for information on interpolations and the MPAS-Atmosphere discrete hydrostatic balance.
- Interpolating sea-surface temperatures (SSTs) and other fields at times other than the initial state time for use by the *atmosphere* core to update these fields during MPAS-Atmosphere integrations.

:::{admonition} MPAS code
:class: note

The *init_atmosphere* core code can be found in `MPAS/src/core_init_atmosphere/`. The main functionality that users access is contained in the source file `MPAS/src/core_init_atmosphere/mpas_init_atm_cases.F`, that contains the routines that drive the idealized and real-data initializations.
:::

### 1.3.2 The atmosphere core

The *atmosphere* core is responsible for integrating the MPAS-Atmosphere state forward in time. The integration is accomplished as depicted in this flow chart given in Figure 1.2. The numbered tasks in the flow chart are explained as follows:

**[Figure 1.2: A top level flow chart for MPAS-Atmosphere. To be added next session.]**

(1) At startup necessary infrastructure is initialized, including I/O, parallel processing capabilities, clocking capabilities, etc.

:::{admonition} MPAS code
:class: note

The code that controls the sequencing of the flow chart items (2) through (5) can be found in `MPAS/src/core_atmosphere/mpas_atm_core.F` and function `atm_core_init`.
:::

(2) There are two modes used to begin an integration - a restart mode (starting from a *restart* file produced by a previous simulation) of starting from an MPAS-Atmosphere initialization file produced be the MPAS *init_atmosphere* core. For more information see the MPAS-Atmosphere Users Guide.

:::{admonition} MPAS code
:class: note

Configuration as a restart run, as opposed to starting from an MPAS-A initialization, is controlled by the `namelist.atmosphere` configuration variable `config_do_restart` in the `&restart` namelist section.
:::

(3) The MPAS-A initialization files only contained the uncoupled atmospheric state (e.g. the horizontal velocity $u$ as opposed to the horizontal momentum $\rho u$) and only fundamental mesh variables. Mesh metrics and other coefficients that can be precomputed are evaluated here. Additionally, the prognostic variables in the MPAS-A solver are the coupled variables (see section 3.1), and here the prognostic state is computed along with any state-dependent diagnostics. In a MPAS-A restart the coupled state would be read from the restart file.

:::{admonition} MPAS code
:class: note

The routine that does the coupling (subroutine `atm_init_coupled_diagnostics`) is found in `MPAS/src/core_atmosphere/dynamics/mpas_atm_time_integration.F`, and diagnostic variables (e.g. vorticity, etc) are computed in subroutine `atm_compute_solve_diagnostics`. These computations are described later in the technical note.
:::

(4) Any aspects of the physics that need initialization is performed, including computing or reading from files any needed data (e.g. lookup tables used in radiation and microphysics) not included in the initialization or restart files.

:::{admonition} MPAS code
:class: note

The MPAS-Atmosphere physics code is located in `MPAS/src/core_atmosphere/physics/`, including any code associated with the init/run/finalize structure of the MPAS-A model.
:::

(5) MPAS-A uses *alarms* to determine when tasks that need to be executed periodically are to be performed during the integration. The bulk of these alarms are used to determine when I/O tasks need to be performed or when physics need to be called.

:::{admonition} MPAS code
:class: note

The code that controls the sequencing of the flow chart items (6) through (10) can be found in `MPAS/src/core_atmosphere/mpas_atm_core.F` and function `atm_core_run`.
:::

(6) MPAS-A will produce stream output files at the initial time for a simulation starting from an MPAS-A initialization (i.e. a file produced by the MPAS-A *init_atmosphere* core). If the simulation is configured as a restart, then these files will not be produced as they will have been produced by the previous simulation that produced the restart file.

(7) At the start of each time step, alarms are checked to see if lateral boundary condition (LBC) data or other data need to be read in and fields updated.

(8) While forcings from most atmospheric physics are updated each MPAS-Atmosphere timestep, some physics forcings (e.g. radiation) are updated less frequently, hence the alarm checks.

(9) This is the main time step for them model, and a high-level overview can be found in section 3.2.

(10) Output (restarts, history and diagnostic output) are performed here.

As described earlier in this chapter, the bulk of this technical note describes the spatial and temporal discretization utilized in the time integration. Atmospheric physics used in this integration is briefly described in Chapter 9. Aspects of the Infrastructure are described in Chapter 10. While some practical aspects of running the MPAS-Atmosphere model are given in the **MPAS code** callouts, a much more detailed and useful guide to running MPAS is found in the MPAS Users Guide:
<https://www2.mmm.ucar.edu/projects/mpas/site/documentation/users_guide.html>.

Additional descriptions of the MPAS systems and instructions for using it are given in the MPAS-Atmosphere tutorials, with tutorial presentations available at
<https://www2.mmm.ucar.edu/projects/mpas/tutorial/Virtual2024/agenda.html>.
