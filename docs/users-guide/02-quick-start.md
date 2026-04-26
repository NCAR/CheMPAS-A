# Chapter 2: CheMPAS-A Quick Start Guide

This chapter provides CheMPAS-A users with a high-level description of the general process of building and running the model. It is meant simply as a brief overview of the process, with more detailed descriptions of each step provided in later chapters.

## Build Process

In general, the build process follows the steps below.

1. Build or locate an MPI implementation (MPICH, OpenMPI, MVAPICH2, MPT, etc; [Section 3.1](03-building.md#31-prerequisites)).
2. Build the serial netCDF library ([Section 3.2.1](03-building.md#321-netcdf)).
3. Build the parallel-netCDF library ([Section 3.2.2](03-building.md#322-parallel-netcdf)).
4. Build the Parallel I/O (PIO) library ([Section 3.2.3](03-building.md#323-pio)).
5. (OPTIONAL) Build the METIS package ([Section 4.1](04-preparing-meshes.md#41-graph-partitioning-with-metis)).
6. Obtain the MPAS-Model source code.
7. Build the *init_atmosphere* and *atmosphere* cores ([Section 3.3](03-building.md#33-compiling-mpas)).

After completing these steps, executable files named `init_atmosphere_model`, `atmosphere_model`, and `build_tables` should have been created in the top-level MPAS-Model directory.

## Running a Basic Global Simulation

Once all three executables have been created, a basic global simulation can be performed using the following steps.

1. Create a run directory.
2. Link the `init_atmosphere_model` and `atmosphere_model` executables to the run directory, as well as physics lookup tables (`src/core_atmosphere/physics/physics_wrf/files/*`).
3. Copy the `namelist.*`, `streams.*`, and `stream_list.*` files to the run directory.
4. Edit the namelist files and the stream files appropriately ([Chapter 7](07-running.md)).
5. (OPTIONAL) Prepare meshes for the simulation ([Chapter 4](04-preparing-meshes.md)).
6. Run `init_atmosphere_model` to create initial conditions ([Section 7.1](07-running.md#71-creating-idealized-ics) or [Section 7.2](07-running.md#72-creating-real-data-ics)).
7. Run `atmosphere_model` to perform model integration ([Section 7.3](07-running.md#73-running-the-model)).
