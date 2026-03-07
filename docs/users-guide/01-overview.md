# Chapter 1: MPAS-Atmosphere Overview

The Model for Prediction Across Scales -- Atmosphere (MPAS-A) is a non-hydrostatic atmosphere model that is part of a family of Earth-system component models collectively known as MPAS. All MPAS models have in common their use of centroidal Voronoi tessellations for their horizontal meshes, which has motivated the development of a common software framework that provides a high-level driver program and infrastructure for providing parallel execution, input and output, and other software infrastructure.

## 1.1 Features

Important features of MPAS-A include:

- Fully-compressible, non-hydrostatic dynamics
- Split-explicit Runge-Kutta time integration
- Exact conservation of dry-air mass and scalar mass
- Positive-definite and monotonic transport options
- Generalized terrain-following height coordinate
- Support for unstructured variable-resolution (horizontal) mesh integrations for the sphere and Cartesian planes
- Support for global and limited-area simulation domains

At present, MPAS-A includes parameterizations of physical processes taken from the Weather Research and Forecasting (WRF) Model[^1]. Specifically, MPAS-A has support for:

- **Radiation:** CAM and RRTMG long-wave and short-wave radiation schemes
- **Land-surface:** NOAH land-surface model
- **Surface-layer:** Monin-Obukhov and MYNN
- **Boundary-layer:** YSU and MYNN PBL schemes
- **Convection:** Kain-Fritsch, Tiedtke, New Tiedtke, and Grell-Freitas convection parameterizations
- **Cloud microphysics:** WSM6, Kessler, and Thompson schemes

[^1]: https://www.mmm.ucar.edu/weather-research-and-forecasting-model

## 1.2 Model Components

MPAS-A is comprised of two main components: the model, which includes atmospheric dynamics and physics; and an initialization component for generating initial conditions for the atmospheric and land-surface state, update files for sea-surface temperature and sea ice, and lateral boundary conditions. Both components (model and initialization) are built as *cores* within the MPAS software framework and make use of the same driver program and software infrastructure. However, each component is compiled as a separate executable.

```
+--------------------------------------+
|        MPAS superstructure           |
+------------------+-------------------+
| init_atmosphere  |    atmosphere     |
|                  |  +------+------+  |
|                  |  |      |      |  |
|                  |  +------+------+  |
+------------------+-------------------+
|      MPAS common infrastructure      |
+--------------------------------------+
```

*Figure 1.1: The initialization and model components of MPAS-A are built as separate cores within the MPAS framework.*

A succinct description of building and running MPAS-A is given in [Chapter 2](02-quick-start.md), the Quick Start Guide. Detailed instructions for building these components are given in [Chapter 3](03-building.md), and the basic steps to create initial conditions and run the MPAS-A model are outlined in [Chapter 7](07-running.md).
