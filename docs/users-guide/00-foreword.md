# CheMPAS-A User's Guide

CheMPAS-A 26.04 (pre-release / beta), based on MPAS-Atmosphere v8.3.1

Last updated: 26 April 2026

---

## Foreword

This user's guide describes the Model for Prediction Across Scales -- Atmosphere (MPAS-A) Version 8.3.0. MPAS-A is the non-hydrostatic atmosphere model built within the MPAS framework. Users guides for other MPAS components, such as MPAS-Ocean, are separate from this guide.

The component models and framework that comprise MPAS are being developed collaboratively between Los Alamos National Laboratory (LANL) and the U.S. National Science Foundation National Center for Atmospheric Research (NSF NCAR). Common functionality required by different MPAS component models, such as parallel input/output, time management, block decomposition, etc., is provided by the MPAS framework, while development of specific component models, referred to in MPAS as *cores*, is handled by the individual development groups. Currently, LANL is responsible for the ocean, land-ice, and sea-ice cores, while NSF NCAR is responsible for the atmospheric core, MPAS-A.

MPAS is very much a collaborative development of both the shared architectures and the component models. There are a number of contributors to the developments leading to the MPAS-A solver, and many of these developments are shared with the ocean core. The C-grid Voronoi discretization is based on critical developments from John Thuburn, Todd Ringler, Bill Skamarock, and Joe Klemp. The mesh generation that enables the MPAS-A development received major contributions from Todd Ringler, Doug Jacobson, Max Gunzburger, and Lili Ju. Significant developments in the transport scheme were accomplished by Bill Skamarock and Almut Gassmann. The atmospheric physics has been taken from the Advanced Research WRF model; these physics have benefitted from the work of hundreds of scientists. On the framework side we leverage a number of outside packages that have received extensive development from a wide community, including netCDF (UCAR/Unidata) and PIO (as used in the Community Earth Systems Model, CESM).

The software developed for MPAS is open source, and it has been copyrighted under a BSD license. The simple copyright statement can be found at the beginning of MPAS source files and the complete copyright statement can be found in this user's guide or in the `LICENSE` file accompanying the source code.

We conclude by noting that this user's guide is a work in progress. We welcome suggestions for improvements to this guide, including additions, corrections, clarifications, etc. Updates to MPAS-A, including the most recent code, user's guide, and test cases, may be found at <https://mpas-dev.github.io/>.

**Contributors to this guide:**
Michael Duda, Laura Fowler, Bill Skamarock, Conrad Roesch, Doug Jacobsen, and Todd Ringler.

*The U.S. National Science Foundation National Center for Atmospheric Research (NSF NCAR) is operated by the University Corporation for Atmospheric Research (UCAR) and is sponsored by the National Science Foundation. Any opinions, findings, conclusions, or recommendations expressed in this publication are those of the authors and do not necessarily reflect the views of the National Science Foundation.*
