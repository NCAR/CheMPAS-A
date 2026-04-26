CheMPAS-A
=========

CheMPAS-A (Chemistry for MPAS - Atmosphere) is an ACOM integration pilot that
couples MUSICA/MICM atmospheric chemistry to MPAS-Atmosphere on its native
unstructured Voronoi mesh. This is CheMPAS-A 26.04, a pre-release currently
in beta testing, based on MPAS-Model v8.3.1.

CheMPAS-A serves as a rapid-prototyping ground for the chemistry coupling:
runtime tracer allocation, MUSICA/MICM state transfer, TUV-x photolysis,
and idealized chemistry test cases.

CheMPAS-A uses calendar versioning (YY.MM), tracked independently from the
MPAS base model version. This documentation includes the imported
MPAS-Atmosphere User's Guide (unmodified from the MPAS v8.3.1 release)
covering build, mesh preparation, I/O, physics, and runtime configuration,
together with the MPAS-Atmosphere Technical Description (verbatim port of
the v8 NCAR Technical Note draft) for the dynamical core, equations, and
spatial discretization.

.. toctree::
   :titlesonly:

   users-guide/index
   technical-description/index
