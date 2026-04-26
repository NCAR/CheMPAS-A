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

.. admonition:: Under construction
   :class: warning

   These docs are an active port and are not yet feature-complete. In
   particular, figures from the MPAS-Atmosphere User's Guide and Technical
   Description still need to be regenerated — placeholders of the form
   ``**[Figure N.M: caption. To be added next session.]**`` mark the
   intended location of each figure in the Technical Description, and
   figure references in the User's Guide (e.g., Figure 9.1, the vertical
   grid schematics in Appendix C) currently render without their source
   images. CheMPAS-specific chemistry chapters covering MUSICA/MICM,
   TUV-x, and runtime tracer allocation are not yet wired into this
   tree; those notes live under ``docs/chempas/`` for now.

.. toctree::
   :titlesonly:

   users-guide/index
   technical-description/index
