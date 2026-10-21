.SUFFIXES: .F .o

.PHONY: all mpas_musica mpas_musica_lib clean

all: mpas_musica_lib

OBJS = \
  mpas_chem_constants.o \
  mpas_chem_timing.o \
  mpas_solar_geometry.o \
  mpas_lightning_nox.o \
  mpas_prescribed_fields.o \
  mpas_tuvx.o \
  mpas_atm_chemistry.o

ifeq ($(CHEMPAS_CHEMISTRY),true)
OBJS += mpas_musica.o mpas_miem.o
endif

mpas_musica: $(OBJS)

mpas_lightning_nox.o: mpas_chem_constants.o
mpas_tuvx.o: mpas_chem_constants.o mpas_prescribed_fields.o mpas_chem_timing.o
mpas_atm_chemistry.o: mpas_chem_constants.o mpas_lightning_nox.o mpas_solar_geometry.o \
  mpas_prescribed_fields.o mpas_tuvx.o \
  $(if $(filter true,$(CHEMPAS_CHEMISTRY)),mpas_musica.o mpas_miem.o)
mpas_musica.o: mpas_chem_constants.o mpas_chem_timing.o
mpas_miem.o: ../../../framework/mpas_derived_types.o ../../../framework/mpas_log.o

mpas_musica_lib: $(OBJS)
	$(RM) ../libchem.a
	ar -cr ../libchem.a $(OBJS)

clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../.. -I../../../framework $(MPAS_ESMF_INC)
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../.. -I../../../framework $(MPAS_ESMF_INC)
endif
