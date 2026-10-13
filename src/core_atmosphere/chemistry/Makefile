.PHONY: all clean

# The host interface and its helpers are required even without MUSICA.
all:
	$(MAKE) -C musica all

clean:
	@# Also remove artifacts left in the parent directory by the old layout.
	$(RM) *.o *.mod *.f90 libchem.a
	$(MAKE) -C musica clean
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i
