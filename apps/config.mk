# Base directory of the HDFIT OpenBLAS repository
OPENBLAS_ROOT = PATH_TO_CUSTOM_OPENBLAS_LIB

# List of HPC applications to be compiled
INSTALL_APPS = qmcpack cp2k nwchem qe remhos hpl

# MPI compilers plus other settings
MPICC = mpicc
MPICXX = mpicxx
MAKE_JOBS = 8

# --------- DO NOT EDIT PAST THIS POINT ---------

SHELL = /bin/bash

hdfit_PATCH = 0001-Integrating-HDFIT.patch

.PRECIOUS: %/.cloned %/.downloaded %/.patched %/.compiled

%/.cloned:
	@echo "Cloning $(@D)..."
	git clone $($(@D)_REPO)
	chmod -R u=rwX,go=rX $(@D)
	@touch $(@D)/.cloned

%/.downloaded:
	@echo "Downloading $(@D)..."
	mkdir -p $(@D)
	cd $(@D) && wget $($(@D)_TARBALL)
	cd $(@D) && tar -zxvf $(shell basename $($(@D)_TARBALL))
	chmod -R u=rwX,go=rX $(@D)
	@touch $(@D)/.downloaded

%/.patched: %/.cloned
	@echo "Patching $(@D)..."
	cd $(@D) && git am ../resources/$(@D)/$(hdfit_PATCH)
	cp -r resources/$(@D)/inputs/* $(@D)/
	@touch $(@D)/.patched
