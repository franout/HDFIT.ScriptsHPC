# Build directories of the HDFIT OpenBLAS and LLTFI distributions - do not add line breaks!
OPENBLAS_ROOT ?= PATH_TO_OPENBLAS_INSTALLATION
LLTFI_ROOT ?= PATH_TO_LLTFI_INSTALLATION

# MPI compilers plus other settings
MPICC ?=  mpicc
MPICXX ?=  mpicxx
CC ?= clang
CXX ?= clang++

MAKE_JOBS = 40

  
# --------- DO NOT EDIT PAST THIS POINT ---------

LLTFI_LINKER = -L$(LLTFI_ROOT)/runtime_lib/ -lllfi-rt
LLTFI_PLUGIN = -Xclang -fpass-plugin=$(LLTFI_ROOT)/llvm_passes/llfi-passes.so -Xclang -load -Xclang $(LLTFI_ROOT)/llvm_passes/llfi-passes.so
LLTFI_SETTINGS = $(LLTFI_PLUGIN) -mllvm -insttype -mllvm -includeinst=fadd -mllvm -includeinst=fmul -mllvm \
                 -includeinst=fsub -mllvm -includeinst=fdiv -mllvm -regloc -mllvm -dstreg

hdfit_PATCH = 0001-Integrating-HDFIT.patch

SHELL = /bin/bash

.PRECIOUS: %/.cloned %/.downloaded %/.touched %/.patched %/.compiled

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

%/.touched: %/.cloned
	@echo "Copying data dependencies for $(@D)..."
	cp -r resources/$(@D)/inputs/* $(@D)/ || :
	@touch $(@D)/.touched

%/.patched: %/.touched 
	@echo "Patching $(@D)..."
	cd $(@D) && git am ../resources/$(@D)/$(hdfit_PATCH)
	@touch $(@D)/.patched
