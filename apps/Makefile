include ../config.mk

# List of HPC applications to be compiled
INSTALL_APPS = qmcpack cp2k nwchem qe remhos hpl

qmcpack_REPO = --branch v3.12.0 https://github.com/QMCPACK/qmcpack qmcpack
cp2k_REPO = --branch v9.1.0 --recursive https://github.com/cp2k/cp2k cp2k
qe_REPO = --branch qe-7.0 https://github.com/QEF/q-e qe
nwchem_REPO = --branch v7.0.2-release https://github.com/nwchemgit/nwchem nwchem
remhos_REPO = https://github.com/CEED/Remhos remhos && cd remhos && git checkout 8d8ae5febbf042272cdcc8319dc9de9a11a8bdbb

hpl_TARBALL = http://www.netlib.org/benchmark/hpl/hpl-2.3.tar.gz

.DEFAULT_GOAL := all
.PHONY: all clean $(INSTALL_APPS)

all: $(INSTALL_APPS)

clean:
	rm -rf $(INSTALL_APPS)

qmcpack: qmcpack/.patched
	@echo "Compiling $(@)..."
	cd deps && make libxml2/.compiled libfftw/.compiled boost/.compiled hdf5/.compiled
	cd $(@)/build && export FFTW_HOME=$(CURDIR)/deps/install && export HDF5_ROOT=$(CURDIR)/deps/install \
	&& cmake -DCMAKE_C_COMPILER=$(MPICC) -DCMAKE_CXX_COMPILER=$(MPICXX) -DQMC_MIXED_PRECISION=0 \
	-DBLAS_openblas_LIBRARY=$(OPENBLAS_ROOT)/libopenblas.so -DLIBXML2_INCLUDE_DIR=$(CURDIR)/deps/install/include/libxml2 \
	-DLIBXML2_LIBRARY=$(CURDIR)/deps/install/lib/libxml2.so -DBoost_INCLUDE_DIR=$(CURDIR)/deps/install/include/ -DENABLE_OMP_TASKLOOP=OFF ../
	cd $(@)/build && make -j $(MAKE_JOBS)

cp2k: cp2k/.patched
	@echo "Compiling $(@)..."
	cd deps && make zlib/.compiled libfftw/.compiled hdf5/.compiled
	cd $(@)/tools/toolchain && export PATH=$(OPENBLAS_ROOT):$(CURDIR)/deps/install/bin:"$(PATH)" \
	&& export LD_LIBRARY_PATH=$(OPENBLAS_ROOT):$(CURDIR)/deps/install/lib:"$(LD_LIBRARY_PATH)" \
	&& ./install_cp2k_toolchain.sh --with-openblas=system --with-cmake=system --with-hdf5=system \
	--with-elpa=no --with-pexsi=no --with-cosma=no --with-sirius=no --with-spfft=no --with-spglib=no && cp install/arch/* ../../arch/
	cd $(@) && source tools/toolchain/install/setup && export LIBRARY_PATH=$(CURDIR)/deps/zlib/zlib-1.2.12:"$(LIBRARY_PATH)" \
	&& make -j $(MAKE_JOBS) ARCH=local VERSION="psmp"

qe: qe/.patched
	@echo "Compiling $(@)..."
	cd $(@) && mkdir -p build && cd build && export PATH=$(OPENBLAS_ROOT):"$(PATH)" \
	&& cmake ../ -DQE_ENABLE_MPI=ON -DQE_FFTW_VENDOR=Internal -DBLA_VENDOR=OpenBLAS && make -j $(MAKE_JOBS)
	cd $(@) && ln -s build/bin/ bin

nwchem: nwchem/.patched
	@echo "Compiling $(@)..."
	cd $(@)/src && export NWCHEM_TOP=$(CURDIR)/$(@) && export NWCHEM_TARGET=LINUX64 && export USE_MPI=y \
	&& export USE_MPIF=y && export NWCHEM_MODULES="all" && export BLASOPT="-L$(OPENBLAS_ROOT) -lopenblas" \
	&& export LAPACK_LIB="-L$(OPENBLAS_ROOT) -lopenblas" && export BLAS_SIZE=4 && make nwchem_config \
	&& make 64_to_32 && make -j $(MAKE_JOBS)

remhos: remhos/.patched
	@echo "Compiling $(@)..."
	cd deps && make mfem/.compiled
	cd $(@) && make MFEM_DIR=$(CURDIR)/deps/mfem -j $(MAKE_JOBS)

hpl/.patched: hpl/.downloaded
	@echo "Patching $(@D)..."
	patch -ruN -p1 -d $(@D)/hpl-2.3 < resources/$(@D)/$(hdfit_PATCH)
	cp -r resources/$(@D)/inputs/* $(@D)/hpl-2.3/
	@touch $(@D)/.patched

hpl: hpl/.patched
	@echo "Compiling $(@)..."
	cd $(@)/hpl-2.3 && export LD_LIBRARY_PATH=$(OPENBLAS_ROOT):"$(LD_LIBRARY_PATH)" \
	&& ./configure LDFLAGS="-L$(OPENBLAS_ROOT)" LIBS="-lopenblas" && make -j $(MAKE_JOBS)
