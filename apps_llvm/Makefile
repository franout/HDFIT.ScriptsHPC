include ../config.mk

# List of HPC applications to be compiled
INSTALL_APPS = qmcpack gromacs seissol milc miniw

qmcpack_REPO = --branch v3.12.0 https://github.com/QMCPACK/qmcpack qmcpack
gromacs_REPO = --branch v2022.3 https://github.com/gromacs/gromacs gromacs
seissol_REPO = --branch v0.9.0 --recursive https://github.com/SeisSol/SeisSol seissol
milc_REPO = https://github.com/milc-qcd/milc_qcd milc && cd milc && git checkout f03f53186e84d2d95cb5c4f4b402808925db685c
miniw_REPO = https://github.com/mrnorman/miniWeather miniw && cd miniw && git checkout 052d87a5543359b1d470cc0b1fcb3cb37edf95fe
yakl_REPO = https://github.com/mrnorman/YAKL && cd YAKL && git checkout a298bdefbf353b2a5c49688960b2dce2c24c7fc4
gadget_REPO = https://gitlab.mpcdf.mpg.de/vrs/gadget4 gadget && cd gadget && git checkout c838e262dc0345ca35698c20fff855f0e5a9f68d

# "1" for mixed precision, "0" for double precision
qmcpack_PREC = -DQMC_MIXED_PRECISION=1
qmcpack_OPTS = -march=haswell
# "off" for mixed precision, "on" for double precision
gromacs_PREC = -DGMX_DOUBLE=off
gromacs_OPTS = -DGMX_SIMD=AVX2_256
# "single" and "double" for single or double precision
seissol_PREC = -DPRECISION=single
seissol_OPTS = -DHOST_ARCH=hsw
# "1" for single precision, "2" for double precision
milc_PREC = PRECISION=1
milc_OPTS = ARCH=hsw
# Comment out variable to compile in double precision
miniw_PREC = -DSINGLE_PREC
miniw_OPTS = -O3 -march=haswell -DSIMD_LEN=4
# Comment out variable to compile in double precision
gadget_PREC = USE_SINGLEPRECISION_INTERNALLY
gadget_OPTS = -O3 -march=haswell

.DEFAULT_GOAL := all
.PHONY: all clean $(INSTALL_APPS)

all: $(INSTALL_APPS)

clean:
	rm -rf $(INSTALL_APPS)

qmcpack: qmcpack/.touched
	@echo "Compiling $(@)..."
	cd ../apps/deps && make libxml2/.compiled libfftw/.compiled boost/.compiled hdf5/.compiled
	cd $(@)/build && export FFTW_HOME=$(CURDIR)/../apps/deps/install && export HDF5_ROOT=$(CURDIR)/../apps/deps/install \
	&& export OMPI_CC=$(CC) && export OMPI_CXX=$(CXX) && export LDFLAGS="$(LLTFI_LINKER)" && cmake -DCMAKE_C_COMPILER=$(MPICC) \
	-DCMAKE_CXX_COMPILER=$(MPICXX) -DCMAKE_C_FLAGS="$(qmcpack_OPTS) -Xclang -no-opaque-pointers $(LLTFI_SETTINGS)" $(qmcpack_PREC) \
	-DCMAKE_CXX_FLAGS="$(qmcpack_OPTS) -Xclang -no-opaque-pointers $(LLTFI_SETTINGS)" -DBLAS_openblas_LIBRARY=$(OPENBLAS_ROOT)/libopenblas.so \
	-DLIBXML2_INCLUDE_DIR=$(CURDIR)/../apps/deps/install/include/libxml2 -DLIBXML2_LIBRARY=$(CURDIR)/../apps/deps/install/lib/libxml2.so \
	-DBoost_INCLUDE_DIR=$(CURDIR)/../apps/deps/install/include/ -DENABLE_OMP_TASKLOOP=OFF ../ && make -j $(MAKE_JOBS)

gromacs: gromacs/.touched
	@echo "Compiling $(@)..."
	cd ../apps/deps && make libfftw/.compiled
	mkdir -p $(@)/build && export OMPI_CC=$(CC) && export OMPI_CXX=$(CXX) && export LDFLAGS="$(LLTFI_LINKER)" \
	&& cd $(@)/build && cmake $(gromacs_OPTS) $(gromacs_PREC) -DCMAKE_C_FLAGS="-Xclang -no-opaque-pointers $(LLTFI_SETTINGS)" \
	-DCMAKE_CXX_FLAGS="-Xclang -no-opaque-pointers $(LLTFI_SETTINGS)" -DCMAKE_PREFIX_PATH=$(CURDIR)/../apps/deps/install \
	-DCMAKE_C_COMPILER=$(MPICC) -DCMAKE_CXX_COMPILER=$(MPICXX) -DGMX_FFT_LIBRARY=fftpack ../ && make -j $(MAKE_JOBS)

seissol: seissol/.touched
	@echo "Compiling $(@)..."
	cd ../apps/deps && make hdf5/.compiled netcdf/.compiled libxsmm/.compiled parmetis/.compiled easi/.compiled eigen/.compiled
	mkdir -p $(@)/build && cd $(@)/build && export PATH=$(CURDIR)/../apps/deps/install/bin:$(PATH) && export OMPI_CC=$(CC) \
	&& export OMPI_CXX=$(CXX) && export LDFLAGS="$(LLTFI_LINKER)" && cmake -DCMAKE_C_COMPILER=$(MPICC) -DCMAKE_CXX_COMPILER=$(MPICXX) \
	-DCOMMTHREAD=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="$(LLTFI_SETTINGS)" -DCMAKE_CXX_FLAGS="$(LLTFI_SETTINGS)" \
	-DCMAKE_PREFIX_PATH=$(CURDIR)/../apps/deps/install -DORDER=4 $(seissol_OPTS) $(seissol_PREC) ../ && make -j $(MAKE_JOBS)

milc: milc/.touched
	@echo "Compiling $(@)..."
	cd ../apps/deps && make libfftw/.compiled
	sed -i -e "s|MPI_Errhandler_create|MPI_Comm_create_errhandler|g" -e "s|MPI_Errhandler_set|MPI_Comm_set_errhandler|g" $(@)/generic/com_mpi.c
	export OMPI_CC=$(CC) && export OMPI_CXX=$(CXX) && export $(milc_PREC) && export $(milc_OPTS) && export MPP=true && export MY_CC=$(MPICC) \
	&& export MY_CXX=$(MPICXX) && export OCFLAGS="$(LLTFI_SETTINGS)" && export OCXXFLAGS="$(LLTFI_SETTINGS)" && export LDFLAGS="$(LLTFI_LINKER)" \
	&& cd $(@)/clover_dynamical && cp ../Makefile ./ && CC=$(MPICC) CXX=$(MPICXX) make su3_rmd su3_hmc \
	&& cd ../ks_imp_rhmc && cp ../Makefile ./ && CC=$(MPICC) CXX=$(MPICXX) make su3_rhmc_hisq su3_rhmd_hisq

miniw: miniw/.touched
	@echo "Compiling $(@)..."
	cd ../apps/deps && make pnetcdf/.compiled
	cd $(@)/cpp && rm -rf YAKL && git clone $(yakl_REPO)
	cd $(@)/cpp/build && export LDFLAGS="$(LLTFI_LINKER)" && export OMPI_CC=$(CC) && export OMPI_CXX=$(CXX) \
	&& cmake -DYAKL_CXX_FLAGS="-I$(CURDIR)/../apps/deps/install/include $(miniw_OPTS) $(miniw_PREC) $(LLTFI_SETTINGS) -DNO_INFORM" \
	-DCMAKE_CXX_COMPILER=$(MPICXX) -DCMAKE_C_COMPILER=$(MPICC) -DLDFLAGS="-L$(CURDIR)/../apps/deps/install/lib -lpnetcdf" \
	-DNX=300 -DNZ=150 -DSIM_TIME=1000 -DOUT_FREQ=100 -DDATA_SPEC=DATA_SPEC_THERMAL .. && make -j $(MAKE_JOBS) && mv mpi mpi_thermal \
	&& cmake -DNX=200 -DNZ=100 -DSIM_TIME=1500 -DOUT_FREQ=150 -DDATA_SPEC=DATA_SPEC_GRAVITY_WAVES .. && make -j $(MAKE_JOBS) && mv mpi mpi_mountain

gadget: gadget/.touched
	@echo "Compiling $(@)..."
	cd ../apps/deps && make gsl/.compiled hdf5/.compiled libfftw/.compiled
	echo $(gadget_PREC) >> $(@)/examples/G2-gassphere/Config.sh
	echo $(gadget_PREC) >> $(@)/examples/G2-galaxy/Config.sh
	printf "GSL_INCL = -I$(CURDIR)/../apps/deps/install/include\nGSL_LIBS = -L$(CURDIR)/../apps/deps/install/lib\n" > $(@)/buildsystem/Makefile.gen.libs
	printf "FFTW_INCL = \$$(GSL_INCL)\nFFTW_LIBS = \$$(GSL_LIBS)\n" >> $(@)/buildsystem/Makefile.gen.libs
	printf "HDF5_INCL = \$$(GSL_INCL)\nHDF5_LIBS = \$$(GSL_LIBS)\n" >> $(@)/buildsystem/Makefile.gen.libs
	printf "CPP = $(MPICXX) -std=c++11\nOPTIMIZE = $(LLTFI_SETTINGS) $(gadget_OPTS)\n" > $(@)/buildsystem/Makefile.comp.gcc
	printf "LINKER = $(MPICXX) $(LLTFI_LINKER)" >> $(@)/buildsystem/Makefile.comp.gcc
	cd $(@) && export OMPI_CC=$(CC) && export OMPI_CXX=$(CXX) && export SYSTYPE="Generic-gcc" && make DIR=examples/G2-galaxy -j $(MAKE_JOBS) \
	&& make DIR=examples/G2-gassphere -j $(MAKE_JOBS)
