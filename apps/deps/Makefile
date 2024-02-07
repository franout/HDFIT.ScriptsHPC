include $(SDENV_APPLICATIONS_DIR)/hpc/config.mk

## franout: most of the libs can be compiled with clang15, 
## 			however, hdf5 and zlig in someway raise errors
##			workaround is to  ,since CC and CXX are env based 
hypre_REPO = https://github.com/hypre-space/hypre hypre && cd hypre && git checkout 14ee602fbf4655715c904fa375165a4565b735bb
mfem_REPO = https://github.com/mfem/mfem.git mfem && cd mfem && git checkout e6c39635c9d21189ff7a9984aa138662059e9b4b
gklib_REPO = https://github.com/KarypisLab/GKlib gklib && cd gklib && git checkout a7f8172703cf6e999dd0710eb279bba513da4fec
metis_REPO = https://github.com/KarypisLab/METIS metis && cd metis && git checkout 94c03a6e2d1860128c2d0675cbbb86ad4f261256
libxsmm_REPO = --branch 1.17 https://github.com/hfp/libxsmm libxsmm
yaml-cpp_REPO = --branch yaml-cpp-0.7.0 https://github.com/jbeder/yaml-cpp yaml-cpp
easi_REPO = --branch v1.1.2 https://github.com/SeisSol/easi easi
impalajit_REPO = https://github.com/uphoffc/ImpalaJIT impalajit && cd impalajit && git checkout 0e4178b7245357a394717cb29f9183b93bafd62f

libxml2_TARBALL = https://gitlab.gnome.org/GNOME/libxml2/-/archive/v2.9.14/libxml2-v2.9.14.tar.gz
libfftw_TARBALL = http://fftw.org/fftw-3.3.10.tar.gz
boost_TARBALL = https://boostorg.jfrog.io/artifactory/main/release/1.79.0/source/boost_1_79_0.tar.gz
hdf5_TARBALL = https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.9/src/hdf5-1.10.9.tar.gz
zlib_TARBALL = https://zlib.net/fossils/zlib-1.2.12.tar.gz
netcdf_TARBALL = https://syncandshare.lrz.de/dl/fiJNAokgbe2vNU66Ru17DAjT/netcdf-4.6.1.tar.gz
pnetcdf_TARBALL = https://parallel-netcdf.github.io/Release/pnetcdf-1.12.3.tar.gz
parmetis_TARBALL = http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
eigen_TARBALL = https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
gsl_TARBALL = https://mirror.ibcp.fr/pub/gnu/gsl/gsl-2.7.1.tar.gz

.DEFAULT_GOAL := libfftw/.compiled
.PHONY: clean

clean:
	rm -rf metis hypre mfem gklib eigen parmetis yaml-cpp impalajit
	rm -rf install libxml2 libfftw boost hdf5 zlib easi libxsmm netcdf pnetcdf

gsl/.compiled: gsl/.downloaded
	@echo "Compiling $(@D)..."
	cd $(@D)/gsl-2.7.1 && ./configure --prefix=$(CURDIR)/install && make -j $(MAKE_JOBS) && make install
	@touch $(@D)/.compiled

eigen/.compiled: eigen/.downloaded
	@echo "Compiling $(@D)..."
	mkdir -p $(@D)/eigen-3.4.0/build && cd $(@D)/eigen-3.4.0/build && cmake -DCMAKE_INSTALL_PREFIX=$(CURDIR)/install ../ \
	&& make install
	@touch $(@D)/.compiled

parmetis/.compiled: parmetis/.downloaded
	@echo "Compiling $(@D)..."
	sed -i 's|#define IDXTYPEWIDTH 32|#define IDXTYPEWIDTH 64|g' $(@D)/parmetis-4.0.3/metis/include/metis.h
	cd $(@D)/parmetis-4.0.3 && make config cc=mpicc cxx=mpicxx prefix=$(CURDIR)/install && make -i install
	cd $(@D)/parmetis-4.0.3 && cp build/Linux-x86_64/libmetis/libmetis.a $(CURDIR)/install/lib
	cd $(@D)/parmetis-4.0.3 && cp metis/include/metis.h $(CURDIR)/install/include
	@touch $(@D)/.compiled

yaml-cpp/.compiled: yaml-cpp/.cloned
	@echo "Compiling $(@D)..."
	## ned to compile with gcc 
	mkdir -p $(@D)/build && cd $(@D)/build && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX=$(CURDIR)/install \
	-DCMAKE_CXX_FLAGS="-Wno-error=maybe-uninitialized" ../ && make -j $(MAKE_JOBS) && make install
	@touch $(@D)/.compiled

impalajit/.compiled: impalajit/.cloned
	@echo "Compiling $(@D)..."
	mkdir -p $(@D)/build && cd $(@D)/build && cmake -DCMAKE_INSTALL_PREFIX=$(CURDIR)/install ../ \
	&& make -j $(MAKE_JOBS) && make install
	@touch $(@D)/.compiled

easi/.compiled: easi/.cloned impalajit/.compiled yaml-cpp/.compiled
	@echo "Compiling $(@D)..."
	cd $(@D) && cmake -DCMAKE_PREFIX_PATH=$(CURDIR)/install -DCMAKE_INSTALL_PREFIX=$(CURDIR)/install \
	-DASAGI=OFF -DIMPALAJIT=ON -DLUA=OFF ./ && make -j $(MAKE_JOBS) && make install
	@touch $(@D)/.compiled

libxsmm/.compiled: libxsmm/.cloned
	@echo "Compiling $(@D)..."
	cd $(@D) && make generator && cp bin/libxsmm_gemm_generator $(CURDIR)/install/bin/
	@touch $(@D)/.compiled

netcdf/.compiled: hdf5/.compiled netcdf/.downloaded
	@echo "Compiling $(@D)..."
	cd $(@D)/netcdf-4.6.1 && export PATH=$(CURDIR)/install/bin:$(PATH) && CFLAGS="-fPIC $(CFLAGS)" \
	CC=h5pcc ./configure --enable-shared=no --enable-logging --prefix=$(CURDIR)/install --disable-dap \
	&& make -j $(MAKE_JOBS) && make install
	@touch $(@D)/.compiled

pnetcdf/.compiled: pnetcdf/.downloaded
	@echo "Compiling $(@D)..."
	cd $(@D)/pnetcdf-1.12.3 && ./configure --prefix=$(CURDIR)/install && make -j $(MAKE_JOBS) && make install
	@touch $(@D)/.compiled

hypre/.compiled: hypre/.cloned
	@echo "Compiling $(@D)..."
	cd $(@D)/src && ./configure --disable-fortran --with-blas-lib="-L$(OPENBLAS_ROOT) -lopenblas" --with-lapack-lib="-L$(OPENBLAS_ROOT) -lopenblas"
	cd $(@D)/src && make -j $(MAKE_JOBS)
	@touch $(@D)/.compiled

gklib/.compiled: gklib/.cloned
	@echo "Compiling $(@D)..."
	cd $(@D) && make config prefix=$(CURDIR)/$(@D)/install && make -j $(MAKE_JOBS) && make install
	@touch $(@D)/.compiled

metis/.compiled: gklib/.compiled metis/.cloned
	@echo "Compiling $(@D)..."
	cd $(@D) && make config gklib_path=$(CURDIR)/gklib/install && make all -i -j $(MAKE_JOBS)
	@touch $(@D)/.compiled

mfem/.compiled: metis/.compiled hypre/.compiled mfem/.cloned
	@echo "Compiling $(@D)..."
	cd $(@D) && make config MFEM_USE_LAPACK=YES MFEM_USE_MPI=YES MFEM_USE_METIS=YES MFEM_USE_METIS_5=YES LAPACK_LIB="-L$(OPENBLAS_ROOT) -lopenblas" \
	METIS_LIB="-L@MFEM_DIR@/../metis/build/libmetis -lmetis -L@MFEM_DIR@/../gklib/install/lib -l:libGKlib.a" \
	METIS_OPT="-I@MFEM_DIR@/../metis/build/xinclude"
	cd $(@D) && make all -j $(MAKE_JOBS)
	@touch $(@D)/.compiled

libxml2/.compiled: libxml2/.downloaded
	@echo "Compiling $(@D)..."
	cd $(@D)/libxml2-v2.9.14 && ./autogen.sh && ./configure --without-python
	cd $(@D)/libxml2-v2.9.14 && make -j $(MAKE_JOBS) &&  make install prefix=$(CURDIR)/install
	@touch $(@D)/.compiled

libfftw/.compiled: libfftw/.downloaded
	@echo "Compiling $(@D)..."
	# Compiling both SP and DP versions of FFTW
	cd $(@D)/fftw-3.3.10 && ./configure --enable-threads --enable-shared --enable-openmp --enable-mpi --prefix=$(CURDIR)/install
	cd $(@D)/fftw-3.3.10 && make -j $(MAKE_JOBS) && make install && make clean
	cd $(@D)/fftw-3.3.10 && ./configure --enable-single --enable-threads --enable-shared --enable-openmp --enable-mpi --prefix=$(CURDIR)/install
	cd $(@D)/fftw-3.3.10 && make -j $(MAKE_JOBS) && make install
	@touch $(@D)/.compiled

boost/.compiled: boost/.downloaded
	@echo "Compiling $(@D)..."
	cd $(@D)/boost_1_79_0 && ./bootstrap.sh --without-libraries=python
	cd $(@D)/boost_1_79_0 && ./b2 -j $(MAKE_JOBS) && ./b2 install --prefix=$(CURDIR)/install
	@touch $(@D)/.compiled

hdf5/.compiled: hdf5/.downloaded
	# broken for clang 15 - needed gcc 
	@echo "Compiling $(@D)..."
	cd $(@D)/hdf5-1.10.9 && ./configure --enable-parallel --prefix=$(CURDIR)/install
	cd $(@D)/hdf5-1.10.9 && make -j $(MAKE_JOBS) && make install
	@touch $(@D)/.compiled

zlib/.compiled: zlib/.downloaded
	# broken for clang 15 - needed gcc
	@echo "Compiling $(@D)..."
	cd $(@D)/zlib-1.2.12 && ./configure
	cd $(@D)/zlib-1.2.12 && make
	@touch $(@D)/.compiled
