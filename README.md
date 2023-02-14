# HDFIT.ScriptsHPC

This repository is part of the Hardware Design Fault Injection Toolkit (HDFIT). HDFIT enables end-to-end fault injection experiments and comprises additionally [HDFIT.NetlistFaultInjector](https://github.com/intel-sandbox/HDFIT.NetlistFaultInjector) and [HDFIT.SystolicArray](https://github.com/intel-sandbox/HDFIT.SystolicArray).

<p align="center" width="100%">
    <img src="HDFIT.png" alt="HDFIT HPC Toolchain" width="80%"/>
</p>

This repository contains the main components of the HDFIT HPC reliability benchmark in order to carry out fault injection experiments on a 
variety of HPC applications, targeting either BLAS GEMM operations (using the proof-of-concept systolic array design implemented in [HDFIT.SystolicArray](https://github.com/intel-sandbox/HDFIT.SystolicArray)) or generic floating-point FPU operations. 

## Directory Structure

The repository is structured in the following directories:

* __apps__: contains code to clone the set of HPC applications supported for BLAS GEMM fault injection, as well as apply patches to them to enable HDFIT. This directory also contains application configurations that can be used to run experiments. Once compiled, the applications can be executed from this location.
* __apps\_llvm__: similar to the above, but focuses on LLVM-based instrumentation for a specific set of HPC applications, targeting generic FPU operations. Application configurations are included, while patching of the source codes is not required.
* __test__: contains scripts to configure and run HDFIT fault injection experiments on the supported HPC applications. The scripts can be used both in a serial context, as well as on distributed HPC clusters for large-scale runs.
* __plot__: contains a series of Python scripts that can be used to process the CSV files produced by HDFIT experiments, in order to generate useful plots and metrics.

For additional details about the components of the HDFIT HPC reliability benchmark, please refer to the README documents in each directory.

## External Dependencies

The main external dependencies of the HPC reliability benchmark are the custom HDFIT OpenBLAS library supporting fault injection, as well as the custom HDFIT version of the LLTFI framework. Before compiling the HPC applications and running experiments, users need to point to the build directories of both, by setting the __OPENBLAS\_ROOT__ and __LLTFI\_ROOT__ variables in the __config.mk__ file. Setting __LLTFI\_ROOT__ is required only for compiling and using the __apps\_llvm__ part of the reliability benchmark - for additional details please refer to the __apps\_llvm__ [README](apps_llvm/README.md) document.

There are other dependencies required to compile the HPC applications and use the Python plotting scripts. These are __make__, __cmake__, __autoconf__, __pkgconf__, __MPI__ and a functional __gcc__ and __gfortran__ toolchain for the former, plus __Python 3__ with the __numpy__, __matplotlib__ and __seaborn__ packages for the latter. The __netCDF4__ Python package is optionally required to perform experiments with the __MiniWeather__ HPC application. In order to compile and use the applications in the __apps\_llvm__ directory, a functional __LLVM__ toolchain (version 15 or above) is additionally required. More details can be found in the README documents in the __apps__, __apps\_llvm__ and __plot__ directories respectively.

## Getting Started

The basic process to run and analyze HDFIT fault injection tests comprises the following steps, considering as an example the GEMM-focused part of the reliability benchmark (__apps__ directory) - the same steps can be applied for LLVM-based FPU fault injection (__apps\_llvm__ directory).
First, the HPC applications need to be compiled:

```
cd apps && make all
```

Then, a test can be run - here we consider as example the CP2K application with the C2H4 input, performing by default 5k fault injection runs:

```
cd cp2k && ../../test/HDFIT_runner.sh CP2K-test-C2H4.env
```

This will eventually produce a __out.C2H4__ directory containing the experiment's results and a CSV summary. It should be noted that the output of each
application run is not printed on the shell, but is directed to separate log files (e.g., __out.C2H4/fi-transient/run10.log__). The CSV summary file
can be further fed into the HDFIT plotting scripts, for example to produce an SDE error curve:

```
cd out.C2H4 && python3 ../../../plot/HDFIT_plot_error_curve.py HDFIT-CP2K-C2H4-29.08.2022-transient.csv
```

This will produce an image file containing the desired plot, as well as display several statistical metrics. Further analysis can be conducted by using the output files resulting from each application run under fault injection.

## License Terms

All original code that is part of the HDFIT HPC reliability benchmark is released under the terms of the GNU Lesser General Public License (LGPL)
version 3 or (at your option) any later version. This includes all files in the __plot__ and __test__ directories of this repository.

The patch files for the individual HPC applications, as well as the associated input configurations, are instead released under the terms of the
respective original licenses. This includes all files under the __apps/resources__ and __apps\_llvm/resources__ directories. A copy of each application's license is included.
