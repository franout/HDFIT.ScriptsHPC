# HDFIT HPC Scripting Infrastructure

The scripts contained in this directory allow to perform and evaluate large-scale fault injection experiments for HPC applications instrumented with HDFIT.

## Running Experiments

The main script to be used in order to launch experiments is __HDFIT\_runner.sh__. This script expects as input argument the path to an __.env__ file 
with environment variables describing the application to be tested, with the associated configuration. Each of the instrumented HDFIT applications 
comes with two distinct testable inputs. the __HDFIT\_test-template.env__ file provides a template for building new test scenarios. 
In the following is an usage example:

```
./HDFIT_runner.sh ../apps/qe/QE-test-AlAs.env &
```

Before running any tests, make sure that HDFIT has been properly configured by editing the __config.mk__ file in the __HDFIT.ScriptsHPC__ root directory, 
as described in the main [README](../README.md) document.

By default, the __HDFIT\_runner.sh__ script performs 5k application runs under fault injection (i.e., transient), storing all relevant raw output 
data in a dedicated directory (e.g., __out.AlAs__ for the example above). Please keep in mind that application output to the standard output and error
channels is also directed to a separate log file for each run, and will not be visible while running an experiment.

The __HDFIT\_runner.sh__ script also computes summary CSV files. Computing the CSV files requires __Python 3__ to be installed with the
__numpy__ package, and is done through the internal __HDFIT\_computeCSV.py__ and __HDFIT\_parsers.py__ scripts. In order to perform experiments
with the __MiniWeather__ application, the __netCDF4__ Python package is additionally required.
The runner's overall behavior can be changed at will by simply editing the script itself. It should be noted that the __HDFIT\_runner\_internal.sh__
script contains most of the actual logic to perform experiments, and it is not meant to be executed directly by users.

The scripts were designed to run on a machine equipped with an 18-core CPU, and hence execute 18 single-core application runs concurrently, 
each mapped to a separate CPU core. This behavior can be tuned with the __FI\_PARRUNS__ variable in each __.env__ file. On top of performing 
fault injection, the scripts will also execute the application once without any corruption, so as to obtain a reference output as well as 
the estimated number of operations (e.g., GEMM or FPU) computed in a single run. The output of this process can be found in the __golden__ 
directory for the specified configuration (e.g., within __out.AlAs__).

The resulting summary CSV files can be analyzed using the Python scripts under the __plot__ directory.

## Additional Functionality

The __HDFIT\_runner.sh__ script can also be used to perform experiments in a distributed HPC cluster environment, on an arbitrary number of nodes. 
In this case, please set __FI\_PARRUNS__ to 1. Additionally, __HDFIT\_runner.sh__ can also be used to perform experiments under permanent faults 
(i.e., stuck-high and stuck-low) if the underlying fault injection model allows for it. To enable this functionality, set the __BLASFI\_MODE__ environment variable within __HDFIT\_runner\_internal.sh__ to __PERMANENT__, and __BLASFI\_CORRUPTION__ to either __STUCKHIGH__ or __STUCKLOW__.

