# Copyright (C) 2022 Intel Corporation
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License, as published
# by the Free Software Foundation; either version 3 of the License,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
#
# SPDX-License-Identifier: LGPL-3.0-or-later

from HDFIT_plot_utils import *
import sys, argparse
import logging as log
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def main():
    log.basicConfig(format="[ %(levelname)s ] %(message)s", level=log.INFO, stream=sys.stdout)

    # Input parsing
    parser = argparse.ArgumentParser(description="HDFIT HW module histogram plotting script")
    parser.add_argument("path", action="store", type=str, help="Path to HDFIT CSV file")
    parser.add_argument("signals", action="store", type=str, help="C++ file with FI signal spec")
    parser.add_argument("-b", action="store", dest="bit_depth", type=HDFIT_intPos, default=sys.maxsize, help="Maximum fault injection bit")
    parser.add_argument("--uoe", action="store_true", dest="exclude_uoe", help="Exclude runs resulting in UOE")
    parser.add_argument("--dmp", action="store_true", dest="dmp_modules", help="Dump module information to a new CSV file")
    parser.add_argument("--rtl", action="store_false", dest="exclude_rtl", help="Include runs containing RTL errors")
    parser.add_argument("--rate", action="store_true", dest="plot_rate", help="Plot failure rate instead of occurrences")
    args = parser.parse_args()
    
    # Getting parsed arguments
    bit_depth    = args.bit_depth
    sig_file     = args.signals
    exclude_uoe  = args.exclude_uoe
    exclude_rtl  = args.exclude_rtl
    dmp_modules  = args.dmp_modules
    plot_rate    = args.plot_rate
    fault_log    = args.path

    # Loading and unpacking data
    dataDict  = loadData(fault_log)
    uuid      = dataDict[HDFIT.assignUUID]
    failure   = dataDict[HDFIT.testFail]
    errors    = dataDict[HDFIT.nrmse]
    modT, mT, sT = loadModuleData(uuid, sig_file)

    # Checking for invalid UUIDs
    if any(u<2 for u in uuid):
        log.warning("One or more invalid UUIDs detected.")

    max_row = len(failure)
    log.info(f"{max_row} faults loaded from the fault log file")

    # Filtering out runs
    ex_runs = filterData(dataDict, bit_depth, exclude_uoe, exclude_rtl)
    # Combining hard program failures with SDEs
    agg_failures = np.asarray([(failure[idx]>0 or errors[idx]>0) and not ex_runs[idx] for idx in range(max_row)])
    
    # Adding dummy failures so that all modules show up in plot
    modules = np.concatenate((modT, mT))
    agg_failures = np.concatenate((agg_failures, np.zeros(mT.shape)))

    # Sorting based on module name
    sortMask     = np.argsort(modules)
    modules      = modules[sortMask]
    agg_failures = agg_failures[sortMask]

    fontsize    = 14
    figsize     = [11, 5]

    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    fig, ax = plt.subplots(figsize=figsize)

    myPlot = sns.barplot(x=modules, y=agg_failures, ax=ax, ci=None, estimator=np.sum if not plot_rate else percent, zorder=2, palette="viridis")
    ax.tick_params(axis='x', rotation=90)
    myPlot.set_xlabel("Module Name", fontsize=fontsize)
    myPlot.set_ylabel("Occurrences" if not plot_rate else FI_METRIC + " [%]", fontsize=fontsize)
    myPlot.set_ylim(bottom=0)
    myPlot.grid(zorder=0)

    module_barplot_file = f"{fault_log[0:-4]}_module_barplot.png"
    plt.savefig(module_barplot_file, bbox_inches='tight')
    log.info(f"Module assignment bar plot saved to {module_barplot_file}")

    if dmp_modules:
        dumpModuleData(fault_log, modT)
        log.info(f"Created a new CSV file containing module data")
 
if __name__ == '__main__':
    sys.exit(main() or 0)
