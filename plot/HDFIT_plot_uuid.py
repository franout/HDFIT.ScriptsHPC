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
    parser = argparse.ArgumentParser(description="HDFIT UUID histogram plotting script")
    parser.add_argument("path", action="store", type=str, help="Path to HDFIT CSV file")
    parser.add_argument("-b", action="store", dest="bit_depth", type=HDFIT_intPos, default=sys.maxsize, help="Maximum fault injection bit")
    parser.add_argument("-sigs", action="store", dest="signals", type=str, default="", help="Use a C++ file to label UUIDs")
    parser.add_argument("-ubins", action="store", dest="uuid_bins", type=HDFIT_intPos, default=64, help="Number of histogram bins")
    parser.add_argument("--uoe", action="store_true", dest="exclude_uoe", help="Exclude runs resulting in UOE")
    parser.add_argument("--rtl", action="store_false", dest="exclude_rtl", help="Include runs containing RTL errors")
    parser.add_argument("--rate", action="store_true", dest="plot_rate", help="Plot failure rate instead of occurrences")
    args = parser.parse_args()
    
    # Getting parsed arguments
    bit_depth    = args.bit_depth
    exclude_uoe  = args.exclude_uoe
    exclude_rtl  = args.exclude_rtl
    signals      = args.signals
    uuid_buckets = args.uuid_bins
    plot_rate    = args.plot_rate
    fault_log    = args.path

    # Loading and unpacking data
    dataDict  = loadData(fault_log)
    uuid      = dataDict[HDFIT.assignUUID]
    failure   = dataDict[HDFIT.testFail]
    errors    = dataDict[HDFIT.nrmse]

    max_row = len(failure)
    log.info(f"{max_row} faults loaded from the fault log file")

    # Checking for invalid UUIDs
    if any(u<2 for u in uuid):
        log.warning("One or more invalid UUIDs detected.")

    # Rounding each UUID to the closest bucket
    uuidStep = int(np.max(uuid) / uuid_buckets)
    # Three possible UUID cases:
    # 1) UUID is greater than uuidStep -> mapped to a UUID bin
    # 2) UUID is lower than uuidStep -> mapped to UUID value 2 (first valid FI signal)
    # 3) UUID is lower than 2 -> mapped to 0, i.e., invalid
    uuid = np.asarray([int(u / uuidStep) * uuidStep if u>=uuidStep else 2 if u>1 else 0 for u in uuid])

    # Filtering out runs
    ex_runs = filterData(dataDict, bit_depth, exclude_uoe, exclude_rtl)
    # Combining hard program failures with SDEs
    agg_failures = np.asarray([(failure[idx]>0 or errors[idx]>0) and not ex_runs[idx] for idx in range(max_row)])
    
    fontsize    = 12
    figsize     = [20, 6]

    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    fig, ax = plt.subplots(figsize=figsize)

    myPlot = sns.barplot(x=uuid, y=agg_failures, ax=ax, ci=None, estimator=np.sum if not plot_rate else percent, zorder=2, palette="viridis")
    ax.tick_params(axis='x', rotation=90)
    # Mapping signal UUID bins to module names, if requested
    if signals != "":
        tickLabels = [int(v.get_text()) for v in ax.get_xticklabels()]
        modules, mT, sT = loadModuleData(tickLabels, signals)
        ax.set_xticklabels([modules[idx] + " - " + str(tickLabels[idx]) for idx in range(len(tickLabels))])
    myPlot.set_xlabel("Assignment UUID", fontsize=fontsize)
    myPlot.set_ylabel("Occurrences" if not plot_rate else FI_METRIC + " [%]", fontsize=fontsize)
    myPlot.set_ylim(bottom=0)
    myPlot.grid(zorder=0)

    uuid_barplot_file = f"{fault_log[0:-4]}_uuid_barplot.png"
    plt.savefig(uuid_barplot_file, bbox_inches='tight')
    log.info(f"UUID assignment bar plot saved to {uuid_barplot_file}")         
 
if __name__ == '__main__':
    sys.exit(main() or 0)
