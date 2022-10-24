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
    parser = argparse.ArgumentParser(description="HDFIT per-bit overhead plotting script")
    parser.add_argument("path", action="store", type=str, help="Path to HDFIT CSV file")
    parser.add_argument("-b", action="store", dest="bit_depth", type=HDFIT_intPos, default=sys.maxsize, help="Maximum fault injection bit")
    parser.add_argument("--uoe", action="store_true", dest="exclude_uoe", help="Exclude runs resulting in UOE")
    parser.add_argument("--rtl", action="store_false", dest="exclude_rtl", help="Include runs containing RTL errors")
    args = parser.parse_args()
    
    # Getting parsed arguments
    bit_depth   = args.bit_depth
    exclude_uoe = args.exclude_uoe
    exclude_rtl = args.exclude_rtl
    fault_log   = args.path

    # Loading and unpacking data
    dataDict  = loadData(fault_log)
    fiBit     = dataDict[HDFIT.fiBit]
    ops       = dataDict[HDFIT.opCnt]
    expOps    = dataDict[HDFIT.opCntExp]

    max_row = len(fiBit)
    log.info(f"{max_row} faults loaded from the fault log file")

    # # Filtering out runs
    ex_runs = filterData(dataDict, bit_depth, exclude_uoe, exclude_rtl)
    ops = np.asarray([ops[idx] for idx in range(max_row) if not ex_runs[idx]]) 
    fiBit = np.asarray([fiBit[idx] for idx in range(max_row) if not ex_runs[idx]])
    medianOps = np.median(expOps)
    ops = (ops - medianOps) * 100 / medianOps

    fontsize    = 12
    figsize     = [20, 6]

    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    fig, ax = plt.subplots(figsize=figsize)

    myPlot = sns.barplot(x=fiBit, y=ops, ax=ax, ci=None, estimator=np.average, zorder=2, palette="viridis")
    myPlot.set_xlabel("Bit Position", fontsize=fontsize)
    myPlot.set_ylabel("Ops Overhead [%]", fontsize=fontsize)
    myPlot.set_ylim(bottom=0)
    myPlot.grid(zorder=0)

    time_barplot_file = f"{fault_log[0:-4]}_time_barplot.png"
    plt.savefig(time_barplot_file, bbox_inches='tight')
    log.info(f"Execution time bar plot saved to {time_barplot_file}")         
 
if __name__ == '__main__':
    sys.exit(main() or 0)
