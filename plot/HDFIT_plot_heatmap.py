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
    parser = argparse.ArgumentParser(description="HDFIT failure heatmap plotting script")
    parser.add_argument("path", action="store", type=str, help="Path to HDFIT CSV file")
    parser.add_argument("-b", action="store", dest="bit_depth", type=HDFIT_intPos, default=64, help="Maximum fault injection bit")
    parser.add_argument("-bbins", action="store", dest="bit_bins", type=HDFIT_intPos, default=8, help="Number of bit position bins")
    parser.add_argument("-cbins", action="store", dest="code_bins", type=HDFIT_intPos, default=8, help="Number of code position bins")
    parser.add_argument("--uoe", action="store_true", dest="exclude_uoe", help="Exclude runs resulting in UOE")
    parser.add_argument("--rtl", action="store_false", dest="exclude_rtl", help="Include runs containing RTL errors")
    args = parser.parse_args()
    
    # Getting parsed arguments
    bit_depth    = args.bit_depth
    exclude_uoe  = args.exclude_uoe
    exclude_rtl  = args.exclude_rtl
    bit_buckets  = args.bit_bins
    code_buckets = args.code_bins
    fault_log    = args.path
   
    # Loading and unpacking data
    dataDict = loadData(fault_log)
    failure  = dataDict[HDFIT.testFail]
    bitPos   = dataDict[HDFIT.fiBit]
    instrPos = dataDict[HDFIT.opFi]
    errors   = dataDict[HDFIT.nrmse]
    maxInstr = np.max(dataDict[HDFIT.opCntExp])

    max_row = len(failure)
    log.info(f"{max_row} faults loaded from the fault log file")

    # Filtering out runs
    ex_runs = filterData(dataDict, bit_depth, exclude_uoe, exclude_rtl)
    # Combining hard program failures with SDEs
    agg_failure = np.asarray([(failure[idx]>0 or errors[idx]>0) and not ex_runs[idx] for idx in range(max_row)])
    
    stepC = np.ceil(maxInstr / code_buckets)
    stepI = np.ceil(bit_depth / bit_buckets)
    cntMat = np.zeros((code_buckets, bit_buckets))
    faultMat = np.zeros((code_buckets, bit_buckets)) + 1

    for idx in range(max_row):
        if bitPos[idx] < bit_depth:
            codeIdx = int(np.floor((instrPos[idx]) / stepC)) if instrPos[idx] < maxInstr else code_buckets - 1
            bitIdx = int(np.floor(bitPos[idx] / stepI))
            faultMat[codeIdx, bitIdx] = faultMat[codeIdx, bitIdx] + 1
            if agg_failure[idx]:
                cntMat[codeIdx, bitIdx] = cntMat[codeIdx, bitIdx] + 1

    fontsize    = 14
    figsize     = [7, 5]

    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    fig, ax = plt.subplots(figsize=figsize)

    ax = sns.heatmap(100 * cntMat / faultMat, linewidths=.5, ax=ax, cbar_kws={'label': FI_METRIC + ' [%]'}, cmap="inferno_r")
    ax.figure.axes[-1].yaxis.label.set_size(fontsize)
    ax.set_xlabel("Bit Position", fontsize=fontsize)
    ax.set_ylabel("Op Position [%]", fontsize=fontsize)

    ax.set_xticklabels([str(int(el * stepI)) for el in ax.get_xticks()])
    ax.set_yticklabels([str(int(100 * el / code_buckets)) for el in ax.get_yticks()])

    error_heatmap_file = f"{fault_log[0:-4]}_failure_rate_heatmap.png"
    plt.savefig(error_heatmap_file, bbox_inches='tight')
    log.info(f"Failure rate heatmap saved to {error_heatmap_file}")         

if __name__ == '__main__':
    sys.exit(main() or 0)
