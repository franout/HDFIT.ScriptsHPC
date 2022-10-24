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


def main():
    log.basicConfig(format="[ %(levelname)s ] %(message)s", level=log.INFO, stream=sys.stdout)

    # Input parsing
    parser = argparse.ArgumentParser(description="HDFIT NRMSE SDE curve plotting script")
    parser.add_argument("path", action="store", type=str, nargs="+", help="Path[s] to HDFIT CSV file[s]")
    parser.add_argument("-b", action="store", dest="bit_depth", type=HDFIT_intPos, default=sys.maxsize, help="Maximum fault injection bit")
    parser.add_argument("-pos", action="store", dest="bit_pos_list", type=HDFIT_intList, default=[], help="Comma-separated list of protection bit positions")
    parser.add_argument("-bins", action="store", dest="num_bins", type=HDFIT_intPos, default=100, help="Number of points to be sampled")
    parser.add_argument("-top", action="store", dest="top", type=HDFIT_intPos, default=100, help="Y axis upper bound")
    parser.add_argument("-rbound", action="store", dest="right_bound", type=int, default=2, help="Order of magnitude of the right bound")
    parser.add_argument("-thresh", action="store", dest="error_thresh", type=int, default=-9, help="Order of magnitude of log cutoff (symlog only)")
    parser.add_argument("--uoe", action="store_true", dest="exclude_uoe", help="Exclude runs resulting in UOE")
    parser.add_argument("--rtlcomp", action="store_true", dest="compare_rtl", help="Compare results with and without RTL errors. Overrides --rtl and -pos")
    parser.add_argument("--rtl", action="store_false", dest="exclude_rtl", help="Include runs containing RTL errors")
    parser.add_argument("--lin", action="store_false", dest="symlog_scale", help="Use linear scale in place of symlog one")
    args = parser.parse_args()
    
    # Getting parsed arguments
    bit_depth    = args.bit_depth
    bit_pos_list = args.bit_pos_list
    num_bins     = args.num_bins
    yTop         = args.top
    error_thresh = args.error_thresh
    right_bound  = args.right_bound
    exclude_uoe  = args.exclude_uoe
    exclude_rtl  = args.exclude_rtl
    compare_rtl  = args.compare_rtl
    symlog_scale = args.symlog_scale
    fault_log    = args.path

    # Additional checks on user input
    if len(bit_pos_list) > 0 and bit_depth == sys.maxsize:
        raise ValueError("Cannot use the -pos argument if -b is not speficied as well.")
    elif any(el > bit_depth for el in bit_pos_list):
        raise ValueError("One or more protection bit positions are higher than precision %d." % bit_depth)
    elif error_thresh > right_bound:
        raise ValueError("Error threshold %d is higher than right bound %d." % (error_thresh, right_bound))

    xLabel       = "NRMSE Error Threshold [%]"
    yLabel       = FI_METRIC + " [%]"
    fontsize     = 14
    figsize      = [7, 5]
    linestyles   = ["-", ":", "-.", "--"]
    markerstyles = ["o", "^", "s", "P", "D"]

    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    fig, ax = plt.subplots(figsize=figsize)
    
    # Constructing bin edges in symlog or linear scale
    if symlog_scale:
        binsLin = np.linspace(0, 10**error_thresh, int(num_bins / 10), endpoint=False)
        binsLog = np.logspace(error_thresh, right_bound, num_bins)
        bins = np.concatenate((binsLin, binsLog))
    else:
        right_bound = 10**right_bound
        bins = np.linspace(0, right_bound, num_bins)

    sIdx = 0
    # Setting bit_pos to 64 means no protection
    bit_pos_list = [bit_depth] + bit_pos_list if not compare_rtl else [False, True]
    for entry in args.path:
        log.info("- File %s" % entry)
        baseName = getFileBaseName(entry)
        for bit_pos in bit_pos_list:
            # Loading data
            dataDict  = loadData(entry)
            # Processing data
            if type(bit_pos) is not bool:
                vals = getDataPartialProtection(dataDict, bit_pos, bins, exclude_uoe, exclude_rtl)
                myLabel = baseName if bit_pos==bit_depth else f"protection [{bit_depth-1}:{bit_pos}]"
                printFailureStatistics(dataDict, bit_pos, bit_depth)
            else:
                vals = getDataPartialProtection(dataDict, bit_depth, bins, exclude_uoe, bit_pos)
                myLabel = baseName if not bit_pos else "hw protection"
                if not bit_pos:
                    printFailureStatistics(dataDict, bit_depth, bit_depth)
                    
            lS = linestyles[sIdx % len(linestyles)]
            mS = markerstyles[sIdx % len(markerstyles)]
            mEvery = 10 + sIdx
            sIdx = sIdx + 1

            ax.plot(bins, vals * 100, linewidth=3, zorder=2, label=myLabel, linestyle=lS, marker=mS, markevery=mEvery, markersize=7)

    ax.set_xlabel(xLabel, fontsize=fontsize)
    ax.set_ylabel(yLabel, fontsize=fontsize)
    ax.set_ylim(bottom=-0.2, top=yTop)
    ax.legend(fontsize=fontsize-2, loc=1)
    ax.grid(zorder=0)
    
    if symlog_scale:
        ax.set_xscale("symlog", linthresh=10**error_thresh)
        ax.set_xlim((bins[0] - 0.1*10**error_thresh, bins[-1]))
        plt.gca().xaxis.set_minor_locator(MinorSymLogLocator(10**error_thresh))
        for l in plt.gca().xaxis.get_ticklabels()[1::2]:
            l.set_visible(False)
    else:
        ax.set_xlim((bins[0] - right_bound*0.01, bins[-1] + right_bound*0.01))
        ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, 3))
    
    error_histogram_file = f"{args.path[0][0:-4]}_error_curve.png" if not compare_rtl else f"{args.path[0][0:-4]}_error_curve_rtlcomp.png"
    plt.savefig(error_histogram_file, bbox_inches='tight')
    log.info(f"Error curve saved to {error_histogram_file}")
 
if __name__ == '__main__':
    sys.exit(main() or 0)
