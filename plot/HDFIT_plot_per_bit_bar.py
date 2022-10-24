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
    parser = argparse.ArgumentParser(description="HDFIT per-bit bar plotting script")
    parser.add_argument("path", action="store", type=str, help="Path to HDFIT CSV file")
    parser.add_argument("-b", action="store", dest="bit_depth", type=HDFIT_intPos, default=64, help="Maximum fault injection bit")
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
    failure   = dataDict[HDFIT.testFail]
    fault     = dataDict[HDFIT.fiBit]
    errors    = dataDict[HDFIT.nrmse]
    
    max_row = len(fault)
    log.info(f"{max_row} faults loaded from the fault log file")

    # Get the bit position distribution of the fault
    bit_list = list(range(bit_depth))
    fault_per_bit = [ (fault == bit).sum() for bit in bit_list ]
    plt.figure(figsize=(20,6))
    p = plt.bar(bit_list, fault_per_bit, align='edge', width=0.8)
    plt.xticks(np.arange(bit_depth-1, -1, -1))
    plt.gca().invert_xaxis()    
    plt.xlabel("Bit Position")
    plt.ylabel("Number of Faults")        
    hist_file = f"{fault_log[0:-4]}_faults_per_bit.png"
    plt.savefig(hist_file, bbox_inches='tight')
    log.info(f"Faults per bit bar plot saved to {hist_file}")
    
    # Filtering out runs
    ex_runs = filterData(dataDict, bit_depth, exclude_uoe, exclude_rtl)
    # Combining hard program failures with SDEs
    agg_failure = np.asarray([(failure[idx]>0 or errors[idx]>0) and not ex_runs[idx] for idx in range(len(failure))])
    failure_per_bit = np.zeros(bit_depth)
    for id in bit_list:
        mask = (fault == id)
        failure_per_bit[id] = np.sum(agg_failure, where=mask)/fault_per_bit[id] if fault_per_bit[id]!=0 else 0
		        
    failure_per_bit_file = f"{fault_log[0:-4]}_failure_rate_per_bit.png"
    plt.figure(figsize=(20,6))
    p = plt.bar(bit_list, failure_per_bit * 100, align='edge', width=0.8)
    plt.xticks(np.arange(bit_depth-1, -1, -1))
    plt.bar_label(p, fmt='%.0f', fontsize=6)    
    plt.gca().invert_xaxis()    
    plt.xlabel("Bit Position")
    plt.ylabel(FI_METRIC + " [%]")
    plt.ylim(0, 105)
    plt.savefig(failure_per_bit_file, bbox_inches='tight')
    log.info(f"Failures per bit bar plot saved to {failure_per_bit_file}")   
 
if __name__ == '__main__':
    sys.exit(main() or 0)
