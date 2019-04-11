#! /usr/bin/env python
from __future__ import print_function
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import sys
import argparse
import os
import pandas as pd
import numpy as np
from scipy.stats import norm
from sklearn import mixture
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

max_k = 10

#argparse
##########################################################################################################################
parser = argparse.ArgumentParser(description="finds potential CNVs using a directory of depth bedfiles from MosDepth")
parser.add_argument("-c", "--counts_matrix", 
    help="CSV file with read counts by region/sample",
    required=True
)
parser.add_argument("-o", "--out_dir", 
    help="directory for output",
    required=True
)
args = parser.parse_args()
##########################################################################################################################
#
#functions 
##########################################################################################################################
def plot_histograms(row, directory, region):
    print(region)
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    ax = row.hist()
    fig = ax.get_figure()
    plotname = os.path.join(directory,region+"hist.png") 
    print(plotname)
    fig.savefig(plotname)

##########################################################################################################################

#main block
##########################################################################################################################
regions = pd.read_csv(args.counts_matrix, index_col=0)
for region, row in regions.iterrows():
    if region[0] == 'b':
        region = region[2:-1]
    #if region != "1_104809_573868":
    #    continue
    #if region != "22_24352143_24386421":
    #    continue
    #if region != "1_87113_398211":
    #    continue

    plot_histograms(row, args.out_dir, region)
    sys.exit()
