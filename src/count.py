#! /usr/bin/env python
from __future__ import print_function
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import pandas as pd
import sys
import argparse
import os
import glob
from scipy.stats import zscore
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pysam

#arg parsing
######################################################################################################################
parser = argparse.ArgumentParser(description="finds potential CNVs using a directory of depth bedfiles from MosDepth")
parser.add_argument("cram_dir")
parser.add_argument("-b", "--bed", 
    help="bedfile of NAHR-prone regions",
    required=True
)
parser.add_argument("-o", "--out_dir", 
    help="directory for output",
    required=True
)

args = parser.parse_args()



#####################################################################################################################

#functions and generators
######################################################################################################################

######################################################################################################################

#main block
######################################################################################################################
regions = []
with open(args.bed, 'r') as regions_file:
    for line in regions_file:
        regions.append(line.split()[:3])

data = pd.DataFrame()
#TODO integrate this into the process for streamlined counting of reads in the samples/regions of interest
for cram in glob.glob(os.path.join(args.cram_dir,"*.cram")):
    sample = os.path.splitext(cram)[0]
    samfile = pysam.AlignmentFile(cram, "rc")
    region_counts = []


    for region in regions:
        iter = samfile.fetch(region[0], int(region[1]), int(region[2]))
        region_counts.append(sum(1 for x in iter if x.mapping_quality > 0 ))
    data[sample] = pd.Series(np.array(region_counts), index=regions)
data.to_csv("readcounts.csv")
