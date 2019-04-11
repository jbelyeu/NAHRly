#! /usr/bin/env python
from __future__ import print_function
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import numpy as np
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
import signal
import time 
from multiprocessing import Pool

#arg parsing
######################################################################################################################
parser = argparse.ArgumentParser(description="finds potential CNVs using a directory of depth bedfiles from MosDepth")
parser.add_argument("cram_dir")
parser.add_argument("-b", "--bed", 
    help="bedfile of NAHR-prone regions",
    required=True
)
parser.add_argument("-c", "--count_file",
    help="filename for output",
    required=True
)
parser.add_argument("-P", "--process_count",
    help="number of processes to use",
    required=False,
    type=int,
    default=1
)

args = parser.parse_args()

#####################################################################################################################
#functions and generators
######################################################################################################################
#this must be global for function access
regions = []
with open(args.bed, 'r') as regions_file:
    for line in regions_file:
        regions.append(line.split()[:3])

def count_sample(cram):
    sample = os.path.split(os.path.splitext(cram)[0])[-1]
    region_counts = []
    for region in regions:
        samfile = pysam.AlignmentFile(cram, "rc")
        iter = samfile.fetch(region[0], int(region[1]), int(region[2]))
        region_counts.append(sum(1 for x in iter if x.mapping_quality > 0 ))
    countseries_pair = sample,pd.Series(np.array(region_counts), index=["_".join(region) for region in regions])
    return countseries_pair
    
######################################################################################################################

#main block
######################################################################################################################
original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
pool=Pool(processes=args.process_count)
signal.signal(signal.SIGINT, original_sigint_handler)
try:
    crams = glob.glob(os.path.join(args.cram_dir,"*.cram"))
    result = pool.map_async(count_sample, (cram for cram in crams))
    while not result.ready():
        time.sleep(1)
except KeyboardInterrupt:
    print("Keyboard interrupt")
    pool.terminate()
counts = result.get()
pool.close()
pool.join()

data = pd.DataFrame()
for sample,count in counts:
    data[sample] = count
data.to_csv(args.count_file)
