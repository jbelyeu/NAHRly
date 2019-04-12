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

def count_sample(counting_arg):
    try:
        cram = counting_arg['cram']
        regions = counting_arg['regions']
        regions_dict = counting_arg['regions_dict']
        sample = os.path.split(os.path.splitext(cram)[0])[-1]
        region_counts = {"_".join(region): 0 for region in regions}
        samfile = pysam.AlignmentFile(cram, "rc")


        for read in samfile:
            if read.is_qcfail or read.is_duplicate or int(read.mapping_quality) <= 0: continue
            chrom = read.reference_name
            
            if chrom in regions_dict:
                pos = min(int(read.reference_start),int(read.reference_end))
                end = max(int(read.reference_start),int(read.reference_end))

                for region_pos, region_end in regions_dict[chrom]:
                    #if the region encompasses any part of the read
                    if (int(region_pos) <= pos <= int(region_end)) or (int(region_pos) <= end <= int(region_end)):
                        region = "_".join([chrom,region_pos,region_end])
                        region_counts[region] += 1
                        
        countseries_pair = [sample,pd.Series(region_counts)]
        return countseries_pair
    except KeyboardInterrupt, e:
        pass
    
######################################################################################################################

#main block
######################################################################################################################
def main():
    regions = []
    regions_dict = {}
    with open(args.bed, 'r') as regions_file:
        for line in regions_file:
            chrom,pos,end = line.strip().split()[:3]

            regions.append([chrom,pos,end])
            if not chrom in regions_dict:
                regions_dict[chrom] = []
            regions_dict[chrom].append([pos,end])

    crams = glob.glob(os.path.join(args.cram_dir,"*.cram"))
    counting_args = [{"regions": regions, "regions_dict": regions_dict, "cram": cram } for cram in crams]
    pool=Pool(processes=args.process_count)
    process = pool.map_async(count_sample, (counting_arg for counting_arg in counting_args))
    
    try:
        counts = process.get(0xFFFF)
    except KeyboardInterrupt:
        print("Keyboard interrupt")
        return
    
    data = pd.DataFrame()
    for sample,count in counts:
        data[sample] = count
    data.to_csv(args.count_file)

if __name__ == "__main__":
    main()
