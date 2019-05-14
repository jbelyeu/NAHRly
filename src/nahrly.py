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
import gzip as gz
from scipy.stats import zscore
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import vcfwriter

#arg parsing
######################################################################################################################
parser = argparse.ArgumentParser(description="finds potential CNVs using a directory of depth bedfiles from MosDepth")
parser.add_argument("bed_dir")
parser.add_argument("-v", "--vcf", help="filename for the ouput in VCF format", required=True)
args = parser.parse_args()
#####################################################################################################################

#functions and generators
######################################################################################################################
def readFile(bedfile):
    with gz.open(bedfile, 'r') as open_bed:
        i = 0
        while True:
            line = open_bed.readline()
            if not line:
                break
            yield float(line.split()[3])

def ext_sname(filename):
    return os.path.basename(filename).split(".")[0]

def depth2CN(region_info):
    region_info["CN"] = np.round(region_info['DP']).astype(int)


    return region_info

######################################################################################################################

#main block
######################################################################################################################
bedfiles = glob.glob(os.path.join(args.bed_dir,"*.regions.bed.gz"))
samples = [ext_sname(x) for x in bedfiles]

regions_bed = bedfiles[0]
regions = []
with gz.open(regions_bed, 'rt') as regions_file:
    for line in regions_file:
        #this shouldn't be necessary, but for now the region pairs suck
        start = int(line.strip().split()[1])
        end = int(line.strip().split()[2])
        length = end-start
        if length < 0:
            temp = start
            start = end
            end = temp
        length = end-start
        if length >= 50:
            regions.append(str("_".join(line.split()[:3])))
data = []
#use generators to read one depth from each file at a time
readers = [readFile(bed) for bed in bedfiles]
#insert the depths into a matrix with samples as cols and regions as rows
for region in regions:
    data.append([next(reader) for reader in readers])
#make that matrix into a pandas dataframe for analysis
depths_matrix = pd.DataFrame(data,index=regions,columns=samples)

#remove regions where all the sample depth is zero
depths_matrix = depths_matrix.loc[(depths_matrix > 0).any(axis='columns')]

#normalize internally by median * 2 (assuming most common value is CN=2)
normalized_depths = depths_matrix
normalized_depths = normalized_depths / normalized_depths.median(axis='rows')
normalized_depths = normalized_depths * 2
normalized_depths.to_csv("internal_norm_depths2.csv")

normalized_depths = normalized_depths


cy_writer = vcfwriter.get_writer(args.vcf,normalized_depths.columns)
for region, row in normalized_depths.iterrows():
    chrom,start,stop = region.split("_")
    start = int(start)
    stop = int(stop)
    region_info = {
        "chrom": chrom,
        "start": start,
        "stop": stop,
        "DP": row.values
    }


    region_info = depth2CN(region_info)
    vcfwriter.write_variant(cy_writer, region_info)
