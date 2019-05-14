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

#arg parsing
######################################################################################################################
parser = argparse.ArgumentParser(description="finds potential CNVs using a directory of depth bedfiles from MosDepth")
parser.add_argument("bed_dir")
parser.add_argument("-z", 
    help="Z-score cutoff (values between -z and z will not be reported)",
    default=4,
    type=float
)
parser.add_argument("-o", "--out_dir", 
    help="directory for output",
    required=True
)
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

def for_peter(df):
    ##here I'm going to output a region for Peter
    #this region should have a mean < 0.8 without external normalization and should be on chromosome 22
    normalized_depths_matrix_22 = normalized_depths_matrix.filter(regex="^22_", axis="rows")
    chr22 = []
    for i,row in normalized_depths_matrix_22.iterrows():
        if row.mean() < 0.75:
            chr22.append(row)
    normalized_depths_matrix_22 = pd.DataFrame(chr22)
    normalized_depths_matrix_22.to_csv("low_depth.csv")

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
depths_matrix.to_csv("depths.csv")

#normalize internally (output for peter), externally
normalized_depths_matrix = depths_matrix
normalized_depths_matrix = normalized_depths_matrix / normalized_depths_matrix.mean(axis='rows')
#for_peter(normalized_depths_matrix)
normalized_depths_matrix.to_csv("internal_norm_depths.csv")
normalized_depths_matrix = normalized_depths_matrix.div(normalized_depths_matrix.median(axis='columns'), axis='rows')
normalized_depths_matrix.to_csv("external_norm_depths.csv")

#find possible deletions and duplications
zscores = normalized_depths_matrix.transpose().apply(zscore, axis=0).transpose()
#possible_dup_regions = normalized_depths_matrix.loc[(zscores > args.z).any(axis=1)]
#possible_del_regions = normalized_depths_matrix.loc[(zscores < -args.z).any(axis=1)]
#
##write out results to a CSV for post-processing
#possible_dup_regions.to_csv("dups.csv")
#possible_del_regions.to_csv("dels.csv")
zscores.to_csv("z.csv")
#
#cluster with kmeans
#from sklearn.cluster import KMeans
#from mpl_toolkits.mplot3d import Axes3D
#import numpy as np
#from Jenks_Breaks import get_jenks_breaks
#
#estimator = KMeans(n_clusters=4)
#pract_region = normalized_depths_matrix.iloc[0].values.reshape(-1,1)
#
#fig = plt.figure('1', figsize=(4, 3))
#ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
#estimator.fit(pract_region)
#labels = estimator.labels_
#
#ax.plot(pract_region, c=labels.astype(np.float), edgecolor='k')
#
#ax.w_xaxis.set_ticklabels([])
#ax.w_yaxis.set_ticklabels([])
#ax.w_zaxis.set_ticklabels([])
#ax.set_title("widget")
#ax.dist = 12
#fig.save("~/public_html/thing.png")
#fig.close()
