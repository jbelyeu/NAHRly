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


#argparse
##########################################################################################################################
parser = argparse.ArgumentParser(description="finds potential CNVs using a directory of depth bedfiles from MosDepth")
parser.add_argument("-d", "--depths_matrix", 
    help="CSV file with depths by region/sample from gnarly",
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
def EMCopyNumber(row):
    depths = row.values.reshape(-1, 1)
    gmm = mixture.GaussianMixture(n_components=5)
    gmm.fit(depths)
    predictions = gmm.predict(depths)
    print(predictions)


##########################################################################################################################

#main block
##########################################################################################################################
regions = pd.read_csv(args.depths_matrix, index_col=0)
for region, row in regions.iterrows():
    region = region[2:-1]
    if region != "22_24352143_24386421":
        continue

    EMCopyNumber(row)
    #binned,centers = EMCopyNumber(row)
#    swarm(binned, centers, args.out_dir, region)
#    print(centers)
#    print(region)
#    for i in range(len(binned)):
#        print (str(i) + ": " + str(binned.iloc[i].count()))
#    print()
#
