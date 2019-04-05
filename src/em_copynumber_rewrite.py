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
    probabilities = gmm.predict_proba(depths)
    return predictions, probabilities


##########################################################################################################################

#main block
##########################################################################################################################
regions = pd.read_csv(args.depths_matrix, index_col=0)
for region, row in regions.iterrows():
    region = region[2:-1]
    if region != "22_24352143_24386421":
        continue

    preds, probs = EMCopyNumber(row)

    print("SAMPLE", "QUAL", "CN", "DEPTH", "PROBS", sep="\t")
    for i,pred in enumerate(preds):
        print(row.index[i], 
                "PASS" if (max(probs[i]) > 0.9) else "FAIL", 
                pred, 
                row.values[i], 
                ",".join(str(x) for x in probs[i]),
                sep="\t"
            )

