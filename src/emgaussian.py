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
def choose_k(bics):
    lowest_bic_idx = -1
    for i,bic in enumerate(bics):
        if (lowest_bic_idx < 0) or bic < bics[lowest_bic_idx]:
            lowest_bic_idx = i
    return lowest_bic_idx

def EMCopyNumber(row):
    depths = row.values.reshape(-1, 1)
    bics = []
    gmms = []
    #create and fit gmms for all k values up to max_k
    for i in range(1,max_k):
        gmm = mixture.GaussianMixture(n_components=i)
        gmm.fit(depths)
        bics.append(float(gmm.bic(depths)))
        gmms.append(gmm)

    #find the number of components that provides the best fit 
    #do this by using the elbow method on the bayesian information-theoretic criterion
    k = choose_k(bics)
    gmm = gmms[k]
   
    #calculate prediction and probabilities from the selected model    
    predictions = gmm.predict(depths)
    probabilities = gmm.predict_proba(depths)
    return predictions, probabilities,k+1

#    #for funzies
#    for i,gmm in enumerate(gmms):
#        predictions = gmm.predict(depths)
#        probabilities = gmm.predict_proba(depths)
#        plot_report(row, predictions, probabilities, i+1, args.out_dir, "22_24352143_24386421")
#    return predictions, probabilities,k+1

def print_report(row, preds, probs, n_clusters):
    print(n_clusters)
    print("SAMPLE", "QUAL", "CN", "DEPTH", "PROBS", sep="\t")
    for i,pred in enumerate(preds):
        print(row.index[i], 
                "PASS" if (max(probs[i]) > 0.9) else "FAIL", 
                pred, 
                row.values[i], 
                ",".join(str(x) for x in probs[i]),
                sep="\t"
            )

def plot_report(data, n_clusters, directory, region):
    print(84)
    if not os.path.exists(directory):
        os.makedirs(directory)
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set(style="whitegrid")
    plt.figure(figsize=(12,9))
    #ax = sns.swarmplot(x=data['cns'],y=data['normdepths'])
    print(data[:10])
    ax = sns.swarmplot(y="normdepths", x=[""]*len(data), hue="cns", data=data)
    

    ax.set_title(region + ": with k=" + str(n_clusters))
    plotname = os.path.join(directory,region+"_"+str(n_clusters)+"_skl.pdf")
    plt.ylabel("Normalized depth")
    plt.xlabel("Predicted copy number")
    plt.savefig(plotname)
    plt.close()

def predict_cn(data, region):
    print(103)
    #data is a dataframe of normalized depth (normdepths), prediction groups (preds), and sample labels (index) for a single region
    #I want to modify it to add a new column with predicted copy number based on the information already in the dataframe
    #to do this, I'll find the largest group and set it as CN=2, then guess the others based on that
    counts = data.groupby("preds")["normdepths"].size()

    largest_cluster_idx = 0
    for i,count in enumerate(counts):
        if count > counts[largest_cluster_idx]:
            largest_cluster_idx = i
    
    #CN=2 is the most common
    copy_numbers = {largest_cluster_idx: 2}

    #the number of copy numbers present is determined by the largest median
    medians = data.groupby("preds")["normdepths"].median()
    largest_median_idx = 0
    for i,median in enumerate(medians):
        if median > medians[largest_median_idx]:
            largest_median_idx = i
    #the maximum copy number in the data for this site is the largest median 
    #divided by the median of the largest cluster (by count of samples in it)
    #I'm also adding a small amount in case the largest cluster is actually at zero which causes a divide by zero error
    max_cn = (2*(round((medians[largest_median_idx] / (medians[largest_cluster_idx]+0.01) ))))
    copy_numbers[largest_median_idx] = max_cn

    #create a list of predicted depths that correspond to the possible copy numbers
    pred_depths = []
    for i in range(int(max_cn) +1):
        pred_depths.append(medians[largest_median_idx] * (i/max_cn))
    
    copy_numbers = {}
    for cn,pred_depth in enumerate(pred_depths):
        min_dist = np.inf
        depth_bestmatch_idx = -1

        for i,median in enumerate(medians):
            if abs(median-pred_depth) < min_dist:
                min_dist = abs(median-pred_depth)
                depth_bestmatch_idx = i
        if depth_bestmatch_idx not in copy_numbers:
            copy_numbers[depth_bestmatch_idx] = cn,min_dist
        else:
            old_cn,old_min_dist = copy_numbers[depth_bestmatch_idx]
            if  min_dist < old_min_dist:
                copy_numbers[depth_bestmatch_idx] = cn,min_dist
    
    #remove unneeded distances
    for k,pair in copy_numbers.items():
        copy_numbers[k] = pair[0]
    
    #add a column to the dataframe with the copy number by mapping it from the cluster
    data['cns'] = data['preds'].map(copy_numbers)
    return data


##########################################################################################################################

#main block
##########################################################################################################################
regions = pd.read_csv(args.depths_matrix, index_col=0)
for region, row in regions.iterrows():
    #region = region[2:-1]
    #if region != "1_104809_573868":
    #    continue
    #if region != "22_24352143_24386421":
    #    continue
    #if region != "1_87113_398211":
    #    continue

    preds, probs, n_clusters = EMCopyNumber(row)
    data = pd.DataFrame(row)
    data.columns = ['normdepths']
    data['preds'] = preds
    
    #print_report(row, preds, probs, n_clusters)
    predict_cn(data, region)
    plot_report(data, n_clusters, args.out_dir, region)
