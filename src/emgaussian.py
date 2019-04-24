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
import operator

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
    if not os.path.exists(directory):
        os.makedirs(directory)
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set(style="whitegrid")
    plt.figure(figsize=(12,9))

    #ax = sns.swarmplot(y="normdepths", x=[""]*len(data), hue="preds", data=data)
    ax = sns.swarmplot(y="normdepths", x=[""]*len(data), hue="cns", data=data)
    
    
    loc_info = region.split("_")
    ax.set_title(loc_info[0]+":"+loc_info[1]+"-"+loc_info[2])
    #ax.set_title(loc_info[0]+":"+loc_info[1]+"-"+loc_info[2] + ": with k=" + str(n_clusters))
    #plotname = os.path.join(directory,region+"_"+str(n_clusters)+"_skl.pdf")
    plotname = os.path.join(directory,region+"_"+str(n_clusters)+"_skl.png")
    plt.ylabel("Normalized depth")
    plt.xlabel("Predicted copy number")
    plt.savefig(plotname)
    plt.close()


def predict_cn(data, region):
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
    # loop through the predicted possible copy numbers to find the closest cluster median
    for cn,pred_depth in enumerate(pred_depths):
        min_dist = np.inf
        depth_bestmatch_idx = -1
        
        # loop through the medians, find smallest distance to a predicted depth's median
        for i,median in enumerate(medians):
            if abs(median-pred_depth) < min_dist:
                min_dist = abs(median-pred_depth)
                depth_bestmatch_idx = i
        
        # if the closest depth to the median is not in the copynumbers list, add it
        if depth_bestmatch_idx not in copy_numbers:
            copy_numbers[depth_bestmatch_idx] = cn,min_dist
        else:
            old_cn,old_min_dist = copy_numbers[depth_bestmatch_idx]
            # TODO fix this. Currently it overwrites if a depth median is closer to the copy number but does nothing with the other cn
            if  min_dist < old_min_dist:
                copy_numbers[depth_bestmatch_idx] = cn,min_dist
    
    #remove unneeded distances
    for k,pair in copy_numbers.items():
        copy_numbers[k] = pair[0]
    
    #add a column to the dataframe with the copy number by mapping it from the cluster
    data['cns'] = data['preds'].map(copy_numbers)
    return data



def simply_predict_cn(data, region):
    medians = data.groupby("preds")["normdepths"].median()
    cluster_medians = []
    for i,median in enumerate(medians):
        cluster_medians.append([i,median])

    #this sort naively solves the problem
    #copy number is assigned to each cluster based on smallest to largest median value
    cluster_medians.sort(key=operator.itemgetter(1))
    copy_numbers = {}

    #to be a bit smarter, we increase the copynumber guess if the distance between copoy numbers is greater than
    #half the distance from cn=0 to cn=2, times 1.5
    
    if len(cluster_medians) <= 1 or True:
        #i is the naive cn prediction, pred_num is the cluster number
        for i,pair in enumerate(cluster_medians):
            pred_num,median = pair
            copy_numbers[pred_num] = i
    else:
        dist_between_cns = -1

        if len(cluster_medians) == 2:
            dist_between_cns = cluster_medians[1][1] - cluster_medians[0][1]
        else:
            dist_between_cns = 0.5*(cluster_medians[2][1] - cluster_medians[0][1])

        stretched_cluster_medians = []
        #i is the naive cn prediction, pred_num is the cluster number
        for i,pair in enumerate(cluster_medians):
            pred_num,median = pair
            cn = i
            if i > 0:
                if median > (cluster_medians[i-1][1] + (dist_between_cns*1.5)):
                    #this means the current cluster is far larger than it should be, so we skipped a cn between them
                    cn += 1
            copy_numbers[pred_num] = cn
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
    #if region != "1_12906129_12927940":
    #    continue
    # missing two clusters in plot
    #if region != "15_45223674_45249700":
    #    continue
    if not region.startswith("1_"):
        continue

    preds, probs, n_clusters = EMCopyNumber(row)
    data = pd.DataFrame(row)
    data.columns = ['normdepths']
    data['preds'] = preds
    
    #print_report(row, preds, probs, n_clusters)
    simply_predict_cn(data, region)

    region_fields = region.split("_")
    interval_len = int(region_fields[2]) - int(region_fields[1])
    if n_clusters > 3 and (20 < interval_len < 1000000):
        plot_report(data, n_clusters, args.out_dir, region)
