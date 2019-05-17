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
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

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


def EMCopyNumber(region_info, max_k=10):
    depths = region_info['DP'].reshape(-1, 1)
    bics = []
    gmms = []
    #create and fit gmms for all k values up to max_k
    for i in range(1,max_k):
        gmm = mixture.GaussianMixture(n_components=i)
        gmm.fit(depths)
        bics.append(float(gmm.bic(depths)))
        gmms.append(gmm)

    #find the number of components that provides the best fit 
    k = choose_k(bics)
    gmm = gmms[k]
   
    #calculate prediction and probabilities from the selected model    
    region_info['PRED'] = gmm.predict(depths)
    region_info['QUAL'] = gmm.predict_proba(depths)
    region_info = predict_cn(region_info)
    return region_info
    

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
    #ax = sns.swarmplot(x=data['cns'],y=data['normdepths'])
    #ax = sns.swarmplot(y="normdepths", x=[""]*len(data), data=data)
    ax = sns.swarmplot(y="normdepths", x=[""]*len(data), hue="cns", data=data)
    

    ax.set_title(region)
    #ax.set_title(region + ": with k=" + str(n_clusters))
    #plotname = os.path.join(directory,region+"_"+str(n_clusters)+"_skl.pdf")
    plotname = os.path.join(directory,region+".pdf")
    plt.ylabel("Normalized depth")
    plt.xlabel("Predicted copy number")
    plt.savefig(plotname)
    plt.close()

def predict_cn(region_info):
    clusters = np.unique(region_info['PRED'])

    depths = {}
    for i,clust in enumerate(region_info['PRED']):
        if clust not in depths: 
            depths[clust] = []
        depths[clust].append(region_info['DP'][i])
    #medians = {i:np.median(cluster) for i,cluster in depths.items()}
    medians = [[k,np.median(depths[k])] for k in depths]
    medians.sort(key=lambda x: x[1])

    map_clust2cn = {}
    cn = 0
    for pred,med in medians:
        map_clust2cn[pred] = cn
        cn +=1
    region_info['CN'] = []
    for i,pred in enumerate(region_info['PRED']):
        region_info['CN'].append(map_clust2cn[pred])
    region_info['CN'] = np.array(region_info['CN'])
    return region_info
    
    
    
    
#    counts = {i:len(cluster) for i,cluster in depths.items()}
#
#    largest_cluster = 0
#    for i,count in counts.items():
#        if count > counts[largest_cluster]:
#            largest_cluster = i
#    
#    #CN=2 is the most common
#    copy_numbers = {largest_cluster: 2}
#
#    #the number of copy numbers present is determined by the highest median
#    highest_median_cluster = 0
#    for i,median in enumerate(medians):
#        if median > medians[highest_median_cluster]:
#            highest_median_cluster = i
#
#    #the maximum copy number in the data for this site is the largest median 
#    #divided by the median of the largest cluster (by count of samples in it)
#    #I'm also adding a small amount in case the largest cluster is actually at zero which causes a divide by zero error
#    max_cn = (2*(round((medians[highest_median_cluster] / (medians[highest_median_cluster]+0.01) ))))
#    copy_numbers[highest_median_cluster] = max_cn
#
#    #create a list of predicted depths that correspond to the possible copy numbers
#    pred_depths = []
#    for i in range(int(max_cn) +1):
#        pred_depths.append(medians[highest_median_cluster] * (i/max_cn))
#    
#    copy_numbers = {}
#    for cn,pred_depth in enumerate(pred_depths):
#        min_dist = np.inf
#        depth_bestmatch_idx = -1
#
#        for i,median in enumerate(medians):
#            if abs(median-pred_depth) < min_dist:
#                min_dist = abs(median-pred_depth)
#                depth_bestmatch_idx = i
#        if depth_bestmatch_idx not in copy_numbers:
#            copy_numbers[depth_bestmatch_idx] = cn,min_dist
#        else:
#            old_cn,old_min_dist = copy_numbers[depth_bestmatch_idx]
#            if  min_dist < old_min_dist:
#                copy_numbers[depth_bestmatch_idx] = cn,min_dist
#    
#    #remove unneeded distances
#    for k,pair in copy_numbers.items():
#        copy_numbers[k] = pair[0]
#    
#    region_info['CN'] = []
#    for clust in region_info['PRED']:
#        region_info['CN'].append(copy_numbers[clust])
#    print(region_info['CN'])
#
#    return region_info


##########################################################################################################################

#main block
##########################################################################################################################
if __name__ == "__main__":
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



