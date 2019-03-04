#! /usr/bin/env python
from __future__ import print_function
import sys
import argparse
import os
import pandas as pd
import numpy as np

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

#MAXIMIZATION
#puts depths in closest bin
def maximization(depths, centers, centers_count):
    binned_dicts = [{} for i in range(0,centers_count)]
    for i in range(len(depths)):
        closest_idx = 0
        dist_to_closest = 0
        for j in range(len(centers)):
            dist_to_curr = abs(depths[i]-centers[j])
            if dist_to_closest == 0 or dist_to_curr < dist_to_closest:
                dist_to_closest = dist_to_curr
                closest_idx = j
        binned_dicts[closest_idx][depths.index[i]] = depths[i]
    
    binned = []
    for i in range(len(binned_dicts)):
        series = pd.DataFrame.from_dict(binned_dicts[i], orient='index')
        binned.append(series)
    return binned

#EXPECTATION
#moves the centers to fit the depths that have been assigned, based on responsibility for those centers
def expectation(depths, binned, centers, eps):
    #adjust CN=2 center
    centers[2] = np.mean(binned[2].values)
    
    #TODO ask about this use of responsibility. What impact does using it only on CN=2 have?
    #Notice that it is used indirectly on the others
    if centers[2] == 0:
        depths_count = len(depths)
        #we exclude the top bin to avoid over-adjusting lambda[2] for extreme outliers.
        for i in range(len(centers)-1):
            #find depth responsibility
            depth_responsibility = binned[i].size / depths_count

            #don't let the coverage at any CN drop below some cutoff eps
            if centers[i] < eps:
                centers[i] = eps

            #scale by depth responsibility
            centers[2] += np.mean(binned[i].values) * (2 / i) * depth_responsibility
    
    #adjust the expected depths of other copy-numbers based on that from CN2
    for i in range(len(centers)):
        centers[i] = centers[2] * (i / 2)

    #expand the range of CN=2
    #this makes sense because if it's more common, there are more opportunities for spread
    cn2_span = centers[2] - centers[1]
    centers[1] -= (cn2_span / 1.5)
    centers[3] += (cn2_span / 1.5)
    
    return centers

#given a list of centers and the list of previous centers, find the largest difference between them
def getMaxChange(centers, last_centers):
    max_change = 0
    for i in range(len(centers)):
        change = abs(centers[i] - last_centers[i])
        if change > max_change:
            max_change = change
    return max_change



#inspired by EMDepth in goleft
#returns a list of integer copy-numbers corresponding to given normalized depths.
#Uses a simple EM to assign depths to copy-number bins with a preference for CN=2.
#Adjusts mean depth after each iteration.
def EMCopyNumber(depths, maxiter=10, eps=0.01, centers_count=5):
    #median of the depths, assumed to be near CN=2
    median = np.median(depths)

    #centers of the copy numbers, should correspond roughly to CN=0, CN=1,CN=2,CN=3,CN>3
    centers = [0 for i in range(centers_count)]

    #for convergence check, store the centers from last iteration
    last_centers = [0 for i in range(centers_count)]
    
    #initialize the CN=2 center 
    centers[2] = median

    #EXPECTATION
    #initialize the centers of bins relative to CN=2
    #this can't be done with the function call to expectation because that's based on bins that haven't been assigned yet
    for i in range(len(centers)):
        if i != 2:
            #multiply the CN=2 depth by 1/2^1.1 gives a slightly scaled value near expected depth for each CN
            #if CN=2 has a depth of 1, this results in 0, 0.55, 1, 1.65, 2.2 being stored
            centers[i] = centers[2]*((i/2)**1.1)
    print(centers)
    
    #put the depths into their bins
    binned = maximization(depths, centers, centers_count)
    
    #store the largest depth change that occurred, starting with large number
    max_change = 10000

    #iterate at most maxiter times 
    #stop early if the max depth change from an iteration is <0.5
    for _ in range(maxiter):
        if max_change < 0.00005:
            break
        #store the last centers
        for i in range(1,len(centers)):
            last_centers[i] = centers[i]

        #MAXIMIZATION
        #puts things in correct bins
        binned = maximization(depths, centers, centers_count)

        #EXPECTATION
        centers = expectation(depths, binned, centers, eps)
        print(centers)

        #update the max_change for convergence testing
        max_change = getMaxChange(centers, last_centers)
    
    return binned

##########################################################################################################################

#main block
##########################################################################################################################
regions = pd.read_csv(args.depths_matrix, index_col=0)
for region, row in regions.iterrows():
    region = region[2:-1]
    if region != "1_1596470_1652900":
        continue
    binned = EMCopyNumber(row)
    print(region)
    for i in range(len(binned)):
        if i != 2 and len(binned[i]) > 0:
            print (str(i) + ": " + str(len(binned[i])))
            print(binned[i])
    print()
