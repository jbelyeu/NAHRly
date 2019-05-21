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
from scipy.stats import ttest_ind as ttest
#import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import vcfwriter

from scipy.spatial import distance_matrix
from scipy.signal import find_peaks, find_peaks_cwt
import seaborn as sns


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

def one_cluster_refine(region_info):
    cn = region_info["CN"]
    zscores = zscore(region_info['DP'])
    cn[zscores > 4] += 1
    cn[zscores < -4] -= 1

    #return region_info
    #step 2.2: using the new groups from 2.1, check each of the values from the next group for better fit 
    #do this again if any points change cn
    found_cns = np.unique(cn)
    cns_changed = True
    while cns_changed:
        cns_changed = False
        for i in range(len(found_cns)-1):
            #size of found_cns changes so have to check
            if len(found_cns) < i+2:
                break
            group1 = region_info['DP'][cn == found_cns[i]]
            group2 = region_info['DP'][cn == found_cns[i+1]]
            #TODO special case for len==1
            gp1_mean,gp1_std = np.mean(group1),np.std(group1)
            gp2_mean,gp2_std = np.mean(group2),np.std(group2)

            zscores_gp1 = []
            zscores_gp2 = []
            for j,dp in enumerate(region_info['DP']):
                if cn[j] == found_cns[i+1]:
                    zscores_gp1.append((dp-gp1_mean) / gp1_std)
                    zscores_gp2.append((dp-gp2_mean) / gp2_std)

                else:
                    zscores_gp1.append(np.inf)
                    zscores_gp2.append(np.inf)
            zscores_gp1 = np.array(zscores_gp1)
            zscores_gp2 = np.array(zscores_gp2)
            
            cn_counts = {}
            for found_cn in found_cns:
                cn_counts[found_cn] = len(cn[cn==found_cn])

            #if the point being tested falls closer by a standarddev to the other cluster, move it
            cn[np.absolute(zscores_gp1) < (np.absolute(zscores_gp2)-1.0)] -= 1
            
            new_found_cns = np.unique(cn)
            if not np.array_equal(new_found_cns, found_cns):
                found_cns = new_found_cns
                cns_changed = True
            else:
                new_cn_counts = {}
                for found_cn in found_cns:
                    new_cn_counts[found_cn] = len(cn[cn==found_cn])
                if new_cn_counts != cn_counts:
                    cns_changed = True

#    found_cns = np.unique(cn)
#    for i in range(len(found_cns)-1):
#        group1 = region_info['DP'][cn == found_cns[i]]
#        group2 = region_info['DP'][cn == found_cns[i+1]]
#        gp1_mean,gp1_std = np.mean(group1),np.std(group1)
#        gp2_mean,gp2_std = np.mean(group2),np.std(group2)
#
#        zscores_gp1 = []
#        zscores_gp2 = []
#        for j,dp in enumerate(region_info['DP']):
#            if cn[j] == found_cns[i+1]:
#                zscores_gp1.append((dp-gp1_mean) / gp1_std)
#                zscores_gp2.append((dp-gp2_mean) / gp2_std)
#
#            else:
#                zscores_gp1.append(np.inf)
#                zscores_gp2.append(np.inf)
#        zscores_gp1 = np.array(zscores_gp1)
#        zscores_gp2 = np.array(zscores_gp2)
#        
#        cn_counts = {}
#        for found_cn in found_cns:
#            cn_counts[found_cn] = len(cn[cn==found_cn])
#        cn[np.absolute(zscores_gp1) < np.absolute(zscores_gp2)] -= 1


    return region_info

def multi_cluster_refine(region_info):
    cluster_depths = []
    mergeable_cns = []
    for i,cn in enumerate(region_info['CN']):
        for _ in range(len(cluster_depths),cn+1):
            cluster_depths.append([])
        cluster_depths[cn].append(region_info['DP'][i])

    mergeable_cns = []
    for i in range(len(cluster_depths)-1):
        if len(cluster_depths[i]) == 0 or len(cluster_depths[i+1]) == 0:
            continue
        
        idxs = (i,i+1) if len(cluster_depths[i]) > len(cluster_depths[i+1]) else (i+1,i)

        small_group = np.array(cluster_depths[idxs[1]])
        large_group = np.array(cluster_depths[idxs[0]] + cluster_depths[idxs[1]])
        lg_mean = np.mean(large_group)
        lg_std = np.std(large_group)

        #zscores = np.array([((x - np.mean(large_group)) / np.std(large_group)) for x in small_group ])

        zscores = np.array([((x - lg_mean) / lg_std) for x in small_group ])
        
        #if the median zscore is less than X from 0, these should be merged
        if np.absolute(np.median(zscores)) < 3.0:
            mergeable_cns.append(idxs)
    
    #merge and return
    return merge_cns(region_info, mergeable_cns)

def merge_cns(region_info, mergeable_cns):
    mergeable_cns = sorted(mergeable_cns, reverse=True) 
    for mergeable_cn in mergeable_cns:
        for i in range(len(region_info['CN'])):
            if region_info['CN'][i] == mergeable_cn[1]:
                region_info['CN'][i] = mergeable_cn[0]
    return region_info


def refine_cn(region_info):
    #return region_info
    #step 1. move CN=1 samples that have been misclassified to CN=0 

    cn = region_info["CN"]
    cn[(cn == 0) & (region_info["DP"] > 0.4)] = 1

    #step 2.1: if there's only one CN called, use z-scores to move outliers by one
    if len(set(region_info['CN'])) == 1:
        region_info = one_cluster_refine(region_info)

    # NOTE: returning early because step 3 makes evaluation much worse.
    return region_info

    #step 3: if two consecutive peaks don't have different means, merge them
    region_info = multi_cluster_refine(region_info)
    return region_info
       
def plot_cns(name,troughs,utroughs,peaks,bins,upeaks,cns,dps,ndps):
    troughs = troughs[:-1]

    fig, axes = plt.subplots(1, 2, figsize=(12, 8))
    fig.suptitle(name)
    axes[0].plot(bins)
    colors = sns.color_palette()
    axes[0].plot(peaks, bins[peaks], "x", color=colors[5])
    axes[0].plot(troughs, bins[troughs], "o", color=colors[4])

    for u in upeaks:
        axes[1].axvline(u, color=colors[5], linewidth=2)
    #for u in utroughs:
    #    axes[1].axvline(u, color=colors[5], linewidth=1)

    df = pd.DataFrame({"dp":ndps, "cn": cns})
    cs = np.array([colors[min(len(colors)-1, d)] for d in cns])
    cs = cs[np.argsort(dps)]

    ax = sns.swarmplot(x="dp", y=['']*len(df), data=df, ax=axes[1], hue="cn")
    #plt.show()
    plt.savefig(os.path.join("/Users/jon/Desktop/nahrpics",name+".png"))

def plot_simple_cns(name,cns,dps):

    fig, axes = plt.subplots(1, 2, figsize=(12, 8))
    fig.suptitle(name)
    colors = sns.color_palette()

    df = pd.DataFrame({"dp":dps, "cn": cns})
    cs = np.array([colors[min(len(colors)-1, d)] for d in cns])
    cs = cs[np.argsort(dps)]

    ax = sns.swarmplot(x="dp", y=['']*len(df), data=df, ax=axes[1], hue="cn")
    plt.show()
 
def cn_probabilities(region_info):
    return region_info

def depth2CN(region_info, method="peaks", plot=True):
    if method == "GMM":
        import emgaussian
        region_info = emgaussian.EMCopyNumber(region_info)
        region_info = refine_cn(region_info)
        if plot:
            plot_simple_cns(region_info["name"],region_info['CN'],region_info['DP'])

        return region_info
    elif method == "PMM":
        import empoisson
        return empoisson.EMCopyNumber(region_info)

    dps = region_info["DP"]
    ndps = 2 * dps / np.maximum(0.4, np.median(dps))
    try:
        imed = np.where(ndps == np.median(ndps))[0][0]
    except IndexError:
        tmp = np.sort(ndps)
        v = tmp[len(tmp)//2]
        imed = np.where(np.abs(ndps - v) < 1e-7)[0][0]

    # NOTE that we could experiment with larger bins at the extremes to catch
    # rare events. or handle that post-hoc
    bins = np.zeros(100)
    #large_scaler = 22.0
    #TODO tune this scaler
    large_scaler = 12.0
    large_scaler = 10.0

    n_samples = len(ndps)
    scale_depth = np.minimum(bins.size - 1, np.round(ndps * (large_scaler if n_samples > 72 else 12.0)).astype(int))

    bins = np.bincount(scale_depth)

    #peak finding can't find peaks at the ends, so adding zero bins to both ends
    bins = np.append(bins,0)
    bins = np.insert(bins,0,0)

    # TODO: prominence of 3 means 3 samples must be in same bin, so won't find
    # de novos. figure how to deal with that post-hoc?
    peaks, peak_info = find_peaks(bins, distance=6, width=0)

    #upeaks (peaks in depth space) are equal to the middle of the peaks divided by the scalar
    upeaks = (peaks - (peak_info['widths']/2.0)) / (large_scaler if n_samples > 72 else 12.0)
    upeaks = peaks / (large_scaler if n_samples > 72 else 12.0)

    if len(upeaks) == 0:
        upeaks = np.array([ndps[imed]])
        peaks = np.array([2])

    # split bins at the troughs (lowest point between 2 peaks)
    troughs = []
    for start, stop in zip(peaks, peaks[1:]):
        # in case of tie, take mean, this splits in a saner place when there is
        # for example a long stretch of 0's with intermittent 1's .
        m = bins[start:stop]
        min_mean, = np.where(m == m.min())
        min_mean = int(0.5 + min_mean.mean())
        troughs.append(start + min_mean) #np.argmin(bins[start:stop]))
    troughs.append(len(bins))
    utroughs = np.array(troughs) / (large_scaler if n_samples > 72 else 12)
    cns = np.zeros(dps.shape[0], dtype=int)

    # the code below splits by trough instead.
    for i, p in enumerate(upeaks):
        if i == 0:
            #set to cn=0 if the depth is less than the first trough
            cns[ndps <= utroughs[0]] = i
        elif i < len(upeaks)-1:
            #set to cn=i where the depth is greater than the prev trough and less than next trough
            cns[(ndps >= utroughs[i-1]) & (ndps <= utroughs[i])] = i
        else:
            #set to i if the depth is greater than the last trough 
            cns[(ndps >= utroughs[i-1]) ] = i

    # re-scale so that CN2 is the value assigned to the median sample.
    cns += (2 - cns[imed])
    cns = np.maximum(0, cns)
    region_info["CN"] = cns
    region_info["NDPS"] = ndps
    raw_cns = list(cns)
    region_info = refine_cn(region_info)
    #region_info = cn_probabilities(region_info)


    if plot:
        plot_cns(region_info["name"],troughs,utroughs,peaks,bins,upeaks,region_info['CN'],dps,ndps)
    
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
normalized_depths = normalized_depths / normalized_depths.median(axis='rows') * 2
variances = normalized_depths.var(axis='rows')

chroms = sorted(set(x.split("_")[0] for x in normalized_depths.index))

cy_writer = vcfwriter.get_writer(args.vcf,normalized_depths.columns, chroms)
for region, row in normalized_depths.iterrows():
    interesting_regions = [
        "1_143009_317718", #test case for merging
        "1_267706_341907", #test case for zscore splitting
        "1_251124_404048", #test case for merging
        "1_121418_235524", #test case for 2.2
        "1_104809_404048", #test case for 2.2
        "1_104809_573868", #test case for 2.2
        "1_13005538_13005884",

    ]
    #if region not in interesting_regions: continue

    chrom,start,stop = region.split("_")
    start = int(start)
    stop = int(stop)
    region_info = {
        "chrom": chrom,
        "start": start,
        "stop": stop,
        "DP": row.values,
        "name": region
    }
    
    methods = ["peaks","GMM","PMM"]
    region_info = depth2CN(region_info, method=methods[1], plot=False)
    vcfwriter.write_variant(cy_writer, region_info)
