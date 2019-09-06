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
import matplotlib.pyplot as plt
import numpy as np
import vcfwriter
import pysam

from sklearn.naive_bayes import GaussianNB
from scipy.spatial import distance_matrix
from scipy.signal import find_peaks
import seaborn as sns


#arg parsing
######################################################################################################################
parser = argparse.ArgumentParser(description="finds potential CNVs using a directory of depth bedfiles from MosDepth")
parser.add_argument("bed_dir")
parser.add_argument("-v", "--vcf", help="filename for the ouput in VCF format", required=True)
parser.add_argument("-s", "--segdups", help="bedfile of segdups", default="../ucsc_superdups_hg19_samechrom_sorted.bed.gz")
parser.add_argument("-r", "--regions", help="bedfile of segdups", default="mosdepth_regions.bed")
args = parser.parse_args()
#####################################################################################################################

#functions and generators
######################################################################################################################
def readFile(bedfile):
    depths = {}
    with gz.open(bedfile, 'rb') as open_bed:
        for line in open_bed:
            fields = line.decode("utf-8").strip().split()

            region = "_".join([str(x) for x in fields[:3]])
            depths[region] = float(fields[3])
    return depths

#TODO this should return a list of the segdups paired, the interstitial regions, and the flanks
def readSegDups(regions_filename):
    region_types = [
        "preflank",
        "segdup1",
        "nahr_region",
        "segdup2",
        "postflank",
    ]

    nahr_regions = {}
    with open(regions_filename, 'r') as regions_file:
        region_info = {
            }
        for line in regions_file:
            fields = line.strip().split()
            region = fields[0] + "_" + fields[1] + "_" + fields[2]
            region_info[region_types[len(region_info)%5]] = region

            if len(region_info) == 5:
                nahr_regions[region_info["nahr_region"]] = region_info
                region_info = {}

    return pd.DataFrame.from_dict(nahr_regions, orient='index')

def getDepths(bedfiles, samples):
    data = []
    for bed in bedfiles:
        sample_depths = readFile(bed)
        data.append(sample_depths)
    return pd.DataFrame(data, index=samples, columns=data[0].keys()).transpose()

def ext_sname(filename):
    return os.path.basename(filename).split(".")[0]

def plot_cns(name,troughs,utroughs,peaks,bins,upeaks,cns,dps,ndps, raw_dp):
    troughs = troughs[:-1]

    fig, axes = plt.subplots(1, 2, figsize=(12, 8))
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(name)
    grid = plt.GridSpec(2, 2)
    dp_hist = fig.add_subplot(grid[0:1, 0:1])
    ogdp_hist = fig.add_subplot(grid[1:2, 0:1])
    dp_swarm   = fig.add_subplot(grid[0:2, 1:2])

    dp_hist.plot(bins)
    colors = sns.color_palette()
    dp_hist.plot(peaks, bins[peaks], "x", color=colors[5])
    dp_hist.plot(troughs, bins[troughs], "o", color=colors[4])

    sns.distplot(raw_dp.values, ax=ogdp_hist, kde=False)


    for u in upeaks:
        dp_swarm.axvline(u, color=colors[5], linewidth=2)

    df = pd.DataFrame({"ndp":ndps, "cn": cns})
    cs = np.array([colors[min(len(colors)-1, d)] for d in cns])
    cs = cs[np.argsort(dps)]

    for i in range(1):
        ax = sns.swarmplot(x="ndp", y=['']*len(df), data=df, ax=dp_swarm, hue="cn")
    #plt.show()
    plt.savefig(os.path.join("/Users/jon/Desktop/tmp/nahrpics",name+".png"))
    #plt.savefig(os.path.join("/Users/jon/Desktop/",name+".png"))
    plt.close()

def cn_probabilities(region_info):
    return region_info

def depth2CN(region_info, plot=True):
    dps = region_info["DP"]

    #scaled+normalized depths are equal to the previous normalized depths times 2 (because diploidy is assumed normal),
    #and divided by either .4 or the median. Don't know why, this is a magic number to me.
    ndps = 2 * dps / np.maximum(0.4, np.median(dps))

    #find the indices of values at or near the median
    try:
        imed = np.where(ndps == np.median(ndps))[0][0]
    except IndexError:
        tmp = np.sort(ndps)
        v = tmp[len(tmp)//2]
        imed = np.where(np.abs(ndps - v) < 1e-7)[0][0]

    # NOTE that we could experiment with larger bins at the extremes to catch
    # rare events. or handle that post-hoc
    bins = np.zeros(100)
    large_scaler = 10.0

    n_samples = len(ndps)
    scale_depth = np.minimum(bins.size - 1, np.round(ndps * (large_scaler if n_samples > 72 else 12.0)).astype(int))

    bins = np.bincount(scale_depth)

    #peak finding can't find peaks at the ends, so adding zero bins to both ends
    bins = np.append(bins,0)
    bins = np.insert(bins,0,0)

    # TODO: prominence of 3 means 3 samples must be in same bin, so won't find
    # de novos. figure how to deal with that post-hoc?

    #distance between peaks should be a little less that half the median
    distance = max((np.mean(scale_depth)/2.0)*0.8, 1)

    peaks, peak_info = find_peaks(bins, distance=distance, prominence=2, rel_height=0.7)

    #upeaks (peaks in depth space) are equal to the middle of the peaks divided by the scalar
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
        thing = np.where(m == m.min())
        min_mean, = np.where(m == m.min())
        min_mean = int(0.5 + min_mean.mean())
        troughs.append(start + min_mean)
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

    # re-scale so that CN2 is the value assigned to the mode CN
    (cn_vals,cn_counts) = np.unique(cns,return_counts=True)
    cn_mode = max(cn_vals[cn_counts == max(cn_counts)])
    
    cns += (2 - cn_mode)
    #if the lowest CN is greater than 0 but depth is very close to 0, the mode is probably at CN=1 instead of CN=2, so scale down
    lowest_cn = min(cns)
    if (lowest_cn > 0) and (np.mean(ndps[cns == lowest_cn]) < 0.2):
        cns -= 1

    cns = np.maximum(0, cns)

    region_info["CN"] = cns
    region_info["NDPS"] = ndps
    region_info['hist_CN'] = cns

    # use Naive Bayes to compute confidence in inferred copy numbers and changeto most probable CNs
    region_info = compute_class_probabilities(region_info)
    region_info['hist_CN'] = np.copy(region_info['CN'])
    region_info['CN'] = region_info['bayes_inferred_copy_number']

    #perform a duphold-esque flank test
    # this will reassign to CN=2 and variant calls that fail #
    # to show a different depth of coverage from flanking regions outside the flanking segdups
    region_info['bayes_CN'] = np.copy(region_info['CN'])
    region_info = flankRefineCN(region_info)
    region_info['GQ'] = region_info['class_probabilities']


    #repeat the naive bayes confidence calculation, this time without changing CN calls
    region_info = compute_class_probabilities(region_info)

    if plot and len(np.unique(region_info['CN'])) > 1:
        plot_cns(region_info["name"],troughs,utroughs,peaks,bins,upeaks,region_info['CN'],dps,ndps, region_info['RAWDP'])

    return region_info

def flankRefineCN(region_info):
    running_sum = 0
    running_count = 0
    for i,cn in enumerate(region_info['CN']):
        if cn == 2: 
            continue

        preflank = region_info["preflank"][i]
        postflank = region_info["postflank"][i]
        segdup1 = region_info["segdup1"][i]
        segdup2 = region_info["segdup2"][i]

        flank_passed = True
        dp = region_info['DP'][i]
        for flank in [preflank,postflank]:
            if flank is not None and flank > 0.0:
                fold_change = (dp/flank)
                if (fold_change > 0.7 and cn < 2) or (fold_change < 1.3 and cn > 2):
                    #print(238)
                    region_info['CN'][i] = 2
                    flank_passed = False
            else:
                region_info['CN'][i] = 2
                flank_passed = False
    return region_info
# https://scikit-learn.org/stable/modules/naive_bayes.html#gaussian-naive-bayes
# https://jakevdp.github.io/PythonDataScienceHandbook/05.05-naive-bayes.html
def compute_class_probabilities(region_info):
    model = GaussianNB()
    X = region_info['DP']
    X = X[..., np.newaxis]
    y = region_info['CN']
    model.fit(X, y)
    region_info['class_probabilities'] = model.predict_proba(X)
    region_info['bayes_inferred_copy_number'] = model.predict(X)
    region_info['class_labels'] = model.classes_
    return region_info
# TODO:
# Down-weigh bad samples with variable coverage when computing p(x) and p(z), 
# where x is normalized read depth and z is copy number. 
# Genome Strip seems to do this: 
# “The model incorporates sample-specific variance terms 
# to model the variation in sequencing depth between samples” 
# https://www.nature.com/articles/ng.3200.pdf



######################################################################################################################

#main block
######################################################################################################################
bedfiles = glob.glob(os.path.join(args.bed_dir,"*.regions.bed.gz"))
samples = [ext_sname(x) for x in bedfiles]

nahr_regions = readSegDups(args.regions)
depths_matrix = getDepths(bedfiles, samples)

#remove regions where all the sample depth is zero
depths_matrix = depths_matrix.loc[(depths_matrix > 0).any(axis='columns')]

#normalize internally by median * 2 (assuming most common value is CN=2)
normalized_depths = depths_matrix
normalized_depths = normalized_depths / normalized_depths.median(axis='rows') * 2
variances = normalized_depths.var(axis='rows')


missing = np.setdiff1d(nahr_regions.index, normalized_depths.index)
missing_idxs = []
for miss in missing:
    #0th index is row indices, all that I need for 1D array
    missing_idxs += list(np.where(nahr_regions.index == miss)[0])
missing_labels = nahr_regions.iloc[missing_idxs].index.tolist()
nahr_regions = nahr_regions[~nahr_regions.index.isin(missing_labels)]

chroms = sorted(set(x.split("_")[0] for x in normalized_depths.index))
cy_writer = vcfwriter.get_writer(args.vcf,normalized_depths.columns, chroms)
tbx = pysam.TabixFile(args.segdups)

count=0
count2=0
for region, row in normalized_depths.loc[nahr_regions.index.tolist(), :].iterrows():
    nahr_region = nahr_regions.loc[region]
    #if region != "2_1538002_1539462":
    #    continue


    chrom,start,stop = region.split("_")
    start = int(start)
    stop = int(stop)

    #for now, don't look at nahr regions with segdup overlaps
    segdup_overlaps = []
    for entry in tbx.fetch("chr"+chrom, start, stop):
        segdup_overlaps.append(entry)
    if len(segdup_overlaps) > 0:
        count+=1
        continue
    count2+=1
 
    region_info = {
        "chrom": chrom,
        "start": start,
        "stop": stop,
        "DP": row.values,
        "name": region,
        "RAWDP": depths_matrix.loc[region],
        "ID": row.index
    }
    npsamples = np.array(samples)

    region_info["preflank"] = normalized_depths.loc[nahr_region["preflank"]] if nahr_region["preflank"] in normalized_depths.index else [None]*len(row)
    region_info["segdup1"] = normalized_depths.loc[nahr_region["segdup1"]] if nahr_region["segdup1"] in normalized_depths.index else [None]*len(row)
    region_info["segdup2"] = normalized_depths.loc[nahr_region["segdup2"]] if nahr_region["segdup2"] in normalized_depths.index else [None]*len(row)
    region_info["postflank"] = normalized_depths.loc[nahr_region["postflank"]] if nahr_region["postflank"] in normalized_depths.index else [None]*len(row)
    region_info = depth2CN(region_info, plot=False)

    #write to file only if there's a variant
    if len(np.unique(region_info['CN'])) > 1 or 2 not in region_info['CN']:
        vcfwriter.write_variant(cy_writer, region_info)
        bed_name = "bedfiles/" + region_info['name']+ ".bed"
        with open(bed_name, 'w') as region_bed:
            print(nahr_region['segdup1'].replace("_","\t"), file=region_bed)
            print(nahr_region['segdup2'].replace("_","\t"), file=region_bed)
        pysam.tabix_index(bed_name, force=True, seq_col=0, start_col=1, end_col=2)

print(count)
print(count2)
