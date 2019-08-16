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

# for data visualization
import random

#arg parsing
######################################################################################################################
parser = argparse.ArgumentParser(description="finds potential CNVs using a directory of depth bedfiles from MosDepth")
parser.add_argument("bed_dir")
parser.add_argument("-v", "--vcf", help="filename for the ouput in VCF format", required=True)
parser.add_argument("-d", "--duplications", help="filename containing duplication coordinates in bed format", required=True)
parser.add_argument("--visualization", help="file to write sample data for visualization", required=True)
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

def getDepths(bedfiles, regions, samples):
    data = []
    for bed in bedfiles:
        sample_depths_filtered = []
        sample_depths = readFile(bed)
        for region in regions:
            sample_depths_filtered.append(sample_depths[region])
        data.append(sample_depths_filtered)
    return pd.DataFrame(data, index=samples, columns=regions).transpose()

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
    plt.show()
    #plt.savefig(os.path.join("/Users/jon/Desktop/nahrpics",name+".png"))
    plt.close()

def cn_probabilities(region_info):
    return region_info

def depth2CN(region_info, plot=True):
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
    large_scaler = 10.0

    n_samples = len(ndps)
    scale_depth = np.minimum(bins.size - 1, np.round(ndps * (large_scaler if n_samples > 72 else 12.0)).astype(int))

    bins = np.bincount(scale_depth)

    #peak finding can't find peaks at the ends, so adding zero bins to both ends
    bins = np.append(bins,0)
    bins = np.insert(bins,0,0)

    # TODO: prominence of 3 means 3 samples must be in same bin, so won't find
    # de novos. figure how to deal with that post-hoc?

    #TODO distance between peaks should be a little less that half the median
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
    cn_mode = cn_vals[cn_counts == max(cn_counts)]


    cns += (2 - cn_mode)
    #if the lowest CN is greater than 0 but depth is very close to 0, the mode is probably at CN=1 instead of CN=2, so scale down
    lowest_cn = min(cns)
    if (lowest_cn > 0) and (np.mean(ndps[cns == lowest_cn]) < 0.2):
        cns -= 1

    cns = np.maximum(0, cns)
    region_info["CN"] = cns
    region_info["NDPS"] = ndps

    if plot and len(np.unique(region_info['CN'])) > 1:
        plot_cns(region_info["name"],troughs,utroughs,peaks,bins,upeaks,cns,dps,ndps, region_info['RAWDP'])

    return region_info

# https://scikit-learn.org/stable/modules/naive_bayes.html#gaussian-naive-bayes
# https://jakevdp.github.io/PythonDataScienceHandbook/05.05-naive-bayes.html
def compute_class_probabilities(region_info):
    from sklearn.naive_bayes import GaussianNB
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

regions = []
with open(args.duplications, 'r') as int_bed:
    for line in int_bed:
        region = line.strip().replace("\t", "_")
        regions.append(region)
depths_matrix = getDepths(bedfiles, regions, samples)

#remove regions where all the sample depth is zero
depths_matrix = depths_matrix.loc[(depths_matrix > 0).any(axis='columns')]

#normalize internally by median * 2 (assuming most common value is CN=2)
normalized_depths = depths_matrix
normalized_depths = normalized_depths / normalized_depths.median(axis='rows') * 2
variances = normalized_depths.var(axis='rows')


chroms = sorted(set(x.split("_")[0] for x in normalized_depths.index))
cy_writer = vcfwriter.get_writer(args.vcf,normalized_depths.columns, chroms)

# for data visualization
regions_info = {}

for region, row in normalized_depths.iterrows():
    interesting_regions = [
        #nice segregators
        "1_1631023_1672080",
        "1_12906129_12927940",
        "1_25594517_25655515",
        "1_152760772_152770178",
        "1_161565123_161599999",
        "1_196740355_196796318",
        "2_1538003_1539461",
        "2_87390047_87999999",
        "2_159712887_159724455",
        "4_4136767_4151983",

        #1-2 realish CNs
        "1_121418_235524",
        "1_251124_404048",
        "1_13005538_13005884",
        "1_13115870_13116215",
        "1_16879369_16996057",
        "1_17064396_17186106",
        "1_22319508_22319647",
        "1_26968746_26972644",
        "1_28543332_28544608",
        "1_104122235_104159277",
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
        "name": region,
        "RAWDP": depths_matrix.loc[region]
    }
    npsamples = np.array(samples)

    region_info = depth2CN(region_info, plot=False)

    # use Naive Bayes to compute confidence in inferred copy numbers
    region_info = compute_class_probabilities(region_info)

    # for data visualization
    regions_info[region_info['name']] = {
        'chromosome': region_info['chrom'],
        'start': region_info['start'],
        'stop': region_info['stop'],
        'possible copy numbers': region_info['class_labels'].tolist(),
        'samples': [
            {
                'normalized read depth': float(normalized_read_depth),
                'random number': random.random(),
                'histogram-inferred copy number': int(histogram_inferred_copy_number),
                'probabilities of possible copy numbers': class_probabilities.tolist(),
                'bayes-inferred copy number': int(bayes_inferred_copy_number)
            }
            for (
                normalized_read_depth,
                histogram_inferred_copy_number,
                class_probabilities,
                bayes_inferred_copy_number
            ) in zip (
                region_info['DP'],
                region_info['CN'],
                region_info['class_probabilities'],
                region_info['bayes_inferred_copy_number']
            )
        ]
    }

    vcfwriter.write_variant(cy_writer, region_info)

# for data visualization
import json
import random
sample_size = 10
random_keys = random.sample(regions_info.keys(), sample_size)
regions_info_sample = {key: regions_info[key] for key in random_keys}
regions_info_sample = json.dumps(regions_info_sample)
regions_info_sample = 'export default {}'.format(regions_info_sample)
with open(args.visualization, 'w') as fp:
    fp.write(regions_info_sample)

