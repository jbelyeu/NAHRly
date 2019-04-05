#! /usr/bin/env python
from __future__ import print_function
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import sys
import argparse
import os
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import Families

#argparse
##########################################################################################################################
parser = argparse.ArgumentParser(description="finds potential CNVs using a directory of depth bedfiles from MosDepth")
#parser.add_argument("-e", "--del_matrix",
#    help="CSV file with putative deletion regions from gnarly",
#    required=True
#)
parser.add_argument("-r", "--regions_matrix", 
    help="CSV file with regions from gnarly",
    required=True
)
#parser.add_argument( "--low_depth_matrix", 
#    help="CSV file with putative duplication regions from gnarly",
#    required=True
#)
parser.add_argument("--zscore_matrix", 
    help="CSV file with Z-scores for each sample for regions with putative CNVs",
    required=True
)
parser.add_argument("-z",
    type=float,
    required=True
)
parser.add_argument("-p", "--pedfile", 
    help="PED file with the relationships of the samples in the set",
    required=True
)
parser.add_argument("-o", "--out_dir", 
    help="directory for output",
    required=True
)
args = parser.parse_args()
##########################################################################################################################


#functions 
##########################################################################################################################
def map_fam(families):
    samples_by_fam = {}
    for family_id in families:
        samples_by_fam[family_id] = list(families[family_id].sample_ids)
    return samples_by_fam

def map_gen(families):
    samples_by_gen = {}
    for family_id in families:
        for sample in families[family_id].sample_ids:
            gen = families[family_id].checkGeneration(sample)
            if gen not in samples_by_gen:
                samples_by_gen[gen] = []
            samples_by_gen[gen].append(sample)
    return samples_by_gen

def hist_count_over_depth(df, directory, region_of_interest=False):
    os.makedirs(directory, exist_ok=True)
    for index, row in df.iterrows():
        region = index
        if sys.version_info[0] >= 3.0:
            region = row.name[2:-1]
        if region_of_interest and region != region_of_interest:
            continue
        row.plot(kind='hist', bins=50)
        plt.xlabel=("Scaled Depth")
        plt.ylabel=("# of Samples")
        plt.savefig(os.path.join(directory,region+".png"))
        plt.close()

def scatter_z_over_depth(df, samples_by_fam,directory, region_of_interest=False):
    os.makedirs(directory, exist_ok=True)
    #TODO find the upper and lower z-score bounds for the row
    for index,row in df.iterrows():
        region = index
        if sys.version_info[0] >= 3.0:
            region = row.name[2:-1]
        if region_of_interest and region != region_of_interest:
            continue        

        for fam in samples_by_fam:
            samples = samples_by_fam[fam]
            depths = row.loc[samples]
            zs = zscores.loc[index].loc[samples]
            plt.plot(depths.values.tolist(), zs.values.tolist(),'o')
            plt.xlabel('Scaled Depth')
            plt.ylabel('Z Score')
            plt.ylim(-10,10)
            plt.yticks(np.arange(-10, 10, 1.0))
            plt.axhline(y=4, linestyle="--")
            plt.axhline(y=-4, linestyle="--")
            title = region.replace("_", ":",1)
            title = title.replace("_", "-",1)
            plt.title(title)

            plotname = os.path.join(directory,region+".png")
        plt.savefig(plotname)
        plt.close()

def getFamID(families, sample):
    family_id = ''
    for fam in families:
        if families[fam].hasSample(sample):
            return fam
    return False

#needs te depths for a region, the zscores for a region, and the families object list
#returns the matrix for splitting the data out, with columns for sample,family,depth,zscore,generation
def summarize_region(df, region_zscores, families):
    labels=["SAMPLE", "FAMILYID","GENERATION","DEPTH","Z-SCORE"]
    summary_matrix = []
    for sample, depth in df.iteritems():
        fam_id = getFamID(families,sample)
        if not fam_id:
        #    print("broken, "+str(sample),file=sys.stderr)
            continue
        gen = families[fam_id].checkGeneration(sample)
        summary_matrix.append([sample,fam_id.split("_")[0],gen,depth,region_zscores.loc[sample]])
    return pd.DataFrame(summary_matrix,columns=labels)


def swarm_z(df,zscores,directory, region_of_interest=False):
    os.makedirs(directory, exist_ok=True)
    print(directory)
    location = "/uufs/chpc.utah.edu/common/home/u1072557/nahr/nahr_cn/temp_data/"
    import seaborn as sns
    sns.set(style="whitegrid")
    for index,row in df.iterrows():
        region = index
        if sys.version_info[0] >= 3.0:
            region = row.name[2:-1]
        if region_of_interest and region != region_of_interest:
            continue        

        upper_depth_cutoff = (row.std()*(args.z)) + row.mean()
        lower_depth_cutoff = (row.std()*(-args.z)) + row.mean()
        region_zs = zscores.loc[index]
        region_summary = summarize_region(row,region_zs,families)
        plt.figure(figsize=(12,9))
        ax = sns.swarmplot(data=region_summary, x="GENERATION", y="DEPTH", hue="FAMILYID", order=["P0", "F1", "F2"])
        ax.axhline(upper_depth_cutoff,linestyle="--", linewidth=.5)
        ax.axhline(lower_depth_cutoff,linestyle="--", linewidth=.5)
        plt.setp(ax.get_legend().get_texts(), fontsize='8')
        plt.setp(ax.get_legend().get_title(), fontsize='10')
        plotname = os.path.join(directory,region+".pdf")
        plt.savefig(plotname)
        plt.close()
        #row.to_csv(location+region+".csv")



##########################################################################################################################


#main block
##########################################################################################################################
#low_depth_regions       = pd.read_csv(args.low_depth_matrix, index_col=0)
regions    = pd.read_csv(args.regions_matrix, index_col=0)
#possible_del_regions    = pd.read_csv(args.del_matrix, index_col=0)
#possible_dup_regions    = pd.read_csv(args.dup_matrix, index_col=0)
zscores                 = pd.read_csv(args.zscore_matrix, index_col=0)
families                = Families.CreateFamilies(args.pedfile, True)
samples_by_fam          = map_fam(families)
samples_by_gen          = map_gen(families)
for fam in families:
    if "1341" in fam:
        print (families[fam].pretty_print())

#plot the regions that might include CNVs
#hist_count_over_depth(possible_del_regions, os.path.join(args.out_dir,"hist","del"))
#hist_count_over_depth(possible_dup_regions, os.path.join(args.out_dir,"hist","dup"))

#plot z over depth
#scatter_z_over_depth(possible_dup_regions, samples_by_fam, os.path.join(args.out_dir,"scatter","del"))
#scatter_z_over_depth(possible_del_regions, samples_by_fam, os.path.join(args.out_dir,"scatter","dup"))

#swarmplot z colored by family
#swarm_z(possible_del_regions, zscores, os.path.join(args.out_dir,"swarm", "del"))
swarm_z(regions, zscores, os.path.join(args.out_dir,"swarm"), "22_24352143_24386421")
#swarm_z(low_depth_regions, zscores, os.path.join(args.out_dir,"swarm", "low_depth"))
