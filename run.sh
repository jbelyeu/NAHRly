#run per-sample in ceph, current example is just one sample
crams="/scratch/ucgd/lustre/UCGD_Datahub/Repository/AnalysisData/2016/A414/16-08-06_WashU-Yandell-CEPH/UGP/Data/PolishedBams"
for filename in $crams/*.cram; do
    sample=$(echo "$filename" | grep -oP "\d+.bam" | grep -oP "\d+")
    echo "mosdepth -n --by ../betwixt_superdups.bed output/$sample $filename"
done
