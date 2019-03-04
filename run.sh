#run per-sample in ceph, current example is just one sample
for filename in /scratch/ucgd/lustre/u1006375/ceph/ceph-bams/*.bam; do
    sample=$(echo "$filename" | grep -oP "\d+.bam" | grep -oP "\d+")
    echo "mosdepth -n --by ../betwixt_superdups.bed output/$sample $filename"
done
