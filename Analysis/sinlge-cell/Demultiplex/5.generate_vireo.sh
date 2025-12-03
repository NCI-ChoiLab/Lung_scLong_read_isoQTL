#!/bin/sh

##sbatch --mem=200g --time=48:00:00 --cpus-per-task=8 run_demuxlet_test.sh

cd /data/Choi_lung/scLongreads/Demuxlet

for BAM in `ls ./BAM/*.bam`
do
NAME=$(basename $BAM _chr.bam)
VCF_file=$(echo "$NAME"_sorted.vcf | sed 's/^.*NCI_/nci/g')
VCF=$(echo ./VCF/$VCF_file)
echo "vireo -t GT -c /data/Choi_lung/scLongreads/B2_percentile/Sample_"$NAME"/cellsnp-lite_bcftools_het -d $VCF \\
    -o /data/Choi_lung/scLongreads/B2_percentile/Sample_"$NAME"/vireo_bcftools_het \\
    --callAmbientRNAs"
done

