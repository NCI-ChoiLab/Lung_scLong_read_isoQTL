#!/bin/sh

cd /data/Choi_lung/scLongreads/Demuxlet

for BAM in `ls ./BAM/*.bam`
do
NAME=$(basename $BAM _chr.bam)
VCF=$(echo ./VCF/"$NAME"_sorted.vcf | sed 's/NCI_/nci/g')
echo "/data/Choi_lung/scLongreads/B2_percentile/Sample_"$NAME"/$NAME.collapsed.sorted.filtered_lite.gff"

done