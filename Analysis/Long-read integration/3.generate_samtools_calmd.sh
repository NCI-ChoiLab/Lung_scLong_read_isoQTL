#!/bin/sh

cd /data/Choi_lung/scLongreads/Demuxlet

for BAM in `ls ./BAM/*.bam`
do
NAME=$(basename $BAM _chr.bam)
VCF=$(echo ./VCF/"$NAME"_sorted.vcf | sed 's/NCI_/nci/g')
echo "samtools calmd /data/Choi_lung/scLongreads/B2_percentile/Sample_"$NAME"/test.sam /data/Choi_lung/lbl/refdata-gex-GRCh38-2020-A/fasta/genome.fa --output-fmt sam > /data/Choi_lung/scLongreads/B2_percentile/Sample_"$NAME"/test.MDtagged.sam"

done