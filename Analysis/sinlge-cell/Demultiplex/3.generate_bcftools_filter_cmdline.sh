#!/bin/sh

##sbatch --mem=200g --time=48:00:00 --cpus-per-task=8 run_demuxlet_test.sh

cd /data/Choi_lung/scLongreads/Demuxlet

for BAM in `ls ./BAM/*.bam`
do
NAME=$(basename $BAM _chr.bam)
VCF=$(echo ./VCF/"$NAME"_sorted.vcf | sed 's/NCI_/nci/g')
echo "bcftools filter -i 'GT=\"0/1\" && QUAL > 51' /data/Choi_lung/scLongreads/B2_percentile/Sample_"$NAME"/bcftools_variants.vcf > /data/Choi_lung/scLongreads/B2_percentile/Sample_"$NAME"/bcftools_variants.hq.vcf"



done

