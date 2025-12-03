#!/bin/sh

##sbatch --mem=200g --time=48:00:00 --cpus-per-task=8 run_demuxlet_test.sh

cd /data/Choi_lung/scLongreads/Demuxlet

for BAM in `ls ./BAM/*.bam`
do
NAME=$(basename $BAM _chr.bam)
VCF=$(echo ./VCF/"$NAME"_sorted.vcf | sed 's/NCI_/nci/g')
echo "cellsnp-lite --cellTAG CB --UMItag XM -s /data/Choi_lung/scLongreads/B2_percentile/Sample_"$NAME"/$NAME.dedup.mapped.bam \\
    -b /data/Choi_lung/scLongreads/B2_percentile/Sample_"$NAME"/"$NAME"_classification_filtered/barcodes1.tsv -R /data/Choi_lung/scLongreads/B2_percentile/Sample_"$NAME"/bcftools_variants.hq.vcf \\
    --minMAF 0.1 --minCount 20 --gzip \\
    -O /data/Choi_lung/scLongreads/B2_percentile/Sample_"$NAME"/cellsnp-lite_bcftools_het"
done

