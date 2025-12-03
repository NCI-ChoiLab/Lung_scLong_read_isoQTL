#!/bin/sh

##sbatch --mem=200g --time=48:00:00 --cpus-per-task=8 run_demuxlet_test.sh

cd /data/Choi_lung/scLongreads/Demuxlet

for BAM in `ls ./BAM/*.bam`
do
NAME=$(basename $BAM _chr.bam)
VCF=$(echo ./VCF/"$NAME"_sorted.vcf | sed 's/NCI_/nci/g')
echo "bcftools mpileup --threads 32 -f /data/Choi_lung/lbl/refdata-gex-GRCh38-2020-A/fasta/genome.fa /data/Choi_lung/scLongreads/B2_percentile/Sample_"$NAME"/$NAME.dedup.mapped.bam \\
    -X pacbio-ccs --indel-size 500 --annotate FORMAT/AD,FORMAT/DP | \\
    bcftools call -f GQ --threads 32 -vm -Ov | \\
    bcftools norm --threads 32 -c s -f /data/Choi_lung/lbl/refdata-gex-GRCh38-2020-A/fasta/genome.fa - | \\
    bcftools filter --threads 32 -i 'QUAL > 20 && INFO/DP > 4' > /data/Choi_lung/scLongreads/B2_percentile/Sample_"$NAME"/bcftools_variants.vcf"
done

