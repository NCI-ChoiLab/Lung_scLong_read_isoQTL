#!/bin/sh

# Download cDNA primers
wget https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/TUTORIAL-DATA-PBMC-single-cell-mini/primers.fasta

# Download cell barcode include list
wget https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/REF-10x_barcodes/3M-february-2018-REVERSE-COMPLEMENTED.txt.gz


for BAM in `ls *.hifi.bam`
do
NAME=$(basename $BAM .hifi.bam)
# remove primer
lima $BAM primers.fasta $NAME.fl.bam --isoseq
# Tag UMI and barcode
isoseq tag $NAME.fl.5p--3p.bam $NAME.flt.bam --design T-12U-16B

isoseq refine $NAME.flt.5p--3p.bam primers.fasta $NAME.fltnc.bam --require-polya
done


ls *.fltnc.bam > fltnc.fofn

# correct cell barcode and sort
samtools sort -t CB fltnc.corrected.bam -o fltnc.corrected.sorted.bam

# Remove the PCR duplication
isoseq groupdedup fltnc.corrected.sorted.bam dedup.bam
# map to GENCODE v32 to keep consistent with short-read mapping
# NOTE: dedup.mapped.bam is shared in GEO
pbmm2 align --preset ISOSEQ --sort dedup.bam ./refdata-gex-GRCh38-2020-A/fasta/genome.fa dedup.mapped.bam
# collapse the transcripts
isoseq collapse dedup.mapped.bam collapsed.gff

pigeon prepare collapsed.gff
# classify transcripts employing SQANTI3
pigeon classify -j 36 collapsed.sorted.gff ./refdata-gex-GRCh38-2020-A/gene/genes.gtf ./refdata-gex-GRCh38-2020-A/fasta/genome.fa --fl abundance.txt
# Filter out transcripts by default parameters
pigeon filter collapsed_classification.txt --isoforms collapsed.sorted.gff
# generate matrx for seurat
pigeon make-seurat --dedup dedup.fasta --group collapse.group.txt -d ./output/ classification.filtered_lite_classification.txt
# saturation analysis
pigeon report --exclude-singletons classification.filtered_lite_classification.txt saturation.txt