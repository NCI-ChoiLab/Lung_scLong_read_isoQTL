# Single-cell long-read RNA sequencing of human lung tissues

This repo includes the code to analyze the single-cell long-read dataset of lung tissues as part of manuscript "Single-cell full-length transcriptome of human lung reveals genetic effects on isoform regulation beyond gene-level expression".

The findings (i.e., isoform structure, isoform expression across lung cell types, and isoQTLs) in our paper could also be found and visulized via [ISOLUTION](https://appshare.cancer.gov/ISOLUTION)

---

## Analysis
### [isoQTL](https://github.com/NCI-ChoiLab/Lung_scLong_read_isoQTL/tree/main/Analysis/isoQTL)
For isoQTL mapping, we utilized negative binomial model by [jaxQTL](https://github.com/mancusolab/jaxqtl) to adress the sparsity of isoform expression. To determine the cell type specificity of isoQTL, [mashr EZ mode](https://github.com/NCI-ChoiLab/Lung_scLong_read_isoQTL/tree/main/Analysis/isoQTL/mashr) was used to harmonize the effect size. For functional annotation, we performed enrichment analysis of isoQTL in RNA-binding protein motifs (downloaded from [ENCODE](https://www.encodeproject.org/metadata/?type=PublicationData&%40id=%2Fpublication-data%2FENCSR456FVU%2F&files.preferred_default=true)). Our [single-cell multiome data](https://www.nature.com/articles/s41467-024-52356-9#Sec33) was also integrated with isoQTL and alternative splicing of isoforms.GWAS integration was performed by coloc and FUSION TWAS. Using our cell-barcode-matched short-read dataset, colocalization of sc-isoQTLs and eQTLs revealed the regulatory mechanism and GWAS contribution of isoQTL independent from eQTL.
### [Long-read integration](https://github.com/NCI-ChoiLab/Lung_scLong_read_isoQTL/tree/main/Analysis/Long-read%20integration)
This pipeline includes single-cell long-read data processing by ISO-seq, isoform integration across 22 single-cell batches, and further QC for protein-coding potential of novel isoforms by TransDecoder.
### [single-cell](https://github.com/NCI-ChoiLab/Lung_scLong_read_isoQTL/tree/main/Analysis/sinlge-cell)
Single cell data were asigned to specific individuals using genotype-based demultiplexing by vireoSNP. The code for single-cell process is all included in this folder. 
Differentially expression isoform (DEI) detection is described in [6.1.Isoform_sugnature_DEI_detect.R](https://github.com/NCI-ChoiLab/Lung_scLong_read_isoQTL/blob/main/Analysis/sinlge-cell/6.1.Isoform_signature_DEI_detect.R), and the benchmarking ofr isoform-level clustering and annotation can be found in [6.3.Clustering_annotation_w_isoform_expr.R](https://github.com/NCI-ChoiLab/Lung_scLong_read_isoQTL/blob/main/Analysis/sinlge-cell/6.3.Clustering_annotation_w_isoform_expr.R).

---

## data
The raw data coule be downloaded on GEO accession number ().
The genotype data used for isoQTL mapping could be found on dbGap under accession: phs004420.
In the data folder, the GTF file for post-QC isoforms identified by long-read sequencing is provided. 
In [mass-spec](https://github.com/NCI-ChoiLab/Lung_scLong_read_isoQTL/tree/main/data/mass-spec) folder, the mass-spec summaries from 4 individuals are included.

---

If you use our single-cell long-read data or ISOLUTION, please cite the following paper:

Li B, Luong T, Sisay E, Yin J, Zhang Z, Vaziripour M, Shin JH, Zhao Y, Byun J, Li Y, Lee CH, O'Neil M, Andresson T, Chang YS, Landi MT, Rothman N, Long E, Lan Q, Amos C, Zhou AX, Zhang T, Lee JG, Shi J, Xia J, Mancuso N, Zhang H, Kim EY, Choi J*. Single-cell full-length transcriptome of human lung reveals genetic effects on isoform regulation beyond gene-level expression. 2025