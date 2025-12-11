# Integrate isoQTL with single-cell multi-omics data of lung tissues
# https://doi.org/10.1038/s41467-024-52356-9

# Overlap isoQTL with cCREs identified by scATAC-seq data

# Bolun

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(sctransform)
library(patchwork)
library(GenomeInfoDb)
library(tidyverse)
library(grid)
library(readr)
library(harmony)
library(GenomicRanges)
library(SeuratData)
source("/data/lib14/R/Rscripts/utilities.R")

# Obtain cCRE called-in matrix
share_seur <- readRDS("/data/Choi_lung/ChiaHan/CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds")
rm(CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB)
peaks_by_ct <- read.csv("/data/Choi_lung/ChiaHan/peaks_by_cell_types_CHL_chranno.csv", header = TRUE, row.names = 1)
peaks_annot <- read.csv("/data/Choi_lung/ChiaHan/Peak_and_motifBreaker_results/peaks_annotation_by_ChIPseeker_CHL_chranno.csv", header = TRUE, row.names = 1)
peak_called_in.list <- strsplit(peaks_by_ct$peak_called_in,",")
celltypes <- as.character(unique(share_seur$CellType))
peak_called_in.list <- lapply(peak_called_in.list, function(x) {
  i <- which(celltypes %in% x)
  idxs <- rep(FALSE, length(celltypes))
  idxs[i] <- TRUE
  return(idxs)
})
peak_called_in <- do.call(rbind, peak_called_in.list)

colnames(peak_called_in) <- celltypes
rownames(peak_called_in) <- paste(peaks_by_ct$seqnames,peaks_by_ct$start, peaks_by_ct$end, sep = "-")

# Load NB significant isoQTLs
NB.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")

final_table <- read.table("/data/lib14/project/scLongread/Tables/NB_eIsoform_sum.txt", sep = "\t",
                          header = TRUE)
length(unique(final_table$variant_id))
head(snp_info)

peaks_annot$peak <- paste(peaks_annot$seqnames, peaks_annot$start, peaks_annot$end, sep = "-")
peaks_annot$peak_called_id <- peaks_by_ct$peak_called_in

# Get the isoQTLs overlapped with cCRE 
leadisoQTL_int_multiome.list <- list()
for (celltype in names(NB.sig.list)) {
  leadisoQTL <- NB.sig.list[[celltype]]
  leadisoQTL$peak <- NA
  leadisoQTL$Celltype <- celltype
  leadisoQTL <- leadisoQTL %>% group_by(phenotype_id) %>% mutate(transcript_name = Search_transcript_name2(phenotype_id),
                                                                 .after = phenotype_id)
  leadisoQTL <- leadisoQTL %>% group_by(variant_id) %>% mutate(chr = paste0("chr", chrom),
                                                               pos = snp_info[variant_id,]$pos,
                                                               ref = snp_info[variant_id,]$ref,
                                                               alt = snp_info[variant_id,]$alt,
                                                               .after = variant_id)
  for(i in 1:nrow(leadisoQTL)){
    cCRE <- peaks_annot$peak[which(peaks_annot$seqnames == leadisoQTL$chr[i] & 
                                     as.numeric(peaks_annot$start) <= as.numeric(leadisoQTL$pos[i]) & 
                                     as.numeric(peaks_annot$end) >= as.numeric(leadisoQTL$pos[i]))]
    if (length(cCRE) == 1) {
      leadisoQTL$peak[i] <- cCRE
    }
  }
  leadisoQTL_subset <- leadisoQTL[which(!is.na(leadisoQTL$peak)),]
  leadisoQTL_subset <- left_join(leadisoQTL_subset, peaks_annot, by = "peak")
  leadisoQTL_int_multiome.list[[celltype]] <- leadisoQTL_subset
  
}
leadisoQTL_int_multiome <- Reduce(rbind, leadisoQTL_int_multiome.list)

# Use plink is faster