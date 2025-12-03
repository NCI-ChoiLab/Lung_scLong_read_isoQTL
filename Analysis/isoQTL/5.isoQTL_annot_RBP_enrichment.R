# Prepare RBP motifs from eCLIP-seq (Van Nostrand et al 2020 Nature)
# For enrichment analysis

# Bolun Li
# Feb 27th, 2025

library(dplyr)
library(ggplot2)
library(stringi)
library(stringr)
source("/data/lib14/R/Rscripts/utilities.R")

setwd("/data/Choi_lung/scLongreads/Enrichment/")
# read downloaded bed files for RBP motifs
files <- list.files("/data/Choi_lung/scLongreads/Enrichment/", pattern = "qced.hg38.bed")
files <- list.files("/data/Choi_lung/scLongreads/Enrichment/", pattern = ".hg38.bed")
colnames <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
bed.list <- list()
bed.list1 <- list()
RBP_names <- NULL
for (file in files) {
  bed <- read.table(file, sep = "\t", header = FALSE)
  colnames(bed) <- colnames
  RBP_names <- c(RBP_names, unique(bed$name))
  bed.list1[[str_split_fixed(file, "\\.", 3)[,1]]] <- bed
}

names(bed.list) <- RBP_names
bed.list <- bed.list[-names(bed.list1)]
# QC with stringent cutoff
bed.qc.list <- lapply(bed.list, function(x){
  x$RBP <- str_split_fixed(x$name, "_", 3)[,1]
  x$cell_line <- str_split_fixed(x$name, "_", 3)[,2]
  x$sample_info <- str_split_fixed(x$name, "_", 3)[,3]
  x$peak <- paste(x$chr, x$start, x$end, sep = "-")
  # x <- subset(x, (signalValue >= 3) & (pValue >=5))
  return(x)
})

K562 <- bed.qc.list[grepl("K562", RBP_names)]
names(K562) <- str_split_fixed(names(K562), "_", 3)[,1]
K562_mtx <- do.call(rbind, K562)
HepG2 <- bed.qc.list[!grepl("K562", RBP_names)]
names(HepG2) <- str_split_fixed(names(HepG2), "_", 3)[,1]
HepG2_mtx <- do.call(rbind, HepG2)
tmp <- sapply(HepG2, function(x)nrow(x))
sum(tmp)
tmp2 <- lapply(HepG2, function(x)(x$peak))
tmp2 <- unique(Reduce(c,tmp2))


RBPs <- unique(str_split_fixed(RBP_names, "_", 3)[,1])
RBPs <- names(K562)
RBP_info <- read.csv("/data/Choi_lung/scLongreads/Enrichment/RBP_annotation.csv", header = TRUE)
rownames(RBP_info) <- RBP_info$RBP.name
RBP_info <- RBP_info[RBPs,c(1:21)]
Exprement_id <- read.table("/data/Choi_lung/scLongreads/Enrichment/ENCSR456FVU_metadata.tsv", header = TRUE, sep = "\t")
Exprement_id <- subset(Exprement_id, File.type == "bed")
Exprement_id <- Exprement_id[c(1:223)*3,]
# write.table(Exprement_id$File.accession, "/data/Choi_lung/scLongreads/Enrichment/IDR.list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# load all common variants in sc-long-read population 
library(snpStats)
fam <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.fam"
bim <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bim"
bed <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bed"
sample <- read.plink(bed, bim, fam)

# Overlap all isoQTLs with RBP motif
NB.nominal.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_sig_list.rds")
isoQTL.list <- lapply(NB.nominal.sig.list, function(x) unique(x$snp))
isoQTLs <- unique(do.call(c, isoQTL.list))
snp_info <- sample$map[isoQTLs,c(1,2,4,5,6)]
colnames(snp_info) <- c("chr", "rsid", "pos", "alt", "ref")
isoQTL_win_RBP_motif_peaks <- lapply(isoQTLs, function(snp){
  cCRE <- K562_mtx$peak[which(K562_mtx$chr == snp_info[snp,]$chr & 
                                as.numeric(K562_mtx$start) <= as.numeric(snp_info[snp,]$pos) & 
                                as.numeric(K562_mtx$end) >= as.numeric(snp_info[snp,]$pos))]
  RBP <- K562_mtx$RBP[which(K562_mtx$chr == snp_info[snp,]$chr & 
                              as.numeric(K562_mtx$start) <= as.numeric(snp_info[snp,]$pos) & 
                              as.numeric(K562_mtx$end) >= as.numeric(snp_info[snp,]$pos))]
  return(paste(cCRE,RBP, sep = "_"))
})
names(isoQTL_win_RBP_motif_peaks) <- isoQTLs
save(list = c("snp_info", "isoQTL_win_RBP_motif_peaks"), file = "/data/Choi_lung/scLongreads/Enrichment/isoQTL_RBP_motif.RData")
save(list = c("K562_mtx", "K562"), file = "/data/Choi_lung/scLongreads/Enrichment/K562_qced.RData")



load("/data/Choi_lung/scLongreads/Enrichment/isoQTL_RBP_motif.RData")
load("/data/Choi_lung/scLongreads/Enrichment/K562_qced.RData")

isoQTLs <- names(isoQTL_win_RBP_motif_peaks)
test.list <- lapply(isoQTLs, function(isoQTL){
  if (length(isoQTL_win_RBP_motif_peaks[[isoQTL]])==0) {
    return(isoQTL_win_RBP_motif_peaks[[isoQTL]])
  }else{
    x <- isoQTL_win_RBP_motif_peaks[[isoQTL]]
    df <- data.frame(RBP_motif = str_split_fixed(x, "_", 2)[,1],
                     RBP = str_split_fixed(x, "_", 2)[,2],
                     snp = isoQTL)
    return(df)
  }
})
test <- do.call(rbind, test.list)


# extract all snps located in RBP motif regions

K562_RBP_region <-K562_mtx[, c(1,2,3,11)]
table(K562_RBP_region$chr)
K562_RBP_region <- K562_RBP_region[!(K562_RBP_region$chr %in% c("chr1_KI270706v1_random","chr14_GL000009v2_random", "chr22_KI270879v1_alt","chr4_GL000008v2_random", "chrUn_KI270742v1")),]
write.table(K562_RBP_region, file = "/data/Choi_lung/scLongreads/Enrichment/K562/K562_RBP_region.txt",
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

# use plink to get all snps within RBP motif regions
arg <- paste0("/usr/local/apps/plink/1.9.0-beta4.4/plink --bfile /data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype --extract ./K562/K562_RBP_region.txt --range --make-bed --out K562_RBP_motif --write-snplist")

system(arg)
RBP_bg_snps <- read.table("/data/Choi_lung/scLongreads/Enrichment/K562_RBP_motif.snplist")
library(snpStats)
fam <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.fam"
bim <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bim"
bed <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bed"
sample <- read.plink(bed, bim, fam)

snp_info <- sample$map[RBP_bg_snps$V1,c(1,2,4,5,6)]
colnames(snp_info) <- c("chr", "rsid", "pos", "alt", "ref")

# generate background snp list for RBP motifs
RBP_bg_snps_peaks.list <- list()
RBP_bg_snps_peaks.list <- lapply(RBP_bg_snps$V1, function(snp){
  cCRE <- K562_mtx$peak[which(K562_mtx$chr == snp_info[snp,]$chr & 
                                as.numeric(K562_mtx$start) <= as.numeric(snp_info[snp,]$pos) & 
                                as.numeric(K562_mtx$end) >= as.numeric(snp_info[snp,]$pos))]
  RBP <- K562_mtx$RBP[which(K562_mtx$chr == snp_info[snp,]$chr & 
                              as.numeric(K562_mtx$start) <= as.numeric(snp_info[snp,]$pos) & 
                              as.numeric(K562_mtx$end) >= as.numeric(snp_info[snp,]$pos))]
  return(paste(cCRE,RBP, sep = "_"))
})
names(RBP_bg_snps_peaks.list) <- RBP_bg_snps$V1
RBP_bg_snps_peaks_df.list <- lapply(RBP_bg_snps$V1, function(snp){
  if (length(RBP_bg_snps_peaks.list[[snp]])==0) {
    return(RBP_bg_snps_peaks.list[[snp]])
  }else{
    x <- RBP_bg_snps_peaks.list[[snp]]
    df <- data.frame(RBP_motif = str_split_fixed(x, "_", 2)[,1],
                     RBP = str_split_fixed(x, "_", 2)[,2],
                     snp = snp)
    return(df)
  }
})
RBP_bg_snps_peaks_df <- do.call(rbind, RBP_bg_snps_peaks_df.list)
RBP_bg_sum <- as.data.frame(table(RBP_bg_snps_peaks_df$RBP))
RBP_isoQTL_sum <- as.data.frame(table(test$RBP))

load("/data/Choi_lung/scLongreads/Enrichment/K562_enrichment_rst.RData")

# load LD pruned SNPs
NB.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")
snp.pruned <- read.table("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/prune.prune.in")

RBP_bg_snps_peaks_df_pruned_cisSNP <- subset(RBP_bg_snps_peaks_df, snp %in% snp.pruned$V1)
test_pruned <- subset(test, snp %in% snp.pruned$V1)
isoQTL_pruned <- isoQTLs[which(isoQTLs%in%snp.pruned$V1)]
functions <- RBP_info[,c(3:20)]
func_category <- colnames(functions)
enrich_tb.list <- list()
for (func in func_category) {
  func_RBPs <- rownames(functions)[(functions[[func]]==1)]
  enriched_tb <- subset(test_pruned, RBP %in% func_RBPs)
  if(nrow(enriched_tb) == 0){
    enrich_tb.list[[func]] <- NULL
  }else{
    enriched_tb$functional <- func
    enrich_tb.list[[func]] <- enriched_tb
  }
  
}
enrich_tb <- Reduce(rbind, enrich_tb.list)

bg.list <- list()
for (func in func_category) {
  func_RBPs <- rownames(functions)[(functions[[func]]==1)]
  bg_tb <- subset(RBP_bg_snps_peaks_df_pruned_cisSNP, RBP %in% func_RBPs)
  if(nrow(bg_tb) == 0){
    bg.list[[func]] <- NULL
  }else{
    bg_tb$functional <- func
    bg.list[[func]] <- bg_tb
  }
  
}
bg_tb <- Reduce(rbind, bg.list)
RBP_bg_sum <- as.data.frame(table(bg_tb$functional))
RBP_isoQTL_sum <- as.data.frame(table(enrich_tb$functional))

RBP_isoQTL_enrich_tb <- left_join(RBP_bg_sum, RBP_isoQTL_sum, by = "Var1")
colnames(RBP_isoQTL_enrich_tb) <- c("Funtional", "background", "isoQTL")
RBP_isoQTL_enrich_tb$non_isoQTL <- (RBP_isoQTL_enrich_tb$background - RBP_isoQTL_enrich_tb$isoQTL)
RBP_isoQTL_enrich_tb[is.na(RBP_isoQTL_enrich_tb)] <- 0

# Test enrichment by fishr's exact test
RBP_isoQTL_enrich_tb <- RBP_isoQTL_enrich_tb %>% group_by(Funtional) %>% mutate(pval = fisher.test(matrix(c((nrow(snp.pruned)-length(isoQTL_pruned)), length(isoQTL_pruned), non_isoQTL, isoQTL),
                                                                                                    nrow = 2, byrow = TRUE))$p.value,
                                                                          logOR = log(fisher.test(matrix(c((nrow(snp.pruned)-length(isoQTL_pruned)), length(isoQTL_pruned), non_isoQTL, isoQTL),
                                                                                                         nrow = 2, byrow = TRUE))$estimate))
RBP_isoQTL_enrich_tb$FDR <- p.adjust(RBP_isoQTL_enrich_tb$pval, method = "BH" )
write.table(RBP_isoQTL_enrich_tb, "/data/Choi_lung/scLongreads/Enrichment/K562/K562_functional_annotation_enrichment_pruned_rst.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
