# To compare isoQTL and eQTL at gene level
# Colocalization between isoQTL and eQTL
# Bolun
# Apr 9th 

library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)
library(future)
library(stringr)
library(stringi)
library(Matrix)
library(forcats)
source("/data/lib14/R/Rscripts/utilities.R")
load("/data/Choi_lung/scLongreads/eisoQTL.RData")
load("/data/Choi_lung/scLongreads/eisoQTL_ver2.RData")
setwd("/data/Choi_lung/scLongreads/eQTL_comparison/")
FT_isoQTL <- read.table("/data/lib14/project/scLongread/Tables/NB_eIsoform_sum.txt", sep = "\t",
                        header = TRUE)

full_sig_list <- readRDS("/vf/users/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/full_sig_list.rds")

for (celltype in names(full_sig_list)) {
  full_sig_list[[celltype]]$Celltype <- celltype
}
FT_eQTL <- Reduce(rbind, full_sig_list)
length(unique(FT_eQTL$phenotype_id))
length(unique(intersect(FT_isoQTL$gene_id, FT_eQTL$phenotype_id)))
shared_eisoGene <- unique(intersect(FT_isoQTL$gene_id, FT_eQTL$phenotype_id))
# write.table(unique(FT_eQTL$Celltype), "eQTL_celltype.txt", row.names = FALSE, quote = FALSE,
#             col.names = FALSE)
# write.table(unique(FT_isoQTL$Celltype), "isoQTL_celltype.txt", row.names = FALSE, quote = FALSE,
#             col.names = FALSE)

495/2229 # 22.2% shared at gene level

CT_match <- read.table("/data/Choi_lung/scLongreads/eQTL_comparison/eQTL_celltype.txt", sep = "\t")
CT_match1 <- read.table("/data/Choi_lung/scLongreads/eQTL_comparison/isoQTL_celltype.txt", sep = "\t")
# FT_eQTL_ct_matched <- subset(FT_eQTL, Celltype %in% CT_match$V1)
# length(unique(FT_eQTL_ct_matched$phenotype_id))

FT_eQTL$Celltype_matching <- FT_eQTL$Celltype
for(i in 1:nrow(CT_match)){
  FT_eQTL$Celltype_matching = replace(FT_eQTL$Celltype_matching, FT_eQTL$Celltype_matching == CT_match$V1[i], CT_match$V2[i])
}
FT_isoQTL$Celltype_matching <- FT_isoQTL$Celltype
for(i in 1:nrow(CT_match1)){
  FT_isoQTL$Celltype_matching = replace(FT_isoQTL$Celltype_matching, FT_isoQTL$Celltype_matching == CT_match1$V1[i], CT_match1$V2[i])
}
FT_isoQTL_shared <- subset(FT_isoQTL, gene_id %in% shared_eisoGene)
FT_eQTL_shared <- subset(FT_eQTL, phenotype_id %in% shared_eisoGene)

FT_isoQTL_shared$gene_ct <- paste(FT_isoQTL_shared$gene_id, FT_isoQTL_shared$Celltype_matching, sep = "-")
FT_eQTL_shared$gene_ct <- paste(FT_eQTL_shared$phenotype_id, FT_eQTL_shared$Celltype_matching, sep = "-")
length(unique(FT_isoQTL_shared$gene_ct))
length(unique(FT_eQTL_shared$gene_ct))
length(unique(intersect(unique(FT_isoQTL_shared$gene_ct),unique(FT_eQTL_shared$gene_ct))))
gene_ct_shared <- intersect(unique(FT_isoQTL_shared$gene_ct),unique(FT_eQTL_shared$gene_ct))
unique(str_split_fixed(gene_ct_shared, "-", 2)[,1]) # 293
isoQTL_unique <- setdiff(unique(FT_isoQTL_shared$gene_ct),unique(FT_eQTL_shared$gene_ct))
unique(str_split_fixed(isoQTL_unique, "-", 2)[,1]) # 307
shared_diff_ct <- setdiff(unique(str_split_fixed(isoQTL_unique, "-", 2)[,1]),
          unique(str_split_fixed(gene_ct_shared, "-", 2)[,1]))


FT_isoQTL_shared_diff_ct <- subset(FT_isoQTL_shared, gene_id %in% shared_diff_ct)
FT_eQTL_shared_diff_ct <- subset(FT_eQTL_shared, phenotype_id %in% shared_diff_ct)
load("/vf/users/Choi_lung/scLongreads/DEI/DEI_ct_specific_eIsofrom_GTEx_compare.RData")
ct_spec_isoforms
tmp <- subset(FT_isoQTL_shared_diff_ct, phenotype_id %in% ct_spec_isoforms)
setdiff(unique(FT_isoQTL_shared_diff_ct$gene_id), unique(tmp$gene_id))

for (celltype in names(DEI_sig.list)) {
  DEI_sig.list[[celltype]]$Celltype <- celltype
}


DEI_sig_allCT <- Reduce(rbind, DEI_sig.list)
length(unique(DEI_sig_allCT$X)) # 19,595
which(FT_isoQTL_shared_diff_ct$transcript_name %in% unique(DEI_sig_allCT$X))

FT_isoQTL_shared_diff_ct$DEI_matching <- paste(FT_isoQTL_shared_diff_ct$transcript_name, 
                                               FT_isoQTL_shared_diff_ct$Celltype, sep = ":")
DEI_sig_allCT$DEI_matching <- paste(DEI_sig_allCT$X, DEI_sig_allCT$Celltype, sep = ":")
# 67 eIsoforms are cell type matching DEI
CT_matching_isoforms <- unique(str_split_fixed(intersect(FT_isoQTL_shared_diff_ct$DEI_matching, DEI_sig_allCT$DEI_matching), ":", 2)[,1])
tmp <- FT_isoQTL_shared_diff_ct[which(FT_isoQTL_shared_diff_ct$transcript_name %in% CT_matching_isoforms),]
DEI_tmp <- subset(DEI_sig_allCT, X %in% unique(str_split_fixed(intersect(FT_isoQTL_shared_diff_ct$DEI_matching, DEI_sig_allCT$DEI_matching), ":", 2)[,1]))
DEI_tmp <- subset(DEI_sig_allCT, DEI_matching %in% intersect(FT_isoQTL_shared_diff_ct$DEI_matching, DEI_sig_allCT$DEI_matching))

FT_isoQTL_shared_same_ct <- subset(FT_isoQTL_shared,  gene_ct %in% gene_ct_shared)
FT_eQTL_shared_same_ct <- subset(FT_eQTL_shared, gene_ct %in% gene_ct_shared)
FT_eQTL_shared_same_ct_for_join <- FT_eQTL_shared_same_ct[!duplicated(FT_eQTL_shared_same_ct$gene_ct),]

FT_isoQTL_shared_same_ct_joined <- left_join(FT_isoQTL_shared_same_ct, FT_eQTL_shared_same_ct_for_join, by = "gene_ct")


isoform_list <- unique(FT_isoQTL_shared_same_ct_joined$phenotype_id.x)
# read wald test isoQTL results
file_paths <- list.files("/data/Choi_lung/scLongreads/jaxqtl/")
idx <- grep("swarm", file_paths)
file_paths <- file_paths[-idx]
file_paths <- file_paths[-c(9,12,25:29,37,39, 40,42,43,46)]

NB.coloc.list <- list()
for (ct in file_paths) {
  files <- list.files(paste0("/data/Choi_lung/scLongreads/jaxqtl/",ct), 
                      pattern = "wald.parquet")
  if(length(files) == 23){
    rst_1k.list <- list()
    for (i in c(1:length(files))) {
      file_path <- paste0("/data/Choi_lung/scLongreads/jaxqtl/", ct, "/", files[i])
      rst <- read_parquet(file_path)
      rst_1k.list[[i]] <- rst
    }
    rst_1k <- do.call(rbind, rst_1k.list)
    rst_sub <- subset(rst_1k, phenotype_id %in% isoform_list)
    NB.coloc.list[[ct]] <- rst_sub
  }else{
    message(ct, " jaxqtl not finished")
    NB.coloc.list[[ct]] <- NULL
  }
}

saveRDS(NB.coloc.list, "/data/Choi_lung/scLongreads/colocalization/isoeQTL_coloc/NB.coloc.list.rds")
NB.coloc.list <- readRDS("/data/Choi_lung/scLongreads/colocalization/isoeQTL_coloc/NB.coloc.list.rds")
nominalfiles <- list.files("/data/Choi_lung/TTL/Files_for_Bolun/Short_read_results_for_RShiny/tensor_nominal/")
celltypes <- str_split_fixed(nominalfiles, pattern = "\\.", 2)[,1]
celltypes <- gsub("_nominal", "",celltypes)

names(nominalfiles) <- celltypes
eQTL.coloc.list <- list()
for (celltype in celltypes) {
  nominalfile_path <- paste0("/data/Choi_lung/TTL/Files_for_Bolun/Short_read_results_for_RShiny/tensor_nominal/",nominalfiles[celltype])
  final_table <- read.table(nominalfile_path, sep = "\t", header = TRUE)
  final_table <- subset(final_table, phenotype_id %in% unique(FT_isoQTL_shared_same_ct_joined$gene_id))
  eQTL.coloc.list[[celltype]] <- final_table
}

saveRDS(eQTL.coloc.list, "/data/Choi_lung/scLongreads/colocalization/isoeQTL_coloc/eQTL.coloc.list.rds")
eQTL.coloc.list <- readRDS("/data/Choi_lung/scLongreads/colocalization/isoeQTL_coloc/eQTL.coloc.list.rds")
eQTL.coloc.list <- eQTL.coloc.list[CT_match$V1[-21]]
names(eQTL.coloc.list) <- CT_match$V2[-21]
unique(FT_isoQTL_shared_same_ct_joined$gene_ct)
df <- str_split_fixed(unique(FT_isoQTL_shared_same_ct_joined$gene_ct), "-",2)
length(unique(df[,2]))



library(coloc)
library(snowfall)
library(foreach)
library(doParallel)
library(parallel)
rst.list <- list()
output.list <- list()
# colocalization of isoQTLs and eQTLs for shared eGene/eIsoform in same cell type
for (i in 1:nrow(df)) {
  gene <- df[i,1]
  celltype <- df[i,2]
  isoforms <- FT_isoQTL_shared_same_ct_joined$phenotype_id.x[which(FT_isoQTL_shared_same_ct_joined$gene_id == gene)]
  isoQTLs <- NB.coloc.list[[celltype]]
  filt_isoQTLs <- subset(isoQTLs, phenotype_id %in% isoforms)
  filt_isoQTLs <- filt_isoQTLs[!is.na(filt_isoQTLs$pval_nominal),] # remove non-converging results
  filt_isoQTLs$rsid = paste(paste0("chr",filt_isoQTLs$chrom), filt_isoQTLs$pos, filt_isoQTLs$a1,filt_isoQTLs$a0,sep = '_')
  filt_isoQTLs = filt_isoQTLs[,c('phenotype_id','rsid','chrom','af','pval_nominal','pos',"slope","slope_se")]
  
  eQTLs <- eQTL.coloc.list[[celltype]]
  filt_eQTLs <- subset(eQTLs, phenotype_id == gene)
  # modify eQTL tables
  filt_eQTLs = filt_eQTLs[,c('phenotype_id','variant_id','phenotype_chr','af','pval_nominal','variant_pos',"slope","slope_se")]
  colnames(filt_eQTLs) <- c('phenotype_id','rsid','chrom','af','pval_nominal','pos',"slope","slope_se")
  
  
  snps <- intersect(filt_isoQTLs$rsid, filt_eQTLs$rsid)
  filt_eQTLs=filt_eQTLs[filt_eQTLs$rsid %in% snps,]
  filt_isoQTLs=filt_isoQTLs[filt_isoQTLs$rsid %in% snps,]
  length(filt_eQTLs$rsid)
  length(unique(filt_isoQTLs$rsid))
  gene_list=as.character(unique(filt_isoQTLs$phenotype_id))
  isoqtl_isoforms=list()
  isoqtl_isoforms=foreach(i=gene_list) %dopar%
    droplevels(filt_isoQTLs[filt_isoQTLs$phenotype_id==i,]) # Dropping unused levels makes a HUGE difference in time and final list size
  names(isoqtl_isoforms)=gene_list
  # change later
  eqtl_isoforms=list()
  eqtl_isoforms=foreach(i=gene_list) %dopar%
    droplevels(filt_eQTLs[na.omit(match(as.character(isoqtl_isoforms[[i]]$rsid),as.character(filt_eQTLs$rsid))),])
  names(eqtl_isoforms)=gene_list

  ct_path <- paste0("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/", celltype, "/samples.txt")
  samples <- read.table(ct_path, header = FALSE)
  sample_size <- nrow(samples)
  
  
  SNP_coloc=list()
  SNP_coloc=foreach(i=(1:length(isoqtl_isoforms)),.packages='coloc') %dopar%
    coloc.abf(dataset1 = list(snp=isoqtl_isoforms[[i]]$rsid,pvalues=isoqtl_isoforms[[i]]$pval_nominal,beta=isoqtl_isoforms[[i]]$slope,varbeta=(isoqtl_isoforms[[i]]$slope_se)^2,position=isoqtl_isoforms[[i]]$pos,N=sample_size,type="quant",MAF=isoqtl_isoforms[[i]]$af),
              dataset2 = list(snp=eqtl_isoforms[[i]]$rsid,pvalues=eqtl_isoforms[[i]]$pval_nominal,beta=eqtl_isoforms[[i]]$slope,varbeta=(eqtl_isoforms[[i]]$slope_se)^2,position=eqtl_isoforms[[i]]$pos,N=sample_size,type="quant",MAF=eqtl_isoforms[[i]]$af))
  names(SNP_coloc)=gene_list
  rst.list[[i]] <- SNP_coloc
  output_genes=list()
  output_genes=foreach(i=gene_list) %dopar%
    capture.output(SNP_coloc[[i]][["summary"]])
  names(output_genes)=gene_list
  output_sum <- lapply(1:length(output_genes), function(x){
    rst <- as.numeric(str_split_fixed(output_genes[[x]][2], " ", 6))
  })
  names(output_sum)=gene_list
  if(length(output_sum) == 1){
    tmp <- Reduce(rbind, output_sum)
    output <- data.frame(nsnps = tmp[1], PP.H0.abf = tmp[2], PP.H1.abf = tmp[3], PP.H2.abf = tmp[4], PP.H3.abf = tmp[5], PP.H4.abf = tmp[6])
    rownames(output) <- gene_list
    output$transcript_id <- rownames(output)
  }else{
    output <- Reduce(rbind, output_sum)
    colnames(output) <- c("nsnps", "PP.H0.abf", "PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")
    rownames(output) <- gene_list
    output <- as.data.frame(output)
    output$transcript_id <- rownames(output)
  }
  output$Celltype <- celltype
  output <- output %>% group_by(transcript_id) %>% mutate(transcript_name = Search_transcript_name2(transcript_id))
  output.list[[i]] <- output
}
names(output.list) <- unique(FT_isoQTL_shared_same_ct_joined$gene_ct)
output_final.list <- lapply(1:700, function(i){
  tmp <- as.data.frame(output.list[[i]])
  tmp$gene_ct <- unique(FT_isoQTL_shared_same_ct_joined$gene_ct)[i]
  return(tmp)
})
output_final <- Reduce(rbind, output_final.list)
save(list = c("rst.list", "output_final.list", "output_final"), file = "/data/Choi_lung/scLongreads/colocalization/isoeQTL_coloc/isoeQTL.coloc.rst.RData")
output_final$transcript_ct <- paste(output_final$transcript_id, output_final$Celltype, sep = "-")
output_final <- subset(output_final, transcript_id %in% unique(FT_isoQTL_shared_same_ct_joined$phenotype_id.x))
FT_isoQTL_shared_same_ct_joined$transcript_ct <- paste(FT_isoQTL_shared_same_ct_joined$phenotype_id.x, FT_isoQTL_shared_same_ct_joined$Celltype.x, sep = "-")
output_final <- subset(output_final, transcript_ct %in% unique(FT_isoQTL_shared_same_ct_joined$transcript_ct))
length(unique(output_final$gene_id[which(output_final$PP.H4.abf > 0.7)])) #194
length(unique(output_final$transcript_id[which(output_final$PP.H4.abf > 0.7)])) #194
length(unique(FT_isoQTL_shared_same_ct_joined$transcript_ct))
FT_isoQTL_shared_same_ct_joined <- left_join(FT_isoQTL_shared_same_ct_joined, output_final[,c(6,12)], by = "transcript_ct")

write.table(FT_isoQTL_shared_same_ct_joined, file = "/data/Choi_lung/scLongreads/eQTL_comparison/isoeQTL_shared_same_ct_coloc_rst.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

FT_isoQTL_shared_same_ct_joined <- read.table("/data/Choi_lung/scLongreads/eQTL_comparison/isoeQTL_shared_same_ct_coloc_rst.txt",
                                              sep = "\t", header = TRUE)
gene = "ENSG00000185250"
celltype = "Multiciliated"
isoQTLs <- NB.coloc.list[[celltype]]
filt_isoQTLs <- subset(isoQTLs, phenotype_id %in% c("ENST00000521072"))
filt_isoQTLs <- filt_isoQTLs[!is.na(filt_isoQTLs$pval_nominal),] # remove non-converging results
filt_isoQTLs$rsid = paste(paste0("chr",filt_isoQTLs$chrom), filt_isoQTLs$pos, filt_isoQTLs$a1,filt_isoQTLs$a0,sep = '_')

filt_isoQTLs = filt_isoQTLs[,c('phenotype_id','rsid','snp','chrom','pos','a1', "a0", 'tss_distance','pval_nominal',"slope","slope_se")]
eQTLs <- eQTL.coloc.list[[celltype]]
filt_eQTLs <- subset(eQTLs, phenotype_id == gene)
# modify eQTL tables
filt_eQTLs = filt_eQTLs[,c('phenotype_id','variant_id','phenotype_chr','variant_pos', "start_distance", 'pval_nominal',"slope","slope_se")]
colnames(filt_eQTLs) <- c('phenotype_id','rsid','chrom','pos',"tss_distance",'pval_nominal',"slope","slope_se")
filt_eQTLs <- left_join(filt_eQTLs, filt_isoQTLs[,c('rsid','snp','a1', "a0")], by = "rsid")
filt_isoQTLs = filt_isoQTLs[,c('phenotype_id','snp','chrom','pos','a1', "a0", 'tss_distance','pval_nominal',"slope","slope_se")]
filt_eQTLs = filt_eQTLs[,c('phenotype_id','snp','chrom','pos','a1', "a0", 'tss_distance','pval_nominal',"slope","slope_se")]
snps <- intersect(filt_isoQTLs$snp, filt_eQTLs$snp)
filt_eQTLs=filt_eQTLs[filt_eQTLs$snp %in% snps,]
filt_isoQTLs=filt_isoQTLs[filt_isoQTLs$snp %in% snps,]
length(filt_eQTLs$snp)
length(unique(filt_isoQTLs$snp))
colnames(filt_eQTLs)[c(1,2,3,5,6)] <- c("gene_id","rsnum", "chr","ref", "alt")
filt_eQTLs$gene_symbol <- "PPIL6"
colnames(filt_isoQTLs)[c(2,3,5,6)] <- c("rsnum", "chr","ref", "alt")
filt_eQTLs$chr <- gsub("chr","", filt_eQTLs$chr)
filt_eQTLs$variant_id <- paste(filt_eQTLs$chr, filt_eQTLs$pos, sep = ":")
filt_eQTLs <- filt_eQTLs[,c(1,11,12,2:10)]

filt_isoQTLs$zscore <- filt_isoQTLs$slope/filt_isoQTLs$slope_se
filt_isoQTLs <- filt_isoQTLs[,c(3,2,4:6,8,11,9,10)]
filt_isoQTLs <- filt_isoQTLs[,c(1,3,4,5,2,6:9)]
colnames(filt_isoQTLs) <- c("chr", "pos","ref","alt", "rsnum","pvalue","zscore","effect", "se")
write.table(filt_eQTLs, file = "/data/Choi_lung/scLongreads/eQTL_comparison/PPIL6_eQTL.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")
write.table(filt_isoQTLs, file = "/data/Choi_lung/scLongreads/eQTL_comparison/PPIL6_isoQTL.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")


# install.packages("devtools")
devtools::install_github("boxiangliu/locuscomparer")

library(locuscomparer)
isoQTL_fn <- filt_isoQTLs[,c("rsnum","pvalue")]
eQTL_fn <- filt_eQTLs[,c("rsnum","pval_nominal")]
gwas <- readRDS("/data/Choi_lung/scLongreads/colocalization/Byun_total_meta_gwas_31loci_1MB.rds")
gwas_loci <- gwas$`6q21`
gwas_fn <- gwas_loci[,c("variant_id","p_value")]
colnames(gwas_fn) <- c("rsid", "pval")
colnames(isoQTL_fn) <- c("rsid", "pval")
colnames(eQTL_fn) <- c("rsid", "pval")
snps_6q21 <- intersect(gwas_fn$rsid, isoQTL_fn$rsid)
gwas_fn <- subset(gwas_fn, rsid %in% snps_6q21)
isoQTL_fn <- subset(isoQTL_fn, rsid %in% snps_6q21)
eQTL_fn <- subset(eQTL_fn, rsid %in% snps_6q21)
write.table(isoQTL_fn, file = "/data/Choi_lung/scLongreads/eQTL_comparison/isoQTL_fn.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(eQTL_fn, file = "/data/Choi_lung/scLongreads/eQTL_comparison/eQTL_fn.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(gwas_fn, file = "/data/Choi_lung/scLongreads/eQTL_comparison/gwas_fn.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
locuscompare(in_fn1 = isoQTL_fn, in_fn2 = eQTL_fn, title = 'isoQTL', title2 = 'eQTL')

