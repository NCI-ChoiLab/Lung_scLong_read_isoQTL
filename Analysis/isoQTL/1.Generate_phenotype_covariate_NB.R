# Generate phenotype files for isoform level QTL analysis using negative binomial regression

# Bolun Li
# Jun 24 2024


.libPaths()
library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)
library(future)
library(stringr)
library(stringi)
library(forcats)
library(Matrix)

# load transcript metadata generated from TALON and integrated with TranDecoder results
TALON_afterqc_orf_secondpass <- readRDS("/data/Choi_lung/scLongreads/Seurat/TALON_afterqc_orf_secondpass_w_known.rds")
# load the final version with cell type annotation
lr <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final_v2_wRNAexpr.RDS")

# prepare splicing QTL analysis
gff = as.data.frame(rtracklayer::import('/data/Choi_lung/scLongreads/TALON_workspace/test10s/NCI_lung_secondpass.isoforms_w_known.gtf'))
table(gff$type)
transcript_gff <- subset(gff, type == "transcript")
rm(gff)
pos <- data.frame(X = transcript_gff$transcript_id,
                  chr = transcript_gff$seqnames, 
                  start = transcript_gff$start, 
                  end = transcript_gff$start,
                  gene_id = transcript_gff$gene_id)
tss <- read.table("/data/Choi_lung/scLongreads/TALON_workspace/test10s/tss.bed", sep = "\t", header = FALSE)
rownames(tss) <- tss$V4
rownames(pos) <- pos$X

library("tidyverse")
`%nin%` <- Negate(`%in%`)
#########
# Covairates
# calculate genotype PCs
gpc <- read.table("/data/Choi_lung/scLongreads/tensorqtl/gtex-pipeline/genotype.maf05_geno01.pruned.pca.evec", 
                  sep = "\t", comment.char = "", header = FALSE)
gpc <- gpc[-1,]
tmp <- gsub("\\s+", ",", gpc)
gpc <- str_split_fixed(tmp, ",", n = 22)
gpc3 <- gpc[,2:4]
gpc3 <- apply(gpc3, 2, as.numeric)
rownames(gpc3) <- str_split_fixed(gpc[,1], ":", 2)[,2]
gpc3 <- t(gpc3)
rownames(gpc3) <- paste0("GPC", 1:3)

# comparing numbers of PFs affecting QTL analysis
generate_covariates <- function(covs_mtx, # covariates info mtx, C covariates * N samples 
                                genotype_pc, # 3 genotype pcs, GPC *N samples
                                peer_mtx, # peer factors N samples * K peer factors
                                K_used,   # Determine how many peer factors used in QTL analysis, could be a vector that does start with 0
                                samples){ # sample list involved into QTL analysis
  covs <- covs_mtx[,samples]
  if (K_used == 0) {
    rst <- rbind(covs, genotype_pc[,samples])
  }else{
    peer_used <- peer_mtx[samples,1:K_used]
    rst <- Reduce(rbind, list(covs, genotype_pc[,samples], t(peer_used)))
  }
  return(rst) # output the list of matrix as covariate input (C+K) covariates * N samples
}

# load batches
lr <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final_v2_wRNAexpr.RDS")
df <- as.data.frame(table(lr$orig.ident, lr$Sample))
df <- (subset(df, Freq != 0))
df <- distinct(df, Var2, .keep_all = TRUE)
df$Var1 <- as.factor(df$Var1)
df$Var1 <- as.numeric(df$Var1)
batch_covs <- t(df[,c(1,2)])
colnames(batch_covs) <- batch_covs[2,]
batch_covs <- batch_covs[-2,]

# load ages
ages <- read.table("/data/Choi_lung/scLongreads/tensorqtl/Age_cov.txt", sep = "\t",row.names = 1, header = TRUE)


jax_path <- "/data/Choi_lung/scLongreads/jaxqtl/"
lr.list.byCT <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr.list.byCT.rds")
celltypes <- levels(lr.list.byCT$`Secretory transitional cells`$Celltype)
# only keep the cell types satisfying the criteria of isoQTL mapping 
# (# of indvs with > 5 cells > 40)
celltypes <- celltypes[-c(9, 13, 25, 35)]
for(celltype in celltypes){
  # pilot analysis
  seur <- lr.list.byCT[[celltype]]
  celltype_for_path <- gsub(" ", "_", celltype)
  cell_number <- table(seur$Sample)
  # filter out indv with <= 5 cells
  ind_filtered <- names(cell_number)[which(cell_number < 5)]
  table(seur$Sample)
  seur <- subset(seur, Sample %nin% ind_filtered)

  # sum counts by indv
  count_mtx <- seur@assays$Isoform@counts
  count_mtx <- count_mtx[TALON_afterqc_orf_secondpass$transcript_name_unique,]
  count_mtx@Dimnames[[1]] <- TALON_afterqc_orf_secondpass$annot_transcript_id
  seur$Sample <- as.character(seur$Sample)
  group <- seur$Sample %>% fct_inorder()
  group_mat <- sparse.model.matrix(~ 0 + group) %>% t
  # Adjust row names to get the correct final row names
  rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
  count_mtx_T <- Matrix::t(count_mtx)
  count_mtx_sum <- group_mat %*% count_mtx_T  
  dim(count_mtx_sum)
  sum(count_mtx_sum)
  count_mtx_sum <- Matrix::t(count_mtx_sum)
  dim(count_mtx_sum)

  # Add isoform information
  pos <- data.frame(X = transcript_gff$transcript_id,
                    chr = transcript_gff$seqnames, 
                    start = transcript_gff$start, 
                    end = transcript_gff$start,
                    gene_id = transcript_gff$gene_id)
  length(unique(pos$X))
  
  rownames(pos) <- pos$X
  pos <- pos[TALON_afterqc_orf_secondpass$annot_transcript_id,]

  non_zero_mtx <- ((count_mtx_sum[idx,]) != 0)
  dim(non_zero_mtx)
  non_zero_transcripts <- rowSums(non_zero_mtx)
  idx_trancripts_.1 <- which(non_zero_transcripts > N*0.2)
  count_mtx_sum = count_mtx_sum[idx,]

  samples <- colnames(count_mtx_sum)
  transcripts <- rownames(count_mtx_sum)
  transcripts <- transcripts[idx_trancripts_.1]
  
  transcripts_ids <- rownames(tss)[which(rownames(tss) %in% transcripts)]
  colnames(tss) <- c("#chr", "start", "end", "trascript_id")
  tss_sub <- tss[transcripts_ids,c(1:4)]
  
  # generate bed file for jaxQTL
  bed <- cbind(tss_sub,count_mtx_sum[rownames(tss_sub),])
  message(nrow(bed), " isofroms in ", celltype)
  dir.create(file.path(jax_path, celltype_for_path), showWarnings = FALSE)
  sample_filter <- cbind(samples,samples)
  colnames(sample_filter) <- c("#FID", "IID")
  write.table(sample_filter, file = paste0(file.path(jax_path, celltype_for_path), "/samples.txt"), quote = F, row.names = F, col.names = TRUE)
  
  # Generate list for jaxQTL mapping
  chrs <- paste("chr", c(1:22,"X"), sep = "")
  
  for (chr in chrs) {
    idx <- which(bed$`#chr` == chr)
    transcripts <- bed$trascript_id[idx]
    write.table(transcripts, paste0(file.path(jax_path, celltype_for_path), "/", chr, "_chunck"), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  bed$`#chr` <- gsub("chr", "", bed$`#chr`)
  write.table(bed, file = paste0(file.path(jax_path, celltype_for_path), "/phenotype.bed"), sep='\t', quote=F, row.names=FALSE, col.names = T, eol = '\n')
  expr<-(bed[,-(1:4)]) 
  tmp <- (expr != 0)
  dim(tmp)

  # Generate expression PCs
  tmp1 <- rowSums(tmp)
  if (length(which(tmp1 <= N*0.2)) == 0) {
    prcompResult<-prcomp(t(expr),center=TRUE,scale.=TRUE) 
  }else{
    expr <- expr[-which(tmp1 <= N*0.2), ]
    prcompResult<-prcomp(t(expr),center=TRUE,scale.=TRUE) 
  }
  PCs<-prcompResult$x
  resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)
  print(resultRunElbow)

  # number of PCs is determined by elbow function
  K_used <- resultRunElbow
  age_used <- ages[,samples]
  
  # load expression PCs
  PFs <- PCs
  
  
  # combine covariates beside expression PCs
  covs_mtx<- rbind(age_used, batch_covs[colnames(age_used)])
  rownames(covs_mtx) <- c("Age", "batch")
  covs_mtx[2,] <- as.numeric(covs_mtx[2,])
  
  covs_final <- generate_covariates(covs_mtx = covs_mtx,
                                    genotype_pc = gpc3,
                                    peer_mtx = PFs,
                                    K_used = K_used,
                                    samples = samples)
  covs_final <- t(covs_final)
  write.table(data.frame("iid"=rownames(covs_final),covs_final), file =  paste0(file.path(jax_path, celltype_for_path), "/covariates_", K_used, "PCs.tsv"), 
              sep = "\t", quote = FALSE,row.names = F)
}

