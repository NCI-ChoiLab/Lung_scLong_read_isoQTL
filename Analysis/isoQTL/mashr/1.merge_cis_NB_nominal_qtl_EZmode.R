# Merge cis nominal QTL results from jaxQTL
# Modified from Natri et al paper
# Bolun Li
# Jul 8th 2024

args <- commandArgs(trailingOnly = TRUE)
celltypes_path <- args[1]
sprintf(celltypes_path)

library(Seurat)
library(dplyr)
library(ggplot2)
library(future)
library(stringr)
library(stringi)
library(forcats)
library(arrow)
library(Matrix)

# read list of cell types 
celltypes_path <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Celltypes.txt"
celltypes <- read.table(file = celltypes_path, sep = "\t")
celltypes <- celltypes$V1

jax_path <- "/data/Choi_lung/scLongreads/mashr/"
dir.create(file.path(jax_path,"output"), showWarnings = FALSE)

for (celltype in celltypes) {
  celltype_for_path <- gsub(" ", "_", celltype)
  qtl_rst_all <- NB.nominal.list[[celltype_for_path]]
  qtl_rst_all$celltype <- celltype_for_path
  dir.create(file.path(jax_path,"output",celltype_for_path), showWarnings = FALSE)
  output_path <- file.path(jax_path,"output",celltype_for_path)
  write.table(qtl_rst_all, file = paste0(output_path,"/qtl_results_all.txt"), sep = "\t",
              quote = FALSE, row.names = FALSE)
}

library(dplyr)
library(tidyr)
library(vroom)
library(collapse)
library(rhdf5)

dirs <- list.dirs(file.path(jax_path,"output"), full.names = TRUE, recursive =FALSE)
dirs <- paste0(dirs, "/qtl_results_all.txt")
dirs <- dirs[file.exists(dirs)]

nCelltypes <- length(dirs)

all_qtl <- vroom(dirs, show_col_types = FALSE, id="celltype",
                 col_select=list(celltype="celltype", feature_id = "phenotype_id",
                                 variant_id = "snp", 
                                 betas = "slope", error = "slope_se", 
                                 pval_nonimal = "pval_nominal"))

all_qtl <- all_qtl %>%
  # drop duplicates that exist if non-biallelic variants were tested
  funique(cols=c("celltype", "feature_id", "variant_id")) %>% 
  fmutate(celltype = basename(dirname(celltype)),
          id = paste(feature_id, variant_id, sep="|"))
table(all_qtl$celltype)

table(all_qtl$celltype)

# calculate z score for EZ mode
all_qtl$zscore <- all_qtl$betas/all_qtl$error
features <- unique(all_qtl$feature_id)
ceiling(length(features)/20)
feature_chunks <- split(features, rep(1:20,ceiling(length(features)/20))[-42940])
beta.list <- list()
for (i in 1:20) {

  betas1 <- subset(all_qtl, feature_id %in% feature_chunks[[i]]) %>% pivot_wider(names_from = celltype, values_from = zscore,
                                                                                 id_cols = id) %>%
    tibble::column_to_rownames(var = "id") %>% qDF()

  message(ncol(betas1), " cell types included")
  tmp <- is.na(betas1)
  tmp <- rowSums(tmp)
  # remove pairs with NA in all cell types
  del_idx <- which(tmp == 33)
  if(length(del_idx) == 0){
    beta.list[[i]] <- betas1
  }else{
    betas1 <- betas1[-del_idx,]
    beta.list[[i]] <- betas1
  }
  
}

betas <- Reduce(rbind, beta.list)
error <- matrix(1,nrow = nrow(betas), ncol = ncol(betas))
error <- as.data.frame(error)
colnames(error) <- colnames(betas)
rownames(error) <- rownames(betas)
length(which(is.na(betas)))




message("Snapshot of merged data...")
n <- ifelse(ncol(betas) > 5, 5, ncol(betas))
message("betas:")
betas[1:5, 1:n]

message("beta standard error:")
error[1:5, 1:n]

# head(betas)

saveRDS(betas, file = "/data/Choi_lung/scLongreads/mashr/all.betas.zscore.rds")
saveRDS(error, file = "/data/Choi_lung/scLongreads/mashr/all.error.zscore.rds")
### Save output
message("Saving merged results as an hdf5...")
setwd("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/")
h5createFile("merge_p02.h5")
ncol <- ifelse(ncol(betas) >= 10, 10, ncol(betas))
nrow <- ifelse(nrow(betas) > 1e4, nrow(betas)/100, nrow(betas)/10)
h5createDataset(file = "merge_p02.h5", dataset = "betas", dims = dim(betas), 
                chunk = c(500, ncol))
h5createDataset(file = "merge_p02.h5", dataset = "error", dims = dim(error), 
                chunk = c(500, ncol))

h5write(as.matrix(betas), "merge_p02.h5", "betas")
h5write(as.matrix(error), "merge_p02.h5", "error")
h5write(rownames(betas), "merge_p02.h5", "rownames")
h5write(colnames(betas), "merge_p02.h5", "colnames")

