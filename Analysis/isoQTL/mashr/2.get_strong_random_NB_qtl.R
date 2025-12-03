# Get the strong and random QTL signals using significant eQTL as strong
# Using Z score
# Bolun Li
# Jul 9th 2024

library(vroom)
library(dplyr)
library(data.table)
library(collapse)
library(rhdf5)


setwd("/data/Choi_lung/scLongreads/mashr/")

betas <-readRDS("all.betas.zscore.rds")
error <-readRDS("all.error.zscore.rds")
min(error)
rownames <- rownames(betas)
colnames <- colnames(betas)
isoQTL_list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_list.rds")
sig_tests <- Reduce(rbind, isoQTL_list)
sig_tests$id <- paste(sig_tests$phenotype_id, sig_tests$variant_id, sep="|")
sig_pairs <- unique(sig_tests$id)
message("Saving ", length(intersect(rownames,sig_pairs)), " strong QTL.")
betas_strong <- betas[intersect(rownames,sig_tests$id), ]
betas_strong <- betas[sig_pairs, ]
error_strong <- error[intersect(rownames,sig_tests$id), ]
error_strong <- error[sig_pairs, ]
save_strong = "strong_all_ez.h5"
h5createFile(save_strong)
ncol <- ifelse(ncol(betas_strong) >= 10, 10, ncol(betas_strong))
nrow <- ifelse(nrow(betas_strong) > 1e4, nrow(betas_strong)/100, nrow(betas_strong)/10)
h5createDataset(file = save_strong, dataset = "betas", dims = dim(betas_strong), 
                chunk = c(1000, ncol))
h5createDataset(file = save_strong, dataset = "error", dims = dim(error_strong), 
                chunk = c(1000, ncol))

h5write(as.matrix(betas_strong), save_strong, "betas")
h5write(as.matrix(error_strong), save_strong, "error")
h5write(rownames(betas_strong), save_strong, "rownames")
h5write(colnames(betas_strong), save_strong, "colnames")


seed <- sample(1:100, 1)
set.seed(seed)
random <- sample(1:nrow(betas), 1000000)
betas_random <- betas[random, ]
error_random <- error[random, ]
rownames(betas_random)[1:5]
save_random <- "random_all_1M_seed18_ez.h5"
h5createFile(save_random)
ncol <- ifelse(ncol(betas_random) >= 10, 10, ncol(betas_random))
nrow <- ifelse(nrow(betas_random) > 1e4, nrow(betas_random)/100, nrow(betas_random)/10)
h5createDataset(file = save_random, dataset = "betas", dims = dim(betas_random), 
                chunk = c(1000, ncol))
h5createDataset(file = save_random, dataset = "error", dims = dim(error_random), 
                chunk = c(1000, ncol))

h5write(as.matrix(betas_random), save_random, "betas")
h5write(as.matrix(error_random), save_random, "error")
h5write(rownames(betas_random), save_random, "rownames")
h5write(colnames(betas_random), save_random, "colnames")

message("Done!")



