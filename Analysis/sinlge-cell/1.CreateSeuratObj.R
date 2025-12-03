# Load PacBio sc-data to Seurat
# Run scrublet for detecting doublets

# Bolun Li
# Jan 11 2024


.libPaths()
library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)
library(future)
library(stringr)
library(stringi)

# Scrublet Python to R
library(reticulate)
scr <- import("scrublet", convert = FALSE)
py_run_string('import scrublet as scr')

file_dir <- "/data/Choi_lung/scLongreads/B2_percentile/"
sample_dir <- list.dirs(file_dir, recursive = FALSE)
name <- list.dirs(file_dir, recursive = FALSE, full.names = FALSE)
name <- str_split_fixed(name, pattern = "_", n = 2)[,2]
name <- name[-23]
dataset_list <- list()
i = 22
for (i in 1:length(name)) {
  filepath = paste(file_dir, paste("/Sample_",name[i],"/", name[i], "_classification_filtered", sep = ""), "/genes_seurat", sep = "")
  filepath_isoforms = paste(file_dir, paste("/Sample_",name[i],"/", name[i], "_classification_filtered", sep = ""), "/isoforms_seurat", sep = "")
  # doubletpath = paste(doublet_dir , name[i], sep = "")
  # Seurat
  seur.data = Read10X(data.dir = filepath)
  seur.isoforms.data = Read10X(data.dir = filepath_isoforms)
  dim(seur.data)
  dim(seur.isoforms.data)
  seur = CreateSeuratObject(counts = seur.data, project = name[i], min.cells=5)
  print(ncol(seur))
  head(colnames(seur))
  seur[["Isoforms"]] <- CreateAssayObject(counts = seur.isoforms.data, min.cells=5)
  # save original barcodes from PacBio
  seur$orig.barcode <- colnames(seur)
  # reverse complementary DNA of barcodes
  CB.reverse <- chartr(old="ATGC", new="TACG", colnames(seur))
  CB.reverse <- stri_reverse(CB.reverse)
  CB.lr <- paste(str_split_fixed(CB.reverse, "-", n =2)[,2], "-1", sep = "")
  # Rename the cells using the reversed barcodes
  seur <- RenameCells(seur, new.names = CB.lr)
  # doublet = read.table(paste(doubletpath, '/predicted_doublets.txt', sep = ""))
  # doublet_score = read.table(paste(doubletpath, '/doublet_scores.txt', sep = ""))
  # seur[['doublet_score']] = doublet_score$V1
  # seur[['doublet']] = doublet$V1
  # seur[['Batch']] = Batch[i]
  # detect doublet by Scrublet
  ## Gene level
  scrub <- scr$Scrublet(t(as.matrix(seur@assays$RNA@counts)))
  doublet_scores = scrub$scrub_doublets()
  doublet_scores = py_to_r(doublet_scores)
  seur$doublet_RNA = doublet_scores[[1]]
  seur$doublet2_RNA = doublet_scores[[2]]
  ## Isoform level
  ptm <- proc.time()
  scrub <- scr$Scrublet(t(as.matrix(seur@assays$Isoforms@counts))) # 40GB
  doublet_scores = scrub$scrub_doublets()
  doublet_scores = py_to_r(doublet_scores)
  seur$doublet_isoform = doublet_scores[[1]]
  seur$doublet2_isoform = doublet_scores[[2]]
  print(proc.time() - ptm)
  # seur <- subset(seur, doublet == 0 & nFeature_RNA < 5000)
  seur[["percent.mt"]] = PercentageFeatureSet(seur, pattern = "^MT-")
  seur[['percent.rb']] = PercentageFeatureSet(seur, pattern = '^RP[SL]')
  dataset_list[[i]] <- seur
}
names(dataset_list) <- name
rm("seur.data", "seur.isoforms.data")
save(list="dataset_list", file = "/data/Choi_lung/scLongreads/Seurat/lr.sc.data.list.21samples.beforeQC.RData")
