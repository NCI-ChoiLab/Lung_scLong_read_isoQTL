# Azimuth annotation

# Bolun Li
# Feb 8 2024


.libPaths()
library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)
library(future)
library(stringr)
library(stringi)
if (!requireNamespace('remotes', quietly = TRUE)) {
  install.packages('remotes')
}
remotes::install_github('satijalab/azimuth', ref = 'master')
rm(list = ls())
gc()

library(Azimuth)

load("/data/Choi_lung/scLongreads/Seurat/lr.sc.data.afterQC.129.RData")
DefaultAssay(lr.seur.afterqc.129)
lr.seur.afterqc.129 <- RunAzimuth(lr.seur.afterqc.129, reference = "lungref")

saveRDS(lr.seur.afterqc.129, file = "/data/Choi_lung/scLongreads/Seurat/lr.sc.data.afterQC.129.azimuth.annotated.rds")
lr.seur.afterqc.129 <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/lr.sc.data.afterQC.129.azimuth.annotated.rds")
table(lr.seur.afterqc.129$predicted.ann_finest_level)
sr.meta <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/annot.meta.rds")
lr.seur.afterqc.129@assays
annot.meta <- left_join(lr.seur.afterqc.129@meta.data, sr.meta ,by = "barcode")

lr.seur.afterqc.129 <- NormalizeData(lr.seur.afterqc.129, verbose = FALSE)
lr.seur.afterqc.129$predicted.ann_level_5
DotPlot(lr.seur.afterqc.129, features = c("EPCAM", "PTPRC", "PECAM1", "COL1A1"), group.by = "predicted.ann_finest_level")
DotPlot(lr.seur.afterqc.129, features = c("EPCAM", "PTPRC", "PECAM1", "COL1A1"), group.by = "predicted.ann_level_3")
DotPlot(lr.seur.afterqc.129, features = c("EPCAM", "PTPRC", "PECAM1", "COL1A1"), group.by = "predicted.ann_level_2")


load("/data/Choi_lung/scLongreads/Seurat/lr.sc.data.afterQC.norm.RData")

