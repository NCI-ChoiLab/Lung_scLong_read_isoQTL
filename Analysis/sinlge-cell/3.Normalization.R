# SCTransfrom and Normalization

# Bolun Li
# Feb 7 2024
.libPaths()
library(Seurat)
library(sctransform)
library(glmGamPoi)
lr.seur.afterqc.129 <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/lr.sc.data.afterQC.129.azimuth.annotated.rds")
lr.data.list <- SplitObject(lr.seur.afterqc.129,split.by = "orig.ident")
# run sctransform
for (i in 1:22) {
  seur <- lr.data.list[[i]]
  seur <- SCTransform(seur, vars.to.regress = "percent.mt", verbose = FALSE)
  lr.data.list[[i]] <- seur
}
lr.seur <- merge(lr.data.list[[1]], y = lr.data.list[2:22], project = "Long-read")
save(list = "lr.seur", file = "/data/Choi_lung/scLongreads/Seurat/lr.sc.data.afterQC.norm.RData")
