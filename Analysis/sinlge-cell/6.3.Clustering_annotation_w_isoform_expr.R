# SCTransfrom and Normalization

# Bolun Li
# Sep 7 2025
.libPaths()
library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)
library(future)
library(reticulate)
use_condaenv("/data/lib14/conda/envs/scvi-env", required = TRUE)
library(SeuratWrappers)
library(sceasy)
library(patchwork)
library(sctransform)
library(glmGamPoi)
options(future.globals.maxSize = 8000 * 1024^2)
source("/data/lib14/R/Rscripts/utilities.R")
lr.seur.afterqc.129 <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/final/lr_final_w_isoform_expr.RDS")
DefaultAssay(lr.seur.afterqc.129)
count_mtx <- lr.seur.afterqc.129@assays$Isoform@counts
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass2$transcript_name_unique,]
lr.seur.afterqc.129[["Curated_isoform"]] <- CreateAssayObject(counts = count_mtx)
DefaultAssay(lr.seur.afterqc.129) <- "Curated_isoform"
lr.seur.afterqc.129[["percent.mt"]] = PercentageFeatureSet(lr.seur.afterqc.129, pattern = "^MT-")
# lr.seur.afterqc.129 <- SCTransform(lr.seur.afterqc.129, assay = "Curated_isoform", vars.to.regress = "percent.mt", verbose = FALSE)
lr.data.list <- SplitObject(lr.seur.afterqc.129,split.by = "orig.ident")
# run sctransform
for (i in 1:22) {
  seur <- lr.data.list[[i]]
  DefaultAssay(seur) <- "Curated_isoform"
  seur <- SCTransform(seur, assay = "Curated_isoform",vars.to.regress = "percent.mt", 
                      verbose = FALSE, return.only.var.genes = FALSE)
  lr.data.list[[i]] <- seur
  message("Finish ", i)
}
# Get all highly variable isoforms
var.features.3.9k <- SelectIntegrationFeatures(object.list = lr.data.list, nfeatures = 3900)
var.features.4k <- SelectIntegrationFeatures(object.list = lr.data.list, nfeatures = 4000)
lr.seur <- merge(lr.data.list[[1]], y = lr.data.list[2:22], project = "Long-read")
FindVariableFeatures(lr.seur, nfeatures = 2000)
VariableFeatures(lr.seur) <- var.features
saveRDS(lr.seur, file = "/data/Choi_lung/scLongreads/Seurat/final/lr_final_w_isoform_expr_norm_sbatch.RDS")

# Load data
lr.seur <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final_w_isoform_expr_norm_sbatch.RDS")
DefaultAssay(lr.seur)
lr.seur <- RunPCA(lr.seur, npcs = 30, verbose = F)

# Normalize isoform expression by SCTransform
# Reduce dimensions using Harmony
lr.seur <- IntegrateLayers(
  object = lr.seur, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = FALSE
)
lr.seur <- FindNeighbors(lr.seur, dims = 1:30, reduction = "harmony")
saveRDS(lr.seur, file = "/data/Choi_lung/scLongreads/Seurat/final/lr_final_w_isoform_expr_norm_sbatch.RDS")
lr.seur <- RunUMAP(lr.seur, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(lr.seur, reduction = "umap.harmony", group.by = "Celltype")
saveRDS(lr.seur, file = "/data/Choi_lung/scLongreads/Seurat/final/lr_final_w_isoform_expr_norm_sbatch.RDS")
lr.seur <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/final/lr_final_w_isoform_expr_norm_sbatch.RDS")

# Performing the calculations on the PC coordinates, like before.
library(bluster)
lr.seur.gene <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/final/lr.sc.data.harmony.RDS")
lr.seur.gene <- subset(lr.seur.gene, cells = colnames(lr.seur))
colnames(lr.seur.gene)[1:5]
colnames(lr.seur)[1:5]
which(colnames(lr.seur.gene) != colnames(lr.seur))
lr.seur.gene$Celltype <- lr.seur$Celltype
length(VariableFeatures(lr.seur))

# Calculate the ct ASW for gene-level clustering using 3000 highly variable genes
sil.approx.gene <- approxSilhouette(lr.seur.gene@reductions$harmony@cell.embeddings, clusters=lr.seur.gene$Celltype)
sil.approx <- approxSilhouette(lr.seur@reductions$harmony@cell.embeddings, clusters=lr.seur$Celltype)
sil.approx.gene$ASW_ct <- (sil.approx.gene$width+1)/2
sil.approx$ASW_ct <- (sil.approx$width+1)/2
plot(sil.approx.gene$ASW_ct, sil.approx$ASW_ct)
hist(sil.approx.gene$width)
t.test(sil.approx.gene$ASW_ct, sil.approx$ASW_ct, paired = TRUE, alternative = "greater")
saveRDS(sil.approx.gene, file = "/data/Choi_lung/scLongreads/Seurat/Gene_level_3kHVG_ASW.rds")

VariableFeatures(lr.seur) -> var.features
var.genes <- TALON_afterqc_orf_secondpass2[var.features,]
var.genes <- unique(var.genes$annot_gene_name)
DefaultAssay(lr.seur.gene)

# Re-define the variable features for clustering 
# using the matching genes of 3000 highly variable isoforms
VariableFeatures(lr.seur.gene) <- var.genes
DefaultAssay(lr.seur.gene)
lr.seur.gene <- RunPCA(lr.seur.gene, npcs = 30, verbose = F)
lr.seur.gene <- IntegrateLayers(
  object = lr.seur.gene, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = FALSE
)
lr.seur.gene <- FindNeighbors(lr.seur.gene, dims = 1:30, reduction = "harmony")
length(VariableFeatures(lr.seur.gene))
sil.approx.gene <- approxSilhouette(lr.seur.gene@reductions$harmony@cell.embeddings, clusters=lr.seur.gene$Celltype)

saveRDS(sil.approx.gene, file = "/data/Choi_lung/scLongreads/Seurat/Gene_level_2532HVG_ASW.rds")

