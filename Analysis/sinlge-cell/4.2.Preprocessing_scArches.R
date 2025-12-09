# convert seurat object into h5ad
# run scArches

library(dior)
library(reticulate)
library(scater)
library(SeuratDisk)
library(loomR)
library(Seurat)
library(stringi)
library(stringr)
load("/data/Choi_lung/scLongreads/Seurat/lr.sc.data.afterQC.norm.RData")

ad <- import("anndata", convert = FALSE)
lr.seur@assays$Isoforms <- NULL
lr.seur@assays$SCT <- NULL
lr.seur@assays$prediction.score.ann_level_1 <- NULL
lr.seur@assays$prediction.score.ann_level_2 <- NULL
lr.seur@assays$prediction.score.ann_level_3 <- NULL
lr.seur@assays$prediction.score.ann_level_4 <- NULL
lr.seur@assays$prediction.score.ann_level_5 <- NULL
lr.seur@assays$prediction.score.ann_finest_level <- NULL
lr.seur@assays$RNA@data <- lr.seur@assays$RNA@counts
DefaultAssay(lr.seur) <- "RNA"

SaveH5Seurat(lr.seur, filename = "/data/Choi_lung/scLongreads/scarches/notebooks/LrHL.h5Seurat")
Convert("/data/Choi_lung/scLongreads/scarches/notebooks/LrHL.h5Seurat", dest = "h5ad")

# convert output of scArches to seurat object
Convert("/data/Choi_lung/scLongreads/Seurat/query_with_refbased_emb_and_anns.h5ad", dest = "h5Seurat")
lr.seur.annot <- LoadH5Seurat("/data/Choi_lung/scLongreads/Seurat/query_with_refbased_emb_and_anns.h5seurat")
ad <- readH5AD("/data/Choi_lung/scLongreads/Seurat/query_with_refbased_emb_and_anns.h5ad")
adata_Seurat <- as.Seurat(ad, counts = "X", data = NULL)
table(adata_Seurat$predicted.ann_finest_level)
adata_Seurat$Celltype_clean <- adata_Seurat$ann_level_5_transferred_label
df <- as.data.frame(str_split_fixed(adata_Seurat$Celltype_clean, pattern = "_", n = 2))
df$V3 <- as.numeric(df$V1)
df$V2[which(is.na(df$V3))] = df$V1[which(is.na(df$V3))]
adata_Seurat$Celltype_clean <- df$V2

write.table(table(adata_Seurat$Celltype_clean, adata_Seurat$predicted.ann_finest_level), file = "/data/Choi_lung/scLongreads/Seurat/Compare_scArches_Azimuth.txt", 
            sep = "\t", quote = FALSE)

################
# Compare the automatic annotation between long-read and short-read

sr.seur.annot <- readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/demuxlet_vireo/az.RDS')
sr.seur.annot@assays$SCT <- NULL
DefaultAssay(sr.seur.annot) <- "RNA"
sr.annot.data.list <- SplitObject(sr.seur.annot, split.by = "orig.ident")
rm(sr.seur.annot)
sr.annot.data.list <- sr.annot.data.list[-c(1,3)]
names(sr.annot.data.list)
names(sr.annot.data.list) <- c("1_NCI_1257", "2_NCI_9_16",  "NCI_17_22","3_NCI_23_29", "4_NCI_30_35", "5_NCI_36_41", "6_NCI_42_47",
                               "7_NCI_48_54", "8_NCI_56_61", "9_NCI_62_68", "10_NCI_69_74", "11_NCI_75_80", 
                               "12_NCI_81_86", "13_NCI_87_92",  "14_NCI_93_98", "15_NCI_99_104",
                               "16_NCI_105_110", "17_NCI_111_116", "18_NCI_117_123_2", "19_NCI_124_129",
                               "20_NCI_130_135", "21_NCI_136_140")
for (i in 1:22) {
  sr.annot.data.list[[i]]$orig.ident <- names(sr.annot.data.list)[i]
}
gc()
sr.seur.annot <- merge(sr.annot.data.list[[1]], y = sr.annot.data.list[2:22], add.cell.ids = names(sr.annot.data.list), project = "short-read")
table(sr.seur.annot$orig.ident)
sr.seur.annot$barcode <- paste0(str_split_fixed(colnames(sr.seur.annot), "-", n=2)[,1], "-1")
rm(sr.annot.data.list)
gc()

library(dplyr)
annot.meta <- left_join(adata_Seurat@meta.data, sr.seur.annot@meta.data,by = "barcode")
write.table(table(annot.meta$azimuth_sorted, annot.meta$predicted.ann_finest_level.y), file = "/data/Choi_lung/scLongreads/Seurat/Compare_sr_lr_Azimuth.txt", 
            sep = "\t", quote = FALSE)

annot.meta$azimuth_sorted <- factor(annot.meta$predicted.ann_finest_level.x, levels = str_sort(unique(annot.meta$predicted.ann_finest_level.x)))

DimPlot(adata_Seurat, reduction = "X_umap", group.by = "Celltype_clean", label = TRUE)+coord_equal()
DimPlot(adata_Seurat, reduction = "X_umap", group.by = "predicted.ann_finest_level", label = TRUE)+coord_equal()
save(list = ls(), file = "/data/Choi_lung/scLongreads/Seurat/lr.sc.data.scArches.annot.129.RData")

# Heat map to show the percentage of agreement
df_per <- data.frame((table(annot.meta$azimuth_sorted, annot.meta$predicted.ann_finest_level.y)))
df_new <- apply(df_per, 2, function(x)x/sum(x))
library(reshape2)
df_new <- melt(df_new)
df_new$Long.read <- df_new$Var1
df_new$Short.read <- df_new$Var2
df_new$Percentage <- df_new$value
ggplot(df_new, aes(Short.read, Long.read, fill= Percentage)) + 
  geom_tile()+ theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))+coord_equal()

