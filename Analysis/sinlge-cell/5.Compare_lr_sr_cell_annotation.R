# Code for comparing final annotation between long-read and short-read data

# Author: Bolun, Elelta
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

lr.final <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/lr_final.RDS")
Epi_celltype_levels <- c("AT1", "AT2", "Alveolar transitional cells",
  "Club", "Goblet", "Secretory transitional cells",
  "Basal", "Multiciliated", "Neuroendocrine")
Immune_celltype_levels <- c("Alveolar macrophages", "Alveolar macrophages CCL3+",
                            "Alveolar macrophages MT+", "Alveolar Mφ proliferating", 
                            "Monocyte-derived Mφ", "Classical monocyte", "Non-classical monocytes",
                            "DC2", "CD4 T cells", "CD8 T cells", "NK T cells", "NK cells","T cell proliferating",
                            "Mast cells", "B cells", "Plasma cell/Plasmacytoid DCs")
Endo_celltype_levels <- c("EC arterial", "EC venous pulmonary",
                          "EC venous systemic", "Lymphatic EC", 
                          "EC aerocyte capillary", "EC general capillary")
Stroma_celltype_levels <- c("Adventitial fibroblasts", "Adventitial/peribronchial fibroblasts",
                            "Alveolar fibroblasts", "Myofibroblast", 
                            "SMC", "Mesothelium")


lr.final$Annotation_recluster <- factor(lr.final$Annotation_recluster, levels = c(Epi_celltype_levels, 
                                                                                  Immune_celltype_levels,
                                                                                  Endo_celltype_levels,
                                                                                  Stroma_celltype_levels))
table(lr.final$Annotation_recluster)
table(lr.final$Annotation_recluster, lr.final$Sample_NCI)
df_stat <- as.data.frame(table(lr.final$Annotation_recluster, lr.final$Sample_NCI))
cellnumber_perSample <- as.data.frame(table(lr.final$Sample_NCI))
colnames(df_stat) <- c("Celltype", "Sample", "Cell number")
df_stat$Cell_num_sum_per_sample <- rep(cellnumber_perSample$Freq, each = 37)

df_stat$`Cell proportion` <- (df_stat$`Cell number`/df_stat$Cell_num_sum_per_sample)*100
sum(df_stat$`Cell proportion`[1:37])
df_stat$Category <- rep(c(rep("Epi", 9), rep("Immune", 16), rep("Endo", 6), rep("Stroma", 6)), 129)
library(ggrepel)
library(ggpubr)
ggboxplot(df_stat, x = "Celltype", y = "Cell proportion", color = "Celltype")+
  facet_wrap(~Category, ncol = 2,scales = "free")+
  RotatedAxis()+NoLegend()
library(tidyverse)

library(gghalves)
categorycolors <- c("coral1","lightslateblue","goldenrod1","lightgray")
ggplot(data = df_stat,aes(x=Celltype, y=`Cell proportion`, fill=Category)) +    
  #geom_half_violin(side = "r", color=NA, alpha=0.5,width=0.5) +    
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width=0.5, linewidth=0.5) +    
  geom_half_point_panel(side = "l", shape=21, size=1, color="white") +  
  scale_fill_manual(values = categorycolors) +
  labs(y="Cell Proportion (%)",x=NULL) +   
  theme_classic() + 
  #facet_wrap(~Category, ncol = 2,scales = "free")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "bottom",axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black")) + NoLegend() 
df_stat_sub <- subset(df_stat, `Cell number` >= 5)
table(df_stat_sub$Celltype)
df <- as.data.frame(table(df_stat_sub$Celltype))
ggplot(df, aes(x=Var1, y=Freq,fill=Var1)) +
  geom_bar(stat="identity")+
  labs(y="Number of individuals (> 5 cells)",x= "Cell type") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  NoLegend()



# load sr-data
stromalU= readRDS('/data/Choi_lung/Elelta/FinalAnnotations/Resolution/FinalResolutions/Approach2Final/UMAPstromapp3_updated.rds')
endoU= readRDS('/data/Choi_lung/Elelta/FinalAnnotations/Resolution/FinalResolutions/Approach2Final/UMAPendoapp3_updated.rds')
epiU = readRDS('/data/Choi_lung/Elelta/FinalAnnotations/Resolution/FinalResolutions/Approach2Final/EpiImmResolutions/UMAPepiapp3.60_updated.rds')
immU = readRDS('/data/Choi_lung/Elelta/FinalAnnotations/Resolution/FinalResolutions/Approach2Final/EpiImmResolutions/UMAPimmapp3.6_updated.rds')
stromalU$Category <- "Stroma"
endoU$Category <- "Endothelial"
epiU$Category <- "Epithelial"
immU$Category <- "Immune"
sr.meta <-  Reduce(rbind, list(stromalU@meta.data,
                               endoU@meta.data,
                               epiU@meta.data,
                               immU@meta.data))
sr.meta$barcode_orig <- str_split_fixed(rownames(sr.meta), "_", n = 2)[,1]
table(stromalU$orig.ident)
nci_samples <-unique(stromalU$orig.ident)
nci_samples <- nci_samples[-1]
NCI_samples <- c("1_NCI_1257", "2_NCI_9_16",  "NCI_17_22","3_NCI_23_29", "4_NCI_30_35", "5_NCI_36_41", "6_NCI_42_47",
                 "7_NCI_48_54", "8_NCI_56_61", "9_NCI_62_68", "10_NCI_69_74", "11_NCI_75_80", 
                 "12_NCI_81_86", "13_NCI_87_92",  "14_NCI_93_98", "15_NCI_99_104",
                 "16_NCI_105_110", "17_NCI_111_116", "18_NCI_117_123_2", "19_NCI_124_129",
                 "20_NCI_130_135", "21_NCI_136_140")
for(i in 1:length(NCI_samples)){
  sr.meta$orig.ident = replace(sr.meta$orig.ident, sr.meta$orig.ident == nci_samples[i], NCI_samples[i])
}
table(sr.meta$orig.ident)
sr.meta$barcode <- paste(sr.meta$orig.ident, sr.meta$barcode_orig, sep = "_")
sr.meta <- subset(sr.meta, orig.ident != "nci1257")
lr.meta <- lr.final@meta.data
combined.meta <- left_join(sr.meta,lr.meta, by = "barcode")
write.table(table(combined.meta$cell_types, combined.meta$Annotation_recluster), file = "/data/lib14/project/scLongread/SR_celltype_LR_celltype.txt", 
            sep = "\t", quote = FALSE)
combined.meta <- combined.meta[-which(is.na(combined.meta$Annotation_recluster)),]
table(combined.meta$cell_types)
library(Metrics)
combined.meta$lr_cell_types <- as.character(combined.meta$Annotation_recluster)
celltype_clean <- read.table("/data/lib14/R/Rscripts/celltype_uniform.txt", header = TRUE, sep = "\t")
for(i in 1:nrow(celltype_clean)){
  combined.meta$lr_cell_types = replace(combined.meta$lr_cell_types, combined.meta$lr_cell_types == celltype_clean$cell_types[i], celltype_clean$lr_cell_types[i])
}
combined.meta$cell_types <- as.character(combined.meta$cell_types)
combined.meta$cell_types[which(combined.meta$cell_types == "CD8+ T cells")] <- "CD8 T cells"
combined.meta$cell_types[which(combined.meta$cell_types == "CD4+ T cells")] <- "CD4 T cells"
table(combined.meta$cell_types)
table(combined.meta$lr_cell_types)
combined.meta$cell_types_sr <- 1
combined.meta$cell_types_lr <- 0
combined.meta$cell_types_lr[which(combined.meta$cell_types == combined.meta$lr_cell_types)] <- 1
compare_df <- combined.meta %>% group_by(cell_types) %>% 
  summarise(Precision = precision(cell_types_sr, cell_types_lr),
            Recall = recall(cell_types_sr, cell_types_lr),
            F1_score = f1(cell_types_sr, cell_types_lr)
            )
compare_df <- compare_df[-which(is.na(compare_df$Precision)),]
table(combined.meta$celltype_level_1)


table(combined.meta$celltype_level_1, combined.meta$Category)
c("Endothelial", "Epithelial", "Immune", "Stroma")

mat <- as.matrix.data.frame(table(combined.meta$celltype_level_1, combined.meta$Category))
colnames(mat) <- c("Endothelial", "Epithelial", "Immune", "Stroma")
rownames(mat) <- c("Endothelial", "Epithelial", "Immune", "Stroma")
precision <- diag(mat) / rowSums(mat)
recall <- (diag(mat) / colSums(mat))
df_catalog <- data.frame(celltype = c("Endothelial", "Epithelial", "Immune", "Stroma"),
                         precision = precision,
                         recall = recall)

precision.list <- list()
recall.list <- list()
for (catalog in c("Endothelial", "Epithelial", "Stroma")) {
  combined.immu <- subset(combined.meta, (Category == catalog) & (Category == celltype_level_1))
  mat.immu <- as.matrix.data.frame(table(combined.immu$lr_cell_types, combined.immu$cell_types))
  colnames(mat.immu) <- names(table(combined.immu$cell_types))
  rownames(mat.immu) <- names(table(combined.immu$lr_cell_types))
  precision.immu <- diag(mat.immu) / rowSums(mat.immu)
  recall.immu <- (diag(mat.immu) / colSums(mat.immu))
  precision.list[[catalog]] <- precision.immu
  recall.list[[catalog]] <-recall.immu
}

mat.immu <- as.matrix.data.frame(table(combined.immu$lr_cell_types, combined.immu$cell_types))
colnames(mat.immu) <- names(table(combined.immu$cell_types))
rownames(mat.immu) <- names(table(combined.immu$lr_cell_types))
precision.immu <- diag(mat.immu) / rowSums(mat.immu)
recall.immu <- (diag(mat.immu) / colSums(mat.immu))
precision.immu
recall.immu
compare_df <- data.frame(celltype = names(Reduce(c,precision.list)),
                         precision = Reduce(c,precision.list),
                         recall = Reduce(c,recall.list))
library(reshape2)
compare_df_long <- melt(compare_df)
ggboxplot(compare_df_long, x = "variable", y = "value", color = "celltype")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(compare_df_long, aes(x=celltype, y = value, fill = variable))+
  geom_bar(position="dodge",stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
df_catalog <- melt(df_catalog)
ggplot(df_catalog, aes(x=variable, y = value, color = celltype))+
  geom_point()+
  scale_y_continuous(limits = c(0.95, 1.05))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


