# Summarize DEI results across lung cell types
# compare DEI with DEG 
# Permutation for FDR control also included in this code

# Author: Bolun

library(stringr)
library(stringi)
library(dplyr)
library(ggplot2)
library(Seurat)

files <- list.files("/data/Choi_lung/Elelta/eQTLProjects/BolunTask/PostIsoformFiltering/WholeGroupCorrected/Outputfeb19/unfiltered/", pattern = ".csv")


celltypes <- str_split_fixed(files, "_", 4)[,4]
celltypes <- str_split_fixed(celltypes, "\\.", 2)[,1]
celltypes <- gsub(" ", "_", celltypes)

DEI.list <- list()
Sig_DEI <- NULL
for (i in 1:37) {
  DEI_rst <- read.csv(paste0("/data/Choi_lung/Elelta/eQTLProjects/BolunTask/PostIsoformFiltering/WholeGroupCorrected/Outputfeb19/unfiltered/",files[i]))
  DEI.list[[celltypes[i]]] <- DEI_rst
  Sig_DEI <- c(length(which(DEI_rst$padj < 0.05 & abs(DEI_rst$log2FoldChange) > 1)), Sig_DEI)
}
names(Sig_DEI) <- celltypes

DEI_sig.list <- list()
DEI_unique.list <- list()
for (i in 1:37) {
  DEI_df <- DEI.list[[i]]
  DEI_sub <- DEI_df[which(DEI_df$padj < 0.05 & DEI_df$log2FoldChange > 1),]
  DEI_sub$Celltype <- celltypes[i]
  DEI_sig.list[[i]] <- DEI_sub
  tmp <- do.call(rbind, DEI.list[-i])
  DEI_other <- unique(tmp$X[which(tmp$padj < 0.05 & tmp$log2FoldChange > 0)])
  DEI_unique <- setdiff(DEI_sub$X, DEI_other)
  DEI_unique.list[[i]] <- DEI_unique
}
names(DEI_sig.list) <- celltypes
names(DEI_unique.list) <- celltypes
save(list = ls(), file = "/data/Choi_lung/scLongreads/DEI/DEI_workspace.RData")
df_DEI <- data.frame(Celltype = celltypes,
                     Number = sapply(DEI_sig.list, nrow))
                    
# Theme for plotting
ct_levels <- c(
  'Secretory transitional cells', 'Club', 'AT2', 'Multiciliated', 'AT1', 'Alveolar transitional cells', 'Goblet', 'Basal', 'Neuroendocrine',  # Epithelial
  'Lymphatic EC', 'EC aerocyte capillary', 'EC arterial', 'EC venous systemic', 'EC venous pulmonary', 'EC general capillary',  # Endothelial
  'CD8 T cells', 'Alveolar macrophages', 'NK T cells', 'NK cells', 'Interstitial macrophages', 'CD4 T cells', 'Classical monocyte', 'Alveolar macrophages MT', 'DC2', 'Plasmacytoid DCs', 
  'T cell proliferating', 'Non_classical monocytes', 'Alveolar macrophages CCL3', 'Mast cells', 'Alveolar Mph proliferating', 'B cells', 'Plasma cells',  # Immune
  'SMC', 'Adventitial fibroblasts', 'Mesothelium', 'Myofibroblast', 'Alveolar fibroblasts'  # Stromal
)
ct_levels <- gsub(" ", "_", levels(lr.isoform.sub$Celltype))
df_DEI$Celltype <- factor(df_DEI$Celltype, levels = ct_levels)
fc <- colorRampPalette(c("#FA877533", "#610C04"))
fc_yellow <- colorRampPalette(c("#8B4000", "#FFD70033"))
fc_purple <- colorRampPalette(c("lightpink","#4B0082"))
fc_blue <- colorRampPalette(c("lightblue", "darkblue"))
fc_yellow(9)
cellcolors <- c(fc_yellow(9)[c(3,2,1,4,8,9,7,5,6)], fc(17)[c(1,17,15,16,14,2,13,4,11,
                                                             3,12,5,
                                                             6,7,9,10,8)],
                fc_purple(6)[c(4,2,6,3,5,1)], fc_blue(5)[c(5,1,2,4,3)])

# Summarize the number of DEI across cell types
ggplot(df_DEI, aes(x = Celltype, y = Number, fill = Celltype)) + 
  geom_col(position = 'dodge') +
  geom_text(aes(label = Number), position = position_dodge(width = 0.9), vjust = -0.25, size = 5) + 
  coord_cartesian(ylim = c(0, max(df_DEI$Number) * 1.1)) + 
  ylab('Number of DEIs') + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10),
        legend.position="top",legend.justification="right") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = cellcolors) + NoLegend()
tmp.list <- lapply(DEI_sig.list, function(x) {x$X})
length(unique(Reduce(c, tmp.list)))

# load DEG results using the same pipeline of DEI
files <- list.files("/data/Choi_lung/Elelta/eQTLProjects/BolunTask/Gene_Level/PhaseII/WholeGroup/Unfiltered/", pattern = ".csv")

files <- files[-c(1,2)]
celltypes <- str_split_fixed(files, "_", 4)[,4]
celltypes <- str_split_fixed(celltypes, "\\.", 2)[,1]
celltypes <- gsub(" ", "_", celltypes)

DEG.list <- list()
Sig_DEG <- NULL
for (i in 1:37) {
  DEI_rst <- read.csv(paste0("/data/Choi_lung/Elelta/eQTLProjects/BolunTask/Gene_Level/PhaseII/WholeGroup/Unfiltered/",files[i]))
  
  DEG.list[[celltypes[i]]] <- DEI_rst
  Sig_DEG <- c(length(which(DEI_rst$padj < 0.05 & DEI_rst$log2FoldChange > 1)), Sig_DEG)
}
names(Sig_DEG) <- celltypes

# load DEG results using the same pipeline of DEI
DEG_sig.list <- list()
DEG_sigpadj.list <- list()
DEG_unique.list <- list()
for (i in 1:37) {
  DEI_df <- DEG.list[[i]]
  DEI_sub <- DEI_df[which(DEI_df$padj < 0.05 & DEI_df$log2FoldChange > 0),]
  DEI_sub$Celltype <- celltypes[i]
  DEG_sigpadj.list[[i]] <- DEI_sub
  DEI_sub <- DEI_df[which(DEI_df$padj < 0.05 & DEI_df$log2FoldChange > 1),]
  DEG_sig.list[[i]] <- DEI_sub
  tmp <- do.call(rbind, DEG.list[-i])
  DEI_other <- unique(tmp$X[which(tmp$padj < 0.05 & tmp$log2FoldChange > 0)])
  DEG_unique <- setdiff(DEI_sub$X, DEI_other)
  DEG_unique.list[[i]] <- DEG_unique
}
names(DEG_sig.list) <- celltypes
names(DEG_sigpadj.list) <- celltypes
names(DEG_unique.list) <- celltypes


# Supplementary figure 5a
# compare top DEI (DEI with max logFC for that gene)
# with logFC of the matching DEG in each cell type
DEI_DEG_comp.list <- list()
celltypes <- gsub(" ", "_", celltypes)
for (celltype in celltypes) {
  AT1_DEI <- DEI_sig.list[[celltype]]
  AT1_DEI$gene <- TALON_afterqc_orf_secondpass[AT1_DEI$X, "annot_gene_name"]
  # AT1_DEI_sig_top
  AT1_DEI_sig_top <- AT1_DEI %>% group_by(gene) %>% mutate(log2FC_DEI_max = max(log2FoldChange))
  AT1_DEI_sig_top <- AT1_DEI_sig_top[which(AT1_DEI_sig_top$log2FoldChange == AT1_DEI_sig_top$log2FC_DEI_max),]
  AT1_DEG <- DEG.list[[celltype]]
  AT1_DEG <- subset(AT1_DEG, log2FoldChange > 0)
  colnames(AT1_DEG)[1] <- "gene"
  AT1_DEI_sig_top_sum <- left_join(AT1_DEI_sig_top, AT1_DEG, by = "gene")
  AT1_DEI_sig_top_sum <- AT1_DEI_sig_top_sum[!is.na(AT1_DEI_sig_top_sum$padj.y),]
  AT1_DEI_sig_top_sum <- subset(AT1_DEI_sig_top_sum, log2FoldChange.y > 0 & log2FoldChange.x > 0)
  plot(AT1_DEI_sig_top_sum$log2FoldChange.x, AT1_DEI_sig_top_sum$log2FoldChange.y)
  AT1_DEI_sig_top_sum$DEI <- FALSE
  AT1_DEI_sig_top_sum$DEI <- (AT1_DEI_sig_top_sum$log2FoldChange.x > 1 & AT1_DEI_sig_top_sum$padj.x < 0.05)
  AT1_DEI_sig_top_sum <- subset(AT1_DEI_sig_top_sum, padj.x < 0.05 & padj.y < 0.05)
  AT1_DEI_sig_top_sum$Celltype <- celltype
  DEI_DEG_comp.list[[celltype]] <- AT1_DEI_sig_top_sum
}

DEI_DEG_comp <- do.call(rbind, DEI_DEG_comp.list)
DEI_DEG_comp$Celltype <- factor(DEI_DEG_comp$Celltype, levels = gsub(" ","_",levels(metadata$Celltype)))
ggplot(DEI_DEG_comp, aes(x = log2FoldChange.y, y = log2FoldChange.x)) + 
  geom_point(aes(color = Celltype), size = .5) + scale_color_manual(values = c(cellcolors))+
  xlab('log2 fold change at gene level') + ylab("log2 fold change of top DEI within each gene")+
  geom_abline(intercept = 0, slope = 1, color = "black")+ facet_wrap(~Celltype, ncol = 8)+
  coord_cartesian(ylim = c(0, max(DEI_DEG_comp$log2FoldChange.x))) + theme_bw() +
  theme(strip.background = element_rect(fill='transparent'),legend.position="none")
ggsave("figures/Fig.S5a.pdf", width = 16, height = 10)


######################
#### permutation #####
######################
# To permute the data and break the association
# The phenotype is shuffled

files <- list.files("/data/Choi_lung/scLongreads/DEI/WholegroupPhaseI/output_permutation/UnfilteredResults/SMC/", pattern = ".csv")


celltypes <- str_split_fixed(files, "_", 4)[,4]
celltypes <- str_split_fixed(celltypes, "\\.", 2)[,1]
celltypes <- gsub(" ", "_", celltypes)

DEI_per.list <- list()
Sig_DEI <- NULL
for (i in 1:1000) {
  DEI_rst <- read.csv(paste0("/data/Choi_lung/scLongreads/DEI/WholegroupPhaseI/output_permutation/UnfilteredResults/SMC/",files[i]))
  DEI_per.list[[i]] <- DEI_rst
  Sig_DEI <- c(length(which(DEI_rst$padj < 0.05)), Sig_DEI)
}
# Sig_Num_list <- list()
Sig_Num_list[["Mast_cells"]] <- Sig_DEI
Sig_Num_list[["SMC"]] <- Sig_DEI
Sig_Num_list[["Myofibroblast"]] <- Sig_DEI
hist(Sig_DEI)
tmp.list <- lapply(DEI_per.list, function(x){
  x <- as.data.frame(x[,c(7)])
  return(x)
})
# names(cellcolors) <- gsub(" ", "_",levels(lr$Celltype))
permutation_100times <- do.call(cbind, tmp.list)
permutation_100times_sig <- rowSums((permutation_100times < 0.05), na.rm = TRUE)

hist(permutation_100times_sig, breaks = 100)
tmp <- data.frame(Sig_times = permutation_100times_sig,
                  Celltype = "Myofibroblast")

# df_per <- data.frame(Sig_times = permutation_100times_sig,
#                      Celltype = "AT2")
df_per <- rbind(df_per, tmp)
library(ggplot2)
library(dplyr)
library(patchwork)

# Visualize potential false positive (ie, claimed as significant more than 50 times in 1000-time permutation)
ggplot(df_per, aes(x = Sig_times, fill = Celltype)) +
  geom_histogram(aes(y = ..density..),position="identity",binwidth = 1) + 
  geom_vline(xintercept = 50, linetype = 2) + 
  scale_fill_manual(values = cellcolors) + facet_wrap(~Celltype, nrow = 1) + 
  ylab("Fraction of tested isoforms") + xlab("# of permutated data where a isoform is falsely claimed as DEI") +
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")
ggsave("/data/lib14/project/scLongread/FigSMethodB.pdf", width = 20,height = 7)
table(df_per$Sig_times)
df_per$FDR_control <- TRUE
df_per$FDR_control[which(df_per$Sig_times < 50)] <- FALSE
table(df_per$Celltype, df_per$FDR_control)

rst <- DEI.list$Myofibroblast
rst <- DEI.list$Mast_cells
rst <- DEI.list$SMC
plot((rst$log2FoldChange), -log10(rst$padj))
rst$Sig_times <- df_per$Sig_times[which(df_per$Celltype == "SMC")]
rst$FP <- df_per$FDR_control[which(df_per$Celltype == "SMC")]
rst$gene_type <- "ns"
rst$gene_type[which(rst$padj < 0.05 & rst$log2FoldChange > 1)] <- "up"
rst$gene_type[which(rst$padj < 0.05 & rst$log2FoldChange < -1)] <- "down"
table(rst$gene_type, rst$FP)
rst$gene_type[which(rst$FP & rst$padj < 0.05 & abs(rst$log2FoldChange) > 1)] <- "FP"
# Add colour, size and alpha (transparency) to volcano plot --------------------
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey", "FP"="red") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1, "FP"=2) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5,"FP"=1)

# Look at the detailed distribution of p value and logFC for permuted data
# Spot check the specific potential false positive genes
rst %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  geom_point(shape = 21, colour ="black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-10, 10, 2)),       
                     limits = c(-10, 10)) +
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")+
  ggtitle("SMC")
ggsave("/data/lib14/project/scLongread/FigSMethodC3_FP_1kperm.pdf", width = 8,height = 7)
