# Figure 3, 4, 5 
# Relevant supplementary figures

# Bolun Li 

library(dplyr)
library(ggplot2)
library(stringr)
library(stringi)
source("/data/lib14/R/Rscripts/utilities.R")
lr <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final_v2_wRNAexpr.RDS")
lr.isoform.sub <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final_w_isoform_expr.RDS")
load("/vf/users/Choi_lung/scLongreads/DEI/DEI_ct_specific_eIsofrom_GTEx_compare.RData")

# Figure 3C
DEI_sig <- do.call(rbind,DEI_sig.list )
DEI_sig$gene_name <-  TALON_afterqc_orf_secondpass2[DEI_sig$X,"annot_gene_name"]
Annotation_gene_list <- unique(c(Epi_sig_list, Imm_sig_list, Endo_sig_list, Stroma_sig_list))
length(which(Annotation_gene_list %in% unique(DEI_sig$gene_name)))
tmp <- subset(DEI_sig, gene_name %in% Annotation_gene_list[(which(Annotation_gene_list %in% unique(DEI_sig$gene_name)))])
tmp <- tmp %>% group_by(gene_name) %>% mutate(rank = rank(-baseMean))
tmp_sub <- subset(tmp, rank == 1)
df_DEI <- data.frame(Celltype = celltypes,
                     Number = sapply(DEI_sig.list, nrow),
                     Type = "Whole group")
df_DEI2 <- data.frame(Celltype = celltypes,
                      Number = sapply(DEI.combined_sig.list, nrow),
                      Type = "Subgroup")
df_DEI <- rbind(df_DEI, df_DEI2)
ct_levels <- gsub(" ", "_", levels(lr.isoform.sub$Celltype))
df_DEI$Celltype <- factor(df_DEI$Celltype, levels = ct_levels)
df_DEI$Type <- factor(df_DEI$Type, levels = c("Subgroup", "Whole group"))
df_DEI$Lable <- df_DEI$Number
df_DEI$Lable[38:74] <- ""
df_DEI$Celltype_type <- paste(df_DEI$Celltype, rep(c("Whole group", "Subgroup"), each = 37), sep = "-")
cor_pal <- c(alpha(cellcolors,0.6),cellcolors)
names(cor_pal) <- paste(rep(levels(df_DEI$Celltype),2), rep(c("Subgroup","Whole group"), each = 37), sep = "-")
df_DEI$Celltype_type <- factor(df_DEI$Celltype_type, levels = paste(rep(levels(df_DEI$Celltype),2), rep(c( "Subgroup", "Whole group"),each = 37), sep = "-"))
df_DEI$Celltype_type <- factor(df_DEI$Celltype_type, levels = paste(rep(levels(df_DEI$Celltype),2), rep(c( "Whole group", "Subgroup"),each = 37), sep = "-"))
df_DEI <- df_DEI[c(38:74, 1:37),]
ggplot(df_DEI, aes(x = Celltype, y = Number, fill = Celltype_type)) + 
  geom_bar(stat="identity",position = "identity") + scale_fill_manual(values = cor_pal) +
  geom_text(aes(label = Lable), vjust = -0.5, size = 5) +
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
        # legend.position="none",
        legend.justification="right") +
  scale_y_continuous(expand = c(0,0))
ggsave("/data/lib14/project/scLongread/Fig3C.pdf", width = 15,height = 7)

legendLabelF1 <- levels(as.factor(df_DEI$Type))
ggplot(df_DEI, aes(x = Celltype, y = Number, fill = Celltype_type)) + 
  geom_bar(aes(alpha = Type),stat="identity",position = "identity") +
  scale_alpha_manual(values=c(1,0.3),
                     name = "",
                     labels = legendLabelF1)+
  geom_text(aes(label = Lable), vjust = -0.5, size = 5) +
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
        # legend.position="none",
        legend.justification="right") +
  scale_y_continuous(expand = c(0,0))
ggsave("/data/lib14/project/scLongread/Fig3C_legend.pdf", width = 15,height = 7)

# Figure 4C
top_beta <- c()
beta_range <- c()
test <- apply(lpairs_pm, 1, function(x){
  min <- min(x)
  max <- max(x)
  beta_range <- abs(abs(max) - abs(min))
  if(min*max > 0){
    if(max > 0){
      x = x/max
      top_beta <- c(top_beta,max)
    }else{
      x = x/min
      top_beta <- c(top_beta,min)
    }
  }else{
    if((max+min) > 0){
      x = x/max
      top_beta <- c(top_beta,max)
    }else{
      x = x/min
      top_beta <- c(top_beta,min)
    }
  }
  return(list((x > 0.5 & x <= 1),top_beta, beta_range))
})
rst_magnitude <- NULL
for (i in 1:length(test)) {
  rst_magnitude <- rbind(rst_magnitude, test[[i]][[1]])
}

rownames(rst_magnitude) <- names(test)
magnitude_sum <- rowSums(rst_magnitude)
rst_sum_pairs <- as.data.frame(t(rst_magnitude))
rst_sum_pairs$Celltype <- rownames(rst_sum_pairs)
rst_sum_pairs$Celltype <- factor(rst_sum_pairs$Celltype, levels = names(isoQTL_sig_list))
rst_sum_pairs <- rst_sum_pairs[order(rst_sum_pairs$Celltype),]
rst_sum_pairs$Categroy <- c(rep("EPI",8),rep("IMMUNE",15), rep("EC", 6), rep("STROMA",4))
df <- data.frame(pairs = names(magnitude_sum),
                 shared_num = magnitude_sum)
df_sum <- as.data.frame(table(df$shared_num))
df_sum$aspect <- "magnitude"
ggplot(df_sum, aes(x=Var1, y=Freq, fill = "#fa8775")) +
  geom_bar(stat="identity") + xlab("Number of shared cell type")+ylab("number of isoQTLs") +
  geom_text(aes(label = Freq),vjust = -.5)+theme_classic()
# Sign
top_beta <- NULL
for (i in 1:length(test)) {
  top_beta <- c(top_beta, test[[i]][[2]])
}
lpairs_sign <- apply(lpairs_pm, 2, function(x) (x*top_beta)>0)

sign_sum <- rowSums(lpairs_sign)
df1 <- data.frame(pairs = names(sign_sum),
                  shared_num = sign_sum)
df_sum1 <- as.data.frame(table(df1$shared_num))
df_sum1$aspect <- "sign"

ggplot(df_sum1, aes(x=Var1, y=Freq, fill = "#0000ff")) +
  geom_bar(stat="identity") + xlab("Number of shared cell type")+ylab("number of isoQTLs") +
  geom_text(aes(label = Freq),vjust = -.5)+theme_classic()+RotatedAxis()+NoLegend()

df_final <- rbind(df_sum, df_sum1)
# df_final$Freq <- round((df_final$Freq/2842)*100, 1)
library(ggbreak)
ggplot(df_final, 
       aes(x=Var1, y=Freq, fill = aspect)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("Number of shared cell type")+ylab("Number of lead isoQTL-eIsoform pairs") +
  geom_text(aes(label = Freq), , vjust=-.1, color="black",
            position = position_dodge(0.9), size=3)+
  scale_y_break(c(400,900))+
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")
ggsave("/data/lib14/project/scLongread/Fig4C.pdf", width = 12,height = 6)

# Figure 4B
NB.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")

df_sum <- data.frame(celltype = names(NB.sig.list), # covar for sec_trans was not generated correctly
                     eIsoform = sapply(NB.sig.list, function(x)nrow(x)))

CT_match <- read.table("/data/Choi_lung/scLongreads/eQTL_comparison/isoQTL_celltype_plotting.txt", sep = "\t")

df_sum$Celltype_matching <- df_sum$celltype
rownames(CT_match) <- CT_match$V1
for(i in 1:nrow(CT_match)){
  df_sum$Celltype_matching = replace(df_sum$Celltype_matching, df_sum$Celltype_matching == CT_match$V1[i], CT_match$V2[i])
}

df_sum$Categroy <- c(rep("EPI",8),rep("IMMUNE",15),rep("EC", 6), rep("STROMA",4))
df_sum$Categroy <- factor(df_sum$Categroy, levels = c("EPI", "IMMUNE", "EC", "STROMA"))
df_sum$Celltype_matching <- factor(df_sum$Celltype_matching, levels = df_sum$Celltype_matching)
ggplot(df_sum, aes(x=Celltype_matching, y=eIsoform, fill = Categroy)) +
  geom_bar(stat="identity") + ylab("Number of eIsoforms")+xlab("Cell type") +
  scale_fill_manual(values = c("#ffd700", "#fa8775","#cd34b5", "#0000ff"))+
  geom_text(aes(label = eIsoform), vjust=-.5) + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 90, hjust = 1, vjust = 0.5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")
ggsave("/data/lib14/project/scLongread/Fig4B.pdf", width = 12,height = 7)

# Figure 4E
df_long$Celltype_matching <- as.character(df_long$celltype)
rownames(CT_match) <- CT_match$V1
for(i in 1:nrow(CT_match)){
  df_long$Celltype_matching = replace(df_long$Celltype_matching, df_long$Celltype_matching == CT_match$V1[i], CT_match$V2[i])
}
df_long$Celltype_matching <- factor(df_long$Celltype_matching, levels = CT_match$V2)
ggplot(df_long, aes(fill=variable, y=value, x=Celltype_matching)) + ylab("Number of cell type sepcific eIsoform")+
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = brewer.pal(n = 8, name = "RdBu")[c(2,7)])+
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 90, hjust = 1, vjust = 0.5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size =16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")
ggsave("/data/lib14/project/scLongread/Fig4E.pdf", width = 8,height = 7)

# Plot isoform structures 
celltypes <- names(DEI_sig.list)
DEI_sig.list <- lapply(celltypes, function(x){
  rst <- DEI_sig.list[[x]]
  rst$Celltype <- x
  return(rst)
})
DEI_sig <- do.call(rbind, DEI_sig.list)
DefaultAssay(lr.isoform.sub)
lr.isoform.sub <- NormalizeData(lr.isoform.sub)
DefaultAssay(lr)
lr <- NormalizeData(lr)
DotPlot(lr, features = "SFTPC", group.by = "Celltype")
DotPlot(lr.isoform.sub, features = "CAV1-202", group.by = "Celltype")
DotPlot(lr.isoform.sub, features = "SCGB1A1-202", group.by = "Celltype")
load("/data/Choi_lung/scLongreads/Seurat/final/lr_isoform_count_ratio_mtx_by_ct_indv.RData")

# Generate isoform ratio matrix
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


rownames(pos) <- pos$X
pos <- pos[TALON_afterqc_orf_secondpass$annot_transcript_id,]
count_mtx <- lr.isoform.sub@assays$Isoform@counts
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass$transcript_name_unique,]

idx <- order(pos$gene_id)
pos <- pos[idx,]
count_mtx <- count_mtx[idx,]

## Group by celltype and indv
lr.isoform.sub$Celltype_chr <- gsub(" ", "_", as.character(lr.isoform.sub$Celltype))
lr.isoform.sub$Celltype_sample <- paste(lr.isoform.sub$Celltype_chr, lr.isoform.sub$Sample, sep = "-")
group <- as.character(lr.isoform.sub$Celltype_sample) %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx)
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
sum(count_mtx_sum)
count_mtx_sum <- Matrix::t(count_mtx_sum)
dim(count_mtx_sum)

# prepare normalized expresssion data
dim(count_mtx_sum)
norm_mtx_sum <- NormalizeData(count_mtx_sum, scale.factor = 10000)
# individual data 
gene_id = "ENSG00000185250" # PPIL6
gene_id = "ENSG00000089127" # OAS1
gene_id = "ENSG00000149021" # SCGB1A1
gene_id = "ENSG00000105974" # CAV1
gene_id = "ENSG00000111671" # SPSB2
gene_id = "ENSG00000186501"


# Fig S11A
gene_id = "ENSG00000185250" # PPIL6
p <- isoform_dist_plot(gene_id, "Multiciliated")
p
ggsave("/data/lib14/project/scLongread/FigS11B.pdf", width = 8, height = 5)
p1 <- plot_top_isoforms_by_ct(gene_id, 1)
p1
ggsave("/data/lib14/project/scLongread/FigS11C.pdf", width = 8, height = 5)

# Fig S13A
gene_id = "ENSG00000089127" # OAS1
p <- isoform_dist_plot(gene_id, "Alveolar_macrophages")
p
ggsave("/data/lib14/project/scLongread/FigS13C.pdf", width = 7, height = 5)
p1 <- plot_top_isoforms_by_ct(gene_id, 4)
p1
ggsave("/data/lib14/project/scLongread/FigS13A.pdf", width = 8, height = 16)

Isoforms = c("SPSB2-TALONT000636472")
Isoforms = c("SCGB1A1-202")
Isoforms = c("CAV1-TALONT003058356", "CAV1-201")
df_long_sub <- subset(df_long, Isoform %in% Isoforms)
df_long_sub <- subset(df_long, Isoform %in% isoform_tested_names)

df_long_sub$Celltype <- factor(df_long_sub$Celltype, levels = gsub(" ", "_",levels(lr.isoform.sub$Celltype)))
df_long_sub$Celltype_matching <- as.character(df_long_sub$Celltype)

for(i in 1:nrow(CT_match)){
  df_long_sub$Celltype_matching = replace(df_long_sub$Celltype_matching, df_long_sub$Celltype_matching == CT_match$V1[i], CT_match$V2[i])
}
df_long_sub$Celltype_matching <- factor(df_long_sub$Celltype_matching, levels = CT_match$V2)
# df_long_sub$Isoform <- factor(df_long_sub$Isoform, levels = names(head(sort(rowSums(count),decreasing = TRUE), 10)))
ggplot(df_long_sub, aes(fill = Celltype_matching, y = value, x = Celltype_matching)) +
  geom_bar(
    position = "dodge", width = 0.8, stat = "summary", fun = "mean",
    color = "black", linewidth = .8
  ) + scale_fill_manual(values = cellcolors) + 
  facet_wrap(~Isoform, nrow = 1, scales = "free") +
  stat_summary(
    fun.data = mean_sdl, geom = "errorbar", color = "black",
    position = position_dodge(0.8), width = 0.2, linewidth = 0.8
  ) + theme(axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
            axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
            axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
            panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent'),
            axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
            legend.title = element_text(size=10),
            legend.text = element_text(size=10),
            legend.position="none",
            axis.text.x = element_text(color = "black", size = 12, angle = 90, vjust = 0.5, hjust=1),
            strip.background = element_blank(), strip.text = element_text(size = 14))+
  geom_point(
    position = position_jitterdodge(0.5, dodge.width = .8),
    alpha = 0.8
  )+ylab("Normalized expression")

ggsave("/data/lib14/project/scLongread/Fig3D.pdf", width = 15,height = 7) # CAV1
ggsave("/data/lib14/project/scLongread/Fig4G.pdf", width = 12,height = 7) # SCGB1A1
ggsave("/data/lib14/project/scLongread/FigS15D.pdf", p,width = 10,height = 7) # OAS1
Isoform_info <- TALON_afterqc_orf_secondpass[,c("annot_gene_id", "transcript_name_unique")]
save(list = c("norm_mtx_sum", "count_mtx_per_final", "Isoform_info", "celltypes", "cellcolors"), file = "/data/Choi_lung/scLongreads/DEI/plots/ISOLUTION.RData")
norm_mtx_sum <- norm_mtx_sum[,]
Isoform_info$abundance <- rowSums(norm_mtx_sum)

tmp <- Isoform_info %>% group_by(annot_gene_id) %>% mutate(order = rank(-abundance))
Isoform_info <- subset(tmp, order < 11)
norm_mtx_sum <- norm_mtx_sum[Isoform_info$transcript_name_unique,]
count_mtx_per_final <- count_mtx_per_final[Isoform_info$transcript_name_unique,]
save(list = c("norm_mtx_sum", "count_mtx_per_final", "Isoform_info", "celltypes", "cellcolors"), file = "/data/Choi_lung/scLongreads/DEI/plots/ISOLUTION.RData")

# Figure 5

# Load library
library(VennDiagram)
GTEx_sGene <- read.table(gzfile("/data/Choi_lung/scLongreads/GTEx/GTEx_Analysis_v8_sQTL/Lung.v8.sgenes.txt.gz"),sep="\t", header = TRUE)
GTEx_sGene$gene_id_new <- str_split_fixed(GTEx_sGene$gene_id, "\\.", 2)[,1]
set1 <- unique(final_table$gene_id)
set2 <- unique(FT_eQTL$phenotype_id)
set3 <- unique(GTEx_sGene$gene_id_new[which(GTEx_sGene$qval < 0.05)])

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Paired")
vennplot_data <- list(sGene=set3, isoGene = set1,sceGene=set2)
library(ggVennDiagram)
ggVennDiagram(vennplot_data)
library("ggvenn")
ggvenn(
  vennplot_data, 
  fill_color = c( "#EFC000FF",  "#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 3,show_percentage = FALSE,
)
ggsave("/data/lib14/project/scLongread/Fig5A.pdf", width = 7,height = 7)
# Chart
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Set 1" , "Set 2 " , "Set 3"),
  filename = '#14_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


# Figure 5B
isoGene_nots <- setdiff(set1,set3)
isoGene_nots_tb <- subset(final_table, gene_id %in% isoGene_nots)
novel_eIsoform_nots <- unique(isoGene_nots_tb$phenotype_id)
# cell type specific
sum(novel_eIsoform_nots %in% ct_spec_isoforms)
unique(subset(isoGene_nots_tb, phenotype_id %in% novel_eIsoform_nots[which(novel_eIsoform_nots %in% ct_spec_isoforms)])$gene_id)
# 504
sum(grepl("TALON", novel_eIsoform_nots))
unique(subset(isoGene_nots_tb, phenotype_id %in% novel_eIsoform_nots[grepl("TALON", novel_eIsoform_nots)])$gene_id)
df <- data.frame(value = c(271,504,308,565),
                 level = c("gene","gene","isoform","isoform"),
                 label= c("271\n(29.4%)", "504\n(54.7%)", "308","565"),
                 reason = c("Novel isoforms", "Cell-type specific", "Novel isoforms", "Cell-type specific"))


ggplot(df, aes(x = reason, y = value, fill = level)) +
  geom_col(position = "dodge") +
  geom_text(
    aes(label = label),
    colour = "black", size = 4.5,
    vjust = 1.5, position = position_dodge(.9)
  )+ylab("Number of genes/isoforms")+ ggtitle("921 novel isoGenes")+
  theme(axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        
        axis.text.x = element_text(color = "black", size = 12),
        strip.background = element_blank(), strip.text = element_text(size = 14))
ggsave("/data/lib14/project/scLongread/Fig5B.pdf", width = 6,height = 7)


# Figure 5C
FT_isoQTL_shared_same_ct_joined <- read.table(file = "/data/Choi_lung/scLongreads/eQTL_comparison/isoeQTL_shared_same_ct_coloc_rst.txt", 
                                              header = TRUE, sep = "\t")


df <- data.frame(transcript_ct = FT_isoQTL_shared_same_ct_joined$transcript_ct,
                 PPH4 = FT_isoQTL_shared_same_ct_joined$PP.H4.abf,
                 Novelty = "Known")
df$Novelty[grepl("TALON", df$transcript_ct)] <- "Novel"
df$Colocalization <- "Non-colocalized"
df$Colocalization[which(df$PPH4 > 0.7)] <- "Colocalized"
library(gghalves)
# categorycolors <- c("coral1","lightslateblue","goldenrod1","lightgray")
ggplot(data = df,aes(x=Novelty, y=PPH4, colour = Novelty)) +    
  geom_boxplot(fill = NA, linewidth = 1, staplewidth = 0.2) +
  geom_jitter(aes(colour = Colocalization),position = position_jitter(seed = 1, width = 0.2)) +
  # scale_fill_manual(values = categorycolors) +
  scale_color_manual(values = c("#B90E0A","#619CFF","black", "#F8766D"))+
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red", linewidth = 1) +
  labs(y="PP H4 of colocalization between \n isoQTL and eQTL in the same cell type",x=NULL) +   
  theme_classic() + 
  #facet_wrap(~Category, ncol = 2,scales = "free")+ 
  theme(legend.position = "bottom",axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black")) + NoLegend() 
ggsave("/data/lib14/project/scLongread/Fig5C.pdf", width = 7,height = 7)


NB.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_list.rds")
Tensor.list <- readRDS("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/isoQTL_list.rds")
AT2_pairs <- paste(NB.list$AT2$phenotype_id, NB.list$AT2$variant_id, sep = ":")
AT2_sig_pairs <- AT2_pairs[which(NB.list$AT2$qval < 0.05)]
which(is.na(NB.list$AT2$pval_beta))
AT2_NB <- read.table("/data/Choi_lung/scLongreads/mashr/output/AT2/qtl_results_all.txt", header = TRUE,
                     sep = "\t")
AT2_NB$pair <- paste(AT2_NB$phenotype_id, AT2_NB$snp, sep = ":")
AT2_tensor <- read.table("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/output/AT2/qtl_results_all.txt", header = TRUE,
                         sep = "\t")
AT2_tensor$pair <- paste(AT2_tensor$phenotype_id, AT2_tensor$variant_id, sep = ":")

AT2_NB_lead <- subset(AT2_NB, pair %in% AT2_pairs)
AT2_tensor_lead <- subset(AT2_tensor, pair %in% AT2_pairs)
rownames(AT2_NB_lead) <- AT2_NB_lead$pair
rownames(AT2_tensor_lead) <- AT2_tensor_lead$pair
AT2_NB_lead <- AT2_NB_lead[AT2_pairs[!is.na(NB.list$AT2$pval_beta)],]
AT2_tensor_lead <- AT2_tensor_lead[AT2_pairs[!is.na(NB.list$AT2$pval_beta)],]
plot(-log10(AT2_tensor_lead$pval_nominal), -log10(AT2_NB_lead$pval_nominal))
abline(0,1)
plot((AT2_tensor_lead$slope), (AT2_NB_lead$slope))
plot((AT2_tensor_lead$slope_se), (AT2_NB_lead$slope_se))
save(list = c("AT2_NB_lead", "AT2_tensor_lead"), file = "/data/Choi_lung/scLongreads/jaxqtl/Model_comparison_AT2.RData")

load("/data/Choi_lung/scLongreads/jaxqtl/Model_comparison_AT2.RData")
df <- data.frame(Tensor = -log10(AT2_tensor_lead$pval_nominal), 
                 NB = -log10(AT2_NB_lead$pval_nominal),
                 pair = AT2_tensor_lead$pair)

rownames(df) <- df$pair
df$Significant <- FALSE
df[AT2_sig_pairs,]$Significant <- TRUE
ggplot(subset(df, Tensor < 10 & NB < 10), aes(x = Tensor, y= NB,colour = Significant)) +
  geom_point(alpha = 0.5)
keep <- grepl(paste(Tensor_sig_iso, collapse="|"), 
              row.names(AT2_tensor_lead))

df$Tensor_sig <- keep
df$If_Significant <- df$Tensor_sig + df$Significant
table(df$If_Significant)
df$Which_sig <- "n.s."
df$Which_sig[(df$Tensor_sig & df$If_Significant == 1)] <- "Tensor Sig only"
df$Which_sig[(df$Significant & df$If_Significant == 1)] <- "NB Sig only"
df$Which_sig[(df$If_Significant == 2)] <- "Sig in both"
table(df$Which_sig)
df <- df[-which.max(df$NB),]
df$Which_sig <- factor(df$Which_sig, levels = rev(c("n.s.", "Tensor Sig only",
                                                    "NB Sig only", "Sig in both")))

# Figure S10
ggplot(df, aes(x = Tensor, y= NB,colour = Which_sig)) +
  geom_point(alpha = 0.3) + ylab("-log(pval_nominal) in NB")  + 
  geom_point(data = subset(df, Which_sig != "n.s."), aes(x = Tensor, y= NB,colour = Which_sig), alpha = 0.7)+
  xlab("-log(pval_nominal) in Tensor") +
  scale_color_manual(values = rev(c("#BEBEBE33","#F8766D","#619CFF", "#C77CFF"))) +
  geom_abline(slope=1, intercept=0) +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),  
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_text(color = "black", size = 16),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10))
ggsave("/data/lib14/project/scLongread/FigS10A.pdf", width = 10,height = 7)
df$Tensor <- AT2_tensor_lead$slope[-which.min(AT2_NB_lead$pval_nominal)]
df$NB <- AT2_NB_lead$slope[-which.min(AT2_NB_lead$pval_nominal)]

ggplot(df, aes(x = Tensor, y= NB,colour = Which_sig)) +
  geom_point(alpha = 0.3) + ylab("Effect size in NB")  + 
  geom_point(data = subset(df, Which_sig != "n.s."), aes(x = Tensor, y= NB,colour = Which_sig), alpha = 0.7)+
  xlab("Effect size in Tensor") +
  scale_color_manual(values = rev(c("#BEBEBE33","#F8766D","#619CFF", "#C77CFF"))) +
  geom_abline(slope=1, intercept=0) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),  
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_text(color = "black", size = 16),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10))
ggsave("/data/lib14/project/scLongread/FigS10B.pdf", width = 7,height = 7)
AT2_cis_NB <- NB.list$AT2
AT2_cis_Tensor <- Tensor.list$AT2
Tensor_sig_iso <- AT2_cis_Tensor$phenotype_id[which(AT2_cis_Tensor$qval < 0.05)]
NB_sig_iso <- AT2_cis_NB$phenotype_id[which(AT2_cis_NB$qval < 0.05)]
exprs <- expr_list$AT2
rownames(exprs) <- expr_list$AT2$trascript_id
rownames(norm_mtx_sum) <- TALON_afterqc_orf_secondpass$annot_transcript_id[idx]
Tensor_eIso_expr <- norm_mtx_sum[Tensor_sig_iso,grepl("AT2", colnames(norm_mtx_sum))]
NB_eIso_expr <- norm_mtx_sum[NB_sig_iso,grepl("AT2", colnames(norm_mtx_sum))]
hist(rowMeans(Tensor_eIso_expr))
hist(rowMeans(NB_eIso_expr))
both_eIso_expr <- norm_mtx_sum[intersect(Tensor_sig_iso, NB_sig_iso),grepl("AT2", colnames(norm_mtx_sum))]
df_expr_both <- data.frame(exprs = rowMeans(both_eIso_expr),
                           isoform = rownames(both_eIso_expr),
                           Model = "Both")
df_expr_tensor <- data.frame(exprs = rowMeans(Tensor_eIso_expr),
                             isoform = rownames(Tensor_eIso_expr),
                             Model = "Tensor")
df_expr_NB <- data.frame(exprs = rowMeans(NB_eIso_expr),
                         isoform = rownames(NB_eIso_expr),
                         Model = "NB")

df_expr_NB <- rbind(df_expr_NB, df_expr_both)
df_expr_tensor <- rbind(df_expr_tensor, df_expr_both)
df_expr <- rbind(df_expr_NB, df_expr_tensor)
df_expr$exprs[which(df_expr$exprs > 1)] <- 1.02

df_expr <- rbind(df_expr_NB, df_expr_tensor)
df_expr$Sig <- df_expr$Model
df_expr$Sig[which(df_expr$isoform %in% intersect(Tensor_sig_iso, NB_sig_iso))] <- "Both"
library(gghalves)

ggplot() +
  geom_half_boxplot(
    data = df_expr,
    aes(x = Model, y = exprs, colour = Model), 
    side = "r", outlier.colour = "black",fill = NA,
    errorbar.draw = FALSE, width=0.7) +
  geom_half_violin(
    data = df_expr %>% filter(Model=="NB"), 
    aes(x = Model, y = exprs, split = Sig, fill = Sig),
    position = "identity", color=NA, alpha=0.5,width=1) +
  geom_half_violin(
    data = df_expr %>% filter(Model=="Tensor"), 
    aes(x = Model, y = exprs, split = Sig, fill = Sig),
    position = "identity", color=NA, alpha=0.5,width=1) +
  geom_half_point(
    data = df_expr,
    mapping = aes(x = Model, y = exprs, colour = Sig, fill=Sig),shape=21, size=1)+
  labs(y="Normalized expression of eIsoform",x=NULL) +   
  scale_fill_manual(values = rev(c("#F8766D","#619CFF", "#C77CFF"))) +
  scale_color_manual(values = rev(c("#F8766D","#619CFF", "#C77CFF"))) +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),  
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_text(color = "black", size = 16),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10), legend.position = "none")
ggsave("/data/lib14/project/scLongread/FigS10C_V3.pdf", width = 10,height = 7)
ggsave("/data/lib14/project/scLongread/FigS10C_legend.pdf", width = 10,height = 7)

# Figure S10C
ggplot(df_expr) +
  geom_histogram(aes(x=exprs,y=..density..,fill=Model),position="identity",alpha=.5,binwidth = 0.02) + 
  geom_density(aes(x=exprs,y=..density.., color=Model), linewidth = 1)+
  scale_fill_manual(values = rev(c("#F8766D","#619CFF", "#C77CFF"))) +
  scale_color_manual(values = rev(c("#F8766D","#619CFF", "#C77CFF"))) +
  ylab("Fraction")+xlab("Normalized expression of eIsoform")+
  # scale_fill_manual(values = c("#6BAED6", "#FC8D59"))+
  theme_bw()+
  theme(axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        strip.text = element_text(size = 16, colour = "black"))

ggsave("/data/lib14/project/scLongread/FigS10C_v2.pdf", width = 10,height = 7)

# Figure 3B
table(lr.isoform.sub$Celltype)
table(lr.isoform.sub$Celltype, lr.isoform.sub$Sample)
df_stat <- as.data.frame(table(lr.isoform.sub$Celltype, lr.isoform.sub$Sample))
cellnumber_perSample <- as.data.frame(table(lr.isoform.sub$Sample))
colnames(df_stat) <- c("Celltype", "Sample", "Cell number")
df_stat$Cell_num_sum_per_sample <- rep(cellnumber_perSample$Freq, each = 37)

df_stat$`Cell proportion` <- (df_stat$`Cell number`/df_stat$Cell_num_sum_per_sample)*100
sum(df_stat$`Cell proportion`[1:37])
df_stat$Category <- rep(c(rep("Epi", 9), rep("Immune", 16), rep("Endo", 6), rep("Stroma", 6)), 129)
library(ggrepel)
library(ggpubr)
library(tidyverse)
library(gghalves)

ggplot(data = df_stat,aes(x=Celltype, y=`Cell proportion`, fill=Celltype)) + 
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width=0.7, linewidth=0.5) +
  geom_half_point_panel(side = "l", shape=21, size=1) +  
  scale_fill_manual(values = cellcolors) +
  labs(y="Cell Proportion per individual (%)",x=NULL) +   
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(color = "black", size = 12, angle = 90, vjust = 0.5, hjust=1),  
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_text(color = "black", size = 16),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10),
        legend.position = "none")
ggsave("/data/lib14/project/scLongread/Fig3B.pdf", width = 12,height = 7)


# Figure 5G
files <- list.files("/data/Choi_lung/scLongreads/TALON_workspace/test10s/", pattern = "_strict.ioe")
files <- files[8:14]
idx <- which(!duplicated(final_table_AS_annot$phenotype_id))
eIsoform_AS_sum <- NULL
for (AS in AS_catalog) {
  final_table_AS_annot[[AS]] <- grepl(AS, final_table_AS_annot$AS_catlog)
  
  eIsoform_AS_sum <- c(eIsoform_AS_sum, sum(final_table_AS_annot[[AS]][idx]))
}

names(eIsoform_AS_sum) <- AS_catalog
df_eisoform <- data.frame(AS_catalog = names(eIsoform_AS_sum), 
                          Number = eIsoform_AS_sum)
ggplot(df_eisoform, aes(x=AS_catalog, y=Number, fill = AS_catalog)) +
  geom_bar(stat="identity") + ylab("Number of eIsoforms have alternative splicing events")+xlab("AS catagory") +
  geom_text(aes(label = Number), vjust=-.5) +
  scale_fill_brewer(palette = "Dark2") + 
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),  
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_text(color = "black", size = 16),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10))

ggsave("/data/lib14/project/scLongread/Fig5G.pdf", width = 7,height = 7)

# Figure 5H
share_seur <- readRDS("/data/Choi_lung/ChiaHan/CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds")
ranges.show <- StringToGRanges("chr1-27319319-27319863")
ranges.show$color <- "orange"

Idents(share_seur) <- "CellType"
CoveragePlot(
  object = share_seur,
  region = "chr1-27319000-27337500",
  region.highlight = ranges.show,
  extend.upstream = 200,
  extend.downstream = 200
)
ggsave("/data/lib14/project/scLongread/Fig5H2.pdf", width = 10,height = 7)
ranges.show <- StringToGRanges(c("chr14-34872798-34875547", "chr14-34743531-34744103"))
ranges.show$color <- c("orange", "darkred")

# Generate scATAC-seq coverage plot
Idents(share_seur) <- "CellType"
CoveragePlot(
  object = share_seur,
  region = "chr14-34743531-35404749",
  region.highlight = ranges.show,
  extend.upstream = 200,
  extend.downstream = 200,
  links = TRUE
)


# Figure S11
NB.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")
NB.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_list.rds")
df_sum <- data.frame(celltype = names(NB.list), # covar for sec_trans was not generated correctly
                     Total_isoform = sapply(NB.list, function(x)nrow(x)),
                     eIsoform = sapply(NB.sig.list, function(x)nrow(x)))

CT_match <- read.table("/data/Choi_lung/scLongreads/eQTL_comparison/isoQTL_celltype_plotting.txt", sep = "\t")

df_sum$Celltype_matching <- df_sum$celltype
rownames(CT_match) <- CT_match$V1
for(i in 1:nrow(CT_match)){
  df_sum$Celltype_matching = replace(df_sum$Celltype_matching, df_sum$Celltype_matching == CT_match$V1[i], CT_match$V2[i])
}

df_sum$Categroy <- c(rep("EPI",8),rep("IMMUNE",15),rep("EC", 6), rep("STROMA",4))
df_sum$Categroy <- factor(df_sum$Categroy, levels = c("EPI", "IMMUNE", "EC", "STROMA"))
df_sum$Celltype_matching <- factor(df_sum$Celltype_matching, levels = df_sum$Celltype_matching)
df_sum$Cell_num <- table(lr$Celltype)[gsub("_", " ", df_sum$celltype)]
cell_num_per_indv <- as.matrix.data.frame(table(lr$Celltype, lr$Sample))
rownames(cell_num_per_indv) <- names(table(lr$Celltype))
rownames(cell_num_per_indv)[16] <- "Non classical monocytes"
cell_num_per_indv[which(cell_num_per_indv <= 5)] <- NA
cell_num_per_indv <- cell_num_per_indv[gsub("_", " ", df_sum$celltype),]
df_sum$Cell_num_per_indv <- rowMeans(cell_num_per_indv, na.rm = TRUE)
ggplot(df_sum, aes(x=Celltype_matching, y=Total_isoform, fill = Categroy)) +
  geom_bar(stat="identity") + ylab("Number of eIsoforms")+xlab("Cell type") +
  scale_fill_manual(values = c("#ffd700", "#fa8775","#cd34b5", "#0000ff"))+
  geom_text(aes(label = Total_isoform), vjust=-.5) + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 90, hjust = 1, vjust = 0.5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")
ggsave("/data/lib14/project/scLongread/FigS11A.pdf", width = 12,height = 7)
library(ggrepel)
library(ggpubr)
df<- df_sum %>%
  arrange(desc(Cell_num)) %>%
  mutate(Celltype_matching = factor(Celltype_matching, Celltype_matching))
df$Cell_num <- as.numeric(df$Cell_num)
df_sum$Cell_num <- as.numeric(df_sum$Cell_num)
ggplot(df_sum,aes(x=Total_isoform,y=eIsoform))+
  geom_point(aes(size = Cell_num, fill = Celltype_matching), alpha = 0.7, shape=21 ) +
  scale_size_area(max_size = 10)+
  stat_smooth(method=lm,formula=y~x) + 
  geom_text_repel(aes(color = Celltype_matching), label=df_sum$Celltype_matching)+
  scale_fill_manual(values = cellcolors)+
  scale_color_manual(values = cellcolors)+
  stat_cor(label.y = 600,aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  ylab("Number of eIsoform")+
  xlab("Number of tested isoforms")+theme_classic()+
  ggtitle("Association between the numbers of eIsoform and tested isoforms")+
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")
ggsave("/data/lib14/project/scLongread/FigS11B.pdf", width = 7,height = 7)
ggplot(df_sum,aes(x=Cell_num,y=eIsoform))+geom_point() +
  stat_smooth(method=lm,formula=y~x) + 
  geom_text_repel(aes(color = Celltype_matching),label=df_sum$Celltype_matching)+
  scale_color_manual(values = cellcolors)+
  stat_cor(label.y = 600,aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  ylab("Number of eIsoform")+
  xlab("Number of cells")+theme_classic()+
  ggtitle("Association between the numbers of eIsoform and cells")+
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")
ggsave("/data/lib14/project/scLongread/FigS11C.pdf", width = 7,height = 7)


files_path <- paste0("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/", names(NB.sig.list),
                     "/samples.txt")
indv_num <- NULL
for (i in 1:length(files_path)) {
  tb <- read.table(files_path[i], sep = "\t", header = FALSE)
  indv_num <- c(indv_num, nrow(tb))
}
df_sum$Indv_num <- indv_num
ggplot(df_sum,aes(x=Indv_num,y=eIsoform))+geom_point() +
  stat_smooth(method=lm,formula=y~x) + 
  geom_text_repel(aes(color = Celltype_matching),label=df_sum$Celltype_matching)+
  scale_color_manual(values = cellcolors)+
  stat_cor(label.y = 600,aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  ylab("Number of eIsoform")+
  xlab("Number of individuals")+theme_classic()+
  ggtitle("Association between the numbers of eIsoform and individuals")+
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")
ggsave("/data/lib14/project/scLongread/FigS11D.pdf", width = 7,height = 7)


eIsoform <- str_split_fixed(rownames(true_sig), "\\|", 2)[,1]
mashr_sig.list <- lapply(ct, function(x){
  idx <- true_sig[,x]
  mashr_sig <- unique(eIsoform[idx])
})
df$eIsoSig_afterMashr <- sapply(mashr_sig.list, function(x) length(x))
ggplot(df, aes(x=Celltype, y=Sig_afterMashr, fill = Celltype)) +
  geom_bar(stat="identity")+scale_fill_manual(values = cellcolors)+ylab("Number of isoQTL after Mashr")+
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 90, hjust = 1, vjust = 0.5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size =16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")
ggplot(df, aes(x=Celltype, y=eIsoSig_afterMashr, fill = Celltype)) +
  geom_bar(stat="identity")+scale_fill_manual(values = cellcolors)+ylab("Number of eIsoform after Mashr")+
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 90, hjust = 1, vjust = 0.5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size =16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none",
        legend.justification="right")
