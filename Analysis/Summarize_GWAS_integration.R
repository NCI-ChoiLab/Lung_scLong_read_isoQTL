# Summarize all colocalization results together
# Summarize TWAS results

# Bolun

.libPaths()
library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)
library(future)
library(stringr)
library(stringi)
library(forcats)
library(Matrix)
library(arrow)
transcripts <- c("HLA-DRB5-TALONT001746251","SECISBP2L-TALONT001050043",
                 "HLA-DRB1-TALONT001751909","CHMP3-203","RSPH4A-TALONT003150695",
                 "PPIL6-207","HLA-A-TALONT001165898","PSMA4-201","HSPA1B-201")
celltypes <- c("AT2", "Multiciliated", "Alveolar_transitional_cells",
               "NK_cells", "Non_classical_monocytes", "Lymphatic_EC")
NB.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")
NB.sig.list <- lapply(names(NB.sig.list), function(x){
  rst <- NB.sig.list[[x]]
  rst$Celltype <- x
  return(rst)
})
isoQTL_sig <- do.call(rbind,NB.sig.list)

isoQTL_sig$sig_iso_ct <- paste(isoQTL_sig$phenotype_id, isoQTL_sig$Celltype, sep = "-")
Byun_total <- read.table("/data/Choi_lung/scLongreads/colocalization/Coloc_Byun_total_lung_PP_sum.txt", header = TRUE)
Byun_luad <- read.table("/data/Choi_lung/scLongreads/colocalization/Byun_coloc/ADE_sig_coloc_table.txt", header = TRUE)
Byun_sqc <- read.table("/data/Choi_lung/scLongreads/colocalization/Byun_coloc/SQC_sig_coloc_table.txt", header = TRUE)
Shi_EAS <- read.table("/data/Choi_lung/scLongreads/colocalization/Coloc_Shi_EAS_PP_sum.txt", header = TRUE)
Shi_meta <- read.table("/data/Choi_lung/scLongreads/colocalization/Coloc_Shi_meta_PP_sum.txt", header = TRUE)
coloc.list <- list(Shi_EAS=Shi_EAS, Shi_meta=Shi_meta, Byun_lung=Byun_total,Byun_luad= Byun_luad[,c(1:9)],Byun_sqc= Byun_sqc[,c(1:9)])

coloc.list <- lapply(1:5,function(i){
  coloc_rst <- coloc.list[[i]]
  coloc_rst$GWAS <- names(coloc.list)[i]
  return(coloc_rst)
})

coloc_rst <- do.call(rbind, coloc.list)
coloc_rst <- subset(coloc_rst, PP.H4.abf > 0.7)

output.fev1 <- readRDS("/data/Choi_lung/scLongreads/colocalization/lung_function/Colocalization_output_table_fev1.rds")
output.ratio <- readRDS("/data/Choi_lung/scLongreads/colocalization/lung_function/Colocalization_output_table.rds")
output.ratio$GWAS <- "FEV1_FVC"
output.fev1$GWAS <- "FEV1"
output.fev1 <- output.fev1[!duplicated(output.fev1$sig_iso_ct),]
output.ratio <- output.ratio[!duplicated(output.ratio$sig_iso_ct),]
coloc_rst$GWAS <- "Lung_cancer"
coloc_rst$isoform_ct <- paste(coloc_rst$transcript_id, coloc_rst$Celltype,sep = "_")
coloc_rst <- coloc_rst[!duplicated(coloc_rst$isoform_ct),]
tmp <- rbind(coloc_rst[,c(1:10)], output.ratio[,c(1:9,12)])
coloc_final <- rbind(tmp, output.fev1[,c(1:9,12)])

coloc_rst_PPH4 <- as.matrix.data.frame(table(coloc_final$Celltype, coloc_final$transcript_name))
rownames <- names(table(coloc_final$Celltype))
rownames(coloc_rst_PPH4) <- rownames
colnames <- names(table(coloc_final$transcript_name))
colnames(coloc_rst_PPH4) <- colnames

transcripts <- unique(coloc_final$transcript_name)
celltypes <- unique(coloc_final$Celltype)
celltype_levels <- levels(lr$Celltype)
celltypes <- intersect(gsub(" ", "_", celltype_levels), celltypes)
coloc_rst_PPH4 <- coloc_rst_PPH4[celltypes,]
p <- pheatmap::pheatmap(coloc_rst_PPH4, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
                   color = c(NA,"yellow", "orange","red"), border_color = NA,cellwidth = 15, cellheight = 15)
ggsave("/data/lib14/project/scLongread/ASHG_coloc_count.pdf", p,width = 15, height = 7)
count_mtx <- lr.isoform.sub@assays$Isoform@counts
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass$transcript_name_unique,]
TALON_afterqc_orf_secondpass$total_counts <- Matrix::rowSums(count_mtx)

TALON_afterqc_orf_secondpass <- TALON_afterqc_orf_secondpass %>% group_by(annot_gene_id) %>%
  mutate(percentage = total_counts/sum(total_counts))
rownames(TALON_afterqc_orf_secondpass) <- TALON_afterqc_orf_secondpass$transcript_name_unique

hm <- matrix(t(TALON_afterqc_orf_secondpass2[colnames,"percentage"]), nrow = 1)
colnames(hm) <- colnames
rownames(hm) <- "Percentage"
p <- pheatmap::pheatmap(hm, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
                        cellwidth = 15, cellheight = 15)
ggsave("/data/lib14/project/scLongread/ASHG_coloc_expr_percentage.pdf", p,width = 15, height = 3.2)

tmp <- data.frame(Freq = rowSums(coloc_rst_PPH4),
                  Celltype = celltypes)

tmp$Celltype <- factor(tmp$Celltype, levels = celltypes)
x.axis.bar.df <- tmp
x.axis.bar.df$Var1 <- factor(x.axis.bar.df$Celltype, levels = celltypes)
x.axis.bar.df$Category <- c(rep("Epithelial",8), rep("Immune", 7), rep("Endothelial", 2))
x.axis.bar.df$Category <- factor(x.axis.bar.df$Category , levels = c("Epithelial", "Immune","Endothelial"))
bp.x <- ggplot(data = x.axis.bar.df, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity", aes(fill = Category)) + theme_classic() +
  theme(axis.title.y = element_blank(), axis.text.y =  element_text(size = 12), 
        axis.ticks.y = element_blank(), axis.text.x = element_text(size = 15), 
        axis.title.x = element_text(size = 20, margin = margin(10,0,0,0)),
        legend.position = "none") +
  scale_fill_manual(values = c("#ffd700","#fa8775","#cd34b5"))
bp.x
bp.x.clean <- bp.x +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x  = element_blank(),
        legend.position="none")
bp.x.clean
ggsave("/data/lib14/project/scLongread/ASHG_coloc_count_bar.pdf", bp.x.clean, width = 7, height = 5)

coloc_final$ID_celltype <- paste(coloc_final$transcript_name, coloc_final$Celltype, sep = "_")
DEI_sig$ID_celltype <- paste(DEI_sig$X, DEI_sig$Celltype, sep = "_")
unique(coloc_final$transcript_id[which(coloc_final$ID_celltype %in% DEI_sig$ID_celltype)])
unique(coloc_final$transcript_id)
df_coloc_sum <- data.frame(transcript_id = unique(coloc_final$transcript_id),
                           transcript_name = unique(coloc_final$transcript_name),
                           percentage = hm[1,]
                           )
df_coloc_sum$abundant <- (df_coloc_sum$percentage > 0.144)
df_coloc_sum$gene_name <- NA
tmp_gene_names <- (TALON_afterqc_orf_secondpass2[df_coloc_sum$transcript_name,"annot_gene_name"])
df_coloc_sum[["gene_name"]] <- tmp_gene_names$annot_gene_name
df_coloc_sum$gene_name[39] <- "HLA-DRB1+HLA-DRB5"
df_coloc_sum$gene_name[15] <- "HLA-DQA1+HLA-DQA2"
df_coloc_sum$replicated <- (df_coloc_sum$gene_name %in% c("SECISBP2L", "RBM23", "MUC1","DBNDD2",
                                                          "MRPL21", "MAP1LC3A", "HLA-DRB5","CD55"))
df_coloc_sum$eQTL_coloc <- (df_coloc_sum$transcript_id %in% unique(eQTL_coloc_sub1$phenotype_id.x))
df_coloc_sum$ct_spec <- (df_coloc_sum$transcript_id %in% ct_spec_isoforms)
df_coloc_sum$novel <- grepl("TALON",df_coloc_sum$transcript_id)


lung_func_gwas <- read.table("Need_clean/Lung_func_gwas_association.txt", header = TRUE, sep = "\t")
lung_func_gwas_fev1 <- read.table("Need_clean/Lung_func_fev1_gwas_association.txt", header = TRUE, sep = "\t")
lung_func_gwas <- lung_func_gwas[,c("STRONGEST.SNP.RISK.ALLELE", "REGION")]
lung_func_gwas <- lung_func_gwas[!duplicated(lung_func_gwas),]
colnames(lung_func_gwas)[1] <- "leadSNP"
lung_func_gwas_fev1 <- lung_func_gwas_fev1[,c("STRONGEST.SNP.RISK.ALLELE", "REGION")]
lung_func_gwas_fev1 <- lung_func_gwas_fev1[!duplicated(lung_func_gwas_fev1),]
colnames(lung_func_gwas_fev1)[1] <- "leadSNP"
output.fev1 <- left_join(output.fev1, lung_func_gwas_fev1, by = "leadSNP")
output.ratio <- left_join(output.ratio, lung_func_gwas, by = "leadSNP")
tmppp1 <- output.fev1[,c("transcript_id", "REGION")]
tmppp1 <- tmppp1[!duplicated(tmppp1),]
tmppp2 <- output.ratio[,c("transcript_id", "REGION")]
tmppp2 <- tmppp2[!duplicated(tmppp2),]
tmppp <- rbind(tmppp2,tmppp1)
tmppp <- tmppp[!duplicated(tmppp),]
df_coloc_sum <- left_join(df_coloc_sum, tmppp, by = "transcript_id")
df_coloc_sum <- df_coloc_sum[!duplicated(df_coloc_sum$transcript_id),]
df_coloc_sum$REGION[2:8] <- c("2p11.2","6q22.1","6p21.33","6p21.32","6p21.33","6q21","6p21.33")
df_coloc_sum$eQTL_coloc[which(!(df_coloc_sum$transcript_id %in% eQTL_coloc$phenotype_id.x))] <- "No_matched"
table(df_coloc_sum$eQTL_coloc)
table(df_coloc_sum$eQTL_coloc, df_coloc_sum$replicated)
df_coloc_sum$eQTL_coloc_replicated <- paste(df_coloc_sum$eQTL_coloc, df_coloc_sum$replicated, sep = "-")
table(df_coloc_sum$eQTL_coloc_replicated, df_coloc_sum$abundant)

# TWAS
colnames(hm) <- hm["transcript_name_unique",]
hm <- (TALON_afterqc_orf_secondpass2[unique(TWAS_tb$transcript_name),"percentage"])
names(hm) <- unique(TWAS_tb$transcript_name)
df_twas_sum <- data.frame(transcript_id = unique(TWAS_tb$ID),
                           transcript_name = unique(TWAS_tb$transcript_name),
                           percentage = hm$percentage
)
df_twas_sum$abundant <- (df_twas_sum$percentage > ninetyth_percentile)

tmp_gene_names <- (TALON_afterqc_orf_secondpass[df_twas_sum$transcript_name,"annot_gene_name"])

df_twas_sum$gene_name <-  tmp_gene_names

# Correct gene names for fusion isoform
rownames(df_twas_sum) <- c(1:74)
df_twas_sum$gene_name[15] <- "HLA-DRB1+HLA-DRB5"
df_twas_sum$gene_name[24] <- "HLA-DRB1+HLA-DRB5"
df_twas_sum$gene_name[11] <- "HLA-DQA1+HLA-DQA2"
df_twas_sum$gene_name[59] <- "ATP5MG+KMT2A"
TWAS_previous <- readxl::read_xlsx("TWAS_sum_previous.xlsx")

TWAS_previous$gene_name <- str_split_fixed(TWAS_previous$Phenptype_id, ":", 2)[,1]
TWAS_previous_iso <- intersect(unique(TWAS_previous$gene_name), unique(df_twas_sum$gene_name))
df_twas_sum$replicated <- (df_twas_sum$gene_name %in% TWAS_previous_iso)
eQTL_coloc_sub1 <- subset(eQTL_coloc, PP.H4.abf > 0.7)
df_twas_sum$eQTL_coloc <- (df_twas_sum$transcript_id %in% unique(eQTL_coloc_sub1$phenotype_id.x))
sum(df_twas_sum$eQTL_coloc)
df_twas_sum$ct_spec <- (df_twas_sum$transcript_id %in% ct_spec_isoforms)
df_twas_sum$novel <- grepl("TALON",df_twas_sum$transcript_id)
tmppppp <- TWAS_tb[,c("ID", "Locus")]
tmppppp <- tmppppp[!duplicated(tmppppp),]
colnames(tmppppp)[1]  <- "transcript_id"
df_twas_sum <- left_join(df_twas_sum, tmppppp, by = "transcript_id")
df_twas_sum <- df_twas_sum[!duplicated(df_twas_sum$transcript_id),]
df_twas_sum$Locus_novel <- "Replicated"
df_twas_sum$Locus_novel[which(df_twas_sum$Locus == "")] <- "New"

Novel_sum <- read.table("Novel_isoform_sum.txt", sep = "\t", header = TRUE)

Novel_sum$pct <- Novel_sum$Novel.isoforms/Novel_sum$Total.number
Novel_sum$pct <- round(Novel_sum$pct * 100, 2)
Novel_sum_long <- melt(Novel_sum[,c(1:3)])

ggplot(Novel_sum_long, aes(x=variable, y=value, fill = variable)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("")+ylab("Number of total/novel isoforms") +
  theme(axis.text.x = element_text(color = "black", size = 12,face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        strip.background = element_blank(), 
        legend.justification="right")+facet_wrap(~Study)
ggsave("LR_studies_novel_sum.pdf", width = 10, height = 7)


TWAS_tb$match_CT_id <- paste(TWAS_tb$transcript_name, TWAS_tb$Celltype, sep = ":")
DEI_sig$Celltype <- str_split_fixed(rownames(DEI_sig), "\\.", 2)[,1]
DEI_sig$match_CT_id <- paste(DEI_sig$X, DEI_sig$Celltype, sep = ":")

ct_specific_TWAS <- TWAS_tb$match_CT_id[which(TWAS_tb$match_CT_id %in% DEI_sig$match_CT_id)]
ct_specific_TWAS <- unique(str_split_fixed(ct_specific_TWAS, ":", 2)[,1])
df_twas_sum$ct_spec_DEI <- (df_twas_sum$transcript_name %in% ct_specific_TWAS)
df_twas_sum$eIsoform <- (df_twas_sum$transcript_id %in% unique(final_table$phenotype_id))
table(df_twas_sum$eIsoform, df_twas_sum$ct_spec_DEI)

twas_isoform_nonHLA_loci <- df_twas_sum$transcript_id[!grepl("6p21", df_twas_sum$Locus)]
TWAS_tb_nonHLA_loci <- subset(TWAS_tb, ID %in% twas_isoform_nonHLA_loci)
table(TWAS_tb_nonHLA_loci$Celltype) # 11 in immune, 9 in epithelial, 5 in endothelial, 3 in stromal
