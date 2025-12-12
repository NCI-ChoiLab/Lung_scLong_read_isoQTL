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

# load transcript metadata generated from TALON and integrated with TranDecoder results
TALON_afterqc_orf_secondpass <- readRDS("/data/Choi_lung/scLongreads/Seurat/TALON_afterqc_orf_secondpass_w_known.rds")
TALON_afterqc_orf_secondpass2 <- readRDS("/data/Choi_lung/scLongreads/Seurat/TALON_afterqc_orf_secondpass_w_known_v3.rds")
table(TALON_afterqc_orf_secondpass$structural_category)
idx1 <- which(TALON_afterqc_orf_secondpass2$structural_category == "fusion")
TALON_afterqc_orf_secondpass2$transcript_catalog <- TALON_afterqc_orf_secondpass2$transcript_novelty
TALON_afterqc_orf_secondpass2$transcript_catalog[idx1] <- "fusion"
table(TALON_afterqc_orf_secondpass2$transcript_catalog)
gtex_compare <- read.csv("/data/lib14/project/scLongread/transcript_compare.txt", 
                         row.names = 1, header = FALSE)

colnames(gtex_compare) <- c("Overlapped/GENCODE_V32", "Overlapped/scLong-read", "Shared/GTEx", "Shared/scLong-read")
gtex_compare$Level <- rownames(gtex_compare)


library(reshape2)
df <- melt(gtex_compare[c(4,5), c(1,2,5)])

ggplot(df, aes(x = variable, y = value, fill = Level)) +
  geom_col(position = "dodge") +
  geom_text(
    aes(label = value),
    colour = "black", size = 3,
    vjust = 1.5, position = position_dodge(.9)
  )+theme_bw()+
  theme(axis.title = element_text(size = 14),        
        axis.text = element_text(colour = 'black',size = 14)
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) 
ggsave("project/scLongread/transcript_compare.pdf", width = 5,height = 6)
lr.isoform.sub <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final_w_isoform_expr.RDS")
table(lr.isoform.sub$Sample)
# lr.isoform.sub$orig.ident <- lr$orig.ident
df_isoform <- as.data.frame.table(table(TALON_afterqc_orf_secondpass2$transcript_catalog))
count_mtx <- lr.isoform.sub@assays$Isoform@counts
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass2$transcript_name_unique,]
lr.isoform.sub$orig.ident <- as.character(lr.isoform.sub$orig.ident)
group <- lr.isoform.sub$orig.ident %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx)
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
sum(count_mtx_sum)
count_mtx_sum <- Matrix::t(count_mtx_sum)
dim(count_mtx_sum)

group <- TALON_afterqc_orf_secondpass2$transcript_catalog %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
count_sum <- group_mat %*% count_mtx_sum  
dim(count_sum)
count_sum <- as.matrix(count_sum)
df_isoform2 <- as.data.frame(t(count_sum))
df_long <- melt(df_isoform2)
df_long$value <- log10(df_long$value + 1)
df_long$variable <- gsub("group", "", df_long$variable)
library(ggrepel)
library(ggpubr)
cat.palette = c( "Known"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "fusion"="goldenrod1", "Genomic"="#969696", "Antisense"="#66C2A4")
df_long$variable <- factor(df_long$variable, levels = c("Known", "ISM", "NIC", 
                                                        "NNC", "fusion","Genomic", "Antisense"))
ggboxplot(df_long, x = "variable", y = "value", color = "variable")+
  scale_color_manual(values = cat.palette) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("/data/lib14/project/scLongread/Fig2B.pdf", device = cairo_pdf, width = 6,height = 8)

df_isoform$Var1 <- factor(df_isoform$Var1, levels = c("Known", "ISM", "NIC", 
                                                      "NNC","fusion", "Genomic", "Antisense"))

ggplot(df_isoform, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = cat.palette) +
  geom_text(aes(label = Freq), vjust=-1) + 
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),  
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_text(color = "black", size = 16),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10))
ggsave("/data/lib14/project/scLongread/Fig2C.pdf", width = 8,height = 8)

count_mtx <- lr.isoform.sub@assays$Isoform@counts
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass2$transcript_name_unique,]
group <- as.character(lr.isoform.sub$Sample) %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx)
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
sum(count_mtx_sum)
count_mtx_sum <- Matrix::t(count_mtx_sum)
dim(count_mtx_sum)
indv_expr_tf <- (count_mtx_sum > 0)
df_density <- data.frame(isoform = rownames(indv_expr_tf),indv_n = rowSums(indv_expr_tf))
df_density$Novelty <- "Known"
df_density$Novelty[grepl("TALON", df_density$isoform)] <- "Novel"
ggplot(df_density, aes(x=indv_n, fill=Novelty)) +
  geom_histogram(binwidth=0.5, alpha=.8, position="identity") + ylab("Number of isoforms")+xlab("Number of individuals")+
  scale_fill_manual(values = c("#6BAED6", "#FC8D59"))+
  facet_wrap(~Novelty, scales = "free")+
  theme_bw()+
  theme(axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        strip.text = element_text(size = 16, colour = "black"))
ggsave("/data/lib14/project/scLongread/Fig2F.pdf",height = 7,width = 15)

# number of novel and known isoforms per gene
tmp <- as.matrix.data.frame(table(TALON_afterqc_orf_secondpass2$annot_gene_id, TALON_afterqc_orf_secondpass2$transcript_catalog))
colnames(tmp) <- names(table(TALON_afterqc_orf_secondpass2$transcript_catalog))
df_niso_per_gene <- data.frame(total_n = rowSums(tmp),
                               known_n = tmp[,"Known"])
df_niso_per_gene$annot_gene_id <- names(table(TALON_afterqc_orf_secondpass2$annot_gene_id))
df_niso_per_gene$novel_n <- df_niso_per_gene$total_n - df_niso_per_gene$known_n
summary(df_niso_per_gene)
quantile(df_niso_per_gene$known_n, 0.95)
quantile(df_niso_per_gene$novel_n, 0.95)
df_niso_per_gene$novel_cata <- df_niso_per_gene$novel_n
df_niso_per_gene$novel_cata[which(df_niso_per_gene$novel_n >= 6 & df_niso_per_gene$novel_n <= 8)] = 6
df_niso_per_gene$novel_cata[which(df_niso_per_gene$novel_n >= 9 & df_niso_per_gene$novel_n <= 11)] = 7
df_niso_per_gene$novel_cata[which(df_niso_per_gene$novel_n >= 12 & df_niso_per_gene$novel_n <= 14)] = 8
df_niso_per_gene$novel_cata[which(df_niso_per_gene$novel_n >= 15 & df_niso_per_gene$novel_n <= 17)] = 9
df_niso_per_gene$novel_cata[which(df_niso_per_gene$novel_n >= 18 & df_niso_per_gene$novel_n <= 24)] = 10
df_niso_per_gene$novel_cata[which(df_niso_per_gene$novel_n >= 25 & df_niso_per_gene$novel_n <= 31)] = 11
df_niso_per_gene$novel_cata[which(df_niso_per_gene$novel_n >= 32 & df_niso_per_gene$novel_n <= 38)] = 12
df_niso_per_gene$novel_cata[which(df_niso_per_gene$novel_n >= 39 & df_niso_per_gene$novel_n <= 45)] = 13
df_niso_per_gene$novel_cata[which(df_niso_per_gene$novel_n >= 46)] = 14
table(df_niso_per_gene$novel_cata)

df_niso_per_gene$known_cata<- df_niso_per_gene$known_n
df_niso_per_gene$known_cata[which(df_niso_per_gene$known_n >=7)]<- 7

mtx_niso_per_gene <- as.matrix.data.frame(table(df_niso_per_gene$known_cata, df_niso_per_gene$novel_cata))
mtx_niso_per_gene <- log2(mtx_niso_per_gene+1)
mtx_niso_per_gene <- t(mtx_niso_per_gene)
rownames(mtx_niso_per_gene) <- c("0","1","2","3","4","5","6-8","9-11","12-14","15-17","18-24","25-31","32-38","39-45","46+")
colnames(mtx_niso_per_gene) <- c("0","1","2","3","4","5","6","7+")
# mtx_niso_per_gene <- mtx_niso_per_gene[rev(c("0","1","2","3","4","5","6-8","9-11","12-14","15-17","18-24","25-31","32-38","39-45","46+")),]
mtx_niso_per_gene <- mtx_niso_per_gene[,rev(c("0","1","2","3","4","5","6","7+"))]
library(ComplexHeatmap)
library(RColorBrewer)

hm<-pheatmap(t(mtx_niso_per_gene), cluster_rows = FALSE, cluster_cols = FALSE,
         row_names_side = c("left"),
         color = (brewer.pal(7,"YlGnBu")), scale = "none",border_color = "white",
         cellwidth = 25,cellheight = 25)

pdf("/data/lib14/project/scLongread/Fig2E.pdf")
print(hm)
dev.off()

gtf <- rtracklayer::import("/data/Choi_lung/scLongreads/TALON_workspace/test10s/NCI_lung_secondpass.isoforms_w_known.gtf")
gtf <- gtf %>% dplyr::as_tibble()
gtf_exon <- gtf %>% filter(type == "exon")
gtf_sum <- gtf_exon %>% group_by(gene_id) %>% summarise(exon_num = length(unique(exon_id)))
cor(df_niso_per_gene$total_n, gtf_sum$exon_num)
cor(df_niso_per_gene$novel_n, gtf_sum$exon_num)
df_niso_per_gene$exon_n <- gtf_sum$exon_num
library(ggpubr)
library(ggbreak)
ggplot(df_niso_per_gene,aes(x=exon_n,y=novel_n))+geom_point() +
  stat_smooth(method=lm,formula=y~x) + 
  stat_cor(label.y = 1100,aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  scale_y_break(c(500,1000)) +
  ylab("Number of novel isoforms")+
  xlab("Number of exons")+theme_classic()+
  ggtitle("Association between numbers of isoforms and exons within gene")+
  NoLegend()


setwd("/data/Choi_lung/scLongreads/DEI/")
df_AS_long <- read.table("/data/Choi_lung/scLongreads/DEI/Fig2F.table.txt", sep = "\t",
                         header = TRUE)
df_AS_long$variable <- factor(df_AS_long$variable, levels = c("Novel", "Known"))
legendLabelF1 <- levels((df_AS_long$variable))
ggplot(df_AS_long, aes(x=AS, y=value, fill = AS, alpha = variable)) +
  geom_bar(position="stack", stat="identity") + ylab("Number of alternative splicing events")+xlab("AS catagory") + 
  scale_alpha_manual(values=c(1,0.5),
                     name = "Novelty",
                     labels = legendLabelF1, guide = "none") +
  geom_text(aes(label = Label), vjust = -1) + 
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
ggsave("/data/lib14/project/scLongread/Fig2G.pdf", width = 8,height = 7)
