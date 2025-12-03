# Summarize the comparison between 2 SMRT cells and 4 SMRT cells
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

file_dir <- "/data/Choi_lung/scLongreads/B2_percentile/indepth_2_samples/"
sample_dir <- list.dirs(file_dir, recursive = FALSE)
name <- list.dirs(file_dir, recursive = FALSE, full.names = FALSE)
name <- str_split_fixed(name, pattern = "_", n = 2)[,2]


talon3 <- read.table("/data/Choi_lung/scLongreads/TALON_workspace/indepth/output3_talon_read_annot.tsv", header = TRUE)


talon1.list <- split(talon1, talon1$dataset)
talon2.list <- split(talon2, talon1$dataset)
talon.list <- c(talon1.list,talon2.list)
talon <- rbind(talon1,talon2)
rm(talon1)
rm(talon1.list)
rm(talon2)
rm(talon2.list)
talon_transcript_id <- unique(talon$annot_transcript_id)
talon_transcript_id_stat <- sapply(talon_transcript_id, function(x){
  df <- subset(talon, annot_transcript_id == x)
  count <- length(unique(df$dataset))
  return(count)
})
talon_transcript_id_stat <- talon %>%
  group_by(annot_transcript_id) %>%
  summarise(unique_dataset = n_distinct(dataset))
table(talon_transcript_id_stat$unique_dataset)
talon_transcript_id_stat_consistent <- subset(talon_transcript_id_stat, unique_dataset == 22)
talon <- readRDS("/data/Choi_lung/scLongreads/TALON_workspace/test10s/talon_read_annot.RDS")
talon <- talon %>%
  group_by(annot_transcript_id) %>%
  mutate(unique_dataset = n_distinct(dataset))
talon_IRF4<- subset(talon, annot_gene_name == "IRF4")

talon$common_detected <- "specific"
talon$common_detected[which(talon$unique_dataset > 19)] <- "common"
mean(talon$read_length)

library(plyr)
cdat <- ddply(talon, "common_detected", summarise, rating.mean=mean(read_length))
cdat
# Overlaid histograms
ggplot(talon, aes(x=read_length, fill=common_detected)) +
  geom_histogram(binwidth=.5, alpha=.5, position="identity") +
  geom_vline(data=cdat, aes(xintercept=rating.mean,  colour=common_detected),
             linetype="dashed", size=1)


table(talon$transcript_novelty)
table(talon$)
tmp <- talon %>%
  group_by(annot_transcript_id) %>%
  summarise(unique_transcripts = count(read_name))

tmp <-talon %>% distinct(annot_transcript_id, transcript_novelty)


####################
i = 2
name[i]
bc_name_tag <- name
bc_name_tag[1] <- "4_NCI_30_35"
Iso.sum.list <- list()
dataset_list <- list()
for (i in 1:length(name)) {
  filepath_isoforms = paste(file_dir, paste("/Sample_",name[i],"/", name[i], "_classification_filtered", sep = ""), "/isoforms_seurat", sep = "")
  # Seurat
  seur.isoforms.data = Read10X(data.dir = filepath_isoforms,gene.column = 1)
  dim(seur.isoforms.data)
  seur = CreateSeuratObject(counts = seur.isoforms.data, project = bc_name_tag[i], min.cells=5)
  # save original barcodes from PacBio
  seur$orig.barcode <- colnames(seur)
  # reverse complementary DNA of barcodes
  CB.reverse <- chartr(old="ATGC", new="TACG", colnames(seur))
  CB.reverse <- stri_reverse(CB.reverse)
  CB.lr <- paste(str_split_fixed(CB.reverse, "-", n =2)[,2], "-1", sep = "")
  # Rename the cells using the reversed barcodes
  seur <- RenameCells(seur, new.names = CB.lr)
  Iso.annot <-as.data.frame(str_split_fixed(rownames(seur), pattern = ":", n = 2))
  colnames(Iso.annot) <- c("isoform", "gene")
  Iso.info <- read.table(paste(file_dir, paste("/Sample_",name[i],"/", name[i], "_classification.filtered_lite_classification.txt", sep = ""),  sep = ""), sep = "\t", header = TRUE)
  Iso.meta <- left_join(Iso.annot, Iso.info, by = "isoform")
  # gffcompar.rst <- str_split_fixed(gffcompare.tracking[,name[i]], "\\|", n = 7)
  talon.rst <- talon3[grepl(bc_name_tag[i], talon3$dataset),]
  colnames(talon.rst) <- c("isoform",colnames(talon.rst)[-1])
  Iso.meta <- left_join(Iso.meta, talon.rst, by = "isoform")
  
  # uniform and merge transcript names
  Iso.meta$transcript_name_unique <- paste(Iso.meta$gene, Iso.meta$annot_transcript_name, sep = "-")
  Iso.meta$transcript_name_unique[which(Iso.meta$transcript_novelty == "Known")] <- Iso.meta$annot_transcript_name[which(Iso.meta$transcript_novelty == "Known")]
  count_mtx <- seur@assays$RNA@layers$counts
  count_mtx <- count_mtx[-which(is.na(Iso.meta$annot_transcript_name)),]
  Iso.meta <- Iso.meta[-which(is.na(Iso.meta$annot_transcript_name)),]
  
  rownames(count_mtx) <- Iso.meta$transcript_name_unique
  dim(count_mtx)
  dim(Iso.meta)
  # use matrix multiplication to deal with sparse matrix
  group <- Iso.meta$transcript_name_unique %>% fct_inorder()
  group_mat <- sparse.model.matrix(~ 0 + group) %>% t
  # Adjust row names to get the correct final row names
  rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
  count_mtx_sum <- group_mat %*% count_mtx  
  Iso.sum <- Iso.meta %>%
    group_by(transcript_novelty)  %>%
    summarise(FL_counts = sum(FL), 
              numbers = n_distinct(annot_transcript_id),
              read_length_mean = mean(read_length),
              read_length_min = min(read_length),
              read_length_max = max(read_length))
  Iso.sum.list[[i]] <- Iso.sum
  colnames(count_mtx_sum) <- colnames(seur)
  seur[["Isoform"]] <- CreateAssayObject(counts = count_mtx_sum)
  seur@assays$Isoform
  seur@assays$Isoform@meta.features <- Iso.meta
  dataset_list[[i]] <- seur
}
saveRDS(dataset_list, file = "/data/Choi_lung/scLongreads/Seurat/lr.list.isoform.expr.indepth.RDS")
saveRDS(Iso.sum.list, file = "/data/Choi_lung/scLongreads/Seurat/lr.list.isoform.expr.sum.indepth.RDS")
Secondpass.pfam.isoforrm <- read.table("/data/Choi_lung/scLongreads/TALON_workspace/test10s/NCI_lung_secondpass.isoforms_w_known")
dataset_list1 <- readRDS("/data/Choi_lung/scLongreads/Seurat/lr.list.isoform.expr.RDS")
iso.info1 <- list()
for (i in 1:22) {
  iso.add <- dataset_list1[[i]]@assays$Isoform@meta.features
  iso.add <- distinct(iso.add, transcript_name_unique, .keep_all = TRUE)
  rownames(iso.add) <- iso.add$transcript_name_unique
  iso.info1[[i]] <- iso.add
}
iso.info1 <- iso.info1[c(16,11)]

Iso.sum.list1 <- readRDS("/data/Choi_lung/scLongreads/Seurat/lr.list.isoform.expr.sum.RDS")

dataset_list <- readRDS("/data/Choi_lung/scLongreads/Seurat/lr.list.isoform.expr.indepth.RDS")
iso.info <- list()
for (i in 1:2) {
  iso.add <- dataset_list[[i]]@assays$Isoform@meta.features
  iso.add <- distinct(iso.add, transcript_name_unique, .keep_all = TRUE)
  rownames(iso.add) <- iso.add$transcript_name_unique
  iso.info[[i]] <- iso.add
}
cols <- colnames(iso.info[[1]])

TALON_afterqc_orf_secondpass <- readRDS("/data/Choi_lung/scLongreads/Seurat/TALON_afterqc_orf_secondpass_w_known_v3.rds")
length(unique(intersect(iso.info[[1]]$annot_transcript_id, TALON_afterqc_orf_secondpass$annot_transcript_id))) 
length(unique(intersect(iso.info1[[1]]$annot_transcript_id, TALON_afterqc_orf_secondpass$annot_transcript_id))) 
length(unique(intersect(iso.info[[2]]$annot_transcript_id, TALON_afterqc_orf_secondpass$annot_transcript_id))) 
length(unique(intersect(iso.info1[[2]]$annot_transcript_id, TALON_afterqc_orf_secondpass$annot_transcript_id)))
load("/data/Choi_lung/scLongreads/Seurat/InDepth_2sample_comparison.RData")


# Improvement of isoform detection at single-cell batch level
df2 = data.frame(Number = c(73880, 83512,  78183, 89569),
                 Batch = rep(c("NCI4", "NCI19"), each = 2),
                 SMRT = rep(c("2 SMART cells","4 SMART cells"),2))
df2$group <- paste(df2$SMRT, df1$Batch, sep = "_")
df2$group <- factor(df2$group , levels = c("2 SMART cells_NCI4", "4 SMART cells_NCI4",
                                           "2 SMART cells_NCI19", "4 SMART cells_NCI19"))
ggplot(df2, aes(x = group, y = Number, fill = SMRT)) + 
  geom_bar(position="stack", stat="identity")+
  ylab('Number of isoforms detected after QC') + xlab("")+
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position="none") +
  scale_y_continuous(expand = c(0,0))
ggsave("/data/lib14/project/scLongread/Fig_realS3B-2.pdf", width = 7,height = 7)




Saturation_NCI4 <- read.table("/data/Choi_lung/scLongreads/B2_percentile/indepth_2_samples/Sample_04_NCI_30_35/4_NCI_30_35_classification.saturation.txt", sep = "\t",header=TRUE)
Saturation_NCI4$batch <-"NCI4"
Saturation_NCI19 <- read.table("/data/Choi_lung/scLongreads/B2_percentile/indepth_2_samples/Sample_19_NCI_124_129/19_NCI_124_129_classification.saturation.txt", sep = "\t",header=TRUE)
Saturation_NCI19$batch <-"NCI19"

df_saturation <- rbind(Saturation_NCI4, Saturation_NCI19)
df_saturation$reads <- df_saturation$reads/1000000
library(ggplot2)
ggplot(df_saturation, aes(x=reads, y=unique_genes_known, colour=batch)) + 
  geom_line(size=0.8) +
  geom_point(size=1) +
  scale_x_continuous(breaks=seq(0,90,5)) +
  ylab("Number of known genes detected") +
  xlab("Number of reads sampled (millions)") +
  geom_vline(xintercept = c(53498213/1000000, 56160854/1000000), color = c("#fb8072","#80b1d3"),
             linetype = "dashed", linewidth = 1)+
  theme_bw() +
  theme(text = element_text(size=12), legend.title=element_blank())#+facet_wrap(~SampleID, nrow = 5)
ggsave("/data/lib14/project/scLongread/Fig_realS3A-1.pdf", width = 10,height = 7)
ggplot(df_saturation, aes(x=reads, y=unique_isoforms_known, colour=batch)) + 
  geom_line(size=0.8) +
  geom_point(size=1) +
  scale_x_continuous(breaks=seq(0,90,5)) +
  ylab("Number of known isoforms detected") +
  xlab("Number of reads sampled (millions)") +
  geom_vline(xintercept = c(53498213/1000000, 56160854/1000000), color = c("#fb8072","#80b1d3"),
             linetype = "dashed", linewidth = 1)+
  theme_bw() +
  theme(text = element_text(size=12), legend.title=element_blank())
ggsave("/data/lib14/project/scLongread/Fig_realS3A-2.pdf", width = 10,height = 7)
ggplot(df_saturation, aes(x=reads, y=unique_isoforms, colour=batch)) + 
  geom_line(size=0.8) +
  geom_point(size=1) +
  scale_x_continuous(breaks=seq(0,90,5)) +
  ylab("Number of isoforms detected") +
  xlab("Number of reads sampled (millions)") +
  geom_vline(xintercept = c(53498213/1000000, 56160854/1000000), color = c("#fb8072","#80b1d3"),
             linetype = "dashed", linewidth = 1)+
  theme_bw() +
  theme(text = element_text(size=12), legend.title=element_blank())
ggsave("/data/lib14/project/scLongread/Fig_realS3A-3.pdf", width = 10,height = 7)

length(iso.info1[[1]]$annot_transcript_name) # 198,280
length(intersect(iso.info1[[1]]$annot_transcript_name, tmp$annot_transcript_name)) # 179,559

tmp <- rbind(iso.info[[1]][,cols], iso.info[[2]][,cols])

tmp <- distinct(tmp, transcript_name_unique, .keep_all = TRUE)
length(unique(tmp$transcript_name_unique))

TALON_afterqc_orf_secondpass <- readRDS("/data/Choi_lung/scLongreads/Seurat/TALON_afterqc_orf_secondpass_w_known_v3.rds")
TALON_afterqc_orf_secondpass2 <- readRDS("/data/Choi_lung/scLongreads/Seurat/TALON_afterqc_orf_secondpass_w_known_v3.rds")
table(TALON_afterqc_orf_secondpass$transcript_catalog)
lr <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final_v2_wRNAexpr.RDS")
NCI4 <- subset(lr, orig.ident == "4_NCI_30_35")
colnames(NCI4) <- NCI4$orig.barcode
NCI4_4SMRT <- dataset_list[[1]]
NCI4_cells <- intersect((NCI4_4SMRT$orig.barcode), NCI4$orig.barcode)
NCI4 <- subset(NCI4, cells = NCI4_cells)
colnames(NCI4_4SMRT) <- NCI4_4SMRT$orig.barcode
NCI4_4SMRT_qc <- subset(NCI4_4SMRT, cells = NCI4_cells)
count_mtx <- NCI4_4SMRT_qc@assays$Isoform@counts
count_mtx <- count_mtx[which(iso.info[[1]]$annot_transcript_id %in% TALON_afterqc_orf_secondpass$annot_transcript_id),]
NCI4_4SMRT_qc[["Curated_isoform"]] <- CreateAssayObject(counts = count_mtx)
DefaultAssay(NCI4_4SMRT_qc) <- "Curated_isoform"
length(which(colnames(NCI4_4SMRT_qc) %in% colnames(NCI4)))
NCI4_4SMRT_qc$Celltype <- NCI4@meta.data[colnames(NCI4_4SMRT_qc),"Celltype"]
NCI4_4SMRT_qc <- NormalizeData(NCI4_4SMRT_qc)
DotPlot(NCI4_4SMRT_qc, features = c("CAV1-TALONT003058356","CAV1-201"), group.by = "Celltype")
DotPlot(NCI4_4SMRT_qc, features = c("CAV1-TALONT003058356"), group.by = "Celltype")


NCI19 <- subset(lr, orig.ident == "19_NCI_124_129")
colnames(NCI19) <- NCI19$orig.barcode
NCI19_4SMRT <- dataset_list[[2]]
NCI19_cells <- intersect((NCI19_4SMRT$orig.barcode), NCI19$orig.barcode)
NCI19 <- subset(NCI19, cells = NCI19_cells)
colnames(NCI19_4SMRT) <- NCI19_4SMRT$orig.barcode
NCI19_4SMRT_qc <- subset(NCI19_4SMRT, cells = NCI19_cells)
count_mtx <- NCI19_4SMRT_qc@assays$Isoform@counts
count_mtx <- count_mtx[which(iso.info[[2]]$annot_transcript_id %in% TALON_afterqc_orf_secondpass$annot_transcript_id),]
NCI19_4SMRT_qc[["Curated_isoform"]] <- CreateAssayObject(counts = count_mtx)
DefaultAssay(NCI19_4SMRT_qc) <- "Curated_isoform"
length(which(colnames(NCI19_4SMRT_qc) %in% colnames(NCI19)))
NCI19_4SMRT_qc$Celltype <- NCI19@meta.data[colnames(NCI19_4SMRT_qc),"Celltype"]
NCI19_4SMRT_qc <- NormalizeData(NCI19_4SMRT_qc)
DotPlot(NCI19_4SMRT_qc, features = c("CAV1-TALONT003058356", "CAV1-201"), group.by = "Celltype")
DotPlot(NCI19_4SMRT_qc, features = c("CAV1-TALONT003058356"), group.by = "Celltype")
lr.isoform.sub <- readRDS("/data/Choi_lung/scLongreads/Seurat/final/lr_final_w_isoform_expr.RDS")
NCI19_2SMRT <- subset(lr.isoform.sub, orig.ident == "19_NCI_124_129")
colnames(NCI19_2SMRT) <- NCI19_2SMRT$orig.barcode
NCI4_2SMRT <- subset(lr.isoform.sub, orig.ident == "4_NCI_30_35")
colnames(NCI4_2SMRT) <- NCI4_2SMRT$orig.barcode
NCI19_2SMRT_qc <- subset(NCI19_2SMRT, cells = NCI19_cells)
NCI4_2SMRT_qc <- subset(NCI4_2SMRT, cells = NCI4_cells)

count_mtx <- NCI4_2SMRT_qc@assays$Isoform@counts
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass2$transcript_name_unique,]
NCI4_2SMRT_qc[["Curated_isoform"]] <- CreateAssayObject(counts = count_mtx)

count_mtx <- NCI19_2SMRT_qc@assays$Isoform@counts
count_mtx <- count_mtx[TALON_afterqc_orf_secondpass2$transcript_name_unique,]
NCI19_2SMRT_qc[["Curated_isoform"]] <- CreateAssayObject(counts = count_mtx)


DefaultAssay(NCI4_2SMRT_qc) <- "Curated_isoform"
DefaultAssay(NCI19_2SMRT_qc) <- "Curated_isoform"
NCI4_2SMRT_qc <- NormalizeData(NCI4_2SMRT_qc)
NCI19_2SMRT_qc <- NormalizeData(NCI19_2SMRT_qc)

################
# Spot check if DEI come from the shallow sequencing depth
idx1 <- colnames(NCI19_4SMRT_qc)
expr_NCI19_2SMRT <- NCI19_2SMRT_qc@assays$Curated_isoform@data
tmp1 <- expr_NCI19_2SMRT["CAV1-TALONT003058356", idx1]
tmp2 <- NCI19_4SMRT_qc@assays$Curated_isoform@data["CAV1-TALONT003058356",]
tmp <- rbind(tmp1, tmp2)
rownames(tmp) <- c("SMRT_cell_2","SMRT_cell_4")
NCI19_4SMRT_qc[["CAV1_novel"]] <- CreateAssayObject(data = tmp)

idx2 <- colnames(NCI4_4SMRT_qc)
expr_NCI4_2SMRT <- NCI4_2SMRT_qc@assays$Curated_isoform@data
tmp1 <- expr_NCI4_2SMRT["CAV1-TALONT003058356", idx2]
tmp2 <- NCI4_4SMRT_qc@assays$Curated_isoform@data["CAV1-TALONT003058356",]
tmp <- rbind(tmp1, tmp2)
rownames(tmp) <- c("SMRT-cell-2","SMRT-cell-4")
NCI4_4SMRT_qc[["CAV1_novel"]] <- CreateAssayObject(data = tmp)
DefaultAssay(NCI4_4SMRT_qc)  <- "CAV1_novel"
DotPlot(NCI4_4SMRT_qc, features = c("SMRT-cell-2", "SMRT-cell-4"), group.by = "Celltype")
ggsave("/data/lib14/project/scLongread/Fig_realS7A-5.pdf", width = 7, height = 7)

DefaultAssay(NCI19_4SMRT_qc)  <- "CAV1_novel"
DotPlot(NCI19_4SMRT_qc, features = c("SMRT-cell-2", "SMRT-cell-4"), group.by = "Celltype")
ggsave("/data/lib14/project/scLongread/Fig_realS7A-6.pdf", width = 7, height = 7)


DotPlot(NCI4_2SMRT_qc, features = c("CAV1-TALONT003058356"), group.by = "Celltype") + ggtitle("NCI4 2 SMRT cells")
ggsave("/data/lib14/project/scLongread/Fig_realS7A-1.pdf", width = 7, height = 7)
DotPlot(NCI4_4SMRT_qc, features = c("CAV1-TALONT003058356"), group.by = "Celltype")+ ggtitle("NCI4 4 SMRT cells")
ggsave("/data/lib14/project/scLongread/Fig_realS7A-2.pdf", width = 7, height = 7)
DotPlot(NCI19_2SMRT_qc, features = c("CAV1-TALONT003058356"), group.by = "Celltype")+ ggtitle("NCI19 2 SMRT cells")
ggsave("/data/lib14/project/scLongread/Fig_realS7A-3.pdf", width = 7, height = 7)
DotPlot(NCI19_4SMRT_qc, features = c("CAV1-TALONT003058356"), group.by = "Celltype")+ ggtitle("NCI19 4 SMRT cells")
ggsave("/data/lib14/project/scLongread/Fig_realS7A-4.pdf", width = 7, height = 7)

DotPlot(NCI4_2SMRT_qc, features = c("CAV1-201"), group.by = "Celltype")
DotPlot(NCI4_4SMRT_qc, features = c("CAV1-201"), group.by = "Celltype")
DotPlot(NCI19_2SMRT_qc, features = c("CAV1-201"), group.by = "Celltype")
DotPlot(NCI19_4SMRT_qc, features = c("CAV1-201"), group.by = "Celltype")
plot(NCI19_2SMRT_qc@assays$Curated_isoform@data["CAV1-TALONT003058356",NCI19_cells], NCI19_4SMRT_qc@assays$Curated_isoform@data["CAV1-TALONT003058356",NCI19_cells])

df_detail_compare <- data.frame(NCI19_2SMRT = NCI19_2SMRT_qc@assays$Curated_isoform@data["CAV1-TALONT003058356",NCI19_cells],
                                NCI19_4SMRT = NCI19_4SMRT_qc@assays$Curated_isoform@data["CAV1-TALONT003058356",NCI19_cells],
                                Celltype = NCI19_4SMRT_qc@meta.data[NCI19_cells,"Celltype"])

table(df_detail_compare$Celltype[which(df_detail_compare$NCI19_2SMRT == 0 & df_detail_compare$NCI19_4SMRT > 0)])
ggplot(df_detail_compare, aes(x = NCI19_2SMRT, y = NCI19_4SMRT, color = Celltype, shape = Celltype)) + 
  geom_point(size = 3) + scale_shape_manual(values=c(1:25,33:44)) + scale_color_manual(values = cellcolors) + 
  geom_abline(slope=1, intercept = 0) + 
  xlab("Normalized expression of CAV1-TALONT003058356 in sample NCI19 with 2 SMRT cells") +
  ylab("Normalized expression of CAV1-TALONT003058356 in sample NCI19 with 4 SMRT cells") +
  theme(axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x = element_text(color = "black", size = 12),
        strip.background = element_blank(), strip.text = element_text(size = 14)) + coord_equal()
ggsave("/data/lib14/project/scLongread/Fig_realS7B-1.pdf", width = 14, height = 10)


df_detail_compare <- data.frame(NCI4_2SMRT = NCI4_2SMRT_qc@assays$Curated_isoform@data["CAV1-TALONT003058356",NCI4_cells],
                                NCI4_4SMRT = NCI4_4SMRT_qc@assays$Curated_isoform@data["CAV1-TALONT003058356",NCI4_cells],
                                Celltype = NCI4_4SMRT_qc@meta.data[NCI4_cells,"Celltype"])

table(df_detail_compare$Celltype[which(df_detail_compare$NCI4_2SMRT == 0 & df_detail_compare$NCI4_4SMRT > 0)])
ggplot(df_detail_compare, aes(x = NCI4_2SMRT, y = NCI4_4SMRT, color = Celltype, shape = Celltype)) + 
  geom_point(size = 3) + scale_shape_manual(values=c(1:25,33:44)) + scale_color_manual(values = cellcolors) + 
  geom_abline(slope=1, intercept = 0) + 
  xlab("Normalized expression of CAV1-TALONT003058356 in sample NCI4 with 2 SMRT cells") +
  ylab("Normalized expression of CAV1-TALONT003058356 in sample NCI4 with 4 SMRT cells") +
  theme(axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x = element_text(color = "black", size = 12),
        strip.background = element_blank(), strip.text = element_text(size = 14)) + coord_equal()
ggsave("/data/lib14/project/scLongread/Fig_realS7B-2.pdf", width = 14, height = 10)

#################
# Improvement by increasing sequencing depth at single-cell level
median(NCI4_2SMRT_qc@meta.data[NCI4_cells,"nFeature_Curated_isoform"])
median(NCI4_4SMRT_qc@meta.data[NCI4_cells,"nFeature_Curated_isoform"])
median(NCI19_2SMRT_qc@meta.data[NCI19_cells,"nFeature_Curated_isoform"])
median(NCI19_4SMRT_qc@meta.data[NCI19_cells,"nFeature_Curated_isoform"])

df_sum_nIsoform <- data.frame(nIsoform = c(NCI4_2SMRT_qc@meta.data[NCI4_cells,"nFeature_Curated_isoform"], 
                                           NCI4_4SMRT_qc@meta.data[NCI4_cells,"nFeature_Curated_isoform"],
                                           NCI19_2SMRT_qc@meta.data[NCI19_cells,"nFeature_Curated_isoform"],
                                           NCI19_4SMRT_qc@meta.data[NCI19_cells,"nFeature_Curated_isoform"]),
                              nSMRT_cell = c(rep(c("2 SMRT cells","4 SMRT cells"),each = length(NCI4_cells)),
                                             rep(c("2 SMRT cells","4 SMRT cells"),each = length(NCI19_cells))),
                              Batch = c(rep("NCI4", 2*length(NCI4_cells)),rep("NCI19", 2*length(NCI19_cells))))

library(ggpubr)
ggboxplot(df_sum_nIsoform, x = "Batch", y ="nIsoform",
          fill = "nSMRT_cell", palette = c("#fb8072","#80b1d3"),font.x=c(9),font.y=c(9)) + 
  font("legend.title",size=8)+font("legend.text",size=8)

ggsave("/data/lib14/project/scLongread/Fig_realS3C.pdf", width = 7, height = 7)

###############################################################
# Compare expression level of isoforms at individual level

colnames(NCI4_4SMRT_qc) <- paste(NCI4_4SMRT_qc$orig.barcode, NCI4_4SMRT_qc$orig.ident, sep = "_")
NCI4_4SMRT_qc$Sample <- NCI4@meta.data[NCI4_4SMRT_qc$orig.barcode,"Sample"]
table(NCI4_4SMRT_qc$Sample)
colnames(NCI19_4SMRT_qc) <- paste(NCI19_4SMRT_qc$orig.barcode, NCI19_4SMRT_qc$orig.ident, sep = "_")
NCI19_4SMRT_qc$Sample <- NCI19@meta.data[NCI19_4SMRT_qc$orig.barcode,"Sample"]
NCI_4SMRT <- merge(NCI4_4SMRT_qc, NCI19_4SMRT_qc)
count_mtx <- NCI_4SMRT@assays$Curated_isoform@counts
NCI_4SMRT$Celltype_chr <- gsub(" ", "_", as.character(NCI_4SMRT$Celltype))
NCI_4SMRT$Celltype_sample <- paste(NCI_4SMRT$Celltype_chr, NCI_4SMRT$Sample, sep = "-")
group <- as.character(NCI_4SMRT$Celltype_sample) %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx)
count_mtx_sum <- group_mat %*% count_mtx_T  
dim(count_mtx_sum)
sum(count_mtx_sum)
count_mtx_sum <- Matrix::t(count_mtx_sum)
dim(count_mtx_sum)
count_mtx_sum_4SMRT <- count_mtx_sum


colnames(NCI4_2SMRT_qc) <- paste(NCI4_2SMRT_qc$orig.barcode, NCI4_2SMRT_qc$orig.ident, sep = "_")
colnames(NCI19_2SMRT_qc) <- paste(NCI19_2SMRT_qc$orig.barcode, NCI19_2SMRT_qc$orig.ident, sep = "_")

NCI_2SMRT <- merge(NCI4_2SMRT_qc, NCI19_2SMRT_qc)
count_mtx <- NCI_2SMRT@assays$Curated_isoform@counts
NCI_2SMRT$Celltype_chr <- gsub(" ", "_", as.character(NCI_2SMRT$Celltype))
NCI_2SMRT$Celltype_sample <- paste(NCI_2SMRT$Celltype_chr, NCI_2SMRT$Sample, sep = "-")
group <- as.character(NCI_2SMRT$Celltype_sample) %>% fct_inorder()
group_mat <- sparse.model.matrix(~ 0 + group) %>% t
# Adjust row names to get the correct final row names
rownames(group_mat) <- rownames(group_mat) %>% str_extract("(?<=^group).+")
count_mtx_T <- Matrix::t(count_mtx)
count_mtx_sum_2SMRT <- group_mat %*% count_mtx_T  
dim(count_mtx_sum_2SMRT)
sum(count_mtx_sum_2SMRT)
count_mtx_sum_2SMRT <- Matrix::t(count_mtx_sum_2SMRT)
dim(count_mtx_sum_2SMRT)

celltypes <- gsub(" ", "_", levels(lr$Celltype))

Isoform_tested <- intersect(rownames(NCI_4SMRT), rownames(TALON_afterqc_orf_secondpass2)[which(TALON_afterqc_orf_secondpass2$annot_transcript_id %in% subsetlist)])
Isoform_tested <- (intersect(Isoform_tested, rownames(NCI_4SMRT)[which(rowSums(count_mtx_sum_4SMRT) > 0)]))
dim(count_mtx_sum_2SMRT)
count_mtx_sum_2SMRT <- NormalizeData(count_mtx_sum_2SMRT)
dim(count_mtx_sum_4SMRT)
count_mtx_sum_4SMRT <- NormalizeData(count_mtx_sum_4SMRT)
rm(lr)
rm(lr.isoform.sub)
save(list= ls(), file = "/data/Choi_lung/scLongreads/Seurat/Compare_isoform_expr_indv_level.RData")

# Correlation of expression between 2 SMRT cells and 4 SMRT cells at individual level
cor_rst_mtx <- matrix(NA, nrow = length(Isoform_tested), ncol = length(celltypes))
for (i in 1:length(celltypes)) {
  celltype = celltypes[i]
  col_idx <- colnames(count_mtx_sum_2SMRT)[grep(celltype, colnames(count_mtx_sum_2SMRT))]
  rst <- c()
  for (j in Isoform_tested) {
    if(sum(count_mtx_sum_2SMRT[j,col_idx]) == 0 | sum(count_mtx_sum_4SMRT[j,col_idx]) == 0){
      corrst <- NA
      rst <- c(rst,corrst)
    }else{
      corrst <- cor(count_mtx_sum_2SMRT[j,col_idx], count_mtx_sum_4SMRT[j,col_idx])
      # message(j)
      rst <- c(rst,corrst)
    }

  }
  cor_rst_mtx[,i] <- rst
  message("finish ", celltype)
}