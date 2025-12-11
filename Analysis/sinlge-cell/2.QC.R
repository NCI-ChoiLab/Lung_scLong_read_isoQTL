# Quality control

# Bolun Li
# Jan 17 2024


.libPaths()
library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)
library(future)
library(stringr)
library(stringi)

# load data list of long-read sequencing
load("/data/Choi_lung/scLongreads/Seurat/lr.sc.data.list.21samples.beforeQC.RData")
lr.data.list <- dataset_list
names(lr.data.list)
rm(dataset_list)

seur <- lr.data.list$`2_NCI_9_16`
vireo_path <- paste("/data/Choi_lung/scLongreads/CS035134/Sample_",sample,"/v39/vireo_bcftools_het/donor_ids.tsv", sep = "")
vireo.rst <- read.table(vireo_path, sep = "\t", header = TRUE)
vireo.rst$orig.barcode <- paste0(vireo.rst$cell, "-1")
meta.old <- seur@meta.data
meta.new <- left_join(meta.old, vireo.rst, "orig.barcode")
seur <- AddMetaData(seur, meta.new)

library(ggrepel)
library(ggpubr)
ggboxplot(meta.new, x = "donor_id", y = "nFeature_RNA", 
          color = "donor_id", palette = "jco",
          add = "jitter") + 
  geom_hline(yintercept=200,linetype=1,color="darkgreen") + RotatedAxis()


# Add vireo results 
for (sample in names(lr.data.list)) {
  seur <- lr.data.list[[sample]]
  vireo_path <- paste("/data/Choi_lung/scLongreads/B2_percentile/Sample_",sample,"/vireo_bcftools_het/donor_ids.tsv", sep = "")
  vireo.rst <- read.table(vireo_path, sep = "\t", header = TRUE)
  vireo.rst$orig.barcode <- paste0(vireo.rst$cell, "-1")
  meta.old <- seur@meta.data
  meta.new <- left_join(meta.old, vireo.rst, "orig.barcode")
  seur <- AddMetaData(seur, meta.new)
  lr.data.list[[sample]] <- seur
}
save(list = "lr.data.list", file = "/data/Choi_lung/scLongreads/Seurat/lr.sc.data.list.22samples.wVireo.beforeQC.RData")

load("/data/Choi_lung/scLongreads/Seurat/lr.sc.data.list.22samples.wVireo.beforeQC.RData")
seur <- lr.data.list$`3_NCI_23_29`
seur$doublet_RNA
seur$doublet2_RNA = TRUE
seur$doublet2_RNA[which(seur$doublet_RNA < 0.12)] = FALSE
lr.data.list$`3_NCI_23_29` <- seur
lr.seur <- merge(lr.data.list[[1]], y = lr.data.list[2:22], add.cell.ids = names(lr.data.list), project = "Long-read")
sr.meta <- readRDS(file = "/data/Choi_lung/scLongreads/Seurat/annot.meta.rds")
lr.seur$barcode <- colnames(lr.seur)
annot.meta <- left_join(lr.seur@meta.data, sr.meta ,by = "barcode")

# QC criteria
# remove the low QC cells
lr.seur.afterqc <- subset(lr.seur, nFeature_RNA > 200 & percent.mt <25)
ncol(lr.seur.afterqc)
annot.meta.afterqc <- subset(annot.meta, nFeature_RNA.x > 200 & percent.mt.x <25)

# remove doublet based on Scrublet and Vireo
annot.meta.afterqc <- subset(annot.meta.afterqc, doublet2_RNA == FALSE & donor_id.x != "doublet")
length(which(is.na(annot.meta.afterqc$Is_doublet)))
table(annot.meta.afterqc$Is_doublet)
write.table(table(annot.meta.afterqc$donor_id.x, annot.meta.afterqc$Sample), file = "/data/Choi_lung/scLongreads/Seurat/CellNumber_across_individuals_afterqc.txt", 
            sep = "\t", quote = FALSE)
table(annot.meta.afterqc$donor_id.x)
annot.meta.afterqc.unassigned <- subset(annot.meta.afterqc, donor_id.x == "unassigned")
table(annot.meta.afterqc.unassigned$best_doublet)
hist(annot.meta.afterqc.unassigned$prob_doublet.x)
hist(annot.meta.afterqc.unassigned$prob_max.x)
max(annot.meta.afterqc$nFeature_RNA.x)


lr.seur.afterqc <- subset(lr.seur.afterqc, doublet2_RNA == FALSE & donor_id != "doublet")
lr.seur.afterqc$Sample <- lr.seur.afterqc$donor_id
lr.seur.afterqc$Sample[which(lr.seur.afterqc$prob_max > lr.seur.afterqc$prob_doublet)] <- lr.seur.afterqc$best_singlet[which(lr.seur.afterqc$prob_max > lr.seur.afterqc$prob_doublet)]
table(lr.seur.afterqc$Sample)

# remove cells w/o assigned individual
lr.seur.afterqc <- subset(lr.seur.afterqc, Sample != "unassigned")

length(unique(lr.seur.afterqc$Sample))
names(lr.data.list)
write.table(table(lr.seur.afterqc$donor_id, lr.seur.afterqc$orig.ident), file = "/data/Choi_lung/scLongreads/Seurat/CellNumber_across_individuals_samples_afterqc2.txt", 
            sep = "\t", quote = FALSE)
lr.seur.afterqc$Ind_Samp <- paste(lr.seur.afterqc$Sample, lr.seur.afterqc$orig.ident, sep = "+")
df <- as.data.frame(table(lr.seur.afterqc$Ind_Samp))
df$Sample <- str_split_fixed(df$Var1, pattern = "\\+", n = 2)[,2]
ggplot(df, aes(x= Var1, y = Freq, color = Sample,fill = Sample)) +
  geom_bar(stat="identity",width = 0.5)+
  geom_text(aes(label=Freq), vjust=1.6, color="black", size=2)+
  theme_minimal() + theme(axis.text.x=element_text(angle = -90, hjust = 0))
ncol(lr.seur.afterqc)
lr.seur.afterqc$barcode <- colnames(lr.seur.afterqc)
annot.meta.afterqc <- left_join(lr.seur.afterqc@meta.data, sr.meta ,by = "barcode")
length(which(is.na(annot.meta.afterqc$Sample.y)))

# Summarize how many cells also captured in long-read
meta.sum <- annot.meta.afterqc %>% group_by(Ind_Samp) %>% summarise(Gene.median = median(nFeature_RNA.x), 
                                                                     Gene.min = min(nFeature_RNA.x),
                                                                     Gene.max = max(nFeature_RNA.x),
                                                                     Iso.median = median(nFeature_Isoforms),
                                                                     Iso.min = min(nFeature_Isoforms),
                                                                     Iso.max = max(nFeature_Isoforms),
                                                                    srGene.median = median(nFeature_RNA.y, na.rm = TRUE), 
                                                                    srGene.min = min(nFeature_RNA.y, na.rm = TRUE),
                                                                    srGene.max = max(nFeature_RNA.y, na.rm = TRUE))
hist(meta.sum$Gene.median)
which(df$Freq < 1000)

save(list = "lr.seur.afterqc", file = "/data/Choi_lung/scLongreads/Seurat/lr.sc.data.afterQC.RData")

# remove one ever-smoker individual
lr.seur.afterqc.129 <- subset(lr.seur.afterqc, Sample != "SC845484_PC62791_G05")
save(list = "lr.seur.afterqc.129", file = "/data/Choi_lung/scLongreads/Seurat/lr.sc.data.afterQC.129.RData")



