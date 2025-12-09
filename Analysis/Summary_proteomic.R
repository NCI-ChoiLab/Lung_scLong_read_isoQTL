# Summarize the validated isoforms at the peptide level
# Integrate 4 samples 

# Author: Bolun

library("seqinr")
library(stringi)
library(stringr)
library(dplyr)
fastafile<- read.fasta(file = "Need_clean/NCI_lung_eIsoforms.fasta", 
                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
cds_fastafile<- read.fasta(file = "/data/Choi_lung/scLongreads/TALON_workspace/test10s/transcripts.fasta.transdecoder.cds", 
                       seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
NB.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_sig_list.rds")
eIsoforms <- unique(do.call(rbind, isoQTL_sig_list)$phenotype_id)
pro_tested <- unique(str_split_fixed(names(fastafile), "\\.", 2)[,1])
sum(grepl("TALONT", pro_tested))
sum(grepl("TALONT", subsetlist))
18697/length(subsetlist)
NCI134 <- read.table("Need_clean/NCI134.txt", header = TRUE, sep = "\t", comment.char = "")
NCI134 <- subset(NCI134, Protein.FDR.Confidence..Combined == "High" )
NCI135 <- read.table("Need_clean/NCI135.txt", header = TRUE, sep = "\t", comment.char = "")
NCI135 <- subset(NCI135, Protein.FDR.Confidence..Combined == "High" )
NCI138 <- read.table("Need_clean/NCI138.txt", header = TRUE, sep = "\t", comment.char = "")
NCI138 <- subset(NCI138, Protein.FDR.Confidence..Combined == "High" )
NCI139 <- read.table("Need_clean/NCI139.txt", header = TRUE, sep = "\t", comment.char = "")
NCI139 <- subset(NCI139, Protein.FDR.Confidence..Combined == "High" )
all_val_pro <- unique(Reduce(union, list(NCI134$Accession, NCI135$Accession, NCI138$Accession, NCI139$Accession)))

pro_detection <- data.frame(Accession = all_val_pro,
                            NCI134 = (all_val_pro %in% NCI134$Accession),
                            NCI135 = (all_val_pro %in% NCI135$Accession),
                            NCI138 = (all_val_pro %in% NCI138$Accession),
                            NCI139 = (all_val_pro %in% NCI139$Accession))
idx <- which(rowSums(pro_detection[,c(2:5)]) >= 2)
pro_detection$MoreThanOne <- rowSums(pro_detection[,c(2:5)]) >= 2
pro_detection$MoreThanTwo <- rowSums(pro_detection[,c(2:5)]) >= 3
pro_detection$All <- rowSums(pro_detection[,c(2:5)]) == 4
pro_detection$transcript_id <- str_split_fixed(pro_detection$Accession, pattern = "\\.", 2)[,1]
pro_detection <- pro_detection %>% group_by(transcript_id) %>% mutate(transcript_name = Search_transcript_name2(transcript_id))
sum(pro_detection$MoreThanOne)/length(subsetlist)

pro_detection <- left_join(pro_detection, NCI134[,c("Accession","X..Unique.Peptides")], by = "Accession")
colnames(pro_detection)[11] <- "Num_unique_peptides_NCI134"
pro_detection <- left_join(pro_detection, NCI135[,c("Accession","X..Unique.Peptides")], by = "Accession")
colnames(pro_detection)[12] <- "Num_unique_peptides_NCI135"
pro_detection <- left_join(pro_detection, NCI138[,c("Accession","X..Unique.Peptides")], by = "Accession")
colnames(pro_detection)[13] <- "Num_unique_peptides_NCI138"
pro_detection <- left_join(pro_detection, NCI139[,c("Accession","X..Unique.Peptides")], by = "Accession")
colnames(pro_detection)[14] <- "Num_unique_peptides_NCI139"
pro_detection <- pro_detection %>% group_by(Accession) %>% mutate(Max_Num_unique_peptides = max(c(Num_unique_peptides_NCI134,
                                                                                 Num_unique_peptides_NCI135,
                                                                                 Num_unique_peptides_NCI138,
                                                                                 Num_unique_peptides_NCI139), na.rm = TRUE))
pro_detection$Novelty <- "known"
pro_detection$Novelty[grepl("TALONT", pro_detection$Accession)] <- "Novel"
table(pro_detection$Novelty)
pro_detection_only_one_unique_peptide <- subset(pro_detection, Max_Num_unique_peptides < 2)
pro_detection_hc <- subset(pro_detection, Max_Num_unique_peptides >= 2 | MoreThanOne == TRUE)
eIsoform_pro_val <- unique(eIsoforms[which(eIsoforms %in% pro_detection_hc$transcript_id)])
sum(grepl("TALONT", eIsoforms))
sum(grepl("TALONT", eIsoform_pro_val))
table(pro_detection_only_one_unique_peptide$Novelty)
cds_fastafile <- cds_fastafile[pro_detection$Accession]

# Novel isoform
sum(grepl("TALONT", pro_detection_hc$Accession))/18697 # 8.7%
sum(grepl("TALONT", pro_detection_hc$Accession[pro_detection_hc$MoreThanOne]))/18697 # 8.0%
sum(grepl("TALONT", pro_detection_hc$Accession[pro_detection_hc$All]))/18697 # 5.6%

# Known isoform
sum(!grepl("TALONT", pro_detection_hc$Accession))/(length(subsetlist) - 18697)  # 22.8%
sum(!grepl("TALONT", pro_detection_hc$Accession[pro_detection_hc$MoreThanOne]))/(length(subsetlist) - 18697)  # 21.5%
sum(!grepl("TALONT", pro_detection_hc$Accession[pro_detection_hc$All]))/(length(subsetlist) - 18697)  # 16.6%



df_overall <- data.frame(pct = c(8.7, 8.0, 5.6, 22.8, 21.5, 16.6), 
                         Novelty = rep(c("Novel", "Known"), each = 3),
                         detection = rep(c("Any", "More than one", "All"), 2))
df_overall$detection <- factor(df_overall$detection, levels = c("Any", "More than one", "All"))
ggplot(subset(df_overall, detection != "More than one"), aes(x = detection, y = pct, fill = Novelty)) +
  geom_col(position = "dodge") +
  geom_text(
    aes(label = pct),
    colour = "black", size = 4.5,
    vjust = 1.5, position = position_dodge(.9)
  )+ylab("Percentage of validated isoforms")+ ggtitle("Validation of 42,824 isoforms")+xlab("")+
  theme(axis.text.y = element_text(color = "black", size = 12, angle = 0, face = "plain"),
        
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        
        axis.text.x = element_text(color = "black", size = 12),
        strip.background = element_blank(), strip.text = element_text(size = 14))

ggsave("/data/lib14/project/scLongread/Fig2G.pdf", width = 7,height = 7)

# External mass-spec validation
# Glinos et al 2022
GTEx_pro_val <- read.table("/data/lib14/R/Rscripts/Glinos_TS3.txt", sep = "\t", header = TRUE)
GTEx_pro_val_sub <- subset(GTEx_pro_val, Mass_Spec == "Yes")
gff_combined_sub <- subset(gff_combined, class_code == "=")
length(which(gff_combined_sub$cmp_ref %in% GTEx_pro_val_sub$Transcript))
validated_overlaps <- unique(GTEx_pro_val_sub$Transcript[which(GTEx_pro_val_sub$Transcript %in% gff_combined_sub$cmp_ref)])
View(GTEx_pro_val_sub[which(GTEx_pro_val_sub$Transcript %in% validated_overlaps),])
length(validated_overlaps) # 756
gff_combined_sub_proval <- subset(gff_combined_sub, cmp_ref %in% GTEx_pro_val_sub$Transcript)
length(which(gff_combined_sub_proval$transcript_id %in% pro_detection$transcript_id))

novel_isoform_val_bytwo <- unique(pro_detection_hc$transcript_id)[which(unique(pro_detection_hc$transcript_id) %in% gff_combined_sub_proval$transcript_id)]

lungfun_coloc_proval <- c("ENST00000344691","TALONT001781255","TALONT003093589","ENST00000296930","ENST00000616870","ENST00000338631", "ENST00000362034","ENST00000367063","TALONT001769911")
sapply(lungfun_coloc_proval, function(x) Search_transcript_name2(x))
lungcancer_TWAS_proval <- unique(TWAS_tb$ID[which(TWAS_tb$ID %in% pro_detection_hc$transcript_id)])
sapply(lungcancer_TWAS_proval, function(x) Search_transcript_name2(x))
