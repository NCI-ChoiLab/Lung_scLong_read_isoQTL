# Plotting transcripts

library(magrittr)
library(dplyr)
library(ggtranscript)
library(ggplot2)
library(rtracklayer)
library(geomtextpath)
source("/data/lib14/R/Rscripts/utilities.R")
table(TALON_afterqc_orf_secondpass$structural_category)
idx1 <- which(TALON_afterqc_orf_secondpass2$structural_category == "fusion")
TALON_afterqc_orf_secondpass2$transcript_catalog <- TALON_afterqc_orf_secondpass2$transcript_novelty
TALON_afterqc_orf_secondpass2$transcript_catalog[idx1] <- "fusion"
table(TALON_afterqc_orf_secondpass2$transcript_catalog)
cat.palette = c( "Known"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "fusion"="goldenrod1", "Genomic"="#969696", "Antisense"="#66C2A4")
NB.nominal.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_sig_list.rds")
AT2_nominal <- NB.nominal.sig.list$AT2
Multiciliated_nominal <- NB.nominal.sig.list$Multiciliated
gtf <- rtracklayer::import("/data/Choi_lung/scLongreads/TALON_workspace/test10s/NCI_lung_secondpass.isoforms_w_known.gtf")
# gtf <- rtracklayer::import("/data/Choi_lung/lbl/refdata-gex-GRCh38-2020-A/genes/genes.gtf")

class(gtf)


gtf <- gtf %>% dplyr::as_tibble()

class(gtf)


# filter your gtf for the gene of interest, here "OAS1"
gene_of_interest <- "ENSG00000111671" #SPSB2
gene_of_interest <- "ENSG00000164265"
gene_of_interest <- "ENSG00000089127"
gene_of_interest <- "ENSG00000186501" # TMEM222
gene_of_interest <- "ENSG00000196735" # HLA-DQA1
gene_of_interest <- "ENSG00000185250" # PPIL6
SPSB2_annotation_from_gtf <- gtf %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_id == gene_of_interest
  ) %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_name
    #transcript_biotype
  )

# extract the required annotation columns
SPSB2_annotation_from_gtf <- SPSB2_annotation_from_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_name
    #transcript_biotype
  )

SPSB2_annotation_from_gtf %>% head()


# extract exons
SPSB2_exons <- SPSB2_annotation_from_gtf %>% dplyr::filter(type == "exon")


SPSB2_exons$transcript_name_unique <- SPSB2_exons$transcript_name
SPSB2_exons$transcript_name_unique[grepl("TALON", SPSB2_exons$transcript_name)] <- paste("PPIL6", SPSB2_exons$transcript_name[grepl("TALON", SPSB2_exons$transcript_name)], sep = "-")
SPSB2_exons$transcript_name_unique[grepl("TALON", SPSB2_exons$transcript_name)] <- paste("SPSB2", SPSB2_exons$transcript_name[grepl("TALON", SPSB2_exons$transcript_name)], sep = "-")
SPSB2_exons$transcript_name_unique[grepl("TALON", SPSB2_exons$transcript_name)] <- paste("TMEM222", SPSB2_exons$transcript_name[grepl("TALON", SPSB2_exons$transcript_name)], sep = "-")
SPSB2_exons$transcript_name_unique[grepl("TALON", SPSB2_exons$transcript_name)] <- paste("HLA-DQA1", SPSB2_exons$transcript_name[grepl("TALON", SPSB2_exons$transcript_name)], sep = "-")
# SPSB2_exons$transcript_name_unique[grepl("TALON", SPSB2_exons$transcript_name)] <- paste("SCGB3A2", SPSB2_exons$transcript_name[grepl("TALON", SPSB2_exons$transcript_name)], sep = "-")
SPSB2_exons$transcript_name_unique[grepl("TALON", SPSB2_exons$transcript_name)] <- paste("OAS1", SPSB2_exons$transcript_name[grepl("TALON", SPSB2_exons$transcript_name)], sep = "-")
colnames(SPSB2_exons)[7] <- "annot_transcript_name"
SPSB2_exons <- left_join(SPSB2_exons, TALON_afterqc_orf_secondpass2[,c("annot_transcript_name", "transcript_catalog")], by = "annot_transcript_name")
SPSB2_exons <- SPSB2_exons %>% dplyr::filter(transcript_name %in% c("SPSB2-205", "TALONT000636440", "TALONT000636472", "SPSB2-204"))
SPSB2_exons <- SPSB2_exons %>% filter(transcript_name %in% c("TMEM222-211", "TALONT001958422"))
# SPSB2_exons <- SPSB2_exons %>% filter(transcript_name %in% c("OAS1-201", "OAS1-202", "OAS1-203"))
SPSB2_exons <- SPSB2_exons %>% filter(transcript_name %in% c("HLA-DQA1-203", "TALONT001781255"))
SPSB2_base_plot <-  pknox1_rescaled_exons %>% 
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = transcript_catalog)
  ) +
  geom_intron(
    data = pknox1_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 300
  )

SPSB2_base_plot <-  SPSB2_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = annot_transcript_name
  )) +
  geom_range(
    aes(fill = transcript_catalog)
  ) +
  geom_intron(
    data = to_intron(SPSB2_exons, "annot_transcript_name"),
    aes(strand = strand), 
    arrow.min.intron.length = 3500
  )

SPSB2_base_plot + theme_bw() + 
  ggtitle("SPSB2 gene") +  #theme_bw() + 
  xlab("chr12") + ylab("Isoforms") + 
  scale_fill_manual(values = cat.palette) +
  geom_vline(xintercept = c(6872998), linetype = 1,
             colour = "darkgreen")+
  theme(legend.position = "top", legend.spacing.x  =  unit(0.01, 'cm'), legend.key.size = unit(0.25, 'cm')) 

ggsave("/data/lib14/project/scLongread/FigS11A.pdf", width = 10,height = 10)

# # only for HLA-DQA1 with longer intron

# # extract exons

# pknox1_rescaled <- shorten_gaps(
#   exons = SPSB2_exons, 
#   introns = to_intron(SPSB2_exons, "transcript_name"), 
#   group_var = "transcript_name"
# )

# # shorten_gaps() returns exons and introns all in one data.frame()
# # let's split these for plotting 
# pknox1_rescaled_exons <- pknox1_rescaled %>% dplyr::filter(type == "exon") 
# pknox1_rescaled_introns <- pknox1_rescaled %>% dplyr::filter(type == "intron") 
