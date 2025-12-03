# New version of violin plot for publication
library(wesanderson)


source("/data/lib14/R/Rscripts/utilities.R")
NB.nominal.sig.list <- readRDS("/data/Choi_lung/scLongreads/jaxqtl/isoQTL_NB_nominal_sig_list.rds")
load("/data/Choi_lung/scLongreads/eisoQTL_ver3.RData")

# prepare genotype data
library(snpStats)
fam <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.fam"
bim <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bim"
bed <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/Genotype/genotype.bed"
sample <- read.plink(bed, bim, fam)
snp_ls <- unique(Nominal_combined$snp)
genotype_mtx <- sample$genotypes[,snp_ls]
snps_mtx <- genotype_mtx@.Data
genotypes <- snps_mtx
snp_info <- sample$map[snp_ls,c(1,2,4,5,6)]
colnames(snp_info) <- c("chr", "rsid", "pos", "alt", "ref")

#
celltypes <- names(NB.nominal.sig.list)
celltype = "DC2"
counts <- count_list[[celltype]][,-c(1:4)]
rownames(counts) <- count_list[[celltype]]$trascript_id
exprs <- expr_list[[celltype]][,-c(1:4)]
rownames(exprs) <- expr_list[[celltype]]$trascript_id
samples <- colnames(counts)
rs = "rs1835471"
transcript = "ENST00000368452"

# count matrix
datainput1 <- data.frame(snp = genotypes[samples,rs],
                         isoform = as.numeric(counts[transcript,]))
datainput2 <- data.frame(snp = genotypes[samples,rs],
                         isoform = as.numeric(exprs[transcript,]))
datainput1$snp_ra <- "0|0"
datainput1$snp_ra[which(datainput1$snp == "01")] <- "1|1"
datainput1$snp_ra[which(datainput1$snp == "02")] <- "0|1"
df <- as.vector(table(datainput1$snp_ra))
datainput1$snp_ra <- paste0("0|0\n", "n=", df[1])
datainput1$snp_ra[which(datainput1$snp == "01")] <- paste0("1|1\n", "n=", df[3])
datainput1$snp_ra[which(datainput1$snp == "02")] <- paste0("0|1\n", "n=", df[2])
sum_res <- datainput1 %>% #
  group_by(snp_ra) %>% 
  summarise(iso_mean = mean(isoform))


p1 <- ggplot(datainput1, aes(x=snp_ra, y=isoform, color=snp_ra)) +
  geom_violin(trim=FALSE,linewidth = 0.75)+
  geom_boxplot(aes(middle=mean(isoform)),
               width=0.1, linewidth = 0.75)+
  geom_smooth(mapping = aes(x = snp_ra, y = isoform, group = 1),formula = y~x, color = "gray",
              method='lm', size = 1, se =TRUE,fill = alpha("gray", .5) ) +
  labs(title="",x=rs, y = paste0("Count of ", transcript))+
  scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1")) +
  ggtitle(paste0(rs, ": ", paste(snp_info[rs,1], snp_info[rs,3], snp_info[rs,5], snp_info[rs,4], "b38", sep = "_")))+
  theme(axis.text.x = element_text(color = "black", size = 12,face = "plain"),
        axis.text.y = element_text(color = "black", size = 12,face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.position="none")

print(p1)
datainput2$snp_ra <- "0|0"
datainput2$snp_ra[which(datainput2$snp == "01")] <- "1|1"
datainput2$snp_ra[which(datainput2$snp == "02")] <- "0|1"
df <- as.vector(table(datainput2$snp_ra))
datainput2$snp_ra <- paste0("0|0\n", "n=", df[1])
datainput2$snp_ra[which(datainput2$snp == "01")] <- paste0("1|1\n", "n=", df[3])
datainput2$snp_ra[which(datainput2$snp == "02")] <- paste0("0|1\n", "n=", df[2])
p2 <- ggplot(datainput2, aes(x=snp_ra, y=isoform, fill=snp_ra)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x=rs, y = paste0("Normalized expression of ", transcript))+
  ggtitle(paste0(rs, ": ", paste(snp_info[rs,1], snp_info[rs,3], snp_info[rs,5], snp_info[rs,4], "b38", sep = "_")))+
  theme(axis.text.x = element_text(color = "black", size = 12,face = "plain"),
        axis.text.y = element_text(color = "black", size = 12,face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1))
Isoform_info <- TALON_afterqc_orf_secondpass2[, c("annot_gene_id", "annot_transcript_id", "transcript_name_unique")]
save(list = c("Nominal_combined", "expr_list", "count_list", "genotypes", "snp_info","Isoform_info"),
     file = "/data/Choi_lung/scLongreads/eisoQTL_ver4.RData")
isoQTL_plot_pub <- function(celltype = "Multiciliated",
                            rs = "rs1835471",
                            transcript = "ENST00000368452",
                            return_count = FALSE
                            ){
  transcript_name <- Search_transcript_name2(transcript)
  counts <- count_list[[celltype]][,-c(1:4)]
  rownames(counts) <- count_list[[celltype]]$trascript_id
  exprs <- expr_list[[celltype]][,-c(1:4)]
  rownames(exprs) <- expr_list[[celltype]]$trascript_id
  samples <- colnames(counts)
  
  # count matrix
  datainput1 <- data.frame(snp = genotypes[samples,rs],
                           isoform = as.numeric(counts[transcript,]))
  datainput2 <- data.frame(snp = genotypes[samples,rs],
                           isoform = as.numeric(exprs[transcript,]))
  datainput1$snp_ra <- "0|0"
  datainput1$snp_ra[which(datainput1$snp == "01")] <- "1|1"
  datainput1$snp_ra[which(datainput1$snp == "02")] <- "0|1"
  df <- as.vector(table(datainput1$snp_ra))
  datainput1$snp_ra <- paste0("0|0\n", "n=", df[1])
  datainput1$snp_ra[which(datainput1$snp == "01")] <- paste0("1|1\n", "n=", df[3])
  datainput1$snp_ra[which(datainput1$snp == "02")] <- paste0("0|1\n", "n=", df[2])
  p1 <- ggplot(datainput1, aes(x=snp_ra, y=isoform, color=snp_ra)) +
    geom_violin(trim=FALSE,linewidth = 0.75)+
    geom_boxplot(aes(middle=mean(isoform)),
                 width=0.1, linewidth = 0.75)+
    geom_smooth(mapping = aes(x = snp_ra, y = isoform, group = 1),formula = y~x, color = "gray",
                method='lm', size = 1, se =TRUE,fill = alpha("gray", .5) ) +
    labs(title="",x=rs, y = paste0("Count of ", transcript))+
    scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1")) +
    ggtitle(paste0(celltype, "\n",transcript_name ,  rs, ": ", paste(snp_info[rs,1], snp_info[rs,3], snp_info[rs,5], snp_info[rs,4], "b38", sep = "_")))+
    theme(axis.text.x = element_text(color = "black", size = 12,face = "plain"),
          axis.text.y = element_text(color = "black", size = 12,face = "plain"),  
          axis.title.x = element_text(color = "black", size = 16, face = "plain"),
          axis.title.y = element_text(color = "black", size = 16, face = "plain"),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent'),
          axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
          legend.position="none")
  datainput2$snp_ra <- "0|0"
  datainput2$snp_ra[which(datainput2$snp == "01")] <- "1|1"
  datainput2$snp_ra[which(datainput2$snp == "02")] <- "0|1"
  df <- as.vector(table(datainput2$snp_ra))
  datainput2$snp_ra <- paste0("0|0\n", "n=", df[1])
  datainput2$snp_ra[which(datainput2$snp == "01")] <- paste0("1|1\n", "n=", df[3])
  datainput2$snp_ra[which(datainput2$snp == "02")] <- paste0("0|1\n", "n=", df[2])
  p2 <- ggplot(datainput2, aes(x=snp_ra, y=isoform, color=snp_ra)) +
    geom_violin(trim=FALSE,linewidth = 0.75)+
    geom_boxplot(aes(middle=mean(isoform)),
                 width=0.1, linewidth = 0.75)+
    geom_smooth(mapping = aes(x = snp_ra, y = isoform, group = 1),formula = y~x, color = "gray",
                method='lm', size = 1, se =TRUE,fill = alpha("gray", .5) ) +
    labs(title="",x=rs, y = paste0("Normalized expression of ", transcript))+
    scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1")) +
    ggtitle(paste0(celltype, "\n",transcript_name , ": ",rs, ": ", paste(snp_info[rs,1], snp_info[rs,3], snp_info[rs,5], snp_info[rs,4], "b38", sep = "_")))+
    theme(axis.text.x = element_text(color = "black", size = 12,face = "plain"),
          axis.text.y = element_text(color = "black", size = 12,face = "plain"),  
          axis.title.x = element_text(color = "black", size = 16, face = "plain"),
          axis.title.y = element_text(color = "black", size = 16, face = "plain"),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent'),
          axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
          legend.position="none")
  
  if(return_count){
    return(list(p1,p2))
  }else{
    return(p2)
  }

}
p1 <- isoQTL_plot_pub(celltype = "Alveolar_macrophages", rs = "rs896962", transcript = "ENST00000449291")
p2 <- isoQTL_plot_pub(celltype = "Multiciliated", rs = "rs896962", transcript = "ENST00000449291")
p3 <- isoQTL_plot_pub(celltype = "Club", rs = "rs7118835", transcript = "ENST00000534397")
p4 <- isoQTL_plot_pub(celltype = "Goblet", rs = "rs7118835", transcript = "ENST00000534397")
print(p1)
print(p2)
ggsave("/data/Choi_lung/scLongreads/DEI/plots/Fig3H1.pdf", p1 ,device = cairo_pdf, width = 6, height = 6)
ggsave("/data/Choi_lung/scLongreads/DEI/plots/Fig3H2.pdf", p2 ,device = cairo_pdf, width = 6, height = 6)
ggsave("/data/Choi_lung/scLongreads/DEI/plots/Fig3G1.pdf", p3 ,device = cairo_pdf, width = 6, height = 6)
ggsave("/data/Choi_lung/scLongreads/DEI/plots/Fig3G2.pdf", p4 ,device = cairo_pdf, width = 6, height = 6)


p5 <- isoQTL_plot_pub(celltype = "Alveolar_macrophages", rs = "rs10774671", transcript = "ENST00000202917")
p6 <- isoQTL_plot_pub(celltype = "Alveolar_macrophages", rs = "rs10774671", transcript = "ENST00000452357")
p7 <- isoQTL_plot_pub(celltype = "Alveolar_macrophages", rs = "rs10774671", transcript = "ENST00000445409")
ggsave("/data/lib14/project/scLongread/Fig5F2.pdf", p5 ,device = cairo_pdf, width = 6, height = 6)
ggsave("/data/lib14/project/scLongread/Fig5F3.pdf", p6 ,device = cairo_pdf, width = 6, height = 6)
ggsave("/data/lib14/project/scLongread/Fig5F4.pdf", p7 ,device = cairo_pdf, width = 6, height = 6)



# plot additional isoforms for SPSB2

p8 <- isoQTL_plot_pub(celltype = "Multiciliated", rs = "rs11064437", transcript = "ENST00000524270") # SPSB2-205
p9 <- isoQTL_plot_pub(celltype = "Multiciliated", rs = "rs11064437", transcript = "TALONT000636472") # SPSB2 novel
p10 <- isoQTL_plot_pub(celltype = "Multiciliated", rs = "rs11064437", transcript = "ENST00000523102") # SPSB2-204

ggsave("/data/lib14/project/scLongread/FigS8A.pdf", p8 ,device = cairo_pdf, width = 6, height = 6)
ggsave("/data/lib14/project/scLongread/FigS8B.pdf", p9 ,device = cairo_pdf, width = 6, height = 6)
ggsave("/data/lib14/project/scLongread/FigS8C.pdf", p10 ,device = cairo_pdf, width = 6, height = 6)

isoQTL_plot_pub(celltype = "Multiciliated", rs = "rs12528822", transcript = "ENST00000521072") 
ggsave("/data/lib14/project/scLongread/FigS11D.pdf", device = cairo_pdf, width = 6, height = 6)
