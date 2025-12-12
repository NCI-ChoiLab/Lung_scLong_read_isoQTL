setwd("/Users/lib14/Desktop/scLongread/")


library(locuscomparer)
isoQTL_fn = 'PPIL6_isoQTL_fn.tsv'
eQTL_fn = 'PPIL6_eQTL_fn.tsv'
gwas_fn = 'GWAS_6q21_byun.tsv'

locuscompare(in_fn1 = eQTL_fn, in_fn2 = isoQTL_fn, title = 'eQTL', title2 = 'isoQTL',
             population = "EAS", snp = "rs7169304")
library(ggplot2)

p1 <- locuscompare(in_fn1 = gwas_fn, in_fn2 = isoQTL_fn, title = 'GWAS', title2 = 'isoQTL',
             population = "EAS", snp = "rs12528822")
p2 <- locuscompare(in_fn1 = gwas_fn, in_fn2 = eQTL_fn, title = 'GWAS', title2 = 'eQTL',
             population = "EAS", snp = "rs12528822")
p3 <- locuscompare(in_fn1 = isoQTL_fn, in_fn2 = gwas_fn, title = 'isoQTL', title2 = 'GWAS',
                   population = "EAS", snp = "rs12528822")
ggsave("./figures/Fig6E1.pdf", p1, width = 12, height = 8)
ggsave("./figures/Fig6E2.pdf", p2, width = 12, height = 8)
ggsave("./figures/Fig6E3.pdf", p3, width = 12, height = 8)

# Summarize CCVs of isoQTls and eQTLs
mytoken = "149dc339070c"
PPIL6_isoQTL <- read.table("PPIL6_isoQTL.txt", header = TRUE, sep = "\t")
PPIL6_eQTL <- read.table("PPIL6_eQTL.txt", header = TRUE, sep = "\t")
PPIL6_eQTL_tmp <- read.table("PPIL6_eQTL.txt", header = TRUE, sep = "\t")
PPIL6_eQTL <- read.table("", header = TRUE, sep = "\t")
test_tmp <- LDmatrix(head(PPIL6_isoQTL$rsnum[order(PPIL6_isoQTL$pvalue)] , n =100), pop = "EAS",token = mytoken, genome_build = "grch38")

df_test_tmp <- test_tmp[,c("RS_number","rs12528822")]
colnames(df_test_tmp)[1] <- "rsnum"
rs_list <- intersect(colnames(genotypes), PPIL6_isoQTL$rsnum)
PPIL6_isoQTL <- left_join(PPIL6_isoQTL, df_test_tmp, by = "rsnum")

# isoQTL lead LD
tb.list <- lapply(PPIL6_isoQTL$rsnum, function(x){
  tryCatch({
    tb <- LDpair("rs12528822", x, pop = "EAS", token = mytoken, genome_build = "grch38")
  }, error=function(e){})
})

df_tmp <- as.data.frame(matrix(as.character(NA), nrow = 1, ncol = 18))

colnames(df_tmp) <- colnames(tb.list[[1]])
tb.list.tmp <- lapply(tb.list, function(x){
  if(is.null(x)){return(df_tmp)}else{return(x[,c(1:18)])}
})

# eQTL lead
library(LDlinkR)
tb.list <- lapply(PPIL6_eQTL$rsnum, function(x){
  tryCatch({
    tb <- LDpair("rs736830", x, pop = "EAS", token = mytoken, genome_build = "grch38")
  }, error=function(e){})
})

df_tmp <- as.data.frame(matrix(as.character(NA), nrow = 1, ncol = 18))

colnames(df_tmp) <- colnames(tb.list[[1]])
tb.list.tmp <- lapply(tb.list, function(x){
  if(is.null(x)){return(df_tmp)}else{return(x[,c(1:18)])}
})
library(data.table)
PPIL6_isoQTL_LD <- rbindlist(tb.list, fill = TRUE)
PPIL6_isoQTL_LD_R2 <- PPIL6_isoQTL_LD[,c("var2", "r2")]
colnames(PPIL6_isoQTL_LD_R2)[1] <- "rsnum"
PPIL6_isoQTL_LD_Dprime <- PPIL6_isoQTL_LD[,c("var2", "d_prime")]
colnames(PPIL6_isoQTL_LD_Dprime)[1] <- "rsnum"
PPIL6_isoQTL <- left_join(PPIL6_isoQTL, PPIL6_isoQTL_LD_R2, by = "rsnum")
PPIL6_isoQTL <- left_join(PPIL6_isoQTL, PPIL6_isoQTL_LD_Dprime, by = "rsnum")

PPIL6_isoQTL_LD <- rbindlist(tb.list, fill = TRUE)
PPIL6_eQTL_LD_R2 <- PPIL6_isoQTL_LD[,c("var2", "r2")]
colnames(PPIL6_eQTL_LD_R2)[1] <- "rsnum"
PPIL6_eQTL <- left_join(PPIL6_eQTL, PPIL6_eQTL_LD_R2, by = "rsnum")

# CCVs in high LD with lead eQTL/isoQTL
PPIL6_isoQTL_lead <- subset(PPIL6_isoQTL, r2 > 0.8 & pvalue < 0.0001)
PPIL6_eQTL_lead <- subset(PPIL6_eQTL, r2 > 0.8 & pval_nominal < 0.0001)
PPIL6_eQTL_lead_CADD <- subset(PPIL6_CADD, ID %in% PPIL6_eQTL_lead$rsnum)
PPIL6_isoQTL_lead_CADD <- subset(PPIL6_CADD, ID %in% PPIL6_isoQTL_lead$rsnum)
PHRED_eQTL <- PPIL6_eQTL_lead_CADD$PHRED[!duplicated(PPIL6_eQTL_lead_CADD$ID)]
PHRED_isoQTL <- PPIL6_isoQTL_lead_CADD$PHRED[!duplicated(PPIL6_isoQTL_lead_CADD$ID)]
mean(PHRED_isoQTL)
mean(PHRED_eQTL)
wilcox.test(PHRED_isoQTL, PHRED_eQTL )
library(readxl)
CADD_category <- read_xlsx("CADD_annotation.xlsx")
idx <- grep("splicing",CADD_category$`Functional type`)
PPIL6_eQTL_lead_CADD_splicing_mtx <- PPIL6_eQTL_lead_CADD[!duplicated(PPIL6_eQTL_lead_CADD$ID),idx[-c(1:2)]]
PPIL6_isoQTL_lead_CADD_splicing_mtx <- PPIL6_isoQTL_lead_CADD[!duplicated(PPIL6_isoQTL_lead_CADD$ID),idx[-c(1:2)]]
PPIL6_eQTL_lead_CADD_splicing_mtx_nonzero <- (PPIL6_eQTL_lead_CADD_splicing_mtx > 0 | PPIL6_eQTL_lead_CADD_splicing_mtx < 0)
PPIL6_isoQTL_lead_CADD_splicing_mtx_nonzero <- (PPIL6_isoQTL_lead_CADD_splicing_mtx > 0 | PPIL6_isoQTL_lead_CADD_splicing_mtx < 0)
rowSums(PPIL6_eQTL_lead_CADD_splicing_mtx_nonzero, na.rm = TRUE)
sum(rowSums(PPIL6_eQTL_lead_CADD_splicing_mtx_nonzero, na.rm = TRUE) > 0 )
rowSums(PPIL6_isoQTL_lead_CADD_splicing_mtx_nonzero, na.rm = TRUE)
sum(rowSums(PPIL6_isoQTL_lead_CADD_splicing_mtx_nonzero, na.rm = TRUE) > 0 )

PPIL6_eQTL_lead_CADD$PHRED[!duplicated(PPIL6_eQTL_lead_CADD$ID)]
