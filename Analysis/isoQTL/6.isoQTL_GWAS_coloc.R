library(dplyr)
library(ggplot2)
library(stringr)
library(stringi)

# read GWAS stats downloaded from GWAS_catalog
# https://www.ebi.ac.uk/gwas/publications/36914875
lung_fun_gwas <- read.table(gzfile("/data/Choi_lung/scLongreads/colocalization/lung_function/GCST90244094_buildGRCh37.tsv.gz"), header = TRUE,
                            sep = "\t")
lung_fun_assoc <- read.table("/data/Choi_lung/scLongreads/colocalization/lung_function/Lung_func_gwas_association.txt",
                             sep = "\t", header = TRUE)

lung_fun_gwas$snp <- paste(lung_fun_gwas$variant_id, lung_fun_gwas$other_allele, lung_fun_gwas$effect_allele, sep = "-")

beds <- lung_fun_gwas[,c(2,3)]
beds$end <- beds$base_pair_location + 1 
beds$chromosome <- paste0("chr", beds$chromosome)
beds$snp <- paste(lung_fun_gwas$variant_id, lung_fun_gwas$other_allele, lung_fun_gwas$effect_allele, sep = "-")
options(scipen = 999)
write.table(beds, file = "/data/Choi_lung/scLongreads/colocalization/lung_function/lung_fun_hg19.bed", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
length(unique(beds$snp))

# read SNP info after liftover
hg38_snps <- read.table("/data/Choi_lung/scLongreads/colocalization/lung_function/lung_fun_hg38.bed",
                        header = FALSE,sep = "\t")
length(unique(hg38_snps$V4))
colnames(hg38_snps) <- c("chr", "pos", "end", "snp")

# clean format to match isoQTL
hg38_snps$effect_allele <- str_split_fixed(hg38_snps$snp, "-", 3)[,3]
hg38_snps$rsid <- str_split_fixed(hg38_snps$snp, "-", 3)[,1]
hg38_snps$STRONGEST.SNP.RISK.ALLELE <- paste(hg38_snps$rsid, 
                                             hg38_snps$effect_allele, sep = "-")
lung_fun_gwas <- left_join(lung_fun_gwas, hg38_snps[,c(1,2,4)], by = "snp")
lung_fun_gwas$STRONGEST.SNP.RISK.ALLELE <- paste(lung_fun_gwas$variant_id, 
                                                 lung_fun_gwas$effect_allele, sep = "-")
lung_fun_assoc <- left_join(lung_fun_assoc, hg38_snps, by = "STRONGEST.SNP.RISK.ALLELE")
lung_fun_assoc <- left_join(lung_fun_assoc, lung_fun_gwas[,c(10:13)], by = "STRONGEST.SNP.RISK.ALLELE")

lung_fun_assoc <- lung_fun_assoc[!duplicated(lung_fun_assoc$STRONGEST.SNP.RISK.ALLELE),]
window_size <- 1000000
lead_snps <- lung_fun_assoc$STRONGEST.SNP.RISK.ALLELE
upstream <- as.numeric(lung_fun_assoc$pos) - window_size/2
downstream <- as.numeric(lung_fun_assoc$pos) + window_size/2
chr <-lung_fun_assoc$chr
gwas_loci_total <- lapply(1:length(lead_snps), function(x){
  rst <- lung_fun_gwas[which((lung_fun_gwas$chr == chr[x]) & (lung_fun_gwas$pos > upstream[x]) & (lung_fun_gwas$pos < downstream[x]) ),]
  return(rst)
})


names(gwas_loci_total) <- lead_snps


saveRDS(gwas_loci_total, "/data/Choi_lung/scLongreads/colocalization/Shrine_et_al_lung_func_406_lead_1MB.rds")
gwas_loci_total <- readRDS("/data/Choi_lung/scLongreads/colocalization/Shrine_et_al_lung_func_406_lead_1MB.rds")
gwas_loci_total <- lapply(gwas_loci_total, function(x){
  idx <- which(paste0("chr", x$chromosome) == x$chr)
  x <- x[idx,]
  return(x)
})
names(gwas_loci_total)
names(gwas_loci_total) -> lead_snps
gwas_loci_total$`rs2058914-A`$STRONGEST.SNP.RISK.ALLELE
i = 0
for (snp in lead_snps) {
  if (snp %in% gwas_loci_total[[snp]]$STRONGEST.SNP.RISK.ALLELE) {
    i= i + 1
    print(i)
  }
}

gwas_loci <- do.call(rbind, gwas_loci_total)

snp_list <- unique(gwas_loci$variant_id)

lSNP_win1MB.list <- lapply(NB.nominal.sig.list, function(x){
  x <- subset(x, snp %in% snp_list)
})

library(coloc)
library(snowfall)
library(foreach)
library(doParallel)
library(parallel)

iso_tbt_list <- lapply(lSNP_win1MB.list, function(x) unique(x$phenotype_id))

# extract isoQTL stats of SNPs w/ GWAS stats
NB.coloc.list <- list()
for (snp in lead_snps) {
  locus <- gwas_loci_total[[snp]]
  if (nrow(locus) > 0) {
    chr <- unique(locus$chromosome)
    snp_list <- unique(locus$variant_id)
    lSNP_win1MB.list <- lapply(NB.nominal.sig.list, function(x){
      x <- subset(x, snp %in% snp_list)
    })
    df <- data.frame(Number_tobetested = sapply(lSNP_win1MB.list, function(x) nrow(x)),
                     Celltype = names(lSNP_win1MB.list))
    celltypes <- df$Celltype[which(df$Number_tobetested != 0)]
    
    ct.list <- list()
    for(ct in celltypes){
      file_path <- paste0("/data/Choi_lung/scLongreads/jaxqtl/",ct, "/chr", chr,"_chunck_jaxqtl_nb.cis_qtl_pairs.", chr, ".wald.parquet")
      rst <- read_parquet(file_path)
      rst_sub <- subset(rst, snp %in% snp_list)
      ct.list[[ct]] <- rst_sub
    }
    NB.coloc.list[[snp]] <- ct.list
  }else{
    NB.coloc.list[[snp]] <- NULL
  }
  message("Finish SNP: ", snp)
}

saveRDS(NB.coloc.list, "/data/Choi_lung/scLongreads/colocalization/Shrine_isoQTL_NB_nominal_wald_lungfun_gwas_overlap_list.rds")
library(coloc)
library(snowfall)
library(foreach)
library(doParallel)
library(parallel)

df <- data.frame(Number_tobetested = sapply(NB.coloc.list, function(x) length(x)),
                 snp = names(NB.coloc.list))
snps_tobetested <- df$snp[which(df$Number_tobetested != 0)]
NB.coloc.list <- NB.coloc.list[snps_tobetested]
rst.list <- list()
output.list <- list()
for(i in 1:length(NB.coloc.list)){
  # Prepare dataset of isoQTL and GWAS summary stats for coloc
  # run each cell type 
  if(!is.null(NB.coloc.list[[i]])){
    filt_eqtl.list <- NB.coloc.list[[i]]
    gwas_loci <- gwas_loci_total[[names(NB.coloc.list)[i]]]
    rst.ct.list <- list()
    output.ct.list <- list()
    for(j in 1:length(filt_eqtl.list)){
      #######
      # for each cell type within specific locus
      #######
      filt_eqtl <- filt_eqtl.list[[j]]
      ct <- names(filt_eqtl.list)[j]
      snp_list <- intersect(gwas_loci$variant_id, filt_eqtl$snp)
      if(length(snp_list) == 0){
        rst.list[[i]] <- NULL
        output.list[[i]] <- NULL
        message("no overlaps in ", names(NB.coloc.list)[i])
      }else{

        filt_eqtl <- filt_eqtl[!is.na(filt_eqtl$pval_nominal),] # remove non-converging results
        filt_eqtl$rsid = paste(filt_eqtl$snp,filt_eqtl$a0, filt_eqtl$a1, sep = '-')
        filt_eqtl = filt_eqtl[,c('phenotype_id','rsid','chrom','af','pval_nominal','pos',"slope","slope_se")]
        
        # subset gwas data
        filt_gwas <- gwas_loci[gwas_loci$variant_id %in% snp_list,]
        
        filt_gwas$SNP_ID = paste(filt_gwas$variant_id,filt_gwas$effect_allele, filt_gwas$other_allele,sep = '-')
        filt_gwas = filt_gwas[,c('chromosome','SNP_ID','pos','Zscore', 'p_value', "N", "effect_allele_frequency")]
        filt_gwas <- filt_gwas[!duplicated(filt_gwas$SNP_ID),]
        colnames(filt_gwas)[2] = 'rsid'
        maf_eqtl <- merge(filt_eqtl, filt_gwas, by="rsid", all=FALSE)
        # snp_alt_list <- intersect(filt_gwas$rsid , unique(filt_eqtl$rsid))
        
        ##filtering snps to only include snps in the eQTL list and vice versa
        filt_gwas=filt_gwas[filt_gwas$rsid %in% maf_eqtl$rsid,]
        filt_eqtl=maf_eqtl[maf_eqtl$rsid %in% filt_gwas$rsid,]
        length(filt_gwas$rsid)
        length(unique(filt_eqtl$rsid))
        gene_list=as.character(unique(filt_eqtl$phenotype_id))
        eqtl_genes=list()
        eqtl_genes=foreach(i=gene_list) %dopar%
          droplevels(filt_eqtl[filt_eqtl$phenotype_id==i,]) 
        names(eqtl_genes)=gene_list
        
        ##note in the drop levels, it needs to be col 6 since that is where the rsIDs are for the eqtl data, otherwise you get an empty snp gene list.
        gwas_genes=list()
        gwas_genes=foreach(i=gene_list) %dopar%
          droplevels(filt_gwas[na.omit(match(as.character(eqtl_genes[[i]][,1]),as.character(filt_gwas$rsid))),]) 
        names(gwas_genes)=gene_list
        
        
        # Add the sample size of isoQTL
        
        ct_path <- paste0("/data/Choi_lung/scLongreads/tensorqtl/isoform_level/", ct, "/samples.txt")
        samples <- read.table(ct_path, header = FALSE)
        sample_size <- nrow(samples)
        
        SNP_coloc=list()
        SNP_coloc=foreach(i=(1:length(eqtl_genes)),.packages='coloc') %dopar%
          coloc.abf(dataset1 = list(snp=eqtl_genes[[i]]$rsid,pvalues=eqtl_genes[[i]]$pval_nominal,beta=eqtl_genes[[i]]$slope,varbeta=(eqtl_genes[[i]]$slope_se)^2,position=eqtl_genes[[i]]$pos.x,N=sample_size,type="quant",MAF=eqtl_genes[[i]]$af),
                    dataset2 = list(snp=gwas_genes[[i]]$rsid,pvalues=gwas_genes[[i]]$p_value,position=gwas_genes[[i]]$pos,type="quant", MAF=gwas_genes[[i]]$effect_allele_frequency, N=gwas_genes[[i]]$N))
        names(SNP_coloc)=gene_list
        rst.ct.list[[ct]] <- SNP_coloc
        output_genes=list()
        output_genes=foreach(i=gene_list) %dopar%
          capture.output(SNP_coloc[[i]][["summary"]])
        names(output_genes)=gene_list
        output_sum <- lapply(1:length(output_genes), function(x){
          rst <- as.numeric(str_split_fixed(output_genes[[x]][2], " ", 6))
        })
        names(output_sum)=gene_list
        if(length(output_sum) == 1){
          tmp <- Reduce(rbind, output_sum)
          output <- data.frame(nsnps = tmp[1], PP.H0.abf = tmp[2], PP.H1.abf = tmp[3], PP.H2.abf = tmp[4], PP.H3.abf = tmp[5], PP.H4.abf = tmp[6])
          rownames(output) <- gene_list
          output$transcript_id <- rownames(output)
        }else{
          output <- Reduce(rbind, output_sum)
          colnames(output) <- c("nsnps", "PP.H0.abf", "PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")
          rownames(output) <- gene_list
          output <- as.data.frame(output)
          output$transcript_id <- rownames(output)
        }
        output$Celltype <- ct
        output <- output %>% group_by(transcript_id) %>% mutate(transcript_name = Search_transcript_name2(transcript_id))
        #######
        
        
        
        ######
        ## get output list for each cell type within locus
        output.ct.list[[ct]] <- output
      }
    }
    #####
    rst.list[[i]] <- rst.ct.list
    output.list[[i]] <- output.ct.list
    
    
  }else{
    rst.list[[i]] <- NULL
    output.list[[i]] <- NULL
  }
}
names(rst.list)  <- names(NB.coloc.list)
names(output.list)  <- names(NB.coloc.list)

output.list1 <- lapply(output.list, function(x) do.call(rbind, x))
output.list.final <- lapply(names(output.list), function(x) {
  rst <- output.list1[[x]]
  rst$leadSNP <- x
  return(rst)
})
output.final <- do.call(rbind, output.list.final)
output.final.sig <- output.final[which(output.final$PP.H4.abf > 0.7),]
saveRDS(rst.list, file = "/data/Choi_lung/scLongreads/colocalization/lung_function/Colocalization_rst_list_lung_function.rds")
saveRDS(output.list, file = "/data/Choi_lung/scLongreads/colocalization/lung_function/Colocalization_output_list_lung_function.rds")
saveRDS(output.final.sig, file = "/data/Choi_lung/scLongreads/colocalization/lung_function/Colocalization_output_table.rds")

