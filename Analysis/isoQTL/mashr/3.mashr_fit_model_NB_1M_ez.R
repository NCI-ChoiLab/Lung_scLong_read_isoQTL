

library(rhdf5)
library(mashr)
library(ashr)
library(data.table)
library(tidyverse)
library(argparse)
library(matrixcalc)
source("/data/lib14/R/Rscripts/utilities.R")

setwd("/data/Choi_lung/scLongreads/mashr/")
data.strong <- fun_h5_2_mashr("strong_all_ez.h5")
data.random <- fun_h5_2_mashr("random_all_1M_seed18_ez.h5", max.missing = 0)

message("Total number of strong associations: ", 
        nrow(data.strong$Bhat))
message("Total number of random associations: ", 
        nrow(data.random$Bhat))

# Estimate the correlation structure in the null tests from the random data
message("Estimating covariates...")

Vhat <- estimate_null_correlation_simple(data.random)

random <- mash_set_data(data.random$Bhat,
                        data.random$Shat, 
                        V=Vhat) 

strong <- mash_set_data(data.strong$Bhat, 
                        data.strong$Shat, 
                        V=Vhat) 
if(ncol(data.strong$Bhat) < 5){
  n <- ncol(data.strong$Bhat)
} else{
  n <- 5
}

covar.pca <- cov_pca(strong, n)
covar.ed <- cov_ed(strong, covar.pca) # Now we use the strong tests to set up data-driven covariances.
covar.c <- cov_canonical(random)

message("Fitting mashr...")
mashFit <- mash(random, Ulist = c(covar.c, covar.ed), outputlevel = 1) 

saveRDS(mashFit, "mashr_all_1M_fit_ez.rds")

