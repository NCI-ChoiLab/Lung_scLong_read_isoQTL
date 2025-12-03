# Type I error control
# summarize pval beta on permuted data

# Bolun

library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")

qqplotdata_permutation <- function(logpvector, logpermutatedp){
  o = sort(logpvector,decreasing=T)
  e = sort(logpermutatedp,decreasing=T)
  qqdata <- data.frame(o,e)
  qqdata$o <- round(qqdata$o,3)
  qqdata$e <- round(qqdata$e,3)
  keepU <- which(!duplicated(qqdata))
  qqdata <- qqdata[keepU,]
  
  N <- length(logpvector) ## number of p-values
  ## create the confidence intervals
  qqdata$c975 <- NA
  qqdata$c025 <- NA
  
  ## the jth order statistic from a
  ## uniform(0,1) sample
  ## has a beta(j,n-j+1) distribution
  ## (Casella & Berger, 2002,
  ## 2nd edition, pg 230, Duxbury)
  
  for(i in 1:length(keepU)){
    j <- keepU[i]
    qqdata$c975[i] <- -log10(qbeta(0.975,j,N-j+1))
    qqdata$c025[i] <- -log10(qbeta(0.025,j,N-j+1))
  }
  return(qqdata)
}

qqplotdata <- function(logpvector){
  o = sort(logpvector,decreasing=T)
  e = -log10(ppoints(length(o)))       
  qqdata <- data.frame(o,e)
  qqdata$o <- round(qqdata$o,3)
  qqdata$e <- round(qqdata$e,3)
  keepU <- which(!duplicated(qqdata))
  qqdata <- qqdata[keepU,]
  
  N <- length(logpvector) ## number of p-values
  ## create the confidence intervals
  qqdata$c975 <- NA
  qqdata$c025 <- NA
  
  ## the jth order statistic from a
  ## uniform(0,1) sample
  ## has a beta(j,n-j+1) distribution
  ## (Casella & Berger, 2002,
  ## 2nd edition, pg 230, Duxbury)
  
  for(i in 1:length(keepU)){
    j <- keepU[i]
    qqdata$c975[i] <- -log10(qbeta(0.975,j,N-j+1))
    qqdata$c025[i] <- -log10(qbeta(0.025,j,N-j+1))
  }
  return(qqdata)
}

input_path <- "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/output"
output_path <- "/data/lib14/project/scLongread/"

for (celltype in celltypes) {
  celltype_for_path <- gsub(" ", "_", celltype)
  message("Celltype: ", celltype_for_path)
  qtl_path <- file.path(input_path,celltype_for_path)
  all_qtl <- read.table(file = paste0(qtl_path,"/qtl_results_all.txt"), 
                        sep = "\t", header = TRUE)
  rownames(all_qtl) <- paste(all_qtl$phenotype_id, all_qtl$variant_id, sep = "-")
  permutated_qtl <- read.table(file = paste0(qtl_path,"/permutation_qtl_results_all.txt"), 
                        sep = "\t", header = TRUE)
  rownames(permutated_qtl) <- paste(permutated_qtl$phenotype_id, permutated_qtl$variant_id, sep = "-")
  rowname <- intersect(rownames(all_qtl), rownames(permutated_qtl))
  permutated_qtl <- permutated_qtl[rowname,]
  set.seed(12345)
  idx <- sample(1:length(all_qtl$pval_nominal), 10000000)
  pvals <- all_qtl$pval_nominal[idx]
  pvals_permutated <- permutated_qtl$pval_nominal[idx]
  legendcol <- character(0)
  ycol <- "log10P"
  plotdata <- qqplotdata_permutation(-log10(pvals), -log10(pvals_permutated))
  fx <- plotdata$e
  fy <- plotdata$o
  fcol <- rep("#E41A1C",length(plotdata$o))
  opt =  list(break.top = 10,
              top.size = 0.125)
  png(filename = paste0(output_path,"QQplot_permutation_vs_oberved_", celltype_for_path, ".png"), width = 8, height = 8, units = "in",res=300)
  xlim <- c(0,max(fx,na.rm=T))
  ylim <- c(0,max(fy,na.rm=T))
  maxY <- max(fy,na.rm=T)
  print("okkkk2")
  par(mar=c(5.1,5.1,4.1,1.1))
  print("okkkk3")
  
  lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
  lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
  lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
  lab2 <- lab2[lab2 > max(lab1)]
  
  # resulting range of top scale in bottom scale units
  top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
  top.data = max(lab2)-opt$break.top
  
  # function to rescale the top part
  rescale = function(y) { opt$break.top+(y-opt$break.top)/(top.data/top.range)}
  rescaled.y = rescale(fy[fy>opt$break.top])
  plot(0,0,
       ylim=c(min(fy),opt$break.top*(1+opt$top.size)),xlim=xlim,axes=FALSE,
       xlab=expression(plain(Permutation)~~group("(",-log[10]*italic(P),")")),
       ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
       cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
       main=opt$maintitle,pch=19)
  
  # add points
  points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)
  
  # identify line & add axis break
  lines(xlim,xlim,col="black",lty = 2)
  axis(1,cex.axis=1.5,cex.lab=1.5)
  par(las=1)
  axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
  par(las=0)
  box()
  par(las=0)
  points(fx[fy>opt$break.top],rescaled.y,cex=1,col=fcol[fy>opt$break.top],pch=19)
  par(las=1)
  axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
  axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
  axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
  lines(range(fx),c(opt$break.top,opt$break.top),col = "grey",lty = 6)
  title(paste0("The expected p-values under the permutation (x-axis) \n against the observed p-values (y-axis) in ",celltype))
  dev.off()
  # null vs observed p
  plotdata <- qqplotdata(-log10(pvals))
  fx <- plotdata$e
  fy <- plotdata$o
  fcol <- rep("#E41A1C",length(plotdata$o))
  opt =  list(break.top = 10,
              top.size = 0.125)
  png(filename = paste0(output_path,"QQplot_null_vs_oberved_", celltype_for_path, ".png"), width = 8, height = 8, units = "in",res=300)
  xlim <- c(0,max(fx,na.rm=T))
  ylim <- c(0,max(fy,na.rm=T))
  maxY <- max(fy,na.rm=T)
  print("okkkk2")
  par(mar=c(5.1,5.1,4.1,1.1))
  print("okkkk3")
  
  lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
  lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
  lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
  lab2 <- lab2[lab2 > max(lab1)]
  
  # resulting range of top scale in bottom scale units
  top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
  top.data = max(lab2)-opt$break.top
  
  # function to rescale the top part
  rescale = function(y) { opt$break.top+(y-opt$break.top)/(top.data/top.range)}
  rescaled.y = rescale(fy[fy>opt$break.top])
  plot(0,0,
       ylim=c(min(fy),opt$break.top*(1+opt$top.size)),xlim=xlim,axes=FALSE,
       xlab=expression(plain(Null)~~group("(",-log[10]*italic(P),")")),
       ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
       cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
       main=opt$maintitle,pch=19)
  
  # add points
  points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)
  
  # identify line & add axis break
  lines(xlim,xlim,col="black",lty = 2)
  axis(1,cex.axis=1.5,cex.lab=1.5)
  par(las=1)
  axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
  par(las=0)
  box()
  par(las=0)
  points(fx[fy>opt$break.top],rescaled.y,cex=1,col=fcol[fy>opt$break.top],pch=19)
  par(las=1)
  axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
  axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
  axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
  lines(range(fx),c(opt$break.top,opt$break.top),col = "grey",lty = 6)
  title(paste0("The expected p-values under the null hypothesis (x-axis) \n against the observed p-values (y-axis) in ",celltype))
  dev.off()
  
  # null vs permutation
  plotdata <- qqplotdata(-log10(pvals_permutated))
  fx <- plotdata$e
  fy <- plotdata$o
  fcol <- rep("#E41A1C",length(plotdata$o))
  opt =  list(break.top = 10,
              top.size = 0.125)
  png(filename = paste0(output_path,"QQplot_null_vs_permutation_", celltype_for_path, ".png"), width = 8, height = 8, units = "in",res=300)
  xlim <- c(0,max(fx,na.rm=T))
  ylim <- c(0,max(fy,na.rm=T))
  maxY <- max(fy,na.rm=T)
  print("okkkk2")
  par(mar=c(5.1,5.1,4.1,1.1))
  print("okkkk3")
  
  lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
  lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
  lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
  lab2 <- lab2[lab2 > max(lab1)]
  
  # resulting range of top scale in bottom scale units
  top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
  top.data = max(lab2)-opt$break.top
  
  # function to rescale the top part
  rescale = function(y) { opt$break.top+(y-opt$break.top)/(top.data/top.range)}
  rescaled.y = rescale(fy[fy>opt$break.top])
  plot(0,0,
       ylim=c(min(fy),opt$break.top*(1+opt$top.size)),xlim=xlim,axes=FALSE,
       xlab=expression(plain(Null)~~group("(",-log[10]*italic(P),")")),
       ylab=expression(plain(Permutated)~~group("(",-log[10]*italic(P),")")),
       cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
       main=opt$maintitle,pch=19)
  
  # add points
  points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)
  
  # identify line & add axis break
  lines(xlim,xlim,col="black",lty = 2)
  axis(1,cex.axis=1.5,cex.lab=1.5)
  par(las=1)
  axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
  par(las=0)
  box()
  par(las=0)
  points(fx[fy>opt$break.top],rescaled.y,cex=1,col=fcol[fy>opt$break.top],pch=19)
  par(las=1)
  axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
  axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
  axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
  lines(range(fx),c(opt$break.top,opt$break.top),col = "grey",lty = 6)
  title(paste0("The expected p-values under the null hypothesis (x-axis) \n against the permutated p-values (y-axis) in ",celltype))
  dev.off()
}

