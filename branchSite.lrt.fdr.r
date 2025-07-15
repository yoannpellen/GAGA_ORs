# Written for the analysis of multiple testing of many branches of a tree (or set of trees) using the branch-site model
# Based on the practice in the SELECTOME pipeline (group of Marc Robinson-Rechavi)
# Uses the qvalue algorithm: Storey and Tibshirani PNAS 2003 "Statistical significance for genome-wide experiments"
# implemented in the qvalue package: http://genomics.princeton.edu/storeylab/qvalue/

# Input: csv file with the following fields: lnull,lalt  (could have additional fields such as ID, which will be ignored but printed to output)
# lnull is the lnL under the null model and lalt is the likelihood under the alternative model
# FDR correction set at 0.05
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("qvalue")
library(qvalue)

args<-commandArgs();
if(length(args)!=5) {
  cat("USAGE: R --slave --args lnlFile outFile < branchSite.lrt.fdr.r\n");
  q();
}
#cat("args = "); cat(args); cat("\n");
lnlFile <- args[4];
outFile <- args[5];

lnl = read.table(lnlFile,header=TRUE,sep=",")
lnl$deltal2 = 2*(lnl$lalt - lnl$lnull) * (lnl$lalt - lnl$lnull > 0)
lnl$pval = 1-pchisq(lnl$deltal2,df=1)

qval = qvalue(lnl$pval, pi0.meth="bootstrap", fdr.level=0.05, robust=TRUE)
lnl$qval = qval$qvalues
lnlSorted <- lnl[order(lnl$pval),]

write.table(lnlSorted,outFile,sep=",",row.names=FALSE,quote=FALSE)
