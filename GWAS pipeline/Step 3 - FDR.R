require(qvalue)
require(plyr)

rm(list=ls())
ori.wd = getwd()
 
  nuclear.df = read.table("Nuclear SNP.txt", header = TRUE)
  epistasis.df = read.table("Epistasis SNP.txt", header = TRUE)
  
  fdr.nuclear = qvalue(p = nuclear.df$Pvalue)
  pdf("Nuclear.pdf")
  plot(fdr.nuclear)
  dev.off()
  pdf("Nuclear Density Histogram.pdf")
  print(hist(fdr.nuclear))
  dev.off()
  fdr.n = as.data.frame(cbind(fdr.nuclear$pvalues, fdr.nuclear$qvalues, fdr.nuclear$lfdr))
  colnames(fdr.n) = c("Pvalue", "Qvalue","local FDR")
  nuclear = cbind(nuclear.df, fdr.n)
  
  fdr.epistasis = qvalue(p = epistasis.df$Pvalue)
  pdf("Mitonuclear.pdf")
  plot(fdr.epistasis)
  dev.off()
  pdf("Mitonuclear Density Histogram.pdf")
  print(hist(fdr.epistasis))
  dev.off()
  fdr.e = as.data.frame(cbind(fdr.epistasis$pvalues, fdr.epistasis$qvalues, fdr.epistasis$lfdr))
  colnames(fdr.e) = c("Pvalue", "Qvalue","local FDR")
  epistasis = cbind(epistasis.df, fdr.e)
  
  write.table(nuclear, "Q-values Nuclear.txt")
  write.table(epistasis, "Q-values Epistasis.txt")
  
  source('Codes for GWAS/Manhattan Plot.R')
  mp.epistasis = manhattan.plot(epistasis.df,epistasis.df.threshold)
  ggsave(paste0("Epistasis.pdf"), width = 15, height = 5)
  dev.off()
  mp.nuclear = manhattan.plot(nuclear.df,nuclear.df.threshold)
  ggsave(paste0("Nuclear.pdf"), width = 15, height = 5)
  dev.off()
  
  