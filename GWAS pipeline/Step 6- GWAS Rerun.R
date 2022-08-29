rm(list=ls())
ori.wd = getwd()

covariates <-as.matrix(fread('Codes for GWAS/Covariates for mt-n mapping.csv'))
cov <-read.table('Codes for GWAS 2.0/Additional covariates for mt-n mapping.csv', header=T, colClasses=c("factor","factor","factor","factor"))
cov <- sapply(cov,as.numeric)-1
covariates = cbind(cov, covariates)
colnames(covariates)= NULL

write.table(covariates, 'Covariates for mt-n mapping.csv', row.names = FALSE)
####At this point, manually move the newly created file to Codes for GWAS 2.0


rm(list=ls())
ori.wd = getwd()

#load in required packages
require(plyr)
#Load in the function to graph Manhattan Plot, make sure it is in the same project directory
source('Codes for GWAS 2.0/Manhattan Plot.R')

#Read in your phenotype data. Your phenotype data should contain Nuclear, Mitotypes ('X273614N' and 'YPS606')
pheno1 = read.table("Mito1.txt", header = TRUE)
pheno2 = read.table("Mito2.txt", header = TRUE)
pheno3 = read.table("Mito3.txt", header = TRUE)

all = rbind(pheno1,pheno2, pheno3)
#Read in files containing all possible strains in both recombinant collection.
rc1 = read.table("Text files/All possible strains in RC1.txt",header = TRUE)
rc2 = read.table("Text files/All possible strains in RC2.txt",header = TRUE)
rc3 = read.table("Text files/All possible strains in RC3.txt",header = TRUE)
rcs = rbind(rc1,rc2,rc3)

#Assign strains without phenotypes NA
pheno = join(rcs,all, by = c("Nuclear","Mito"))
pheno$Petitefrequency = 100*pheno$Petite/(pheno$Petite+pheno$Grande)

nuclear.df.threshold = 5.7
epistasis.df.threshold = 3

  source("Codes for GWAS 2.0/Load in SNP tables, covariates, SNP maps, and annotations.R")
  source("Codes for GWAS 2.0/GWAS function for binomial generalized linear model.R")

  GWAS(pheno,'Petite','Grande')

  setwd(ori.wd)
  write.table(epistasis.df, "New Epistasis SNP.txt")
  write.table(nuclear.df, "New Nuclear SNP.txt")

  mp.epistasis = manhattan.plot(epistasis.df,epistasis.df.threshold)
  ggsave(paste0("New Epistasis.pdf"), width = 15, height = 5)
  dev.off()
  mp.nuclear = manhattan.plot(nuclear.df,nuclear.df.threshold)
  ggsave(paste0("New Nuclear.pdf"), width = 15, height = 5)
  dev.off()
