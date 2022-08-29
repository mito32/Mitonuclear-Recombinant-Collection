require(plyr)
require(stringi)
require(ggplot2)
require(tidyr)

rm(list=ls())
ori.wd = getwd()

source("Codes for GWAS/Load in SNP tables, covariates, SNP maps, and annotations.R")

annotation$Distance = ifelse(annotation$CONSEQUENCE =="upstream_gene_variant", as.numeric(gsub("\\D", "", annotation$Extra)),
                             ifelse(annotation$CONSEQUENCE =="downstream_gene_variant", -1 * as.numeric(gsub("\\D", "", annotation$Extra)),0 ))
#Read in files containing all possible strains in both recombinant collection.
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

  
  nuclear.df = read.table("Q-values Nuclear.txt", header = TRUE)
  
  nuclear.df.sig = subset(nuclear.df, Qvalue < 0.001)
  final.nuclear.df = join(nuclear.df.sig, annotation)
  nuclear.candidates = aggregate(neg_logP ~ GENE + CHROM + POS + CONSEQUENCE + Distance, FUN = mean, data = final.nuclear.df)
  nuclear.candidates = subset(nuclear.candidates, !(CONSEQUENCE %in% c( "downstream_gene_variant","synonymous_variant")))
  nuclear.candidates = subset(nuclear.candidates, nuclear.candidates$Distance < 500)
 
  write.table(nuclear.candidates, "Nuclear Candidates.txt", row.names = FALSE)
  
  
  nuclear.candidates$ID = paste0(nuclear.candidates$CHROM,":",nuclear.candidates$POS)
  nuclear.candidates$Gene.ID = paste0(nuclear.candidates$GENE,":",nuclear.candidates$CHROM,":",nuclear.candidates$POS)
  nuclear.candidate = nuclear.candidates[,(ncol(nuclear.candidates)-1):ncol(nuclear.candidates)]
  nuclear.candidate = unique(nuclear.candidate[,(1:2)])
  reduced.map = subset(map,ID %in% unique(nuclear.candidate$ID))
  reduced.map$location = row.names(reduced.map)
  reduced.map = join(reduced.map,nuclear.candidate)
  #SNP table should have been imported at this stage
  reduced.snp = as.data.frame(geno.working[,(as.integer(reduced.map$location))])
  colnames(reduced.snp) = as.character(reduced.map$Gene.ID)
  #pheno$Petitefrequency = pheno$Petitefrequency*100
  final <- cbind(pheno,reduced.snp)
  
  
  hits = reduced.map$Gene.ID
  
  consolidate = data.frame()
  for (hit in hits){
    working.df = cbind(hit,final$Petitefrequency, final[hit])
    working.df = na.omit(working.df)
    colnames(working.df) = c("SNP",'Pheno','Variant')
    tempdf = aggregate(data = working.df, working.df$Pheno ~ ., FUN = mean)
    a = as.data.frame(count(working.df[c(1:181),], vars = "Variant"))
    a$MAF = ifelse(a$freq == min(a$freq), "Minor", "Major")
    a$Allele.frequency = as.numeric(round(a$freq/sum(a$freq),2))
    a$freq = NULL
    
    
    tempdf = join(tempdf,a)
    colnames(tempdf) = c("SNP",'Variant','Average.Pheno', 'MAF','Allele.Frequency')
    
    minor = subset(tempdf, MAF =="Minor")
    major = subset(tempdf, MAF =="Major")
    
    maf = minor$Allele.Frequency[1]
    
    tempdf$GENE = (unlist(strsplit(hit,":")))[1]
    tempdf$CHROM = (unlist(strsplit(hit,":")))[2]
    tempdf$POS = (unlist(strsplit(hit,":")))[3]
    
    delta = tempdf[, c("GENE", "CHROM", "POS", "SNP", "Variant","MAF","Allele.Frequency","Average.Pheno")]
    delta$Effect.Size = abs((major$Average.Pheno -minor$Average.Pheno)/2)*maf
    
    consolidate = rbind(consolidate,delta)
  }
  
  
  consolidate.nuclear = join(consolidate,nuclear.candidates[c(1:6)])
  nuclear = subset(nuclear.df, ID %in% nuclear.candidate$ID, c(1,2,4,5,7))
  consolidate.nuclear = join(consolidate.nuclear,nuclear)
  
