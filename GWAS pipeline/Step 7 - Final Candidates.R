rm(list=ls())
ori.wd = getwd()


require(qvalue)
require(plyr)

source("Codes for GWAS/Load in SNP tables, covariates, SNP maps, and annotations.R")
annotation$Distance = ifelse(annotation$CONSEQUENCE =="upstream_gene_variant", as.numeric(gsub("\\D", "", annotation$Extra)),
                             ifelse(annotation$CONSEQUENCE =="downstream_gene_variant", -1 * as.numeric(gsub("\\D", "", annotation$Extra)),0 ))

#Read in files containing all possible strains in both recombinant collection.
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

### Epistasis Candidates
  
epistasis.df = read.table("New Epistasis SNP.txt", header = TRUE)

fdr.epistasis = qvalue(p = epistasis.df$Pvalue)
pdf("Mitonuclear 2.0.pdf")
plot(fdr.epistasis)
dev.off()
pdf("Mitonuclear Density Histogram 2.0.pdf")
print(hist(fdr.epistasis))
dev.off()
fdr.e = as.data.frame(cbind(fdr.epistasis$pvalues, fdr.epistasis$qvalues, fdr.epistasis$lfdr))
colnames(fdr.e) = c("Pvalue", "Qvalue","local FDR")
epistasis = cbind(epistasis.df, fdr.e)

write.table(epistasis, "Q-values Epistasis 2.0.txt")

epistasis.df.sig = subset(epistasis, Qvalue < 0.05)
final.epistasis.df = join(epistasis.df.sig, annotation)
epistasis.candidates = aggregate(neg_logP ~ GENE + CHROM + POS + CONSEQUENCE + Distance, FUN = mean, data = final.epistasis.df)
epistasis.candidates = subset(epistasis.candidates, !(CONSEQUENCE %in% c( "downstream_gene_variant","synonymous_variant")))
epistasis.candidates = subset(epistasis.candidates, epistasis.candidates$Distance < 500)

write.table(epistasis.candidates, "Mitonuclear Candidates.txt", row.names = FALSE)

epistasis.candidates$ID = paste0(epistasis.candidates$CHROM,":",epistasis.candidates$POS)
epistasis.candidates$Gene.ID = paste0(epistasis.candidates$GENE,"-",epistasis.candidates$CHROM,"-",epistasis.candidates$POS)
epistasis.candidate = epistasis.candidates[,(ncol(epistasis.candidates)-1):ncol(epistasis.candidates)]
epistasis.candidate = unique(epistasis.candidate[,(1:2)])
reduced.map = subset(map,ID %in% unique(epistasis.candidate$ID))
reduced.map$location = row.names(reduced.map)
reduced.map = join(reduced.map,epistasis.candidate)
#SNP table should have been imported at this stage
reduced.snp = as.data.frame(geno.working[,(as.integer(reduced.map$location))])
colnames(reduced.snp) = as.character(reduced.map$Gene.ID)
#pheno$Petitefrequency = pheno$Petitefrequency*100
final <- cbind(pheno,reduced.snp)
hits = reduced.map$Gene.ID
consolidate = data.frame()
for (hit in hits){
  working.df = cbind(hit,final$Mito,final$Petitefrequency, final[hit])
  working.df = na.omit(working.df)
  colnames(working.df) = c("GeneID",'Mitotype','Pheno','Hit')
  tempdf = aggregate(data = working.df, working.df$Pheno ~ ., FUN = mean)
  
  change = data.frame()
  for (snp in unique(tempdf$Hit)){
    SNP.table = subset(tempdf, Hit == snp)
    mito1.2 = subset(SNP.table, Mitotype %in% c("X273614N","YPS606"))
    mito1.2 = data.frame(Gene = mito1.2$GeneID[1], Allele =mito1.2$Hit[1], Mitotype.1 = "X273614N", Mitotype.2 = "YPS606" ,Change = mito1.2$`working.df$Pheno`[1]-mito1.2$`working.df$Pheno`[2])
    tempDF = mito1.2
    change = rbind(change,tempDF)
  }
  
  delta = data.frame()
  for (mito in unique(change$Mitotype.1)){
    SNP.table = subset(change, Mitotype.1 == mito)
    SNP.table$Delta.Change = abs(SNP.table$Change[1]-SNP.table$Change[2])
    delta = rbind(delta,SNP.table)
  }
  
  
  delta$GENE = (unlist(strsplit(hit,"-")))[1]
  delta$CHROM = (unlist(strsplit(hit,"-")))[2]
  delta$POS = (unlist(strsplit(hit,"-")))[3]
  consolidate = rbind(consolidate,delta)
}
consolidate.epistasis = join(consolidate,epistasis.candidates[c(1:6)])

#### Nuclear Candidates


nuclear.df = read.table("New Nuclear SNP.txt", header = TRUE)

fdr.nuclear = qvalue(p = nuclear.df$Pvalue)
pdf("Nuclear 2.0.pdf")
plot(fdr.nuclear)
dev.off()
pdf("Nuclear Density Histogram 2.0.pdf")
print(hist(fdr.nuclear))
dev.off()
fdr.n = as.data.frame(cbind(fdr.nuclear$pvalues, fdr.nuclear$qvalues, fdr.nuclear$lfdr))
colnames(fdr.n) = c("Pvalue", "Qvalue","local FDR")
nuclear = cbind(nuclear.df, fdr.n)

write.table(nuclear, "Q-values Nuclear 2.0.txt")


nuclear.df.sig = subset(nuclear, Qvalue < 0.001)
final.nuclear.df = join(nuclear.df.sig, annotation)
nuclear.candidates = aggregate(neg_logP ~ GENE + CHROM + POS + CONSEQUENCE + Distance, FUN = mean, data = final.nuclear.df)
nuclear.candidates = subset(nuclear.candidates, !(CONSEQUENCE %in% c( "downstream_gene_variant","synonymous_variant")))
nuclear.candidates = subset(nuclear.candidates, nuclear.candidates$Distance < 500)

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



