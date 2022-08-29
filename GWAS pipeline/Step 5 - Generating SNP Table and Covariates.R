require(plyr)
require(stringi)
require(data.table)

rm(list=ls())
ori.wd = getwd()


old.location = paste0(ori.wd,"/Codes for GWAS")
oldfiles = list.files(old.location)

geno.working<-as.matrix(fread('Codes for GWAS/Final SNP Tables for mt-n mapping.csv', header = FALSE))
covariates <-as.matrix(fread('Codes for GWAS/Covariates for mt-n mapping.csv'))
map = read.table('Codes for GWAS/Final SNP Map.txt', header = TRUE)

geno = cbind(map, t(geno.working))

nuclear.df = read.table("Nuclear Candidates.txt", header = TRUE)
nuclear = subset(nuclear.df, GENE == "MIP1")
nuclear = nuclear[!duplicated(nuclear[3]), ]
nuclear = nuclear[with(nuclear,order(-Effect.Size)),]
nuclear = nuclear[1:4,]
nuclear$ID = paste0(nuclear$CHROM,":",nuclear$POS)

new.geno = subset(geno, !(ID %in% nuclear$ID),c(4:382))
new.geno = t(new.geno)
new.map = subset(map, !(ID %in% nuclear$ID))

cov = as.data.frame(t(subset(geno, ID %in% nuclear$ID, c(4:382))))

rownames(cov) = NULL
colnames(cov) = paste0("V",c(23:26))

dir.create('Codes for GWAS 2.0')
do.call(file.remove, list(list.files("Codes for GWAS 2.0", full.names = TRUE)))
new.location = 'Codes for GWAS 2.0'
file.copy(file.path(old.location,oldfiles), new.location)
setwd(new.location)
write.table(cov, "Additional covariates for mt-n mapping.csv", row.names = FALSE)
write.table(new.geno, "Final SNP Tables for mt-n mapping.csv", row.names = FALSE)
write.table(new.map, "Final SNP map.txt", row.names = FALSE) 

setwd(ori.wd)

