#Load in required packages
require(data.table)

#Set project directory to codes' location

setwd("Codes for GWAS")

#Load the genotype 
geno.working<-as.matrix(fread('Final SNP Tables for mt-n mapping.csv', header = FALSE))


#load in covariates
covariates <-as.matrix(fread('Covariates for mt-n mapping.csv'))

#Load in SNP Map
map = read.table('Final SNP Map.txt', header = TRUE)
map$ID = paste0(map$CHROM,":",map$POS)

#Load in Annotation File
annotation = read.table('Annotations.txt', header = TRUE)

#Reset working directory
setwd(ori.wd)
