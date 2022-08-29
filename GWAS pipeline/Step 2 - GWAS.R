
#Get project directory
ori.wd = getwd()

# Load in all required files. This will take a minute or two
source("Codes for GWAS/Load in SNP tables, covariates, SNP maps, and annotations.R")
source("Codes for GWAS/GWAS function for binomial generalized linear model.R")

#Create empty dataframe to store generated data
nuclear.df = data.frame()
epistasis.df = data.frame()


#GWAS function requires you to fill in the formula as followed:
## GWAS(yourdata,yourfirstvariable, yoursecondvariable)

GWAS(pheno,'Petite','Grande')


write.table(nuclear.df, "Nuclear SNP.txt")
write.table(epistasis.df, "Epistasis SNP.txt")

