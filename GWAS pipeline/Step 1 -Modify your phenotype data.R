rm(list=ls())

#This R code will help you remove/add strains to the phenotype datafile. Missing strains will get NA values.
##Strains that do not have genotypic data will be removed. 

#Get project directory
ori.wd = getwd()

#load in required packages
require(plyr)

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

write.table(pheno, "Petite frequency phenotype.txt")

