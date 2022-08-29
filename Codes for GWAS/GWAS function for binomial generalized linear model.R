
#Load the phenotype, genotype and covariate files
GWAS = function(yourdata,yourfirstvariable, yoursecondvariable){
  
  first = yourdata[[yourfirstvariable]]
  second = yourdata[[yoursecondvariable]]
  mito = as.matrix(yourdata$Mito)
  # need to go through and eliminate any of the individuals that don't have a phenotype
  # need to remove the correspond genotype files also
  # This is done because there is a chance that a nonmonomorphic marker will become monomorphic because of missing phenotype data
  #'remove' is a vector with TRUE if the pheno vector has a missing value at that row
  remove <-is.na(first)
  remove <-is.na(second)
  #convert all the genotypes for that individual into NA
  geno.working[remove == TRUE]<-NA
  covariates[remove ==TRUE] <- NA
  mito[remove ==TRUE] <-NA
  #need to convert the genotype matrix into factors (othewise it might assume something of an additive model)
  #nSNP is the number of SNPs determined by the number of colums in the matrix 'geno'
  # geno.working =na.omit(geno.working)
  # first = na.omit(first)
  # second = na.omit(second)
  # covariates = na.omit(covariates)
  # mito =na.omit(mito)
  
  source("Codes for GWAS/Load covariates.R")
  mito = as.factor(mito)
  
  nSNP<-ncol(geno.working)
  
  #preallocate a pvalue matrix with -1 values
  n.pvalue <-matrix(-1,nSNP)
  m.pvalue <- matrix(-1, nSNP)
  mn.pvalue <- matrix(-1,nSNP)
  #loop through and test each of the SNPs and store the P-value from the linear model	
  for (i in 1:nSNP){
    if (i %in% c(2495,4991,7486,9982,12477,14972,17468,19964,22458,24955)){
      print(paste0(round(100*i/nSNP, digits = 0),"% Completed!"))
      #convert to factors
      geno.temp<-factor(geno.working[,i])
      #use table and sum to determine how many factor levels have at least 2 entries
      n.allele<-matrix(table(geno.temp))
      notsingle<-sum(n.allele >1)
      
      #If only a single type of marker is present in more than 1 individual then set the pvalue to be NA
      #Otherwise run the general linear model and store the pvalue
      
      if (notsingle < 2) {
        n.pvalue[i]<- NA
        m.pvalue[i]<- NA
        mn.pvalue[i]<- NA
      }else{tempresults<-glm(cbind(first,second) ~ cov1+ cov2+cov3+cov4+cov5+cov6+cov7+cov8+cov9+cov10+cov11+cov12+cov13+cov14+cov15+cov16+cov17+cov18+cov19+cov20+cov21+cov22+geno.temp*mito, family = binomial)
      anovatable<-anova(tempresults, test = "Chi", dispersion = tempresults[[10]]/tempresults[[16]])
      n.pvalue[i]<-anovatable['geno.temp',][[5]]
      m.pvalue[i]<-anovatable['mito',][[5]]
      mn.pvalue[i]<-anovatable['geno.temp:mito',][[5]]
      }
    }else{
      geno.temp<-factor(geno.working[,i])
    n.allele<-matrix(table(geno.temp))
    notsingle<-sum(n.allele >1)
    if (notsingle < 2) {
      n.pvalue[i] <- NA
      m.pvalue[i]<- NA
      mn.pvalue[i] <- NA
    }else{tempresults<-glm(cbind(first,second) ~cov1+ cov2+ cov3+cov4+cov5+cov6+cov7+cov8+cov9+cov10+cov11+cov12+cov13+cov14+cov15+cov16+cov17+cov18+cov19+cov20+cov21+cov22+ geno.temp*mito, family = binomial)
    anovatable<-anova(tempresults, test = "Chi", dispersion = tempresults[[10]]/tempresults[[16]])
    n.pvalue[i]<-anovatable['geno.temp',][[5]]
    m.pvalue[i]<-anovatable['mito',][[5]]
    mn.pvalue[i]<-anovatable['geno.temp:mito',][[5]]
    }
  }
}
  nuclear = as.data.frame(n.pvalue[,1])
  mito = as.data.frame(m.pvalue[,1])
  epistasis = as.data.frame(mn.pvalue[,1])
  
  #convert the P-value to log
  n.logresults = (-log10(n.pvalue[,1]))
  m.logresults = (-log10(m.pvalue[,1]))
  mn.logresults = (-log10(mn.pvalue[,1]))
  
  #Combine SNP map with pvalues
  nuclear.df <<- cbind(map,nuclear,n.logresults)
  epistasis.df <<- cbind(map,epistasis, mn.logresults)
  mito.df <<-cbind(map,nuclear,m.logresults)
  colnames(nuclear.df)<<-c("CHROM", "POS","ID", "Pvalue", "neg_logP")
  colnames(mito.df) <<-c("CHROM", "POS","ID", "Pvalue", "neg_logP")
  colnames(epistasis.df)<<-c("CHROM", "POS","ID", "Pvalue", "neg_logP")
  
  print("GWAS Finished!")
}
