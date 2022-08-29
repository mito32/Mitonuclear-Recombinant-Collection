
allele.extraction = function(yourdata, yourphenocolumn, youroutput){
  require(plyr)
  yourdata$ID = paste0(yourdata$CHROM,":",yourdata$POS)
  yourdata$Gene.ID = paste0(yourdata$GENE,"-",yourdata$POS)
  yourdata = yourdata[,(ncol(yourdata)-1):ncol(yourdata)]
  yourdata = unique( yourdata[,(1:2)])
  reduced.map = subset(map,ID %in% unique(yourdata$ID))
  reduced.map$location = row.names(reduced.map)
  reduced.map = join(reduced.map,yourdata)
  #SNP table should have been imported at this stage
  reduced.snp = as.data.frame(geno.working[,(as.integer(reduced.map$location))])
  colnames(reduced.snp) = as.character(reduced.map$Gene.ID)
  final <<- cbind(pheno,reduced.snp)
  
  
  hits = reduced.map$Gene.ID
  
  consolidate = data.frame()
  for (hit in hits){
    working.df = cbind(hit,final$Mito,final[[yourphenocolumn]], final[hit])
    working.df = na.omit(working.df)
    colnames(working.df) = c("GeneID",'Mitotype','Pheno','Hit')
    tempdf = aggregate(data = working.df, working.df$Pheno ~ ., FUN = mean)
    maxdf =aggregate(data = tempdf, tempdf$`working.df$Pheno` ~ tempdf$Mitotype, FUN = 'max')
    colnames(maxdf) = c('Mitotype','Max.Pheno')
    colnames(tempdf) = c('GeneID', 'Mitotype', 'Allele','Mean.Pheno')
    tempdf = merge(tempdf,maxdf)
    tempdf$Gene = (unlist(strsplit(hit,"-")))[1]
    tempdf$Position = (unlist(strsplit(hit,"-")))[2]
    consolidate = rbind(consolidate,tempdf)
  }
    assign(youroutput,consolidate,envir = .GlobalEnv)  
}

print("Read the following instructions to use the function correctly: The function format is : allele.extraction (yourdata, yourphenocolumn, youroutput) 1. The first element that should be included is the name of the dataframe containing your list of significant candidate genes. 2. The second element is the name of your phenotype columm in the original pheno dataframe in Step 2. 3. The last element is the name of your output dataframe.")
