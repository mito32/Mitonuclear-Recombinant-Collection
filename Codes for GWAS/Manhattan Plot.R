require(ggplot2)
palette = c("#666666", "#984ea3", "#377eb8", "#000075", "#66a61e", "#808000", "#800000", "#a65628", "#a6761d", "#ff7f00", "#e6ab02", "#1b9e77","#008080", "#4daf4a", "#7570b3", "#e41a1c")

manhattan.plot = function(yourdata, yourthreshold){
  ##Omit NA values
  yourdata = na.omit(yourdata)
  ##Extract chromosome numbers
  yourdata$Chr.number = as.numeric(gsub("[^[:digit:]]*", "", yourdata$CHROM))
  ## To put all chromosomes in one, the POS needs to be continuous without overlap.
  ### Maximum POS of each chromosome is extracted
  df.max = aggregate(POS~Chr.number, data = yourdata, FUN = max)
  ### The maximum POS value of the following chromosome is the sum of maximum POS of the previous chromosomes, with the first chromosome having the value of 0. 
  df.max$max = cumsum(df.max$POS)- df.max$POS
  df.max$POS <- NULL 
  ##New column containing the max value is then added to yourdata file (merged by Chromosome number)
  yourdata = merge(yourdata,df.max, by = "Chr.number")
  ##The new POS value is calculated by adding the max value to the old POS
  yourdata$newPOS = yourdata$POS + yourdata$max
  
  ### The new POS values are used to specify where to label the chromosome numbers on the x-axis 
  df.mean = aggregate(newPOS~Chr.number, data = yourdata, FUN = mean)
  x.ticks = df.mean$newPOS
  x.ticks.labels = as.roman(unique((yourdata$Chr.number)))
  max.log = max(yourdata$neg_logP)
  
  print(ggplot(data =yourdata, aes(y =neg_logP, x= newPOS))+
          geom_point(aes(colour = as.factor(Chr.number)))+
          scale_color_manual(values= palette)+
          scale_y_continuous(breaks = c(0:max.log)*1,name = "-log10(P)")+
          scale_x_continuous(name = "Chromosome", breaks = x.ticks, labels = x.ticks.labels) +
          geom_hline(yintercept= yourthreshold, linetype="dashed", color = "black")+
          theme(text = element_text(size=20, color = 'black'),
                axis.title.x = element_text(vjust=-1.5),
                axis.ticks = element_blank(),
                plot.margin = margin(0, -10, 15, 10),
                panel.background = element_rect(fill = "white"), legend.position = "none"))
  
}

