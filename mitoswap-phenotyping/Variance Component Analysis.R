require(ggplot2)
require(reshape2)
require(lme4)
require(plyr)
require(Rmisc)
require(ggbeeswarm)
require(gridExtra)
require(grid)
require(tidyverse)

ori.wd = getwd()

df = read.table("Processed Summary.txt", header = T, na.strings = NA)
data = data.frame()
for (condition in sort(unique(df$Metadata_Plate))){
  subcondition = subset(df, Metadata_Plate == condition)
  for (individual in sort(unique(subcondition$Strain))){
    substrain = subset(subcondition, Strain == individual)
    for (batch in sort(unique(substrain$Metadata_Date))){
      subdate = subset(substrain, Metadata_Date == batch)
      n = plyr::count(subdate, "Metadata_Well")
      outlier = subset(n, freq <5 | freq > 50) ##remove well with low worms or too many worms
      ##remove outliers
      tempdf = subset(subdate, worm_length_um < mean(subdate$worm_length_um)+ sd(subdate$worm_length_um) & worm_length_um > mean(subdate$worm_length_um)- sd(subdate$worm_length_um))
      if (nrow(outlier)==0){
        tempdf = tempdf
        data = rbind(data, tempdf)
      }else{
        tempdf = subset(tempdf, !(Metadata_Well %in% outlier$Metadata_Well))
        data = rbind(data, tempdf)
      }
    }
  }
}

df = data[c(14,8:10,2:3,7)]
colnames(data) = c("Condition","Strain","Nuclear","Mito","Replicate","Batch","WORMLENGTH")

conditions = sort(unique(df$Condition))
variance.component= data.frame()
for(condition in conditions){
  
  subcondition = subset(df, Condition == condition)
  FM = lmer(WORMLENGTH ~ 1 + (1|Batch) +(1|Nuclear) + (1|Mito) + (1|Nuclear:Mito), data = subcondition) 
  
  vDF = as.data.frame(VarCorr(FM))
  
  proportion_variance = round((100*vDF$vcov[c(1,2,3,5)] / sum(vDF$vcov[c(1,2,3,5)])), digits = 1)
  
  tempdf = data.frame(Condition = condition,
                      H2_Nuclear = proportion_variance[3],
                      H2_Mitochondrial = proportion_variance[2],
                      H2_Mitonuclear = proportion_variance[1],
                      Environment = proportion_variance[4],
                      Variance = sum(vDF$vcov[c(1,2,3,5)]),
                      H2 = 100-proportion_variance[4],
                      Vnuc_Vg = round((100*vDF$vcov[3] / sum(vDF$vcov[c(1,2,3)])), digits = 1),
                      Vmt_Vg = round((100*vDF$vcov[2] / sum(vDF$vcov[c(1,2,3)])), digits = 1),
                      Vmitonuc_Vg = round((100*vDF$vcov[1] / sum(vDF$vcov[c(1,2,3)])), digits = 1))
  variance.component = rbind(variance.component, tempdf)
  
}


proportion = variance.component[c(1,8:10,7)]
colnames(proportion) = c("Condition","Nuclear", "Mitochondrial","Epistasis", "H2")

proportion = reshape2::melt(proportion, id=c("Condition","H2"))

proportion[['variable']] = factor(proportion[['variable']], levels = rev(unique(proportion$variable)))

proportion$size = 2+4*(proportion$H2-min(proportion$H2))/(max(proportion$H2)-min(proportion$H2))

colors = c("#c66b3d","#c4a35a","#26495c" )
ggplot(proportion, aes(x= 2,y=value,fill=variable,width = 2))+
  geom_bar(stat='identity')+
  #geom_text(aes(label = value))+
  facet_wrap(~Condition, ncol=3) +
  coord_polar("y")+
  scale_fill_manual(values = colors, name ="Component", labels=c("Mt-n Epistasis", "Mitochondrial", "Nuclear"))+
  theme(panel.grid.major = element_line(colour = "gray", linewidth =0.3), 
        axis.text.x = element_blank(),axis.title.x=element_blank(),
        axis.text.y = element_blank(),axis.title.y =element_blank(),axis.ticks = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white"))

df$Status = ifelse(df$Nuclear==df$Mito, "Matched","Mismatched")

write.table(df, "FinalDataframe.txt", row.names = FALSE)

