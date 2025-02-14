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

# Read in tables, with stringsAsFactors = T
# StrainCodes.txt modified to remove Originals, leaving only Synthetics
df = read.table("Processed Summary.txt", header = T, na.strings = NA)
df = subset(df, n>=5 & n <75)
#df = subset(df, !(Nuclear == "ECA2367" & Condition == "Copper"))
df = subset(df, !(Nuclear == "ECA2367" & Metadata_Date =="2012021"))
df = subset(df, Strain != "QG5050")
a = plyr::count(df, c("Strain","Condition"))


#REMOVE OUTLIER WELLS - Wells with median worm length within 2 standard deviations, and number of worms within 2 SDs, among replicates of that strain in that condition *within the same batch* 

# This removes 567 wells, leaving 8505

newdata = data.frame()
for (condition in sort(unique(df$Condition))){
  subcondition = subset(df, Condition == condition)
  for (individual in sort(unique(subcondition$Strain))){
    substrain = subset(subcondition, Strain == individual)
    for (batch in sort(unique(substrain$Metadata_Date))){
      subdate = subset(substrain, Metadata_Date == batch)
      tempdf = subset(subdate, median_wormlength_um < mean(subdate$median_wormlength_um)+ 2*sd(subdate$median_wormlength_um) & median_wormlength_um > mean(subdate$median_wormlength_um)- 2*sd(subdate$median_wormlength_um))
      tempdf = subset(tempdf, n < mean(subdate$n)+ 2*sd(subdate$n) & n > mean(subdate$n)- 2*sd(subdate$n))
      if (nrow(tempdf) < 1) {
        next  # skip this subgroup if there are not enough observations
      }
      tempdf$Replicate = c(1:nrow(tempdf))
      newdata = rbind(newdata, tempdf)
    }
  }
}

mydata = newdata[c(27,21,22,23,28,14,7,13)]
colnames(mydata) = c("Condition","Strain","Nuclear","Mito","Replicate","Batch","WORMLENGTH","n")

ggplot(mydata, aes(x= Condition,y=WORMLENGTH))+xlab("CONDITION")+
  geom_beeswarm(cex=0.65)+scale_y_continuous(breaks = seq(0,1300, by = 50))+
  theme(panel.grid.major = element_line(colour = "gray", linewidth =0.3), 
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white"))

conditions = sort(unique(mydata$Condition))
correction = data.frame()
variance.component= data.frame()
for(condition in conditions){
  
  subcondition = subset(mydata, Condition == condition)
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



# What all these columns represent:
# H2_Nuclear: proportion of phenotypic variance (after exclusion of batch effects) attributable to nuclear genotype
# H2_Mitochondrial: proportion of phenotypic variance (after exclusion of batch effects) attributable to mitochondrial genotype
# H2_Mitonuclear: proportion of phenotypic variance (after exclusion of batch effects) attributable to mitonuclear interactions
# Environment: proportion of phenotypic variance (after exclusion of batch effects) attributable to environmental effects
# The four above add up to 1.
# Variance: phenotypic variance (after exclusion of batch effects)
# H2: heritability, i.e., proportion of phenotypic variance (after exclusion of batch effects) attributable to genetic effects
# Vnuc_Vg: Fraction of genetic variance attributable to nuclear genotype
# Vmt_Vg: Fraction of genetic variance attributable to mitochondrial genotype
# Vmitonuc_Vg: Fraction of genetic variance attributable to mitonuclear interactions


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

df = mydata
df$Status = ifelse(df$Nuclear==df$Mito, "Matched","Mismatched")

write.table(df, "FinalDataframe.txt", row.names = FALSE)