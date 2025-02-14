require(ggplot2)
require(reshape2)
require(lme4)
require(plyr)
require(Rmisc)
require(ggbeeswarm)
require(gridExtra)
require(grid)
require(tidyverse)
library(cowplot)


df = read.table("FinalDataframe.txt", header = T, na.strings = NA)

working = aggregate(WORMLENGTH ~ Condition + Nuclear + Mito + Status, data = df, mean)

mean_diff = data.frame()
for(condition in unique(working$Condition)){
  subDF = subset(working, Condition == condition)
  cond.= subDF$Condition[1]
  nuclears = sort(unique(subDF$Nuclear))
  mito = subDF$Mito
  a = data.frame()
  for (nuclear in nuclears){
    
    sub_df = subset(subDF, subDF$Nuclear == nuclear)
    sub_df = sub_df[order(sub_df$Status),]
    nat = sub_df$Status
    c= sub_df$Condition[1]
    nuc = sub_df$Nuclear[1]
    mit = sub_df$Mito
    
    ori_mean = sub_df$WORMLENGTH[1]
    diff = 1+(sub_df$WORMLENGTH-ori_mean)/ori_mean
 
    tempDF = data.frame(Condition = c, Nuclear=nuc, Mito=mit, WORMLENGTH = sub_df$WORMLENGTH, Diff = diff, Native = nat)
    a = rbind(a, tempDF)
  }
  
  mean_diff = rbind(mean_diff , a)}


summary.data = data.frame()

for (condition in conditions) {
  subcondition = subset(df, Condition == condition)
  nuclears = sort(unique(subcondition$Nuclear))
  for (nuclear in nuclears) {
    subnuclear = subset(subcondition, Nuclear == nuclear)
    synthetic = subset(subnuclear, Mito != nuclear)
    original = subset(subnuclear, Mito == nuclear)
    
    # Skip analysis if "original" is empty
    if (nrow(original) == 0) {
      next  # Skip this iteration and move to the next one
    }
    mitos = sort(unique(synthetic$Mito))
    for (mito in mitos) {
      submito = subset(subnuclear, Mito == mito)
      work = rbind(original, submito)
      
      num_batches <- length(unique(work$Batch))
      
      if (num_batches > 1) {
        FM <- lmer(WORMLENGTH ~ Mito + (1|Batch), data = work)
        FM1 <- lmer(WORMLENGTH ~ (1|Batch), data = work)
        a = anova(FM, FM1)
        p = a$`Pr(>Chisq)`[2]
      } else {
        FM <- lm(WORMLENGTH ~ Mito, data = work)  # Use lm() if Batch has only 1 level
        FM1 <- lm(WORMLENGTH ~ 1, data = work)
        a = anova(FM, FM1)
        p = a$`Pr(>F)`[2]
      }
      
      tempDF = data.frame(Condition = condition, Nuclear = nuclear, Mito = mito, Pvalue = p)
      summary.data = rbind(summary.data, tempDF)
    }
  }
}



all = join(mean_diff,summary.data)
all$Significance = ifelse(all$Pvalue> 0.05,"NS",
                               ifelse(all$Diff> 1,"black","red"))


#####FOR GRAPHING
tempDF = all
tempDF$Native = "Matched"
tempDF$Diff = 1

final = rbind(tempDF,all)
conditions = sort(unique(final$Condition))

for (condition in unique(final$Condition)) {
  subcondition = subset(final, Condition == condition)
  plot.condition <- list()
  for (nuclear in unique(subcondition$Nuclear)) {  # Loop through each Nuclear
    subnuclear = subset(subcondition, Nuclear == nuclear)  
    title = paste(nuclear)
    p <- ggplot(data = subnuclear, aes(x = Native, y = Diff, group = Mito)) +
      geom_line(aes(color = Significance)) +
      geom_point(fill = 'black', shape = 21, size = 2, stroke = 0) +
      scale_y_continuous(limits = c(0.8, 1.15), breaks = seq(0.8, 1.15, by = 0.1), name = expression(paste(Delta, "Size"))) +
      scale_color_manual(name = NULL, values = c('NS' = 'gray90', 'red' = 'red', 'black' = 'black')) +
      xlab("Mitonuclear Combination") +
      ggtitle(title) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            strip.background = element_blank(),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            legend.justification = c(1, 0.5),
            legend.background = element_blank(),
            panel.background = element_rect(fill = "white"))
    plot.condition[[length(plot.condition) + 1]] <- p
  }
  # Create a combined plot layout with uniform sizes using cowplot
  combined_plot <- plot_grid(plotlist = plot.condition, ncol = 3, align = "v")
  # Save all plots for this condition into one PDF file
  ggsave(paste0(condition, ".pdf"), plot = combined_plot, width = 7.5 , height = 5)
}
 