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
require(ggpmisc)
require(ggpubr)

df = read.table("FinalDataframe.txt", header = T, na.strings = NA)
working = aggregate(WORMLENGTH ~ Condition + Nuclear + Mito + Status, data = df, mean)

mito_diff = data.frame()
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
    diff = abs(sub_df$WORMLENGTH-ori_mean)
    
    tempDF = data.frame(Condition = c, Nuclear=nuc, Mito=mit, WORMLENGTH = sub_df$WORMLENGTH, Diff = diff, Native = nat)
    a = rbind(a, tempDF)
  }
  
  mito_diff = rbind(mito_diff , a)}

nuclear_diff = data.frame()
for(condition in unique(working$Condition)){
  subDF = subset(working, Condition == condition)
  cond.= subDF$Condition[1]
  mitos = sort(unique(subDF$Mito))
  nuclear = subDF$Nuclear
  b = data.frame()
  for (mito in mitos){
    
    sub_df = subset(subDF, subDF$Mito == mito)
    sub_df = sub_df[order(sub_df$Status),]
    nat = sub_df$Status
    c= sub_df$Condition[1]
    mit = sub_df$Mito[1]
    nuc = sub_df$Nuclear
    
    ori_mean = sub_df$WORMLENGTH[1]
    diff = abs(sub_df$WORMLENGTH-ori_mean)
    
    tempDF = data.frame(Condition = c, Nuclear=nuc, Mito=mit, WORMLENGTH = sub_df$WORMLENGTH, Diff = diff, Native = nat)
    b = rbind(b, tempDF)
  }
  
  nuclear_diff = rbind(nuclear_diff , b)}

nuclear.dist = subset(genetic_distance, Type == "nuclear.diff")
mito.dist = subset(genetic_distance, Type == "mito.diff")

library(ggplot2)
library(cowplot)
library(ggpubr)

# -------------------------------
# Nuclear
# -------------------------------

# Get global y-axis limits for Diff in nuclear
nuclear_ymin <- min(nuclear_diff$Diff, na.rm = TRUE)
nuclear_ymax <- max(nuclear_diff$Diff, na.rm = TRUE)

correlation <- data.frame()
nuclear_plot_list <- list()

for (condition in unique(nuclear_diff$Condition)) {
  subcondition <- subset(nuclear_diff, Condition == condition)
  working.df <- merge(subcondition, nuclear.dist)
  working.df <- subset(working.df, Distance > 0)
  
  p <- ggplot(working.df, aes(x = Distance, y = Diff)) +
    geom_point(size = 0.75) +
    ggtitle(condition) +
    scale_y_continuous(name = "", limits = c(nuclear_ymin, nuclear_ymax)) +
    scale_x_continuous(name = "") +
    geom_smooth(method = "lm", se = FALSE, colour = 'red') +
    theme_light() +
    stat_cor(method = "pearson", cor.coef.name = "r")
  
  cor_result <- cor.test(x = working.df$Distance, y = working.df$Diff)
  summary.nuclear <- data.frame(
    condition = condition,
    R = cor_result$estimate,
    p.value = cor_result$p.value,
    Genome = "Nuclear"
  )
  correlation <- rbind(summary.nuclear, correlation)
  
  nuclear_plot_list[[condition]] <- p
}

nuclear_combined_plot <- plot_grid(plotlist = nuclear_plot_list, ncol = 2)


# -------------------------------
# Mitochondrial
# -------------------------------

# Get global y-axis limits for Diff in mitochondrial
mito_ymin <- min(mito_diff$Diff, na.rm = TRUE)
mito_ymax <- max(mito_diff$Diff, na.rm = TRUE)

plot_list <- list()

for (condition in unique(mito_diff$Condition)) {
  subcondition <- subset(mito_diff, Condition == condition)
  working.df <- merge(subcondition, mito.dist)
  working.df <- subset(working.df, Distance > 0)
  
  p <- ggplot(working.df, aes(x = Distance, y = Diff)) +
    geom_point(size = 0.75) +
    ggtitle(condition) +
    scale_y_continuous(name = "", limits = c(mito_ymin, mito_ymax)) +
    scale_x_continuous(name = "") +
    geom_smooth(method = "lm", se = FALSE, colour = 'red') +
    theme_light() +
    stat_cor(method = "pearson", cor.coef.name = "r")
  
  cor_result <- cor.test(x = working.df$Distance, y = working.df$Diff)
  summary.mito <- data.frame(
    condition = condition,
    R = cor_result$estimate,
    p.value = cor_result$p.value,
    Genome = "Mito"
  )
  correlation <- rbind(correlation, summary.mito)
  
  plot_list[[condition]] <- p
}

combined_plot <- plot_grid(plotlist = plot_list, ncol = 2)
