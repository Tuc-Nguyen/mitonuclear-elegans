require(ggplot2)
require(reshape2)
require(lme4)
require(plyr)
require(Rmisc)
require(ggbeeswarm)
require(gridExtra)
require(grid)
require(tidyverse)



df = read.table("FinalDataframe.txt", header = T, na.strings = NA)

conditions <- sort(unique(df$Condition))
final= data.frame()
for (condition in conditions) {
  subCondition <- subset(df, Condition == condition)
  second <- data.frame()
  for (i in unique(subCondition$Nuclear)) {
    first <- data.frame()
    for (j in unique(subCondition$Nuclear)) {
      tryCatch({
        if (i == j) {
          print("Nada!")
        } else {
          specific.combo <- subset(subCondition, Nuclear %in% c(i, j) & Mito %in% c(i, j))
          title <- paste("Mito-Nuclear Epistasis Analysis between", i, "and", j, "pair combination in", condition)
          
          a = aggregate( WORMLENGTH ~ Nuclear+Mito+Status,FUN = mean,specific.combo)
          print(title)
          FM <- lmer(WORMLENGTH ~ Nuclear + Mito + Nuclear *Mito + (1|Batch), data = specific.combo)
          FM2 <- lmer(WORMLENGTH ~ Nuclear + Mito + (1|Batch), data = specific.combo)
          
          a <- anova(FM,FM2)
          p = a$`Pr(>Chisq)`[2]
          
          tempDF <- data.frame(Condition = condition, Background.1 = i, Background.2 = j, Epistatic = p)
          tempDF$Significance <- ifelse(tempDF$Epistatic < 0.05, "S", "NS")
          first <- rbind(first, tempDF)
        }
      }, error = function(e) {
        cat("Error occurred for condition:", condition, "\n")
        return(NULL)
      })
    }
    second <- rbind(second, first)
  }
  final <- rbind(final, second)
}


library(ggplot2)
library(gridExtra)
library(cowplot)

# Subset significant data
significant <- subset(final, Significance == "S")

# Create a new column that combines sorted Background.1 and Background.2 to represent unique pairs
significant$UniquePair <- apply(significant[, c("Background.1", "Background.2")], 1, function(x) paste(sort(x), collapse = "_"))
significant$UniqueConditionPair <- paste(significant$Condition, significant$UniquePair, sep = "_")
# Subset the significant dataframe to keep only unique pairs within each condition
significant <- significant[!duplicated(significant$UniqueConditionPair), ]


# Get unique conditions and create directory
conditions <- sort(unique(significant$Condition))
dir.create("2025 Epistatic Plot")
setwd("2025 Epistatic Plot")

for (condition in conditions) {
  subCondition <- subset(significant, Condition == condition)
  # Create an empty list to store plots for the current condition
  plot_list <- list()
  # Generate plots for all combinations of Background.1 and Background.2
  for (i in 1:nrow(subCondition)) {
  strain1 <- subCondition$Background.1[i]
  strain2 <- subCondition$Background.2[i]
  specific.combo <- subset(df, Condition == condition & Nuclear %in% c(strain1, strain2) & Mito %in% c(strain1, strain2))
  title <- paste( strain1, "and", strain2, "pairings")
  print(title)
        
  # Calculate means and standard errors for each combination of Nuclear and Mito
  s.c <- specific.combo %>%
    group_by(Nuclear, Mito, Status) %>%
    summarise(
      x = mean(WORMLENGTH, na.rm = TRUE),
      se = sd(WORMLENGTH, na.rm = TRUE) / sqrt(n())
    ) %>%
    ungroup()
  
  p <- ggplot(s.c) +
    aes(x = as.factor(Mito), color = as.factor(Nuclear), group = as.factor(Nuclear), shape = Status, y = x) +
    ggtitle(title) +
    xlab("Mitochondrial Haplotype") +
    ylab(expression(paste("Worm Length (", mu, "m)"))) +  # y-axis label with micron symbol
    geom_errorbar(aes(ymax = x + se, ymin = x - se), width = .05) +
    scale_color_manual(name = "Nuclear", values = c('red3', 'black')) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          panel.grid.major = element_line(colour = "gray", size = 0.3),
          strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = c(1, 0.5), 
          legend.justification = c(1, 0.5),
          legend.background = element_blank(),
          panel.background = element_rect(fill = "white")) +
    stat_summary(fun = mean, geom = "point", size = 4) +
    stat_summary(fun = mean, geom = "line")
  # Append plot to list
  plot_list[[length(plot_list) + 1]] <- p
        
        
      }
  # If there are plots to display, save them to a single PDF
  if (length(plot_list) > 0) {
    # Determine the number of rows and columns based on the number of plots
    ncol <- min(3, length(plot_list))  # Set max columns to 3, or less if fewer plots
    nrow <- ceiling(length(plot_list) / ncol)  # Determine number of rows needed
    
    # Adjust PDF size dynamically
    pdf_width <- ncol * 5  # Width per plot
    pdf_height <- nrow * 5  # Height per plot
    
    # Create a combined plot layout with uniform sizes using cowplot
    combined_plot <- plot_grid(plotlist = plot_list, ncol = ncol, align = "v", rel_widths = rep(1, ncol), rel_heights = rep(1, nrow))
    
    # Save all plots for this condition into one PDF file
    ggsave(filename = paste0(condition, ".pdf"), plot = combined_plot, width = pdf_width, height = pdf_height)
  }
  print("And... That's a wrap!")
}

nodes <- as.data.frame(c(1,2,3,4,5))
nodes$id <- paste0("s",seq.int(nrow(nodes)))
nodes$strain <- nodes$`c(1,2,3,4,5)`
nodes$strainname <- c("ECA2602","ECA2546","ECA1493","ECA1298","ECA1229")
nodes[1] <- NULL

require(igraph)
p <- par(mfrow=c(2,3))
for(condition in unique(final$Condition)){
  i = subset(final, Condition ==condition)
  i$from <- nodes$id[match(i$Background.1,nodes$strainname)] 
  i$to <- nodes$id[match(i$Background.2,nodes$strainname)] 
  i$weight <- i$Significance
  i$type <- i$Significance
  i[c(1:5)]<- NULL
  links<- na.omit(i)
  
  net <- graph_from_data_frame(d=links, vertices=nodes, directed = F)
  l <- layout_in_circle(net)
  title = paste(condition)
  #pdf(title) 
  plot<- plot(net,vertex.label=NA,edge.curved=0,vertex.color=c("#3A86AD","#05828B","#9B1A2D","#1B6831","#F9A003"),
              edge.lty = ifelse(links$weight =="NS", "dotted","solid"),
              edge.color= ifelse(links$weight =="NS", "gray90","gray40")
              ,layout=layout_in_circle, main = title)
  
  #dev.off()
}
par(p)

