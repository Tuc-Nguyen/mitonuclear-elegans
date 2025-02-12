library(vcfR)
library(ape)
library(phytools)
library(readxl)

TableS2 <- read_excel("TableS2.xlsx")
functional = subset(TableS2, Variant_Type %in% c("Missense","NonCoding"))

all <- read.vcfR("FINAL.annotated.vcf", verbose = FALSE)
df = cbind(as.data.frame(all@fix)[c(1,2,4,5)],extract.gt(all, element = 'GT', as.numeric = F))
df = subset(df, POS %in% functional$POS)

genomat <- as.matrix(t(df[,5:544]))

allele.counts <- apply(genomat, 2, function(col) length(unique(col)))

AlleleCount <- data.frame(
  SNP = names(allele.counts),  # Column names from genomat
  Counts = allele.counts   # Corresponding unique counts
)

MultiCount = subset(AlleleCount, Counts > 3)
MultiCount = subset(AlleleCount, Counts == 3)
AlleleCount = subset(AlleleCount, Counts == 2)


convert_genotype <- function(genotype) {
  if (is.na(genotype)) {
    return(NA)
  } else if (genotype == "0/0") {
    return(0)
  } else if (genotype == "1/1") {
    return(1)
  } else if (genotype == "2/2") {
    return(2)
  } else if (genotype == "3/3") {
    return(3)
  } else {
    return(NA)
  }
}

tree <- read.tree(file = "FINAL.treefile")

genomat = genomat[rownames(genomat) %in% tree$tip.label,]
genomat <- apply(genomat, c(1, 2), convert_genotype)
geno.matrix = genomat[,(colnames(genomat) %in% MultiCount$SNP)]
genomat = genomat[,(colnames(genomat) %in% AlleleCount$SNP)]

# make a temp stochastic mapping
make.simmap(tree, x = genomat[,15]) -> testsimmap

# extract the transition matrix
fixedq <- testsimmap$Q

# the output of make.simmap has these parts:
# edge: 731x2, naming each edge by its two nodes
# edge.length: 731x1, branchlengths
# Nnode: 363, number of internal nodes
# tip.label: 363, names of tips
# root.edge, 1 
# maps: a list with 673 elements, one per edge. 
# It's one or two elements, indicating how much of the length of that branch spent in each state, 0 and 1 if both, or just one if one.
# So names(x$maps[[integer]]) returns the state along that branch if no change
# names(x$maps[[integer]]) returns the state along that 
# each element of maps - corresponding element of edge.length gives the time spent in that state
# With parsimony-like reconstruction, this gives the branches that have changes:
# lx <- NULL
# for(i in 1:673){
	# lx[i] <- length(x$maps[[i]])
# }
# x$edge[lx == 2,]


# Infer the branches for each mutation, under the assumption that each site changes very rarely.
# This is parsimony, but we're doing it using stochastic character mapping because R packages don't have a parsimony option. 

# This takes a long time to run. I've saved the results (AllTheChanges) so don't run it again.


AllTheChanges <- matrix(nrow = 727, ncol = 419, data = 0)

for(snp in which(apply(genomat, 2, sum) >0)){
	snpchange <- make.simmap(tree, x = genomat[,snp], Q = 0.001*fixedq)
	for(i in 1:727){
	AllTheChanges[i,snp] <- length(snpchange$maps[[i]])-1
	}
	print(paste("finished",snp))
	}

colnames(AllTheChanges) = colnames(genomat)
write.table(AllTheChanges, "Allthechanges.csv", row.names = FALSE)

####for multiallelic 

# extract the transition matrix
multichange <- matrix(nrow = 727, ncol = 43, data = 0)
make.simmap(tree, x = geno.matrix[,3]) -> testsimmap

# extract the transition matrix
fixedq <- testsimmap$Q

for(snp in which(apply(geno.matrix, 2, sum) >0)){
  snpchange <- make.simmap(tree, x = geno.matrix[,snp], Q = 0.001*fixedq)
  for(i in 1:727){
    multichange[i,snp] <- length(snpchange$maps[[i]])-1
  }
  print(paste("finished",snp))
}

colnames(multichange) = colnames(geno.matrix)
write.table(multichange, "Multichanges.csv", row.names = FALSE)


multi <- matrix(nrow = 727, ncol = 3, data = 0)
make.simmap(tree, x = geno.matrix[,3]) -> testsimmap

# extract the transition matrix
fixedq <- testsimmap$Q

for(snp in which(apply(geno.matrix, 2, sum) >0)){
  snpchange <- make.simmap(tree, x = geno.matrix[,snp], Q = 0.001*fixedq)
  for(i in 1:727){
    multi[i,snp] <- length(snpchange$maps[[i]])-1
  }
  print(paste("finished",snp))
}

colnames(multi) = colnames(geno.matrix)
write.table(multi, "Fourchanges.csv", row.names = FALSE)

#allthechanges = cbind(AllTheChanges, multichange, multi)

###After running all SNPs, saved AllTheChanges file and reload it to save time
#library(readxl)
allthechanges <- read_excel("Allthechanges.xlsx")

functional$ID = paste0("MtDNA_",functional$POS)

# This method doesn't work with missing data. Of the 488 sites, 13 have some NAs, so these 13 sites get thrown out. 	

# RESULTS
# 1. There are lots of sites here that require multiple independent mutations to fit on the tree:
table(apply(allthechanges, 2, sum))

# 0   1   2   3   4   5   6   9  13 
# 13 365  57  19   4   3   1   1   1 
# That is, of 488 sites examined, 365 require a single mutation on the tree, but 86 require recurrent mutations. 
# One of them, site MtDNA_2042, requires 13 changes!!! 
# Also, of 595 mutations on the tree, 230 are due to the 86 recurrent-mutation sites. 
# Worth looking into these recurrent sites a bit-- what are their characteristics?


columns_with_changes <- colnames(allthechanges)[which(apply(allthechanges,2,sum) > 1)]
changes_count <- data.frame(Change_Count = apply(allthechanges[, columns_with_changes, drop = FALSE], 2, sum))
changes_count$SNP = row.names(changes_count)

recurrent = subset(functional, ID %in% colnames(allthechanges)[which(apply(allthechanges,2,sum) > 1)])
recurrent$Change_Count <- changes_count$Change_Count[match(recurrent$ID, changes_count$SNP)]

# 2. The number of mutations per branch ranges from 0 to 10.many branches are 0. The internal branches contains 35 
# mutations, which means the remaining 560 mutations are found at the tips of the tree. 

table(apply(allthechanges, 1, sum))


# 0   1   2   3   4   5   6   7   8   9  10 
# 424 175  57  30  18   8   8   3   2   1   1


library(dplyr)
library(ggtree)
library(ggplot2)

#3. Diagnostic changes for each mitogroup

# A.
colnames(allthechanges)[which(allthechanges[37,] ==1)]
# "MtDNA_382"	              nduo-6   
# "MtDNA_911"               s-rRNA  
# "MtDNA_10407"             l-rRNA 
# "MtDNA_12906" 406V>406M   nduo-5
# "MtDNA_10003" 119L>119M   ctc-2
# B. 
colnames(allthechanges)[which(allthechanges[39,] ==1)]
# "MtDNA_227"     39V>39I   nduo-6
# "MtDNA_2476"    238K>238N nduo-1
# "MtDNA_9634"              tRNA-Gly  
# "MtDNA_10569"             l-rRNA    
# C. 
colnames(allthechanges)[which(allthechanges[437,] ==1)]
# "MtDNA_155" 15I>15V       nduo-6  
# "MtDNA_10372"             tRNA-His
# "MtDNA_9242" 
# D.
colnames(allthechanges)[which(allthechanges[262,] ==1)]
# "MtDNA_911"               s-rRNA
# "MtDNA_5668"              tRNA-Leu
# E. 
colnames(allthechanges)[which(allthechanges[483,] ==1)]
# "MtDNA_911"               s-rRNA   
# "MtDNA_3619" 68G>68S      nduo-2 
# "MtDNA_10415"             l-rRNA 
# "MtDNA_11347"             l-rRNA
# F. 
colnames(allthechanges)[which(allthechanges[492,] ==1)]
# "MtDNA_4202" 262T>262I    nduo-2 
# "MtDNA_4249" 278V>278I    nduo-2 
# "MtDNA_9596"              tRNA-Gly
# G. 
colnames(allthechanges)[which(allthechanges[481,] ==1)]
# "MtDNA_883"               tRNA-Glu  
# "MtDNA_1278"              s-rRNA 
# "MtDNA_1643"              tRNA-Ser  
# "MtDNA_4663"      54P>54S ctb-1 
# "MtDNA_11380"     9V>9L   nduo-3
# H. 
colnames(allthechanges)[which(allthechanges[529,] ==1)]
# "MtDNA_1111"              s-rRNA 
# "MtDNA_3994" 193S>193G    nduo-2 
# "MtDNA_9502"              tRNA-Met
# "MtDNA_2279"  173L>173M   nduo-1
# "MtDNA_4284"              tRNA-Ile
# "MtDNA_5668"         	    tRNA-Leu
# A+H. 
colnames(allthechanges)[which(allthechanges[38,] ==1)]
# "MtDNA_3614"  66C>66F     nduo-2 



# So some of the recurrent mutations show up multiple times here-- 911.

edge_matrix = as.data.frame(tree$edge)
#edge_matrix$branch = rownames(edge_matrix)
colnames(edge_matrix)= c("parent", "node")


fin = as.data.frame(allthechanges)
functional$ID = paste0(functional$CHROM, "_", functional$POS)
functional <- functional %>%
  mutate(CDS = case_when(
    grepl("rRNA", ORF) ~ "rRNA",        # If ORF contains "rRNA", assign "rRNA"
    grepl("tRNA", ORF) ~ "tRNA",        # If ORF contains "tRNA", assign "tRNA"
    TRUE ~ ORF                         # Otherwise, keep the original value
  ))


p <- readRDS("final_ggtree.rds")


for (orf in unique(functional$CDS)){
  sub.orf = subset(functional, CDS ==orf)
  
  snps = unique(sub.orf$ID)
  some.changes = fin[colnames(fin) %in% snps]
  some.changes$Sum = rowSums(some.changes)
  
  some.changes$color = ifelse(some.changes$Sum>0, T,F)
  a = as.integer(ncol(some.changes))
  b = some.changes[,c(a-1,a)]
  some.changes = cbind(edge_matrix, some.changes[,c(a-1,a)])
  #some.changes = some.changes[c(3:5)]
  
  colnames(some.changes) = c("parent","node","changes","color")
  class(some.changes$node) <- "integer"
  assign(orf, some.changes)
  
  plot.temp = p
  plot.temp$data <- plot.temp$data %>%
    left_join(some.changes, by = c("parent","node"))
  
  final_plot <- plot.temp + geom_tree(aes(color = as.factor(color), size = ifelse(changes == 0, 0.05, 0.5*changes)),layout = "daylight")+
    geom_tippoint(aes(color = "white"), size = 2) +
    scale_size_continuous(
      range = c(0.05, 1),  # Branch thickness range: Thin for no changes, thicker for more changes
    ) +
    scale_color_manual(values = c("FALSE" = "grey90", "TRUE" = "red"))+
    #ggtitle(orf)+
    theme(plot.title = element_text(hjust = 0.5),
      legend.position = "none")
  
  # Save to PDF
  ggsave(
    filename = paste0(orf, ".pdf"),  # Save file as ORF name
    plot = final_plot,
    device = "pdf",
    width = 15,  # Set width
    height = 15  # Set height
  )
}

# A couple different ways to plot things
# branchlengths scaled by number of functional mutations
# Calculate the sum of changes for each branch
change_sums <- apply(allthechanges, 1, sum)

# Which are the branches that connect all the mitogroups?
internalbranches <-c(37:39,529,260:262,480:483,492,437)

internalbranch_data <- data.frame(
  parent = tree$edge[internalbranches,1],
  node = tree$edge[internalbranches,2], 
  sum_changes = apply(allthechanges[internalbranches, ], 1, sum)  # Sums of changes
)


p$data <- p$data %>%
  left_join(internalbranch_data, by = c("parent","node"))

p = p + geom_tippoint(aes(color = Clade)) + #geom_tiplab(check.overlap = T)+
  scale_color_manual(values = c("A"= '#45919B',"B"='#1686AD',"C"='#777CB2',"D"='#EF8B07', "E"='#1B6831',"F"='#90A04A', "G" = "gray30","H" = "#BC242F"))+
  theme(legend.position = "none")+
  geom_label(aes(x = branch, label = sum_changes),
    size =5, # Adjust label size
    color = "blue" )

