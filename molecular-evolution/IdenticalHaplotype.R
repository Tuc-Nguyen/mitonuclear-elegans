# Load necessary libraries
library(vcfR)
library(ape)
library(phytools)
library(ggtree)
library(treeio)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)



all <- read.vcfR("FINAL.annotated.vcf", verbose = FALSE)
df = cbind(as.data.frame(all@fix)[c(1,2,4,5)],extract.gt(all, element = 'GT', as.numeric = F))
geno = df
genomat <- as.matrix(t(geno[,5:544]))
allele.counts <- apply(genomat, 2, function(col) length(unique(col)))

# Function to count individuals for each allele
count_alleles_multiallelic <- function(genotype_column) {
  # Get the unique alleles (removing NA)
  alleles <- unique(genotype_column[!is.na(genotype_column)])
  
  # Count occurrences of each allele
  allele_counts <- sapply(alleles, function(allele) sum(genotype_column == allele, na.rm = TRUE))
  
  # Return a named vector with allele counts
  return(allele_counts)
}

AlleleCount <- data.frame(
  SNP = names(allele.counts),  # Column names from genomat
  Counts = allele.counts   # Corresponding unique counts
)

# Apply the function to all columns of the genotype matrix
allele_counts_all <- apply(genomat, 2, count_alleles_multiallelic)

# Convert the result to a data frame for better readability
allele_counts_list <- lapply(names(allele_counts_all), function(snp) {
  counts <- allele_counts_all[[snp]]
  data.frame(
    SNP = snp,
    Allele = names(counts),
    Count = counts
  )
})

# Combine all allele counts into a single data frame
allele_frequency <- do.call(rbind, allele_counts_list)
rm(allele_counts_list)
rm(allele_counts_all)

# Display the allele counts
singletons = subset(allele_frequency, Count ==1)
het = subset(allele_frequency, !(Allele %in% c("0/0","1/1","2/2","3/3","4/4","5/5","6/6")))
hetsites = subset(allele_frequency, SNP %in% unique(het$SNP))

#turn all heterozygous values to NA
df[] <- lapply(df, function(col) {
  ifelse(col %in% unique(het$Allele), NA, col)
})

# Initialize an empty data frame to store results
results.1 <- data.frame()
results.2 <- data.frame()
# Loop over all pairs of columns between 5 and 554
for (i in 5:543) {  # Outer loop for the first column
  for (j in (i + 1):544) {  # Inner loop for the second column
    
    # Get the names of the columns being compared
    col1_name <- colnames(df)[i]
    col2_name <- colnames(df)[j]
    
    # Compare the values between the two columns
    diff_count.1 <- sum((df[[i]] != df[[j]]) | is.na(df[[i]]) | is.na(df[[j]]))
    diff_count.2 <- sum(df[[i]] != df[[j]], na.rm = TRUE)
    # Store the result in the result dataframe
    results.1 <- rbind(results.1, data.frame(Strain.1 = col1_name, Strain.2 = col2_name, Diff_Count = diff_count.1))
    results.2 <- rbind(results.2, data.frame(Strain.1 = col1_name, Strain.2 = col2_name, Diff_Count = diff_count.2))
  }
}

write.csv(results.1, "Results.1.csv", row.names = F)
write.csv(results.2, "Results.2.csv", row.names = F)

identicalpairs = subset(results.1, Diff_Count ==0)
a = data.frame( Strain.1 = unique(c(identicalpairs$Strain.1,identicalpairs$Strain.2)), Number.1 = c(1:218))
b = data.frame( Strain.2 = unique(c(identicalpairs$Strain.1,identicalpairs$Strain.2)), Number.2 = c(1:218))
identicalpairs = plyr::join(identicalpairs,a, type = "left")
identicalpairs = plyr::join(identicalpairs,b, type = "left")

assigned = data.frame()
for(i in c(1:218)){
  substrain = subset(identicalpairs, Number.1 == i | Number.2 == i )
  first = min(unique(substrain$Number.1), unique(substrain$Number.2))
  isotype = subset(a, Number.1 == first)[1,1]
  substrain$Isotype = isotype
  assigned = rbind(assigned, substrain)
}
final.isotype = data.frame()
for(iso in unique(assigned$Isotype)){
  subdf = subset(assigned, Isotype == iso)
  tempdf = data.frame(Strain = unique(c(subdf$Strain.1, subdf$Strain.2)), Haplotype = iso)
  final.isotype = rbind(final.isotype, tempdf)
}

final.isotype$Haplotype = ifelse(final.isotype$Haplotype == "ECA250", "MY1",
                                 ifelse(final.isotype$Haplotype =="CB4932", "N2",
                                        ifelse(final.isotype$Haplotype =="ECA1074", "ECA2546", 
                                               ifelse(final.isotype$Haplotype =="ECA1887", "ECA2041",
                                                      ifelse(final.isotype$Haplotype =="JT11398", "PB303",
                                                             ifelse(final.isotype$Haplotype =="JU2106", "MY10", 
                                                                    ifelse(final.isotype$Haplotype =="ECA1843", "ECA2367",
                                                                           ifelse(final.isotype$Haplotype =="EG4347", "PX179", final.isotype$Haplotype))))))))


others = data.frame(Strain = unique(c(results.1$Strain.1, results.1$Strain.2)))
others$Haplotype = others$Strain
others = subset(others, !(Strain %in% final.isotype$Strain))
others = rbind(others, final.isotype)
write.table(others, "FINAL ISOTYPES.txt", row.names = F)