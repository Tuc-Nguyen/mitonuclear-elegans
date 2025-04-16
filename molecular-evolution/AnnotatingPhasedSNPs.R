# Load required packages
library(vcfR)        # For reading and processing VCF files
library(ape)         # Phylogenetics tools
library(phytools)    # Additional phylogenetics tools
library(ggtree)      # Tree visualization
library(treeio)      # Tree import/export tools
library(ggplot2)     # Visualization
library(dplyr)       # Data manipulation
library(tidyr)       # Data reshaping
library(stringr)     # String manipulation
library(readxl)      # Reading Excel files

# Source custom function definitions
source("Functions.R", echo=TRUE)

# Read in annotated VCF data
all <- read.vcfR("FINAL.annotated.vcf", verbose = FALSE)

# Read list of phased variants
PHASE <- read_excel("PHASE.xlsx")

# Extract genotype table and join with VCF fixed fields (CHROM, POS, REF, ALT)
df <- cbind(as.data.frame(all@fix)[c(1,2,4,5)], extract.gt(all, element = 'GT', as.numeric = FALSE))

# Filter VCF to keep only the phased positions
PHASE <- subset(df, POS %in% PHASE$POS)

# Convert genotype matrix (GT format like 0|1) into numeric encodings (0, 1, 2)
genomat <- as.data.frame(apply(PHASE[c(5:544)], c(1, 2), convert_genotype))

# Add back SNP position from rownames, removing "MtDNA_" prefix
genomat$POS <- as.numeric(gsub("MtDNA_", "", rownames(genomat)))

# -----------------------------------------
# Group nearby positions together into potential codons (1-3 bp apart)
# -----------------------------------------
numbers <- sort(genomat$POS)
groups <- integer(length(numbers))
group_id <- 1
i <- 1

while (i <= length(numbers)) {
  groups[i] <- group_id
  j <- i + 1
  while (j <= length(numbers) && (numbers[j] - numbers[i]) <= 2) {
    groups[j] <- group_id
    j <- j + 1
  }
  i <- j
  group_id <- group_id + 1
}

# Assign group information to each position
genomat$Group <- groups

# Manually override one group (likely to correct an edge case at position 11389)
genomat$Group <- ifelse(genomat$POS == 11389, 111, genomat$Group)

# -----------------------------------------
# Count allele combinations per group (codon)
# -----------------------------------------
stats <- data.frame()
for (i in 1:max(genomat$Group)) {
  # Extract genotypes for the group and transpose to get individuals as rows
  codon <- as.data.frame(t(subset(genomat, Group == i, select = c(1:540))))
  codon <- na.omit(codon)
  
  # Concatenate alleles for each strain into single allele string (e.g., 000, 012)
  if (ncol(codon) == 2) {
    codon$Allele <- apply(codon[, 1:2], 1, function(x) paste0(x, collapse = ""))
  } else {
    codon$Allele <- apply(codon[, 1:3], 1, function(x) paste0(x, collapse = ""))
  }
  
  # Count number of times each codon configuration appears
  tempdf <- plyr::count(codon$Allele)
  tempdf$Group <- i
  
  stats <- rbind(stats, tempdf)
}

# Join codon frequencies with group and SNP metadata
phased <- stats
phased <- merge(phased, genomat[, c(541:542)])
phased <- merge(phased, df[, c(1:4)])

# Load reference sequence for codon reconstruction
reference <- read.table("Reference.txt", header = TRUE)

# Load SNP annotations including reference amino acid for each variant
TableS2 <- read_excel("TableS2.xlsx")

# -----------------------------------------
# Define translation function using mitochondrial genetic code (Table 5)
# -----------------------------------------
translate <- function(dna_seq) {
  dna_seq <- toupper(gsub("[^ATGC]", "", dna_seq))  # Remove non-ATGC characters
  dna_seq <- substr(dna_seq, 1, floor(nchar(dna_seq) / 3) * 3)  # Trim to codon multiple
  
  codon_table_5 <- c(
    "TTT" = "F", "TTC" = "F", "TTA" = "L", "TTG" = "L",
    "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
    "TAT" = "Y", "TAC" = "Y", "TAA" = "*", "TAG" = "*",
    "TGT" = "C", "TGC" = "C", "TGA" = "W", "TGG" = "W",
    "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
    "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    "CAT" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
    "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
    "ATT" = "I", "ATC" = "I", "ATA" = "M", "ATG" = "M",
    "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
    "AAT" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
    "AGT" = "S", "AGC" = "S", "AGA" = "S", "AGG" = "S",
    "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
    "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    "GAT" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
    "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
  )
  
  codons <- substring(dna_seq, seq(1, nchar(dna_seq), by = 3), seq(3, nchar(dna_seq), by = 3))
  aa <- sapply(codons, function(codon) codon_table_5[[codon]])
  paste(aa, collapse = "")
}


# Function to generate all possible allele combinations for codons involving two SNP positions
gen2 <- function(subcodon, change) {
  library(dplyr)
  
  # Step 1: Construct a long-format genotype table
  geno <- data.frame()
  for (snp in unique(subcodon$POS)) {
    subsnp <- subset(subcodon, POS == snp)
    subchange <- subset(change, POS == snp)
    
    REF <- subsnp$REF[1]
    
    if (nrow(subchange) == 1) {
      # Only one ALT allele
      ALT <- subsnp$ALT[1]
      tempdf <- data.frame(POS = snp, REF = REF, ALT = ALT, GR = 0, GA = 1)
    } else {
      # Multiple ALT alleles, limit to 3
      alts <- strsplit(subsnp$ALT[1], ",")[[1]]
      alts <- alts[1:min(3, length(alts))]
      tempdf <- data.frame(
        POS = rep(snp, length(alts)),
        REF = rep(REF, length(alts)),
        ALT = alts,
        GR = rep(0, length(alts)),
        GA = seq_along(alts)
      )
    }
    geno <- rbind(geno, tempdf)
  }
  
  # Step 2: Get unique alleles at each position
  allele_options <- geno %>%
    group_by(POS) %>%
    summarise(alleles = list(unique(c(REF, ALT))), .groups = "drop") %>%
    arrange(POS)
  
  if (nrow(allele_options) < 2) stop("Need at least two positions to generate combinations.")
  
  # Extract position labels and allele sets
  pos1 <- allele_options$POS[1]
  pos2 <- allele_options$POS[2]
  alleles_pos1 <- allele_options$alleles[[1]]
  alleles_pos2 <- allele_options$alleles[[2]]
  
  # Step 3: Generate all pairwise combinations
  combos <- expand.grid(
    allele_pos1 = alleles_pos1,
    allele_pos2 = alleles_pos2,
    stringsAsFactors = FALSE
  )
  
  # Create codon.gen based on 0-indexed allele position
  combos$codon.gen <- mapply(function(a1, a2) {
    idx1 <- match(a1, alleles_pos1) - 1
    idx2 <- match(a2, alleles_pos2) - 1
    paste0(idx1, idx2)
  }, combos$allele_pos1, combos$allele_pos2)
  
  # Step 4: Combine with a header showing POS info
  header <- data.frame(POS1 = pos1, POS2 = pos2, codon.gen = NA, stringsAsFactors = FALSE)
  combo_df <- data.frame(
    POS1 = combos$allele_pos1,
    POS2 = combos$allele_pos2,
    codon.gen = combos$codon.gen,
    stringsAsFactors = FALSE
  )
  
  result <- rbind(header, combo_df)
  return(result)
}
# Function to generate all possible codons from three SNP positions
gen3 <- function(subcodon, change) {
  library(dplyr)
  
  # Step 1: Build the expanded SNP-allele dataframe
  geno <- data.frame()
  for (snp in unique(subcodon$POS)) {
    subsnp <- subset(subcodon, POS == snp)
    subchange <- subset(change, POS == snp)
    
    REF <- subsnp$REF[1]
    alts <- if (nrow(subchange) == 1) {
      subsnp$ALT[1]
    } else {
      strsplit(subsnp$ALT[1], ",")[[1]][1:min(3, nrow(subchange))]
    }
    
    tempdf <- data.frame(
      POS = rep(snp, length(alts)),
      REF = rep(REF, length(alts)),
      ALT = alts,
      GR = rep(0, length(alts)),
      GA = seq_along(alts)
    )
    
    geno <- rbind(geno, tempdf)
  }
  
  # Step 2: Summarize unique alleles at each position
  allele_options <- geno %>%
    group_by(POS) %>%
    summarise(alleles = list(unique(c(REF, ALT))), .groups = "drop") %>%
    arrange(POS)
  
  if (nrow(allele_options) < 3) stop("Need at least three SNPs for 3-position codons.")
  
  # Get the allele options and labels for the first 3 SNPs
  pos1 <- allele_options$POS[1]
  pos2 <- allele_options$POS[2]
  pos3 <- allele_options$POS[3]
  a1 <- allele_options$alleles[[1]]
  a2 <- allele_options$alleles[[2]]
  a3 <- allele_options$alleles[[3]]
  
  # Step 3: Create the full codon and its genotype code
  combos <- expand.grid(a1 = a1, a2 = a2, a3 = a3, stringsAsFactors = FALSE)
  combos$codon <- paste0(combos$a1, combos$a2, combos$a3)
  
  # Generate codon.gen as string of 0-indexed allele positions
  combos$codon.gen <- mapply(function(x, y, z) {
    paste0(match(x, a1) - 1, match(y, a2) - 1, match(z, a3) - 1)
  }, combos$a1, combos$a2, combos$a3)
  
  # Step 4: Return the output table
  header <- data.frame(POS1 = pos1, POS2 = pos2, POS3 = pos3, Codon.Geno = NA)
  combo_df <- data.frame(
    POS1 = combos$a1,
    POS2 = combos$a2,
    POS3 = combos$a3,
    Codon.Geno = combos$codon.gen
  )
  
  return(rbind(header, combo_df))
}

# Initialize three data frames to collect processed codon data
biallelic <- data.frame()       # For codons with 2 SNPs and 2 alleles total
multiallelic <- data.frame()    # For codons with 2 SNPs and >2 alleles or 3 SNPs
others <- data.frame()          # (Unused in this code block, placeholder for other conditions)

# Loop over each codon group
for (codon in unique(phased$Group)) {
  # Subset data for the current codon group
  subcodon <- subset(phased, Group == codon)
  
  # Extract non-duplicated rows containing chromosome, group, genotype code, and frequency
  later <- subcodon[c(5, 2, 3, 4)]
  later <- later[!duplicated(later), ]
  colnames(later) <- c("CHROM", "Group", "Codon.Geno", "Frequency")
  
  # Get SNPs and allele information from reference table
  change <- subset(TableS2, POS %in% subcodon$POS)
  n_snp <- length(unique(change$POS))
  n_alleles <- nrow(change)
  
  # ==== CASE 1: Two SNPs with exactly two total alleles ====
  if (n_snp == 2 & n_alleles == 2) {
    pos1 <- min(subcodon$POS)
    pos2 <- max(subcodon$POS)
    geno <- gen2(subcodon, change)  # Generate genotype combinations
    
    if (pos1 + 2 == pos2) {
      # SNPs are at positions 1 and 3 → ref base must be at position 2
      pos3 <- pos2 - 1
      subref <- subset(reference, POS == pos3)
      
      colnames(geno) <- c("POS1", "POS3", "Codon.Geno")
      geno$POS2 <- subref$REF[1]
      geno$POS2[1] <- subref$POS[1]
      
      final <- geno[c(1, 4, 2, 3)]
      final <- plyr::join(final, later, type = "left")
      biallelic <- rbind(biallelic, final)
    } else {
      # SNPs are at adjacent positions, determine codon frame
      pos3 <- pos2 + 1
      pos0 <- pos1 - 1
      subref <- subset(reference, POS %in% c(pos0, pos1, pos2, pos3))
      change <- subset(TableS2, POS %in% c(pos1, pos2))
      
      # Try two possible reading frames to find correct amino acid context
      possible.codon1 <- paste0(subref$REF[1], subref$REF[2], subref$REF[3])
      possible.codon2 <- paste0(subref$REF[2], subref$REF[3], subref$REF[4])
      possible.aa1 <- translate(possible.codon1)
      possible.aa2 <- translate(possible.codon2)
      
      if (possible.aa1 == change$AA_REF[1]) {
        # SNPs are at positions 2 and 3, ref is at position 1
        subref <- subset(subref, POS == pos0)
        colnames(geno) <- c("POS2", "POS3", "Codon.Geno")
        geno$POS1 <- subref$REF[1]
        geno$POS1[1] <- subref$POS[1]
        final <- geno[c(4, 1, 2, 3)]
      } else {
        # SNPs are at positions 1 and 2, ref is at position 3
        subref <- subset(subref, POS == pos3)
        colnames(geno) <- c("POS1", "POS2", "Codon.Geno")
        geno$POS3 <- subref$REF[1]
        geno$POS3[1] <- subref$POS[1]
        final <- geno[c(1, 2, 4, 3)]
      }
      
      final <- plyr::join(final, later, type = "left")
      biallelic <- rbind(biallelic, final)
    }
    
    # ==== CASE 2: Two SNPs with more than two alleles (multiallelic) ====
  } else if (n_snp == 2 & n_alleles > 2) {
    pos1 <- min(subcodon$POS)
    pos2 <- max(subcodon$POS)
    geno <- gen2(subcodon, change)
    
    if (pos1 + 2 == pos2) {
      # SNPs at positions 1 and 3
      pos3 <- pos2 - 1
      subref <- subset(reference, POS == pos3)
      
      colnames(geno) <- c("POS1", "POS3", "Codon.Geno")
      geno$POS2 <- subref$REF[1]
      geno$POS2[1] <- subref$POS[1]
      final <- geno[c(1, 4, 2, 3)]
      final <- plyr::join(final, later, type = "left")
      multiallelic <- rbind(multiallelic, final)
    } else {
      # Need to resolve frame ambiguity
      pos3 <- pos2 + 1
      pos0 <- pos1 - 1
      subref <- subset(reference, POS %in% c(pos0, pos1, pos2, pos3))
      change <- subset(TableS2, POS %in% c(pos1, pos2))
      
      possible.codon1 <- paste0(subref$REF[1], subref$REF[2], subref$REF[3])
      possible.codon2 <- paste0(subref$REF[2], subref$REF[3], subref$REF[4])
      possible.aa1 <- translate(possible.codon1)
      possible.aa2 <- translate(possible.codon2)
      
      if (possible.aa1 == change$AA_REF[1]) {
        subref <- subset(subref, POS == pos0)
        colnames(geno) <- c("POS2", "POS3", "Codon.Geno")
        geno$POS1 <- subref$REF[1]
        geno$POS1[1] <- subref$POS[1]
        final <- geno[c(4, 1, 2, 3)]
      } else {
        subref <- subset(subref, POS == pos3)
        colnames(geno) <- c("POS1", "POS2", "Codon.Geno")
        geno$POS3 <- subref$REF[1]
        geno$POS3[1] <- subref$POS[1]
        final <- geno[c(1, 2, 4, 3)]
      }
      
      final <- plyr::join(final, later, type = "left")
      multiallelic <- rbind(multiallelic, final)
    }
    
    # ==== CASE 3: More than 2 SNPs, handled as multiallelic ====
  } else {
    geno <- gen3(subcodon, change)                    # Generate all codon combinations
    final <- plyr::join(geno, later, type = "left")   # Merge with observed data
    multiallelic <- rbind(multiallelic, final)
  }
}


# ==== For Biallelic Data ====

# Construct codon sequence (GENO) only if all three positions are valid nucleotides
biallelic$GENO <- ifelse(
  biallelic$POS1 %in% c("A", "T", "C", "G") &
    biallelic$POS2 %in% c("A", "T", "C", "G") &
    biallelic$POS3 %in% c("A", "T", "C", "G"),
  paste0(biallelic$POS1, biallelic$POS2, biallelic$POS3),
  NA
)

# Translate codon to amino acid using a row-wise function
biallelic$AA <- apply(biallelic, 1, function(row) {
  codon <- row["GENO"]
  if (!is.na(codon) && nchar(codon) == 3) {
    translate(codon)  # Use your previously defined codon-to-AA function
  } else {
    NA
  }
})


# ==== For Multiallelic Data ====

# Construct GENO for valid codons (3-letter A/T/C/G only)
multiallelic$GENO <- ifelse(
  multiallelic$POS1 %in% c("A", "T", "C", "G") &
    multiallelic$POS2 %in% c("A", "T", "C", "G") &
    multiallelic$POS3 %in% c("A", "T", "C", "G"),
  paste0(multiallelic$POS1, multiallelic$POS2, multiallelic$POS3),
  NA
)

# Translate GENO into amino acids for multiallelic codons
multiallelic$AA <- apply(multiallelic, 1, function(row) {
  codon <- row["GENO"]
  if (!is.na(codon) && nchar(codon) == 3) {
    translate(codon)
  } else {
    NA
  }
})
