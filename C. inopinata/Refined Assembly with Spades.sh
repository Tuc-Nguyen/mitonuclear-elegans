# Step 0: Set up environment by loading required tools
module purge  # Clear any previously loaded modules to ensure a clean environment
module load bwa/intel/0.7.17  # Load BWA for sequence alignment

# Step 1: Index the reference genome
# Index the C. inopinata mtDNA fasta file to prepare it for alignment
bwa index C.inopinata.fasta

# Step 2: Align unmapped reads to the C. inopinata reference genome using BWA
# Perform sequence alignment for each of the FASTQ files using BWA MEM and output SAM files
bwa mem C.inopinata.fasta unmapped.DRR093016.fastq > DRR093016.sam
bwa mem C.inopinata.fasta unmapped.DRR093017.fastq > DRR093017.sam
bwa mem C.inopinata.fasta unmapped.DRR093018.fastq > DRR093018.sam
bwa mem C.inopinata.fasta unmapped.DRR093019.fastq > DRR093019.sam

# Step 3: Convert SAM files to BAM format for further processing
# SAM -> BAM conversion
samtools view -bS DRR093016.sam -o mtDNA.DRR093016.bam
samtools view -bS DRR093017.sam -o mtDNA.DRR093017.bam
samtools view -bS DRR093018.sam -o mtDNA.DRR093018.bam
samtools view -bS DRR093019.sam -o mtDNA.DRR093019.bam

# Step 4: Extract only aligned reads (i.e., those that mapped to C. inopinata mtDNA)
# The -F 4 flag is used to filter out unmapped reads and retain only mapped ones
samtools view -b -F 4 mtDNA.DRR093016.bam -o aligned.DRR093016.bam
samtools view -b -F 4 mtDNA.DRR093017.bam -o aligned.DRR093017.bam
samtools view -b -F 4 mtDNA.DRR093018.bam -o aligned.DRR093018.bam
samtools view -b -F 4 mtDNA.DRR093019.bam -o aligned.DRR093019.bam

# Step 5: Convert aligned BAM files back to FASTQ format
# Convert the aligned reads (mapped reads) to FASTQ for further downstream analysis
bedtools bamtofastq -i aligned.DRR093016.bam -fq mapped.DRR093016.fastq
bedtools bamtofastq -i aligned.DRR093017.bam -fq mapped.DRR093017.fastq
bedtools bamtofastq -i aligned.DRR093018.bam -fq mapped.DRR093018.fastq
bedtools bamtofastq -i aligned.DRR093019.bam -fq mapped.DRR093019.fastq

# Step 6: Convert aligned BAM files to FASTA format to obtain the aligned reads in FASTA format
# Extract FASTA sequences of the aligned reads (optional for downstream use)
samtools fasta -F 0x2 -o mtDNA-DRR093016-reads.fasta aligned.DRR093016.bam
samtools fasta -F 0x2 -o mtDNA-DRR093017-reads.fasta aligned.DRR093017.bam
samtools fasta -F 0x2 -o mtDNA-DRR093018-reads.fasta aligned.DRR093018.bam
samtools fasta -F 0x2 -o mtDNA-DRR093019-reads.fasta aligned.DRR093019.bam

# Step 7: Sort the BAM files to prepare them for indexing and variant calling
# Sorting the aligned BAM files by coordinate
samtools sort aligned.DRR093016.bam -o sorted_DRR093016.bam
samtools sort aligned.DRR093017.bam -o sorted_DRR093017.bam
samtools sort aligned.DRR093018.bam -o sorted_DRR093018.bam
samtools sort aligned.DRR093019.bam -o sorted_DRR093019.bam

# Step 8: Index the sorted BAM files
# Indexing the sorted BAM files for fast access during variant calling
samtools index sorted_DRR093016.bam
samtools index sorted_DRR093017.bam
samtools index sorted_DRR093018.bam
samtools index sorted_DRR093019.bam

# Step 9: Perform variant calling using Samtools mpileup and bcftools
# Generate a consensus FASTQ file for each of the samples
samtools mpileup -uf C.inopinata.fasta sorted_DRR093016.bam | bcftools call -c | vcfutils.pl vcf2fq > DRR093016.consensus.fastq
samtools mpileup -uf C.inopinata.fasta sorted_DRR093017.bam | bcftools call -c | vcfutils.pl vcf2fq > DRR093017.consensus.fastq
samtools mpileup -uf C.inopinata.fasta sorted_DRR093018.bam | bcftools call -c | vcfutils.pl vcf2fq > DRR093018.consensus.fastq
samtools mpileup -uf C.inopinata.fasta sorted_DRR093019.bam | bcftools call -c | vcfutils.pl vcf2fq > DRR093019.consensus.fastq


# Step 10: Perform hybrid assembly using SPAdes
# Using both Illumina and PacBio reads, with a trusted reference contig (C. inopinata fasta)

# The command below combines short-read sequencing (Illumina) and long-read sequencing (PacBio)
spades.py \
   --s1 mapped.DRR093016.fastq \  # Illumina paired-end read 1 (mapped reads)
   --s1 mapped.DRR093017.fastq \  # Illumina paired-end read 2 (mapped reads)
   --s1 mapped.DRR093018.fastq \  # Illumina paired-end read 3 (mapped reads)
   --s1 mapped.DRR093019.fastq \  # Illumina paired-end read 4 (mapped reads)
   --pacbio mapped.DRR093029.fastq \  # PacBio long-read sequencing data
   --trusted-contigs C.inopinata.fasta \  # Trusted contigs to guide the assembly, based on the C. inopinata reference genome
   -o Mapped-Hybrids  # Output directory for the hybrid assembly
