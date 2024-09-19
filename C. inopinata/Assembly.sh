
# Load necessary modules for the analysis
module load minimap2/2.22          # Load Minimap2 for sequence alignment
module load bedtools/intel/2.29.2  # Load Bedtools for handling BAM/FastQ format conversions
module load samtools/intel/1.14    # Load Samtools for working with SAM/BAM files

## Step 1: Index the reference genome ##
# Create an index for the reference genome (PRJDB5687_WS289.fa) to speed up the alignment process.
# REF.mmi will be the index file that minimap2 will use for alignment.
minimap2 -d REF.mmi PRJDB5687_WS289.fa  

## Step 2: Align the reads to the reference genome ##
# Align the reads from 'reads.fastq' to the reference genome index (REF.mmi).
# The output will be in SAM format (alignment.sam), which contains the alignment information.
minimap2 -a REF.mmi reads.fastq > alignment.sam 

## Step 3: Convert SAM to BAM format and extract unmapped reads ##
# Convert the SAM file (alignment.sam) into BAM format (aligned.bam) using samtools. 
# BAM format is binary, more compact, and quicker to process.
samtools view -bS -o aligned.bam alignment.sam

# Extract unmapped reads from the aligned BAM file (aligned.bam) and save them in a new BAM file (unmapped.bam).
# The -f 4 flag in samtools specifies that we are looking for unmapped reads.
samtools view -b -f 4 -o unmapped.bam aligned.bam

# Convert the unmapped BAM reads (unmapped.bam) back to FastQ format (unmapped.fastq).
# This allows you to realign these reads later, possibly to another reference.
bedtools bamtofastq -i unmapped.bam -fq unmapped.fastq

## Step 4: Map unmapped reads to the mitochondrial DNA of C. elegans ##
# Now, create an index for the C. elegans mitochondrial DNA (elegans-MtDNA.fasta).
# REF-elegans.mmi will be the index used for alignment against mitochondrial DNA.
minimap2 -d REF-elegans.mmi elegans-MtDNA.fasta  

# Align the unmapped reads (unmapped.fastq) to the C. elegans mitochondrial DNA index (REF-elegans.mmi).
# The output SAM file (inopinata-MtDNA.sam) will contain the alignment information for reads that potentially map to mtDNA.
minimap2 -a REF-elegans.mmi unmapped.fastq > inopinata-MtDNA.sam

## Step 5: Convert the SAM file to BAM and extract aligned reads ##
# Convert the mitochondrial DNA alignment SAM file (inopinata-MtDNA.sam) to BAM format (inopinata-MtDNA.bam).
samtools view -bS -o inopinata-MtDNA.bam inopinata-MtDNA.sam

# Extract the reads that are aligned (mapped) to the mitochondrial DNA.
# The -F 4 flag ensures that we only keep reads that are mapped.
samtools view -b -F 4 -o inopinata-MtDNA-aligned.bam inopinata-MtDNA.bam

# Convert the aligned mitochondrial reads from BAM format (inopinata-MtDNA-aligned.bam) to FASTA format (mtDNA-reads.fasta).
# This FASTA file will be used for subsequent assembly steps.
samtools fasta -F 0x2 -o mtDNA-reads.fasta inopinata-MtDNA-aligned.bam

## Step 6: Assemble the mitochondrial reads using Flye ##
# Use Flye to assemble the mitochondrial reads (mtDNA-reads.fasta) into a mitochondrial genome.
# This step assumes the reads are from PacBio data (--pacbio-raw), and you can adjust the --genome-size to fit the expected size of the mitochondrial genome.
# The output will be stored in the 'mtDNA-flye-1000' directory.
# The --min-overlap 1000 option sets a minimum overlap of 1000 bases for read assembly, and --threads 10 allows Flye to use 10 CPU threads for faster execution.
python /scratch/tn2220/C.inopinata/Flye/bin/flye --pacbio-raw mtDNA-reads.fasta --out-dir mtDNA-flye-1000 --genome-size 14000 --min-overlap 1000 --threads 10


# Step 7: BLAST search for homology between assembled mtDNA and C. elegans mtDNA
module load blast+/2.13.0  # Load BLAST+ for sequence similarity searches

# Perform BLAST search of assembled mtDNA against C. elegans mtDNA
blastn -query mtDNA-flye-1000/assembly-26kb.fasta -subject elegans-MtDNA.fasta -outfmt 10 -perc_identity 80 -evalue 1e-10 -out alignment.txt 
# Extract nucleotide positions 10386-24268 from the assembled contig to generate trusted contig for further analysis
samtools faidx assembly-26kb.fasta contig_1:10386-24268 > C.inopinata.fasta

