#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=inopinata
#SBATCH --mail-type=END
#SBATCH --mail-user=tn2220@nyu.edu
#SBATCH --output=C.inopinata_%j.out

module purge 
module load minimap2/2.22
module load bedtools/intel/2.29.2
module load samtools/intel/1.14

#minimap2 -d REF.mmi PRJDB5687_WS289.fa  
#minimap2 -a REF.mmi reads.fastq > alignment.sam 

samtools view -bS -o aligned.bam alignment.sam
samtools view -b -f 4 -o unmapped.bam aligned.bam
bedtools bamtofastq -i unmapped.bam -fq unmapped.fastq
