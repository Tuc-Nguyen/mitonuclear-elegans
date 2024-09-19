
##generating fasta files from vcf file. 
bcftools query -l isotype-mtDNA.vcf.gz | while read sample; do
    echo $sample
    bcftools consensus -f WS283.MtDNA.fasta -s "$sample" isotype-mtDNA.vcf.gz | awk -v header=">$sample" 'NR==1{print header; next}1' > FASTA/"$sample.fasta"
done

##combining all fasta files into one and align them using MAFFT to generate a .phy - the format that is compatible with RAxML
cat FASTAs/*.fasta > mtdna.fa
mafft --auto mtDNA.fa > aligned.phy

##generate ML phylogenetic tree with RAxML
# This will identify identical genomes and exclude them and output reduced MSA in a different .phy file
raxmlHPC -m GTRGAMMA -p 12345 -s aligned.phy -# 5 -n RAxML-TEST

# Generate random seed numbers
main_seed=$(shuf -i 1-99999 -n 1)
bootstrap_seed=$(shuf -i 1-99999 -n 1)
# Run RAxML command with random seed numbers on reduced MSA 
raxmlHPC -m GTRGAMMA -p $main_seed -# 1000 -s aligned.phy.reduced -n 1000ReducedBootstrap -T 20 -f a -x $bootstrap_seed
