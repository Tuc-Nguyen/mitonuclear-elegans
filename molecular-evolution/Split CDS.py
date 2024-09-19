from Bio import SeqIO

# Input files
input_fasta = "mtdna.fa"
genes_file = "genes.txt"

# Function to extract sequence based on start and end positions
def extract_sequence(start, end, record):
    return record.seq[start-1:end+3]

# Read the genes file
with open(genes_file) as genes:
    for line in genes:
        start, end, gene_name = line.strip().split()
        start, end = int(start), int(end)
        
        output_file = f"{gene_name}.fasta"
        
        with open(output_file, "w") as outfile:
            # Read the input FASTA file
            for record in SeqIO.parse(input_fasta, "fasta"):
                extracted_seq = extract_sequence(start, end, record)
                record.seq = extracted_seq
                SeqIO.write(record, outfile, "fasta")

print("Sequences have been extracted and saved to individual files.")
