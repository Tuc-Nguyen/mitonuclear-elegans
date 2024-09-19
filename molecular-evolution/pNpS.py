import sys
import os
import warnings
from Bio import SeqIO
from Bio.Data import CodonTable

# Suppress specific Biopython warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

# Get the invertebrate mitochondrial codon table
codon_table = CodonTable.unambiguous_dna_by_id[5]

def translate_codon(codon):
    """
    Translate a codon into an amino acid using the invertebrate mitochondrial genetic code.
    """
    return codon_table.forward_table.get(str(codon).upper(), '*')

def get_codon_substitution_types(codon):
    """
    Determine the number of synonymous and nonsynonymous sites for a given codon.
    """
    syn_sites = 0  # Initialize the count of synonymous sites
    nonsyn_sites = 0  # Initialize the count of nonsynonymous sites
    codon = str(codon).upper()  # Convert the codon to uppercase to ensure uniformity
    original_aa = translate_codon(codon)  # Translate the original codon to its amino acid

    for i in range(3):  # Loop through each position in the codon
        syn_count = 0  # Count of synonymous changes at this position
        nonsyn_count = 0  # Count of nonsynonymous changes at this position
        for nt in "ACGT":  # Loop through each possible nucleotide
            if codon[i] != nt:  # Only consider changes
                if i == 0:
                    mutated_codon = nt + codon[1:]  # Create the mutated codon for the first position
                elif i == 1:
                    mutated_codon = codon[0] + nt + codon[2:]  # Create the mutated codon for the second position
                elif i == 2:
                    mutated_codon = codon[:2] + nt  # Create the mutated codon for the third position

                mutated_aa = translate_codon(mutated_codon)  # Translate the mutated codon
                if mutated_aa == original_aa:  # Check if the mutation is synonymous
                    syn_count += 1  # Increment synonymous count
                else:
                    nonsyn_count += 1  # Increment nonsynonymous count

        syn_sites += syn_count / 3  # Average synonymous sites at this position
        nonsyn_sites += nonsyn_count / 3  # Average nonsynonymous sites at this position

    return syn_sites, nonsyn_sites  # Return the total synonymous and nonsynonymous sites

def calculate_nucleotide_diversity(pairwise_diffs, total_sites):
    """
    Calculate nucleotide diversity (π) from pairwise differences.
    """
    if total_sites == 0:
        return 0  # Avoid division by zero
    return pairwise_diffs / total_sites  # Calculate nucleotide diversity

def process_alignment(alignment_file):
    """
    Process the alignment file to calculate summary statistics.
    """
    gene_name = os.path.basename(alignment_file).split(".")[0]  # Extract gene name from file name
    alignment = list(SeqIO.parse(alignment_file, "fasta"))  # Parse the alignment file in FASTA format
    num_sequences = len(alignment)  # Get the number of sequences in the alignment
    
    seq_synonymous_sites = [0] * num_sequences  # Initialize a list to store synonymous sites count for each sequence
    seq_nonsynonymous_sites = [0] * num_sequences  # Initialize a list to store nonsynonymous sites count for each sequence
    total_synonymous_diffs = 0  # Initialize the total count of synonymous differences
    total_nonsynonymous_diffs = 0  # Initialize the total count of nonsynonymous differences

    for codon_index in range(0, len(alignment[0].seq), 3):  # Loop through codons in steps of 3 nucleotides
        codons = [str(record.seq[codon_index:codon_index+3]) for record in alignment]  # Extract codons at the current index
        
        if any(len(codon) < 3 for codon in codons):  # Skip incomplete codons
            continue

        for seq_index, codon in enumerate(codons):  # Loop through each codon
            syn, nonsyn = get_codon_substitution_types(codon)  # Get synonymous and nonsynonymous sites for the codon
            seq_synonymous_sites[seq_index] += syn  # Accumulate synonymous sites for the sequence
            seq_nonsynonymous_sites[seq_index] += nonsyn  # Accumulate nonsynonymous sites for the sequence

        # Collect unique differences between sequences
        unique_diffs = set()
        for i in range(num_sequences):
            for j in range(i + 1, num_sequences):
                if codons[i] != codons[j]:  # Check if codons are different
                    unique_diffs.add((codons[i], codons[j]))  # Add unique differences

        for codon1, codon2 in unique_diffs:  # Loop through unique differences
            amino_acid1 = translate_codon(codon1)
            amino_acid2 = translate_codon(codon2)
            
            if amino_acid1 == amino_acid2:
                total_synonymous_diffs += count_pairwise_differences(codon1, codon2)  # Count synonymous differences
            else:
                total_nonsynonymous_diffs += count_pairwise_differences(codon1, codon2)  # Count nonsynonymous differences

    avg_synonymous_sites = sum(seq_synonymous_sites) / num_sequences  # Calculate average synonymous sites per sequence
    avg_nonsynonymous_sites = sum(seq_nonsynonymous_sites) / num_sequences  # Calculate average nonsynonymous sites per sequence

    pi_s = calculate_nucleotide_diversity(total_synonymous_diffs, avg_synonymous_sites)  # Calculate nucleotide diversity for synonymous sites
    pi_n = calculate_nucleotide_diversity(total_nonsynonymous_diffs, avg_nonsynonymous_sites)  # Calculate nucleotide diversity for nonsynonymous sites
    
    # Export the count of synonymous and nonsynonymous sites for each sequence to a text file
    sites_output_file = f"{gene_name}_sites_counts.txt"
    with open(sites_output_file, "w") as f_out:
        f_out.write("Sequence_ID\tSynonymous_Sites\tNonsynonymous_Sites\n")
        for i, record in enumerate(alignment):
            f_out.write(f"{record.id}\t{seq_synonymous_sites[i]}\t{seq_nonsynonymous_sites[i]}\n")

    return gene_name, avg_synonymous_sites, avg_nonsynonymous_sites, total_synonymous_diffs, total_nonsynonymous_diffs, pi_s, pi_n
    
def count_pairwise_differences(seq1, seq2):
    """
    Count the number of pairwise differences between two sequences.
    """
    differences = sum(a != b for a, b in zip(seq1, seq2))  # Count differences between corresponding positions
    return differences

def main(alignment_file):
    """
    Main function to process the alignment file and save summary statistics.
    """
    output_file = f"{os.path.basename(alignment_file).split('.')[0]}_summary_statistics.txt"  # Define the output file name
    result = process_alignment(alignment_file)  # Process the alignment file

    with open(output_file, "w") as f_out:  # Write summary statistics to the output file
        f_out.write("Gene\tNumber Synonymous Sites\tNumber NonSyn Sites\tTotal Synonymous Diffs\tTotal NonSyn Diffs\tpiS\tpiN\n")
        f_out.write("\t".join(map(str, result)) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 2:  # Ensure the correct number of arguments is provided
        print("Usage: python piNPiS_summary.py <alignment_file>")  # Print usage information
        sys.exit(1)  # Exit with error code 1
    alignment_file = sys.argv[1]  # Get the alignment file from the command-line arguments
    main(alignment_file)  # Run the main function
