import sys
import os
import warnings
import math
from Bio import SeqIO
from Bio.Data import CodonTable

# Suppress specific Biopython warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

# Get the invertebrate mitochondrial codon table
codon_table = CodonTable.unambiguous_dna_by_id[5]

def translate_codon(codon):
    """ Translate a codon into an amino acid using the invertebrate mitochondrial genetic code. """
    return codon_table.forward_table.get(str(codon).upper(), '*')

def get_codon_substitution_types(codon):
    """ Determine the number of synonymous and nonsynonymous sites for a given codon. """
    syn_sites = 0
    nonsyn_sites = 0
    codon = str(codon).upper()
    original_aa = translate_codon(codon)
    
    for i in range(3):  # Iterate over each nucleotide position in the codon
        for nt in "ACGT":  # Check all possible nucleotide substitutions
            if codon[i] != nt:
                mutated_codon = codon[:i] + nt + codon[i+1:]
                mutated_aa = translate_codon(mutated_codon)
                if mutated_aa == original_aa:
                    syn_sites += 1/3  # Increment synonymous sites
                else:
                    nonsyn_sites += 1/3  # Increment nonsynonymous sites
                    
    return syn_sites, nonsyn_sites

def count_pairwise_differences(seq1, seq2):
    """ Count the number of pairwise differences between two codons. """
    differences = sum(a != b for a, b in zip(seq1, seq2))
    return differences

def calculate_nucleotide_diversity(pairwise_diffs, total_sites):
    """ Calculate nucleotide diversity (π) from pairwise differences. """
    if total_sites == 0:
        return 0
    return pairwise_diffs / total_sites

def calculate_pi_finite(pi, L):
    """ Calculate nucleotide diversity (π) with Jukes-Cantor correction. """
    pi_per_site = pi / L
    if 1 - 4 * pi_per_site / 3 > 0:
        pi_finite = pi_per_site / (1 - 4 * pi_per_site / 3)
        return pi_finite
    else:
        return float('nan')  # Handle cases where correction cannot be applied

def process_alignment(alignment_file):
    print(f"Processing alignment: {alignment_file}")
    gene_name = os.path.basename(alignment_file).split(".")[0]
    alignment = list(SeqIO.parse(alignment_file, "fasta"))
    num_sequences = len(alignment)
    L = len(alignment[0].seq)  # Assume all sequences are the same length
    
    print(f"Number of sequences: {num_sequences}")

    if num_sequences == 0:
        print(f"No sequences found in {alignment_file}")
        return None

    total_synonymous_sites = 0
    total_nonsynonymous_sites = 0
    total_synonymous_diffs = 0
    total_nonsynonymous_diffs = 0

    pairwise_diffs_syn = 0
    pairwise_diffs_nonsyn = 0

    for codon_index in range(0, len(alignment[0].seq), 3):
        codons = [str(record.seq[codon_index:codon_index+3]) for record in alignment]
        
        if any(len(codon) < 3 for codon in codons):
            continue

        unique_codons = set(codons)
        
        # Calculate synonymous and nonsynonymous sites for each unique codon
        syn_sites = 0
        nonsyn_sites = 0
        for codon in unique_codons:
            syn, nonsyn = get_codon_substitution_types(codon)
            syn_sites += syn
            nonsyn_sites += nonsyn

        total_synonymous_sites += syn_sites / len(unique_codons)
        total_nonsynonymous_sites += nonsyn_sites / len(unique_codons)

        # Collect unique pairs of codons for difference calculation
        unique_diffs = set()
        for i in range(num_sequences):
            for j in range(i + 1, num_sequences):
                if codons[i] != codons[j]:
                    unique_diffs.add((codons[i], codons[j]))

        # Calculate synonymous and nonsynonymous differences based on codon pairs
        for codon1, codon2 in unique_diffs:
            amino_acid1 = translate_codon(codon1)
            amino_acid2 = translate_codon(codon2)
            
            if amino_acid1 == amino_acid2:
                total_synonymous_diffs += count_pairwise_differences(codon1, codon2)
                pairwise_diffs_syn += count_pairwise_differences(codon1, codon2)
            else:
                total_nonsynonymous_diffs += count_pairwise_differences(codon1, codon2)
                pairwise_diffs_nonsyn += count_pairwise_differences(codon1, codon2)

    pi_s = calculate_nucleotide_diversity(total_synonymous_diffs, total_synonymous_sites)
    pi_n = calculate_nucleotide_diversity(total_nonsynonymous_diffs, total_nonsynonymous_sites)
    pi_s_finite = calculate_pi_finite(pi_s, L)
    pi_n_finite = calculate_pi_finite(pi_n, L)

    return gene_name, total_synonymous_sites, total_nonsynonymous_sites, total_synonymous_diffs, total_nonsynonymous_diffs, pi_s, pi_n, pi_s_finite, pi_n_finite

def main(alignment_file):
    output_file = f"{os.path.basename(alignment_file).split('.')[0]}_summary_output.txt"
    result = process_alignment(alignment_file)

    if result is None:
        print(f"No valid results for {alignment_file}")
        return

    # Write results to file
    with open(output_file, "w") as f_out:
        f_out.write("Gene\tNumber Synonymous Sites\tNumber NonSyn Sites\tTotal Synonymous Diffs\tTotal NonSyn Diffs\tpiS\tpiN\tpiS_Finite\tpiN_Finite\n")
        f_out.write("\t".join(map(str, result)) + "\n")

    print(f"Results written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python extract_data.py <alignment_file>")
        sys.exit(1)
    alignment_file = sys.argv[1]
    main(alignment_file)
