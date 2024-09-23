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
    """
    Translate a codon into an amino acid using the invertebrate mitochondrial genetic code.
    """
    return codon_table.forward_table.get(str(codon).upper(), '*')

def get_codon_substitution_types(codon):
    """
    Determine the number of synonymous and nonsynonymous sites for a given codon.
    """
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
    """
    Count the number of pairwise differences between two codons.
    """
    differences = sum(a != b for a, b in zip(seq1, seq2))
    return differences

def harmonic_number(n):
    """
    Calculate the harmonic number H(n-1) = 1 + 1/2 + 1/3 + ... + 1/(n-1)
    """
    return sum(1.0 / i for i in range(1, n + 1))

def calculate_theta_finite(S, L, num_sequences):
    """
    Calculate Watterson's theta with Jukes-Cantor correction for finite sites.
    S = number of segregating sites, L = sequence length, num_sequences = number of sequences.
    """
    a = harmonic_number(num_sequences)
    sum_i_squared_inverse = sum(1 / i**2 for i in range(1, num_sequences + 1))
    b = (4 * a / 3) - (3.5 * (a**2 - sum_i_squared_inverse)) / (3 * a)
    
    # Apply the Jukes-Cantor correction for theta_w (Watterson's theta), first calculate per-site theta_w, then multiply by the length of the whole alignment to get theta
    theta_finite = L* (S / L) / (a - b * (S / L))
    return theta_finite

def calculate_pi_finite(pi, L):
    """
    Calculate nucleotide diversity (π) with Jukes-Cantor correction
    pi = nucleotide diversity, L = sequence length.
    """
    pi_per_site = pi / L
    if 1 - 4 * pi_per_site / 3 > 0:
        pi_finite = pi_per_site / (1 - 4 * pi_per_site / 3)
        return pi_finite
    else:
        return float('nan')  # Handle cases where correction cannot be applied

def calculate_nucleotide_diversity(pairwise_diffs, total_sites):
    """
    Calculate nucleotide diversity (π) from pairwise differences.
    """
    if total_sites == 0:
        return 0
    return pairwise_diffs / total_sites

def calculate_tajima_d(pi, theta_w, segregating_sites, num_sequences):
    """
    Calculate Tajima's D using the mean pairwise differences (pi), Watterson's theta (theta_w),
    and the number of segregating sites (S).
    """
    if segregating_sites == 0:
        return float('nan')
    
    # Constants for Tajima's D
    a1 = harmonic_number(num_sequences)  # a1 = harmonic number H(n-1)
    a2 = sum(1.0 / (i ** 2) for i in range(1, num_sequences))  # a2 = sum(1 / i^2)

    b1 = (num_sequences + 1) / (3 * (num_sequences - 1))
    b2 = 2 * (num_sequences ** 2 + num_sequences + 3) / (9 * num_sequences * (num_sequences - 1))

    c1 = b1 - 1 / a1
    c2 = b2 - (num_sequences + 2) / (a1 * num_sequences) + a2 / (a1 ** 2)

    e1 = c1 / a1
    e2 = c2 / (a1 ** 2 + a2)

    # Variance of Tajima's D
    variance = math.sqrt(e1 * segregating_sites + e2 * segregating_sites * (segregating_sites - 1))

    # Tajima's D formula
    return (pi - theta_w) / variance if variance > 0 else float('nan')

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

    total_synonymous_segregating_sites = 0
    total_nonsynonymous_segregating_sites = 0

    total_segregating_sites = 0

    pairwise_diffs_syn = 0
    pairwise_diffs_nonsyn = 0
    total_syn_pairs = 0
    total_nonsyn_pairs = 0

    for codon_index in range(0, len(alignment[0].seq), 3):
        codons = [str(record.seq[codon_index:codon_index+3]) for record in alignment]
        
        if any(len(codon) < 3 for codon in codons):
            continue

        # Get unique codons at this position
        unique_codons = set(codons)
        
        # If more than one unique codon is present, it's a segregating site
        if len(unique_codons) > 1:
            total_segregating_sites += 1  # Count this as a segregating site

            # Check whether the difference is synonymous or non-synonymous
            reference_aa = translate_codon(list(unique_codons)[0])
            is_synonymous = True  # Start assuming it's synonymous
            
            for codon in unique_codons:
                if translate_codon(codon) != reference_aa:
                    is_synonymous = False
                    break
            
            # Count the segregating site
            if is_synonymous:
                total_synonymous_segregating_sites += 1
            else:
                total_nonsynonymous_segregating_sites += 1

        # Calculate synonymous and nonsynonymous sites for each unique codon
        syn_sites = 0
        nonsyn_sites = 0
        for codon in unique_codons:
            syn, nonsyn = get_codon_substitution_types(codon)
            syn_sites += syn
            nonsyn_sites += nonsyn

        # Average the sites across unique codons
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
                total_syn_pairs += 1
            else:
                total_nonsynonymous_diffs += count_pairwise_differences(codon1, codon2)
                pairwise_diffs_nonsyn += count_pairwise_differences(codon1, codon2)
                total_nonsyn_pairs += 1

    # Uncorrected values
    pi_s = calculate_nucleotide_diversity(total_synonymous_diffs, total_synonymous_sites)
    pi_n = calculate_nucleotide_diversity(total_nonsynonymous_diffs, total_nonsynonymous_sites)
    theta_w_syn = total_synonymous_segregating_sites / harmonic_number(num_sequences)
    theta_w_nonsyn = total_nonsynonymous_segregating_sites / harmonic_number(num_sequences)
    tajima_d_syn = calculate_tajima_d(pi_s, theta_w_syn, total_synonymous_segregating_sites, num_sequences)
    tajima_d_nonsyn = calculate_tajima_d(pi_n, theta_w_nonsyn, total_nonsynonymous_segregating_sites, num_sequences)

    # Corrected values
    pi_s_finite = calculate_pi_finite(total_synonymous_diffs, L)
    pi_n_finite = calculate_pi_finite(total_nonsynonymous_diffs, L)
    theta_w_syn_finite = calculate_theta_finite(total_synonymous_segregating_sites, L, num_sequences)
    theta_w_nonsyn_finite = calculate_theta_finite(total_nonsynonymous_segregating_sites, L, num_sequences)
    tajima_d_syn_finite = calculate_tajima_d(pi_s_finite, theta_w_syn_finite, total_synonymous_segregating_sites, num_sequences)
    tajima_d_nonsyn_finite = calculate_tajima_d(pi_n_finite, theta_w_nonsyn_finite, total_nonsynonymous_segregating_sites, num_sequences)

    print(f"Finished processing {gene_name}")
    return gene_name, total_synonymous_sites, total_nonsynonymous_sites, total_synonymous_diffs, total_nonsynonymous_diffs, pi_s, pi_n, total_synonymous_segregating_sites, total_nonsynonymous_segregating_sites, theta_w_syn, tajima_d_syn, theta_w_nonsyn, tajima_d_nonsyn, pi_s_finite, pi_n_finite, theta_w_syn_finite, tajima_d_syn_finite, theta_w_nonsyn_finite, tajima_d_nonsyn_finite

def main(alignment_file):
    output_file = f"{os.path.basename(alignment_file).split('.')[0]}_summary_output.txt"
    result = process_alignment(alignment_file)

    if result is None:
        print(f"No valid results for {alignment_file}")
        return

    with open(output_file, "w") as f_out:
        f_out.write("Gene\tNumber Synonymous Sites\tNumber NonSyn Sites\tTotal Synonymous Diffs\tTotal NonSyn Diffs\tpiS\tpiN\tSynonymous Segregating Sites\tNonSynonymous Segregating Sites\tTheta_w_Syn\tTajima_D_Syn\tTheta_w_NonSyn\tTajima_D_NonSyn\tpiS_Finite\tpiN_Finite\tTheta_w_Syn_Finite\tTajima_D_Syn_Finite\tTheta_w_NonSyn_Finite\tTajima_D_NonSyn_Finite\n")
        f_out.write("\t".join(map(str, result)) + "\n")

    print(f"Results written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python piNPiS_summary.py <alignment_file>")
        sys.exit(1)
    alignment_file = sys.argv[1]
    main(alignment_file)
    
    
