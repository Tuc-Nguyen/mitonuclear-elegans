#!/usr/bin/env python

"""Convert fasta files to VCF files.

This script converts a fasta file and a given fasta reference file to
a file in VCF. A reference genome can be provided in fasta format as
input. If no reference is given, the first sequence in the fasta file
will be used as reference.

"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

ver = "1.0"

parser = argparse.ArgumentParser(prog="FastaToVCF.py",
                                 description="convert fasta to VCF; script version " + ver)
parser.add_argument("fastafile",
                    help="path to fasta file")
parser.add_argument("output",
                    help="name of VCF output file")
parser.add_argument("-r", "--reference",
                    help="path to reference genome in fasta format")
args = parser.parse_args()

def read_fasta(file_path):
    with open(file_path, "r") as handle:
        return list(SeqIO.parse(handle, "fasta"))

def get_aa_change(ref_seq, alt_seq, pos, codon_table):
    codon_pos = (pos - 1) % 3
    ref_codon_start = pos - 1 - codon_pos
    alt_codon_start = ref_codon_start
    ref_codon = ref_seq.seq[ref_codon_start:ref_codon_start + 3]
    alt_codon = alt_seq.seq[alt_codon_start:alt_codon_start + 3]
    
    if len(ref_codon) == 3 and len(alt_codon) == 3:
        ref_aa = ref_codon.translate(table=codon_table)
        alt_aa = alt_codon.translate(table=codon_table)
        if ref_aa == alt_aa:
            return "synonymous", f"{ref_aa}->{alt_aa}"
        else:
            return "nonsynonymous", f"{ref_aa}->{alt_aa}"
    return "synonymous", ""

def write_vcf(sequences, ref_seq, output_file, codon_table):
    sample_ids = [record.id for record in sequences]
    
    with open(output_file, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_ids) + "\n")
        
        for pos, ref_base in enumerate(ref_seq.seq, start=1):
            alt_bases = []
            for record in sequences:
                alt_base = record.seq[pos - 1]
                if alt_base != ref_base:
                    alt_bases.append(alt_base)
            
            if alt_bases:
                alt_bases = list(set(alt_bases))  # Get unique alt bases
                alt_str = ",".join(alt_bases)
                syn_status, aa_change = get_aa_change(ref_seq, record, pos, codon_table)
                info_field = f"{syn_status};AA_CHANGE={aa_change}"
                vcf.write(f"MtDNA\t{pos}\t.\t{ref_base}\t{alt_str}\t.\tPASS\t{info_field}\tGT")
                
                for record in sequences:
                    if record.seq[pos - 1] == ref_base:
                        vcf.write("\t0")
                    else:
                        alt_index = alt_bases.index(record.seq[pos - 1]) + 1
                        vcf.write(f"\t{alt_index}")
                
                vcf.write("\n")

# Read the sequences
sequences = read_fasta(args.fastafile)

# Determine the reference sequence
if args.reference is not None:
    ref_sequences = read_fasta(args.reference)
    ref_seq = ref_sequences[0]  # Assuming the first sequence in the reference file is the reference
else:
    ref_seq = sequences[0]  # Use the first sequence in the fasta file as reference
    sequences = sequences[1:]  # Remove the reference sequence from the list

# Define the mitochondrial invertebrate codon table
codon_table = 5  # Table 5 for invertebrate mitochondria

# Write the VCF file
write_vcf(sequences, ref_seq, args.output, codon_table)
