import pysam
import sys

import numpy as np
from scipy.spatial.distance import pdist, squareform
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
import io


# Define the input BAM and reference genome files
bam_file = sys.argv[1]
reference_fasta = sys.argv[2]
output_vcf = sys.argv[3]

# Open the BAM and reference genome files using pysam
bam = pysam.AlignmentFile(bam_file, "rb")
reference = pysam.FastaFile(reference_fasta)

# Initialize a dictionary to hold SNPs
snps = {}

# Iterate over all reads in the BAM file
for read in bam.fetch():
    # Get the aligned positions and corresponding reference positions
    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if query_pos is None or ref_pos is None:
            continue

        # Get the base from the read and the reference genome
        read_base = read.query_sequence[query_pos]
        ref_base = reference.fetch(read.reference_name, ref_pos, ref_pos+1)

        # Check if there is a SNP (mismatch between read and reference)
        if read_base != ref_base:
            if ref_pos not in snps:
                # Store the first occurrence of a SNP
                snps[ref_pos] = {
                    'ref_base': ref_base,
                    'alt_bases': [read_base],
                    'read_count': 1
                }
            else:
                # Add new alternative allele to the existing SNP
                if read_base not in snps[ref_pos]['alt_bases']:
                    snps[ref_pos]['alt_bases'].append(read_base)
                snps[ref_pos]['read_count'] += 1

# iterate over all reads, tag by 0 or 1 vector based on if reference or alternate at all snp positions
vectors = []
for read in bam.fetch():
    # Initialize a list to store the genotype of the read
    genotype = []
    # Get the aligned positions and corresponding reference positions
    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if ref_pos in snps:
            if query_pos is None or ref_pos is None:
                continue
        # Check if the reference position is a SNP
            # Get the base from the read and the reference genome
            read_base = read.query_sequence[query_pos]
            ref_base = reference.fetch(read.reference_name, ref_pos, ref_pos+1)
            # Check if the read base is the reference base
            if read_base == snps[ref_pos]['ref_base']:
                genotype.append(0)
            elif read_base in snps[ref_pos]['alt_bases']:
                for i, alt_base in enumerate(snps[ref_pos]['alt_bases']):
                    if read_base == alt_base:
                        genotype.append(i+1)
                        break
    vectors.append((str(read.query_name), genotype))

for ( geno, vec) in vectors:
    print(''.join(str(x) for x in vec) , geno)

# Write the SNPs to a VCF file
with open(output_vcf, 'w') as vcf_out:
    # Write VCF header
    vcf_out.write("##fileformat=VCFv4.2\n")
    vcf_out.write(f"##reference={reference_fasta}\n")
    vcf_out.write(f'##FILTER=<ID=PASS,Description="All filters passed">\n')
    vcf_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    # Write each SNP to the VCF
    for pos, snp in sorted(snps.items()):
        chrom = bam.get_reference_name(read.reference_id)
        ref_base = snp['ref_base']
        alt_bases = ','.join(snp['alt_bases'])
        vcf_out.write(f"{chrom}\t{pos+1}\t.\t{ref_base}\t{alt_bases}\t.\t.\t.\n")

# Close files
bam.close()
reference.close()

