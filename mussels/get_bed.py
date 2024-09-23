import pysam

def find_snps_in_region(vcf_file, chrom, start, end):
    """Count SNPs in a region by parsing the VCF file directly."""
    snp_count = 0
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            vcf_chrom = parts[0]
            vcf_pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]

            if vcf_chrom == chrom and start <= vcf_pos < end:
                # Check if it's a SNP (single nucleotide polymorphism)
                if len(ref) == 1 and len(alt) == 1 and ref != alt:
                    snp_count += 1
    return snp_count


def calculate_average_depth(bam, chrom, start, end):
    """Calculate the average depth in a region from BAM file."""
    total_depth = 0
    positions_covered = 0
    for pileupcolumn in bam.pileup(chrom, start, end, truncate=True, min_mapping_quality = 15, stepper = 'all'):
        total_depth += pileupcolumn.n
        positions_covered += 1

    if positions_covered == 0:
        return 0
    return total_depth / positions_covered

def main(bam_file, vcf_file, output_bed):
    # Open BAM file for fetching references and coverage information
    bam = pysam.AlignmentFile(bam_file, "rb")

    with open(output_bed, 'w') as bed_out:
        for chrom in bam.references:
            chrom_length = bam.get_reference_length(chrom)
            if chrom_length< 8000:
                continue
            # Start from 0, increment by 10kb for each step
            start = 0
            while start + 3000 <= chrom_length:
                end = start + 3000

                # Check if the region has > 10 SNPs
                snp_count = find_snps_in_region(vcf_file, chrom, start, end)
                
                # Check if the coverage is between 170 and 230
                print(snp_count, start)
                if snp_count > 10:
                    avg_coverage = calculate_average_depth(bam, chrom, start, end)
                    print(avg_coverage)
                    if 200 < avg_coverage < 250:
                        # Write region to BED file
                        bed_out.write(f"{chrom}\t{start}\t{end}\n")
                
                # Move 10kb apart for the next region
                start += 5000

if __name__ == "__main__":
    import sys
    bam_file = './reads_on_mussel.bam'
    vcf_file = './lofreq.vcf'
    output_bed = './out.bed'
    main(bam_file, vcf_file, output_bed)

