import pysam
import sys
import scipy.stats as stats
from collections import Counter

def analyze_bam_file(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    results = []

    for pileupcolumn in bam.pileup():
        position = pileupcolumn.pos
        base_counts = Counter()
        
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                base_counts[base] += 1
        
        print(base_counts)
        if len(base_counts) < 2 or sum(base_counts.values()) < 10:
            continue

        max_allele, max_count = base_counts.most_common(1)[0]
        second_max_allele, second_max_count = base_counts.most_common(2)[1]

       # p_value = stats.binom_test(second_max_count, max_count, 0.05, alternative='greater')
        p_value = 10
        print(p_value, position, max_allele, max_count, second_max_allele, second_max_count)
        
        if p_value < 0.05:
            result = {
                'position': position,
                'max_allele': max_allele,
                'max_count': max_count,
                'second_max_allele': second_max_allele,
                'second_max_count': second_max_count,
                'p_value': p_value
            }
            results.append(result)

    bam.close()
    return results

def main():
    bam_file = sys.argv[1]
    results = analyze_bam_file(bam_file)
    
    with open("results.txt", "w") as f:
        for result in results:
            f.write(f"Position: {result['position']}, Max Allele: {result['max_allele']} ({result['max_count']}), "
                    f"Second Max Allele: {result['second_max_allele']} ({result['second_max_count']}), "
                    f"P-value: {result['p_value']}\n")

if __name__ == "__main__":
    main()
