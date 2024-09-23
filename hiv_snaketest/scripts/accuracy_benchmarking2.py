import matplotlib.pyplot as plt
import re
import os.path
import seaborn as sns
import pandas as pd
import yaml
import ot
import scipy
import numpy as np
import sys
import pysam

plt.rcParams.update({'font.size': 7.0})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})
plt.rcParams['svg.fonttype'] = 'none'

def parse_usrbintime_output(output):
    """
    Parse the output of 'usrbintime -v' and extract the relevant information into a dictionary.
    """
    result = {}

    # Regular expressions to capture each field
    regex_patterns = {
        "command": re.compile(r"Command being timed:\s+\"(.+?)\""),
        "user_time": re.compile(r"User time \(seconds\):\s+([0-9.]+)"),
        "system_time": re.compile(r"System time \(seconds\):\s+([0-9.]+)"),
        "cpu_percent": re.compile(r"Percent of CPU this job got:\s+([0-9.]+)%"),
        "elapsed_time": re.compile(r"Elapsed \(wall clock\) time \(.*?\):\s+([0-9:.]+)"),
        "max_resident_set_size": re.compile(r"Maximum resident set size \(kbytes\):\s+([0-9]+)"),
        "wall_clock_time": re.compile(r"Elapsed \(wall clock\) time \(.*?\):\s+([0-9:.]+)")
    }

    # Parse each field from the output
    for key, pattern in regex_patterns.items():
        match = pattern.search(output)
        if match:
            result[key] = match.group(1)

    #format wall_clock_time into seconds
    if 'wall_clock_time' in result:
        time = result['wall_clock_time']
        time = time.split(':')
        seconds = 0
        for i in range(len(time)):
            seconds += float(time[i]) * 60**(len(time) - i - 1)
        result['wall_clock_time'] = seconds

    result['cpu_time'] = float(result['user_time']) + float(result['system_time'])

    return result

def distance_hamming_asym(snp_vector1, snp_vector2):
    distance = 0
    for pos in snp_vector1:
        if pos in snp_vector2:
            if snp_vector1[pos] != snp_vector2[pos]:
                distance += 1
    return distance

def af_and_acc(snp_vectors1, truth_vec):
    covered_bases = [set() for _ in range(len(truth_vec))]
    scores = []
    for v1 in snp_vectors1:
        best_ind = -1
        best_score = sys.maxsize
        for i,v2 in enumerate(truth_vec):
            score = distance_hamming_asym(v1, v2)
            if score < best_score:
                best_score = score
                best_ind = i
        for pos in v1:
            covered_bases[best_ind].add(pos)
        scores.append(best_score / len(truth_vec[best_ind]))

    af = [len(cb)/len(truth_vec[0]) for cb in covered_bases]
    af = np.sum(af)/len(af)
    hamming_acc = np.sum(scores)/len(scores)
    return hamming_acc, af


def wasserstein_distance(snp_vec, abundances, true_snp_vector, true_qnames, genomes_to_coverages):

    if len(snp_vec) == 0:
        return 100

    M = []
    w1 = [x/np.sum(abundances) for x in abundances]
    w2 = [genomes_to_coverages[qname] for qname in true_qnames]
    for vec1 in snp_vec:
        list2 = []
        for vec2 in true_snp_vector:
            list2.append(distance_hamming_asym(vec1, vec2))
        M.append(list2)

    
    d = ot.emd2(w1, w2, M)
    return d




def read_vcf(vcf_file):
    """Read VCF file and extract SNP positions and alleles"""
    vcf_reader = pysam.VariantFile(vcf_file)
    snp_data = {}

    for record in vcf_reader:
        # Extract SNP positions and ref/alt alleles
        pos = record.pos
        ref_allele = record.ref
        alt_alleles = record.alts
        snp_data[pos] = {
            'ref': ref_allele,
            'alt': [str(allele) for allele in alt_alleles]
        }

    return snp_data

def generate_snp_vector_for_alignment(alignment, snp_data, chromosome):
    """Generate SNP vector for a single alignment (one genome)."""
    snp_vector = {}

    # Fetch the aligned pairs (reference and read positions)
    aligned_pairs = alignment.get_aligned_pairs(matches_only=True, with_seq=True)

    for read_pos, ref_pos, ref_base in aligned_pairs:
        if ref_pos + 1 in snp_data:
            read_base = alignment.seq[read_pos]
            ref = snp_data[ref_pos+1]['ref']
            alt = snp_data[ref_pos+1]['alt']
            if read_base == ref:
                snp_value = 0  # Reference allele
            elif read_base in alt:
                snp_value = alt.index(read_base) + 1  # Alt allele
            else:
                snp_value = -1  # Unexpected base

            if snp_value != -1:
                snp_vector[ref_pos + 1] = snp_value

    return snp_vector

def generate_snp_vectors_for_genomes(bam_file, vcf_file, chromosome, algo = None):
    """Generate SNP vectors for multiple genomes aligned to a reference."""
    snp_data = read_vcf(vcf_file)
    snp_vectors = []
    qnames = []
    abundances = []

    bam = pysam.AlignmentFile(bam_file, "rb")

    for alignment in bam.fetch(chromosome):
        snp_vector = generate_snp_vector_for_alignment(alignment, snp_data, chromosome)
        snp_vectors.append(snp_vector)
        qnames.append(alignment.query_name)
        query_name = alignment.query_name
        if algo == 'dbghaplo':
            ab = float(query_name.split(',')[-2].split(':')[-1]) / 100
        elif algo == 'rvhaplo':
            ab = float(query_name.split('_')[5])
        elif algo == 'cliqueSNV':
            ab = float(query_name.split('_')[2])
        else:
            ab = None

        abundances.append(ab)

    return (snp_vectors, qnames, abundances)


# Reading a YAML file
with open("config2.yaml", 'r') as file:
    yaml_data = yaml.safe_load(file)

algo = ['rvhaplo', 'cliqueSNV', 'igda', 'dbghaplo']
palette = sns.color_palette("muted")
p = [palette[1], palette[2], palette[3], palette[0]]


results = []
for setup in yaml_data.keys():
    print(setup)
    read_accuracy = 95
    length = 9000
    coverage = 160
    for num_strains in yaml_data[setup]['num_strains']:
        if num_strains > 30:
            break

        setup = 'setup2'
        # Printing the loaded YAML data
        truth_bam = f'results/{setup}/{num_strains}_true.bam'
        truth_vcf = f'results/{setup}/{num_strains}_true.vcf.gz'
        ref_file = "ref.fa"
        #open ref_file get chromosome name
        f = open(ref_file, 'r')
        chromosome = f.readline().strip().split(' ')[0][1:]

        true_snp_vector, true_qnames, _ = generate_snp_vectors_for_genomes(truth_bam, truth_vcf, chromosome)

        for alg in algo:
            genomes_to_coverages = {}
            total_rat = 0
            instance = yaml_data[setup]
            for i in range(num_strains):
                genome = instance['genomes'][i] + '.1'
                ratio = instance['ratio'][i]
                genomes_to_coverages[genome] = ratio
                total_rat += ratio
            for genome in genomes_to_coverages:
                genomes_to_coverages[genome] = genomes_to_coverages[genome] / total_rat

            file = "results/{}/{}_{}_{}_{}_{}/results.bam".format(setup, alg, length, read_accuracy, coverage, num_strains)
            if not os.path.exists(file):
                continue
            time_file = "benchmarks/{}/{}_{}_{}_{}_{}.benchmark.usrbintime".format(setup, alg, length, read_accuracy, coverage, num_strains)
            time_results = parse_usrbintime_output(open(time_file, 'r').read())

            snp_vector, _, abundances = generate_snp_vectors_for_genomes(file, truth_vcf, chromosome, alg)
            hamming_acc, af = af_and_acc(snp_vector, true_snp_vector)
            num_contigs = len(snp_vector)

            if 'igda' in alg:
                d = None
            else:
                d = wasserstein_distance(snp_vector, abundances, true_snp_vector, true_qnames, genomes_to_coverages)
            results.append({'algorithm': alg, 'Wasserstein distance': d, 'Fraction recovered': af, 'Haplotype overestimate': num_contigs - num_strains, 'coverage': coverage, length: length, 'read accuracy': read_accuracy, 'setup': setup, 'Hamming SNP distance': hamming_acc, 'system_time': float(time_results['system_time']), 'Memory (GB)': float(time_results['max_resident_set_size'])/1_000_000, 'user_time': float(time_results['user_time']), 'num_strains': num_strains, 'Wall time (s)': float(time_results['wall_clock_time']), 'CPU time (s)': float(time_results['cpu_time']) })

df = pd.DataFrame(results)

#print as tsv
df.to_csv(f'/home/jshaw/projects/temp/amr/figures/accuracy_benchmarking2_{length}_{read_accuracy}.tsv', sep='\t', index=False)

# plot number of contigs, af, hamming_acc, emd versus coverage for each method (line plot)

ms = 5
fig, axs = plt.subplots(2, 3, figsize=(16/2.54, 8/2.54))
markers = ['s', 'd', 'D' , 'o']
sns.lineplot(data=df, x='num_strains', y='Haplotype overestimate', hue='algorithm', ax=axs[0, 0], marker = 'o', palette=p, markersize=ms, style = 'algorithm', markers = markers)
#axs[0,0].plot([0, 30], [0, 30], 'k--')
sns.lineplot(data=df, x='num_strains', y='Fraction recovered', hue='algorithm', ax=axs[0, 1], marker = 'o', palette=p, markersize=ms, style = 'algorithm', markers = markers)
sns.lineplot(data=df, x='num_strains', y='Hamming SNP distance', hue='algorithm', ax=axs[1, 0], marker='o', palette=p, markersize=ms, style = 'algorithm', markers = markers )
sns.lineplot(data=df, x='num_strains', y='Wasserstein distance', hue='algorithm', ax=axs[1, 1], marker = 'o', palette=p, markersize=ms, style = 'algorithm', markers = markers)
g = sns.lineplot(data=df, x='num_strains', y='Wall time (s)', hue='algorithm', ax=axs[0, 2], marker = 'o', palette=p, markersize=ms, style = 'algorithm', markers = markers)
#set ylim to 10000
g.set(ylim=(0, 3000))
sns.lineplot(data=df, x='num_strains', y='Memory (GB)', hue='algorithm', ax=axs[1, 2], marker = 'o', palette=p, markersize=ms, style = 'algorithm', markers = markers)

#ohrizontal line
axs[0,0].axhline(y=0, color='black', linestyle='--')
axs[0,1].axhline(y=1, color='black', linestyle='--')
axs[1,0].axhline(y=0, color='black', linestyle='--')
axs[1,1].axhline(y=0, color='black', linestyle='--')


for ax in axs.flat:
    #borderless
    sns.despine(ax=ax)
    ax.legend(frameon=False)
    # no legend
    ax.get_legend().remove()
    ax.set_xlabel('Number of HIV strains')

for ax in axs[0, :]:
    ax.set_xlabel('')
    
  

plt.savefig(f'/home/jshaw/projects/temp/amr/figures/accuracy_benchmarking2_{length}_{read_accuracy}.svg')
plt.show()
