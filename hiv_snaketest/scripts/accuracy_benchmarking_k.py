import matplotlib.pyplot as plt
from  matplotlib.ticker import FuncFormatter
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

def mad(x):
    return np.mean(np.abs(x))

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
        return None

    if len(abundances) == 0:
        return None

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
        if algo == None:
            ab = None
        else: 
            ab = float(query_name.split(',')[-2].split(':')[-1]) / 100

        abundances.append(ab)

    return (snp_vectors, qnames, abundances)

# Reading a YAML file
with open("config.yaml", 'r') as file:
    yaml_data = yaml.safe_load(file)

# Printing the loaded YAML data
truth_bam = 'hivs-genome-alns.bam'
truth_vcf = 'out.vcf.gz'
ref_file = "ref.fa"
#open ref_file get chromosome name
f = open(ref_file, 'r')
chromosome = f.readline().strip().split(' ')[0][1:]

true_snp_vector, true_qnames, _ = generate_snp_vectors_for_genomes(truth_bam, truth_vcf, chromosome)
#algo = ['dbghaplo', 'rvhaplo', 'cliqueSNV', 'igda']
palette = sns.color_palette("muted")
p = [palette[1], palette[2], palette[3], palette[0]]
#algo = ['rvhaplo', 'cliqueSNV', 'igda', 'dbghaplo']
algo = ['dbghaplo']
#colors = ['orange', 'green', 'red', 'blue']
#hue_order = ['igda', 'rvhaplo', 'cliqueSNV', 'dbghaplo']
#order = ['igda', 'rvhaplo', 'cliqueSNV', 'dbghaplo']
#color_order = ['orange', 'green', 'red', 'blue']


results = []
for setup in yaml_data.keys():
    print(setup)
    read_accuracy = 95
    length = 9000
    for coverage in yaml_data[setup]['coverages']:
        for alg in algo:
            for k in range(10,31,3):
                genomes_to_coverages = {}
                total_rat = 0
                instance = yaml_data[setup]
                for i in range(len(instance['ratio'])):
                    genome = instance['genomes'][i] + '.1'
                    ratio = instance['ratio'][i]
                    genomes_to_coverages[genome] = ratio
                    total_rat += ratio
                for genome in genomes_to_coverages:
                    genomes_to_coverages[genome] = genomes_to_coverages[genome] / total_rat

                file = "results-k/setup1/{}_{}_{}_results-{}.bam".format(length, read_accuracy, coverage, k)
                if not os.path.exists(file):
                    print(f'{file} does not exist')
                    continue
                time_file = "benchmarks/{}/{}_{}_{}_{}.benchmark.usrbintime".format(setup, alg, length, read_accuracy, coverage)
                time_results = parse_usrbintime_output(open(time_file, 'r').read())
                if 'dbghaplo' in alg:
                    lofreq_file = "benchmarks/{}/{}_{}_{}.lofreq.benchmark.usrbintime".format(setup, length, read_accuracy, coverage)
                    lofreq_results = parse_usrbintime_output(open(lofreq_file, 'r').read())
                    time_results['user_time'] = float(time_results['user_time']) + float(lofreq_results['user_time'])
                    time_results['system_time'] = float(time_results['system_time']) + float(lofreq_results['system_time'])
                    time_results['wall_clock_time'] = float(time_results['wall_clock_time']) + float(lofreq_results['wall_clock_time'])
                    time_results['max_resident_set_size'] = max(float(time_results['max_resident_set_size']), float(lofreq_results['max_resident_set_size']))

                snp_vector, _, abundances = generate_snp_vectors_for_genomes(file, truth_vcf, chromosome, alg)
                print(file)
                if 'igda' in alg:
                    d = None
                else:
                    d = wasserstein_distance(snp_vector, abundances, true_snp_vector, true_qnames, genomes_to_coverages)
                hamming_acc, af = af_and_acc(snp_vector, true_snp_vector)
                num_contigs = len(snp_vector)
                results.append({'k': k, 'algorithm': 'devider', "Earth mover's distance (average)": d, 'Fraction recovered (average)': af, 'Haplotype error (average)': num_contigs - 7, 'coverage': coverage, length: length, 'read accuracy': read_accuracy, 'setup': setup, 'Hamming SNP error (average)': hamming_acc * 100, 'system_time': float(time_results['system_time']), 'Memory (GB)': float(time_results['max_resident_set_size'])/1_000_000, 'user_time': float(time_results['user_time']), 'Wall time (s)': float(time_results['wall_clock_time']), 'CPU time (s)': float(time_results['cpu_time'])})

df = pd.DataFrame(results)
df.to_csv('../figures/accuracy_benchmark2.tsv', sep='\t')
#print average over fraction recovered
print(df.groupby(['algorithm'])['Fraction recovered (average)'].mean())
print(df.groupby(['algorithm'])['Haplotype error (average)'].mean())
print(df.groupby(['algorithm'])['Hamming SNP error (average)'].mean())
print(df.groupby(['algorithm'])["Earth mover's distance (average)"].mean())

# plot number of contigs, af, hamming_acc, emd versus coverage for each method (line plot)
#small markers
alpha = 1.0
ms= 6
markers = ['s', 'd', 'D' , 'o']

# get mean number of haps for each method and k over coverage


fig, axs = plt.subplots(2, 2, figsize=(12/2.54, 8/2.54))
sns.lineplot(data=df, x='k', y='Haplotype error (average)', ax=axs[0, 0], markers = markers, alpha=alpha, palette=p, style='algorithm' , errorbar=None, estimator=mad)
sns.lineplot(data=df, x='k', y='Fraction recovered (average)', ax=axs[0, 1], markers =markers, alpha = alpha, palette=p, style='algorithm', errorbar=None)
sns.lineplot(data=df, x='k', y='Hamming SNP error (average)', ax=axs[1, 0], markers =markers, alpha=alpha, palette=p, style='algorithm', errorbar=None)
sns.lineplot(data=df, x='k', y="Earth mover's distance (average)",  ax=axs[1, 1], markers =markers, alpha=alpha, palette=p, style='algorithm', errorbar=None)

plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda x, _: int(x)))

for ax in axs.flat:
    #borderless
    sns.despine(ax=ax)
    # no legend
    ax.get_legend().remove()
    ax.set_xlabel('k-mer length')

for ax in axs[0, :]:
    ax.set_xlabel('')

#ohrizontal line
axs[0,0].axhline(y=0, color='black', linestyle='--')
axs[0,0].axhline(y=0.92, color = palette[3], label = 'igda')
axs[0,0].axhline(y=7 - 5.46, color = palette[1], label = 'rvhaplo')
axs[0,0].axhline(y=7 - 5.38, color = palette[2], label = 'cliqueSNV')


axs[0,1].axhline(y=1, color='black', linestyle='--')
axs[0,1].axhline(y=0.86, color = palette[3], label = 'igda')
axs[0,1].axhline(y=0.76, color = palette[1], label = 'rvhaplo')
axs[0,1].axhline(y=0.74, color = palette[2], label = 'cliqueSNV')


axs[1,0].axhline(y=0, color='black', linestyle='--')
axs[1,0].axhline(y=0, color=palette[3],  label = 'igda')
axs[1,0].axhline(y=0.84, color=palette[1], label = 'rvhaplo')
axs[1,0].axhline(y=0.996, color=palette[2], label = 'cliqueSNV')

axs[1,1].axhline(y=0, color='black', linestyle='--')
axs[1,1].axhline(y=5.89, color=palette[1], label = 'rvhaplo')
axs[1,1].axhline(y=5.94, color=palette[2], label = 'cliqueSNV')



#legend at top
handles, labels = axs[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=4, frameon=False)
plt.tight_layout()
plt.savefig(f'/home/jshaw/projects/temp/amr/figures/accuracy_benchmarking_{length}_{read_accuracy}.svg')
    
plt.show()
