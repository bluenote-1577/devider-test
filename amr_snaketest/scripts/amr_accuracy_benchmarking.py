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
first_plot = True
second_plot = True

plt.rcParams.update({'font.size': 7.0})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})
plt.rcParams['svg.fonttype'] = 'none'
sns.set_palette("muted")


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
    }

    # Parse each field from the output
    for key, pattern in regex_patterns.items():
        match = pattern.search(output)
        if match:
            result[key] = match.group(1)

    return result

def distance_hamming_asym(snp_vector1, snp_vector2):
    distance = 0
    for pos in snp_vector1:
        if pos in snp_vector2:
            if snp_vector1[pos] != snp_vector2[pos]:
                distance += 1
    return distance

def af_and_acc(snp_vectors1, truth_vec, true_qnames, qnames_to_coverage, qnames_to_ani, arg):
    covered_bases = [set() for _ in range(len(truth_vec))]
    scores = []
    qnames_to_abundance = {}

    total_cov = sum([c for qname, c in qnames_to_coverage.items()])
    for qname, c in qnames_to_coverage.items():
        qnames_to_abundance[qname] = c / total_cov

    for v1 in snp_vectors1:
        best_ind = -1
        best_score = sys.maxsize
        for i,v2 in enumerate(truth_vec):
            score = distance_hamming_asym(v1, v2)
            if score < best_score:
                best_score = score
                best_ind = i
        for pos in v1:
            if pos in truth_vec[best_ind]:
                covered_bases[best_ind].add(pos)
        scores.append(best_score / len(truth_vec[best_ind]))

    af = [len(cb)/len(truth_vec[0]) for cb in covered_bases]
    af = np.sum(af)/len(af)
    hamming_acc = np.sum(scores)/len(scores)

    afs_by_stats = []
    for i in range(len(truth_vec)):
        qname = true_qnames[i]
        coverage = qnames_to_coverage[qname]
        if (qname, arg) in qnames_to_ani:
            ani = float(qnames_to_ani[(qname, arg)])
        else:
            ani = 100
        af_temp = len(covered_bases[i])/len(truth_vec[i])
        abund = qnames_to_abundance[qname]
        afs_by_stats.append((coverage, ani, af_temp, abund))

    return hamming_acc, af, afs_by_stats 

def wasserstein_distance(snp_vec, true_snp_vector):
    M = []
    for vec1 in snp_vec:
        list2 = []
        for vec2 in true_snp_vector:
            list2.append(distance_hamming_asym(vec1, vec2))
        M.append(list2)
    w1 = np.ones(len(snp_vec)) / len(snp_vec)
    w2 = np.ones(len(true_snp_vector)) / len(true_snp_vector)

    d = ot.emd2(w1, w2, M)
    return d


def parse_skani_results(file):
    results = {}
    with open(file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            gn1 = line[-2]
            gn2 = line[-1]
            ani = line[2]
            results[(gn1, gn2)] = ani
            results[(gn2, gn1)] = ani
    return results



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
    abundances = []

    bam = pysam.AlignmentFile(bam_file, "rb")

    for alignment in bam.fetch(chromosome):
        # get query name
        query_name = alignment.query_name
        snp_vector = generate_snp_vector_for_alignment(alignment, snp_data, chromosome)
        snp_vectors.append((snp_vector, query_name))
        if algo == 'dbghaplo':
            ab = float(query_name.split(',')[-2].split(':')[-1]) / 100
        if algo == 'rvhaplo':
            ab = float(query_name.split('_')[5]) / 100
        if algo == 'cliqueSNV':
            ab = float(query_name.split('_')[2]) / 100
        else:
            ab = None

        abundances.append(ab)

    if algo is None:
        return snp_vectors
    else:
        return snp_vectors, abundances

# Reading a YAML file
with open("auto_config_megares.yaml", 'r') as file:
    yaml_data = yaml.safe_load(file)


qnames_to_ani = parse_skani_results("skani_megares.tsv")
results = []
results_stratified_by_genome = []
algo = ['dbghaplo', 'rvhaplo', 'cliqueSNV', 'igda']
for instance in yaml_data['instances']:
    print(instance)
    arg = instance['arg']
    genomes_to_coverages = {}
    for i in range(len(instance['genomes'])):
        genome = instance['genomes'][i]
        coverage = instance['coverages'][i]
        genomes_to_coverages[genome] = coverage


    # Printing the loaded YAML data
    truth_bam = f'data_megares_v1/{arg}_true.bam'
    truth_vcf = f'data_megares_v1/{arg}_true.vcf.gz'
    ref_file = f"data_megares_v1/{arg}.fasta"
    #open ref_file get chromosome name
    f = open(ref_file, 'r')
    chromosome = f.readline().strip().split(' ')[0][1:]

    true_snp_vector_qname = generate_snp_vectors_for_genomes(truth_bam, truth_vcf, chromosome)
    true_snp_vector = [x[0] for x in true_snp_vector_qname]
    true_qnames = [x[1] for x in true_snp_vector_qname]

    for alg in algo:
        if 'dbghaplo' in alg:
            file = "results/{}_{}/results.bam".format(arg, alg)
        else:
            file = "results_megares_v1/{}_{}/results.bam".format(arg, alg)
        if not os.path.exists(file):
            continue
        time_file = "benchmarks/{}_{}.benchmark.usrbintime".format(arg, alg)
        time_results = parse_usrbintime_output(open(time_file, 'r').read())

        snp_vector_qname = generate_snp_vectors_for_genomes(file, truth_vcf, chromosome)
        #d = wasserstein_distance(snp_vector, true_snp_vector)
        d = 0
        snp_vector = [x[0] for x in snp_vector_qname]
        qnames = [x[1] for x in snp_vector_qname]
        hamming_acc, af, af_by_stats = af_and_acc(snp_vector, true_snp_vector, true_qnames, genomes_to_coverages, qnames_to_ani, arg)
        for coverage, ani, af_temp, abund in af_by_stats:
            results_stratified_by_genome.append({'algorithm': alg, 'Coverage': coverage, 'ANI': ani, 'AF': af_temp, 'Abundance': abund})
        num_contigs_diff = len(snp_vector) - len(true_snp_vector)
        results.append({'algorithm': alg, 'Wasserstein distance': d, 'Fraction recovered': af, 'Haplotype\noverestimation': num_contigs_diff, 'SNP Hamming error': hamming_acc, 'System time': float(time_results['system_time']), 'Memory (GB)': float(time_results['max_resident_set_size'])/1_000_000, 'User time': float(time_results['user_time'])})

df = pd.DataFrame(results)
print(df)

# plot number of contigs, af, hamming_acc, emd versus coverage for each method (line plot)

if first_plot:
    fig, axs = plt.subplots(1, 3, figsize=(16/2.54, 4/2.54))
    #sns.boxplot(x='algorithm', y='num_contigs_diff', data=df, ax=axs[0,0], hue='algorithm')
    #sns.boxplot(x='algorithm', y='hamming_acc', data=df, ax=axs[1,0], hue='algorithm')
    #sns.boxplot(x='algorithm', y='af', data=df, ax=axs[0,1], hue='algorithm')
    #sns.boxplot(x='algorithm', y='emd', data=df, ax=axs[1,1],  hue='algorithm')
    #sns.boxplot(x='algorithm', y='system_time', data=df, ax=axs[0,2], hue='algorithm')
    #sns.boxplot(x='algorithm', y='max_rss', data=df, ax=axs[1,2], hue='algorithm')


    sns.boxplot(x='algorithm', y='Haplotype\noverestimation', data=df, ax=axs[0], hue = 'algorithm' )
    sns.boxplot(x='algorithm', y='SNP Hamming error', data=df, ax=axs[1], hue = 'algorithm' )
    sns.boxplot(x='algorithm', y='Fraction recovered', data=df, ax=axs[2], hue = 'algorithm' )
    #sns.boxplot(x='algorithm', y='Wasserstein distance', data=df, ax=axs[1,1], )
    #sns.boxplot(x='algorithm', y='System time', data=df, ax=axs[0,2], )
    #sns.boxplot(x='algorithm', y='Memory (GB)', data=df, ax=axs[1,2], )

    for ax in axs.flat:
        #borderless
        sns.despine(ax=ax)
        ax.legend(frameon=False)
        ax.set_xlabel('')
        #rotate xticks
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)

    ## add horizontal line to [0,0]
    axs[0].axhline(0, color='black', linestyle='--')

    plt.tight_layout()
    plt.savefig('../figures/amr_nospike_boxes.svg')
    plt.show()

strat_df = pd.DataFrame(results_stratified_by_genome)
print(strat_df['ANI'])
#bin into 20 bins for coverage and ani
strat_df['Coverage'] = pd.cut(strat_df['Coverage'], bins=5)
strat_df['Coverage'] = strat_df['Coverage'].apply(lambda x: x.mid)
strat_df["% identity to ref."] = pd.cut(strat_df['ANI'], bins=10)
strat_df['% identity to ref.'] = strat_df['% identity to ref.'].apply(lambda x: x.mid)
strat_df['Abundance'] = pd.cut(strat_df['Abundance'], bins=5)
strat_df['Abundance'] = strat_df['Abundance'].apply(lambda x: x.mid)
strat_df['Fraction recovered'] = strat_df.groupby(['algorithm', 'Coverage', 'ANI'])['AF'].transform('mean')

# set ANI dtype to float
strat_df['% identity to ref.'] = strat_df['% identity to ref.'].astype(float)
# filter 98
#strat_df = strat_df[strat_df['ANI'] > 98]

print(strat_df)
fig, axs = plt.subplots(1, 3, figsize=(16/2.54, 4/2.54))
eb = 'se'
markers = ['o', 's', '^', 'D']


sns.lineplot(x='Coverage', y='Fraction recovered', data=strat_df, ax=axs[0], hue='algorithm', errorbar=eb, marker = 'o')
sns.lineplot(x='% identity to ref.', y='Fraction recovered', data=strat_df[strat_df['% identity to ref.'] > 97.5], ax=axs[1], hue='algorithm', errorbar=eb, marker = 'o')
sns.lineplot(x='Abundance', y='Fraction recovered', data=strat_df, ax=axs[2], hue='algorithm', errorbar=eb, marker = 'o')

for ax in axs.flat:
    #borderless
    sns.despine(ax=ax)
    #remove legend
    ax.get_legend().remove()

#shared legend up top
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=4, frameon=False)
plt.savefig('../figures/amr_nospike_AF.svg')

plt.tight_layout()
plt.show()
