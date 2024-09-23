import yaml
import random

def load_top_genomes(tsv_file, top_n=50):
    """Load the top N genomes from the TSV file."""
    genomes = []
    
    with open(tsv_file, 'r') as f:
        # Skip the header
        next(f)
        
        for line in f:
            if len(genomes) >= top_n:
                break
            # Split the line by tab and take the 6th column (the Query_name column)
            query_name = line.split('\t')[0].split()[0]  # Extract the genome ID (before any spaces)
            genomes.append(query_name)
    
    return genomes

def generate_yaml(config_file, genomes):
    """Generate a YAML configuration file with the top genomes."""
    config_data = {
        'setup1': {
            'ref_genome': 'sfu_hivs/OR483991.1.fasta',
            'input_folder': 'sfu_hivs',
            'genomes': genomes,
            'ratio': [round(random.uniform(0.1, 0.9), 2) for _ in range(len(genomes))],  # Random ratio between 0.1 and 0.9
            'coverages': [10, 20, 40, 80, 160, 320],  # Cyclic pattern for coverages
            'lengths': [3000, 9000],
            'read_identities': [98, 95]
        }
    }
    
    # Truncate the coverages list to match the number of genomes
    config_data['setup1']['coverages'] = config_data['setup1']['coverages'][:len(genomes)]
    
    with open(config_file, 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False)

# Example usage
tsv_file = '483991-vs-all.tsv'
config_file = 'config2.yaml'

# Load the top 50 genomes
genomes = load_top_genomes(tsv_file, top_n=75)

# Generate the YAML configuration
generate_yaml(config_file, genomes)

print(f"YAML configuration written to {config_file}")

