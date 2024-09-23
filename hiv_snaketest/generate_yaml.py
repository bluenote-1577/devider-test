import random
import yaml

# Function to read the input file and process the lines
def read_clusters(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return [line.strip().split() for line in lines]

# Function to select 2-10 random entries and coverages
def select_random_entries(entries):
    if len(entries) < 10:
        return
    num_entries = random.randint(2, 10)
    selected_entries = random.sample(entries, num_entries)
    coverages = [random.randint(50, 1000) for _ in range(num_entries)]
    return selected_entries, coverages

# Function to create the YAML structure
def create_yaml_structure(lines):
    instances = []
    for i, entries in enumerate(lines):
        ret = select_random_entries(entries)
        if not ret:
            continue
        selected_entries, coverages  = ret
        argname = selected_entries[0].split('|')[-1]
        instance = {
            'arg': argname,
            'genomes': selected_entries,
            'coverages': coverages
            'readlength': readlengths
        }
        instances.append(instance)
    return {'instances': instances}

# Function to write the YAML structure to a file
def write_yaml(file_path, data):
    with open(file_path, 'w') as file:
        yaml.dump(data, file, default_flow_style=False)

# Main function to read input, process data, and write output
def main():
    output_file = 'auto_config.yaml'
    
    lines = read_clusters(clusters_file)
    yaml_structure = create_yaml_structure(lines)
    write_yaml(output_file, yaml_structure)

if __name__ == "__main__":
    main()

