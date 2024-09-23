import yaml
import subprocess
import os

# Load the YAML file
with open('./auto_config.yaml', 'r') as file:
    data = yaml.safe_load(file)

# Function to extract sequences from ref.fa based on genome IDs
def extract_sequences(genomes, output_fasta):
    with open("data/card.fasta", "r") as ref_file, open(output_fasta, "w") as out_fasta:
        write = False
        for line in ref_file:
            if any(genome.split('|')[1] in line for genome in genomes):
                write = True
                out_fasta.write(line)
                continue
            elif line.startswith(">"):
                write = False
            if write:
                out_fasta.write(line)

# Iterate through each `arg` in the YAML data
for entry in data['instances']:
    arg = entry['arg']
    genomes = entry['genomes']
    output_fasta = f"{arg}-gn.fa"
    
    # Extract the sequences from ref.fa
    extract_sequences(genomes, output_fasta)
    
    # Run minimap2
    output_sam = f"{arg}-gn.sam"
    ref2 = "./dereplicated_card_0.90.fa"  # Replace with your actual reference file
    subprocess.run(["minimap2", "-ax", "map-ont", ref2, output_fasta, "-o", output_sam])
    
    print(f"Completed processing for {arg}. Output saved to {output_sam}")

print("All tasks completed.")
