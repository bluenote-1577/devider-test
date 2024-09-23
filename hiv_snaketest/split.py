from Bio import SeqIO
import sys
import os

# Function to split a FASTA file by each record
def split_fasta_by_record(input_fasta, output_dir):
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Parse the input FASTA file and write each record to a separate file
    for record in SeqIO.parse(input_fasta, "fasta"):
        # Create a filename based on the sequence ID
        output_filename = f"{record.id}.fasta"
        output_path = os.path.join(output_dir, output_filename)

        # Write the record to a new FASTA file
        with open(output_path, "w") as output_file:
            SeqIO.write(record, output_file, "fasta")

        print(f"Written {output_filename}")

# Example usage
input_fasta = sys.argv[1]
output_dir = sys.argv[2]
split_fasta_by_record(input_fasta, output_dir)

