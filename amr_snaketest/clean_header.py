import re
from Bio import SeqIO

def sanitize_header1(header):
    # Define the regex pattern
    pattern = re.compile(r'^[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*$')
    
    # If the header matches the pattern, return it as is
    if pattern.match(header):
        return header

    # Replace invalid characters with a temporary character, e.g., '_'
    temp_header = re.sub(r'[^0-9A-Za-z!#$%&+./:;?@^_|~-]', 'xx', header)
    
    # Ensure the first character is valid
    if not re.match(r'^[0-9A-Za-z!#$%&+./:;?@^_|~-]', temp_header[0]):
        temp_header = '_' + temp_header[1:]
    
    return temp_header



def sanitize_header2(header):
    #replace () with -
        
    #return header.split("|")[-1]
    head = header
    head = head.replace("(", "xx")
    head = head.replace(")", "xx")
    head = head.replace('-', "xx")
    head = head.replace("|", "xx")
    return head

def rename_fasta_headers(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            record.id = sanitize_header1(record.id)
            record.id = sanitize_header2(record.id)
            record.description = ""
            SeqIO.write(record, outfile, "fasta")

input_file = "../megares_database_v3.00.fasta"
output_file = "./clean_megares.fasta"

rename_fasta_headers(input_file, output_file)

print(f"Headers have been sanitized and written to {output_file}")

