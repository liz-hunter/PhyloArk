import os
import sys
import re

def load_accession_taxid_map(map_file):
    accession_taxid = {}
    with open(map_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                full_accession = parts[1]  # e.g., GCA_036940625.1
                taxid = parts[2]
                accession_taxid[full_accession] = taxid
    return accession_taxid

def modify_fasta_headers(fasta_file, accession_taxid_map):
    # Extract accession from filename using regex
    filename = os.path.basename(fasta_file)
    match = re.match(r'^(GCA|GCF)_\d+\.\d+', filename)
    if not match:
        print(f"Could not extract accession from filename: {filename}")
        return
    accession = match.group(0)
    taxid = accession_taxid_map.get(accession)
    if not taxid:
        print(f"TaxID not found for accession {accession} in {fasta_file}")
        return

    temp_file = fasta_file + '.tmp'
    with open(fasta_file, 'r') as infile, open(temp_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                header = line[1:].strip()
                new_header = f">{accession}|kraken:taxid|{taxid} {header}\n"
                outfile.write(new_header)
            else:
                outfile.write(line)
    os.replace(temp_file, fasta_file)

def main(fasta_dir, map_file):
    accession_taxid_map = load_accession_taxid_map(map_file)
    for root, dirs, files in os.walk(fasta_dir):
        for file in files:
            if file.endswith('.fna'):
                fasta_path = os.path.join(root, file)
                modify_fasta_headers(fasta_path, accession_taxid_map)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python modify_fasta_headers.py <fasta_directory> <accession2taxid.map>")
        sys.exit(1)
    fasta_directory = sys.argv[1]
    accession2taxid_map = sys.argv[2]
    main(fasta_directory, accession2taxid_map)
