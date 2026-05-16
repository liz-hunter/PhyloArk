import csv

# Input files
accessions_file = "accessions.txt"
assembly_summary_file = "all_assembly.txt"
output_file = "accession2taxid.map"
missing_output_file = "missing_accessions.txt"

# Load your list of accessions
with open(accessions_file, "r") as f:
    accessions = set(line.strip() for line in f if line.strip())

# Map accession.version -> (base accession, taxid)
accession_to_info = {}

# Parse the assembly_summary file
with open(assembly_summary_file, "r") as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        if row[0].startswith("#") or len(row) < 6:
            continue  # skip comments and incomplete lines
        assembly_accession = row[0]  # ex: GCA_000001405.28
        taxid = row[5]              # taxid column (6th column)

        if assembly_accession in accessions:
            base_accession = assembly_accession.split('.')[0]  # remove version
            accession_to_info[assembly_accession] = (base_accession, taxid)

# Write output in Kraken2 accession2taxid format
with open(output_file, "w") as out:
    for full_acc, (base_acc, taxid) in accession_to_info.items():
        out.write(f"{base_acc}\t{full_acc}\t{taxid}\t0\n")

# Identify and write missing accessions
matched_accessions = set(accession_to_info.keys())
missing_accessions = accessions - matched_accessions

with open(missing_output_file, "w") as f:
    for acc in sorted(missing_accessions):
        f.write(f"{acc}\n")

print(f"Done! Wrote {len(accession_to_info)} entries to {output_file}")
print(f"{len(missing_accessions)} accessions from {accessions_file} were not found in {assembly_summary_file}")
print(f"Missing accessions written to {missing_output_file}")

