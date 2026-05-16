#!/usr/bin/env bash

accessions_file="$1"
filenames_list="$2"
fna_dir="$3"
output_dir="$4"

# Create output dir if it doesn't exist
mkdir -p "$output_dir"

# Step 1: Count how many times each accession appears
# Output: accession<TAB>count
sort "$accessions_file" | uniq -c | awk '{print $2 "\t" $1}' > accession_counts.tsv

# Step 2: Read all filenames into an array
mapfile -t filenames < "$filenames_list"

# Step 3: Match accessions and create symlinks
while IFS=$'\t' read -r accession count; do
    match=""
    for fname in "${filenames[@]}"; do
        # Extract accession from filename
        acc=$(echo "$fname" | grep -oE '^GC[AF]_[0-9]+\.[0-9]+')
        if [[ "$acc" == "$accession" ]]; then
            match="$fname"
            break
        fi
    done

    if [[ -n "$match" && -f "$fna_dir/$match" ]]; then
        for i in $(seq 1 "$count"); do
            ln -s "$(realpath "$fna_dir/$match")" "$output_dir/${accession}_dup${i}.fna"
        done
    else
        echo "Warning: No file found for accession [$accession] (match: [$match])" >&2
    fi
done < accession_counts.tsv