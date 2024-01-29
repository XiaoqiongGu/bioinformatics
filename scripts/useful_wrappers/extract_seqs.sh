#!/bin/bash

# Path to your FASTA file
fasta_file="allamgenomesmax10ref.fasta"

# Path to your list of IDs
id_file="allamabumax10ref.top10.txt"

#
output_file="allamgenomes.top10.fasta"

# Read each line in the ID file
while IFS= read -r line
do
    # Use grep to find the line with the ID and the next line (the sequence)
    grep -A 1 -w "$line" "$fasta_file" >> "$output_file"
done < "$id_file"

