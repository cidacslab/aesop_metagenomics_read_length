#!/bin/bash
:<<DOC
Author: Pablo Viana
Created: 2023/06/21

Script used to create a multifasta file from an accession list and individual genome files
DOC

# create alias to echo command to log time at each call
echo() {
    command echo "[$(date +"%Y-%m-%dT%H:%M:%S%z")]: $@"
}
# exit when any command fails
set -e
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command ended with exit code $?." >&2' EXIT

# Function parameters
tsv_file="$1"
genome_folder="$2"
output_file="$3"

count=0  # Variable to track the number of accessions used

# Get the accession numbers from the CSV file (excluding the first line/header), shuffle the list, and save it to a temporary file
awk -F '\t' 'NR>1 {print $1}' "$tsv_file" > accession_file.txt

# Create an empty output file
>"$output_file"

# Iterate over each shuffled accession number
while IFS= read -r accession; do
  # Generate the filename by appending the accession number with the appropriate extension
  filename="${accession}.fasta"

  # Path to the genome file
  genome_file="${genome_folder}/${filename}"

  # Check if the genome file exists
  if [ -f "$genome_file" ]; then
    echo "Copying file to output multifasta: $filename"
    # Remove empty lines and lines with only spaces or tabs from the genome file and copy it to the output multifasta file
    awk 'NF && !/^[[:space:]]*$/' "$genome_file" >> "$output_file"
    count=$((count+1))  # Increment the count
  else
    echo "Genome file '$genome_file' does not exist."
  fi
done < accession_file.txt

# Remove the temporary file
echo "Removing intermediate file: rm accession_file.txt"
rm accession_file.txt

echo "Multifasta file '$output_file' created successfully with '$count' accession numbers"


