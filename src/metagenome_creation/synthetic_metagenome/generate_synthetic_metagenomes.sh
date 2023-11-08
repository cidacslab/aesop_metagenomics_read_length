#!/bin/bash
:<<DOC
Author: Pablo Viana
Created: 2023/06/21

Script used to create metagenomes for the designated taxonomic composition using different read size models.
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

# Absolute path of project files
project_path="/mnt/c/Users/pablo/Documents/github/taxonomy_analysis"
# Path to the CSV files containing the accession numbers composition
accession_path="${project_path}/results/new_mocks/mocks_composition"
# Path to the folder containing the genome files
genome_path="${project_path}/results/genomes"
# Output multifasta file
output_path="${project_path}/results/new_mocks"
multigenome_file="genomes_output.fasta"

# creates output path if it doesnt exist
mkdir -p $output_path

# Source the script containing the create_multifasta function
create_multifasta_script="${project_path}/src/synthetic_metagenome/generate_multifasta_genomes_from_accessions.sh"

# Command to run the InSilicoSeq synthetic metagenome generator
iss_command="docker run -v ${output_path}:/mnt/data hadrieng/insilicoseq iss generate --cpus 4 --genomes /mnt/data/${multigenome_file}"

# Iterate over all files in the accession folder and its subfolders
for tsv_file in $(find "$accession_path" -type f -name "*.tsv"); do
  echo "Creating fasta for dataset: $tsv_file"  
  # Create the multifasta genome file
  echo "$create_multifasta_script $tsv_file $genome_path ${output_path}/$multigenome_file"
  $create_multifasta_script "$tsv_file" "$genome_path" "${output_path}/$multigenome_file"
  # Capture the return value
  accession_count=$?
  
  # Get the filename without the path and extension
  filename=$(basename "$tsv_file" .tsv)
  abundance_file="mocks_composition/${filename}.tsv"
  # Define a dictionary of models
  declare -A models
  models=( ["novaseq"]="150_reads" )
  #["miseq"]="300_reads"
  #models=( ["ref_read75.npz"]="ref_75_reads")

  # Iterate through the list of models
  for model in "${!models[@]}"; do
    # create the filename for the metagenome
    mock_filename="${filename}_${models[$model]}"
    # call the iss software
    echo "$iss_command --abundance_file /mnt/data/${abundance_file} -m /mnt/data/$model -o /mnt/data/${mock_filename}"
    $iss_command --abundance_file /mnt/data/${abundance_file} -m $model -o /mnt/data/${mock_filename}
    echo "Synthetic metagenome created successfully: ${mock_filename}"
    echo ""
  done

done