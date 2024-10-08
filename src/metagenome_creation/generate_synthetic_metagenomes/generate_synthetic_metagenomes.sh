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
project_path=""
# Path to the CSV files containing the accession numbers composition
accession_path="results/mocks_throat_based/composition2"
# Path to the folder containing the genome files
genome_path="data/genomes"
# Output multifasta file
output_path="results/mocks_throat_based/mock_metagenomes"
multigenome_file="data/genomes_output.fasta"

# creates output path if it doesnt exist
mkdir -p $output_path

# Source the script containing the create_multifasta function
create_multifasta_script="src/metagenome_creation/4-generate_synthetic_metagenomes/generate_multifasta_genomes_from_accessions.sh"

# Command to run the InSilicoSeq synthetic metagenome generator
iss_command="iss generate --seed 1234 --cpus 25 -t metagenomics --genomes ${multigenome_file}"

# Iterate over all files in the accession folder and its subfolders
for tsv_file in $(find "$accession_path" -type f -name "*.tsv"); do
  
  echo "Creating fasta for dataset: $tsv_file"  
  # Create the multifasta genome file
  echo "$create_multifasta_script $tsv_file $genome_path $multigenome_file"
  $create_multifasta_script "$tsv_file" "$genome_path" "$multigenome_file"
  
  # Get the filename without the path and extension
  filename=$(basename "$tsv_file" .tsv)

  # Define a dictionary of models
  declare -A models
  # read_sizes=( ["75"]="75bp_reads"  ["150"]="150bp_reads"  ["300"]="300bp_reads")
  models["-m novaseq"]="150bp_reads"
  models["-m miseq"]="300bp_reads"
  models["--mode basic"]="75bp_reads"
  #models=( ["ref_read75.npz"]="ref_75_reads")

  # Iterate through the list of models
  for model in "${!models[@]}"; do
    # create the filename for the metagenome
    mock_filename="${output_path}/${filename}_${models[$model]}"
    # call the iss software
    echo "$iss_command --abundance_file $tsv_file $model -n 2m -o $mock_filename"
    $iss_command --abundance_file $tsv_file $model -n 2m -o $mock_filename
    echo "Synthetic metagenome created successfully: ${mock_filename}"
    echo ""
  done  
done

echo "Finished generating synthetic metagenome script!"