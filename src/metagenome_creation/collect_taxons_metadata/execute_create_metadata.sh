#!/bin/bash

# Specify the path to the folder containing the files
input_path="results/mock_composition"

# Specify the path to the new script to run on each file
metadata_script="./src/collect_taxons_metadata/create_metadata_for_taxons.sh"

# Specifiy the path to the output folder
output_path="results/mock_composition_metadata"
mkdir -p $output_path

# Iterate over each file in the folder
for file in $(find "$input_path" -type f -name "*.csv"); do
    filename=$(basename $file ".csv")
    output_file="${output_path}/${filename}_metadata.csv"
    # Execute the new script on the file
    bash -c "$metadata_script $file $output_file"
done
