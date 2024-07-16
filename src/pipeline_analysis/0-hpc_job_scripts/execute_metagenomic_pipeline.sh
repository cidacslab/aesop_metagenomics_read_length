#!/bin/bash
:<<DOC
Author: Pablo Viana
Created: 2023/03/16

Template script used to run a script over the biome metagenomic samples.

params $1 - Number os parallel processes to be executed
DOC

# create alias to echo command to log time at each call
echo() {
    command echo "B_PID: $BASHPID [$(date +"%Y-%m-%dT%H:%M:%S%z")]: $@"
}
# exit when any command fails
set -e
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command ended with exit code $?." >&2' EXIT

# Start job profile
start=$(date +%s.%N)
echo "Started running job!"

#number of parallel processes
num_processes=$1
# name of the run folder
run_name=$2


################################################################################
############################### ATTENTION !!!!! ################################
################################################################################
################### FOR EACH ANALYSIS FILL THESE INFORMATION ###################
################################################################################

# run_name="rs01"
dataset_name="aesop_${run_name}"

# old_dataset_path="/scratch/pablo.viana/aesop/pipeline_v2/dataset_${run_name}"
raw_data_path="/home/work/aesop/github/aesop_metagenomics_read_length/results/mocks_${run_name}/mock_metagenomes"
base_dataset_path="/home/work/aesop/github/aesop_metagenomics_read_length/results/mocks_${run_name}/pipeline_outputs"

# Kraken2 database
kraken2_database="/dev/shm/k2_pluspfp_20240605"

# Location of src folder in the github directory
repository_src="/home/work/aesop/github/aesop_metagenomics_read_length/src/pipeline_analysis"

# Script to execute the tasks
custom_script="$repository_src/0-hpc_job_scripts/execute_custom_script.sh"


################################################################################
###################################  FASTP  ####################################
################################################################################
input_suffix="_R1.fastq"
params=("$num_processes"
        "$repository_src/1-analysis_pipeline/1-quality_control-fastp_filters.sh"
        "$dataset_name"
        "$input_suffix"
        "$raw_data_path"
        "$base_dataset_path/1-fastp_output")

# $custom_script "${params[@]}"

# echo "Tar gziping report files: tar -czf ${dataset_name}_fastp_filters_reports.tar.gz *.html *.json"
# tar -czf "${dataset_name}_fastp_filters_reports.tar.gz" *.html *.json

# rm -vf *.html *.json


################################################################################
##################################  KRAKEN2  ###################################
################################################################################
input_suffix="_1.fastq"
params=("$num_processes"
        "$repository_src/1-analysis_pipeline/2-taxonomic_annotation-kraken2.sh"
        "$dataset_name"
        "$input_suffix"
        "$base_dataset_path/1-fastp_output"
        "$base_dataset_path/2-kraken_results"
        "$kraken2_database")

# $custom_script "${params[@]}"


################################################################################
##################################  BRACKEN  ###################################
################################################################################

declare -A input_suffixes=( ["_75bp_reads.kreport"]="75" ["_150bp_reads.kreport"]="150" ["_300bp_reads.kreport"]="300" )

for input_suffix in "${!input_suffixes[@]}"; do
  read_length=${input_suffixes[$input_suffix]}

  params=("$num_processes"
          "$repository_src/1-analysis_pipeline/3-taxonomic_annotation-bracken.sh"
          "$dataset_name"
          "$input_suffix"
          "$base_dataset_path/2-kraken_results"
          "$base_dataset_path/3-bracken_results"
          "$kraken2_database"
          "$read_length")

  $custom_script "${params[@]}"

  mv ${dataset_name}_3-taxonomic_annotation-bracken_logs.tar.gz \
    ${dataset_name}_3-taxonomic_annotation-bracken${input_suffix}_logs.tar.gz
done


################################################################################
################################################################################

mv *.tar.gz $base_dataset_path

# echo ""
# df
# du -hd 4 /scratch/pablo.viana
# find /scratch/pablo.viana 

#  Finish pipeline profile
finish=$(date +%s.%N)
runtime=$(awk -v a=$finish -v b=$start 'BEGIN{printf "%.3f", (a-b)/60}')
echo ""
echo "Total elapsed time: ${runtime} min."
echo "B_PID: $BASHPID [$(date +"%Y-%m-%dT%H:%M:%S%z")]: Finished complete pipeline!"