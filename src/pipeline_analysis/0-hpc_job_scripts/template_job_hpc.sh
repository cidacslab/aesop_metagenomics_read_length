#!/bin/bash
############################################################
#       BEGIN Job                                          #
############################################################
#SBATCH --job-name='AESOP Job'              # Job name
#SBATCH --partition=cpu_iterativo           # CPU batch queue
#SBATCH --nodes=1                           # Maxinum amount of nodes
#SBATCH --cpus-per-task=40                  # Maxinum amount of cores
#SBATCH --mem=1024GB                        # Maxinum amount of memory
#SBATCH --time=999:00:00                    # Time limit hrs:min:sec
#SBATCH --output=aesop_%j.log               # Standard output log
#SBATCH --error=aesop_%j.err                # Standard error log
###########################################################
:<<DOC
Author: Pablo Viana
Created: 2023/03/16

Template script used to start a SLURM batch job on biome metagenomic samples.
DOC

# Number os parallel processes to be executed
num_processes=3

repository_src="/home/work/aesop/github/aesop_metagenomics_read_length/src/pipeline_analysis"
# repository_src="/home/pablo.viana/metagenomics_src"

# Script that call the pipeline for each dataset
script_for_datasets="$repository_src/0-hpc_job_scripts/execute_script_for_datasets.sh"

# Pipeline script to be executed
script="$repository_src/0-hpc_job_scripts/execute_metagenomic_pipeline.sh"
# script="/home/pablo.viana/metagenomics_src/0-hpc_job_scripts/execute_report_scripts.sh"

echo "Execute script: $script"

# All datasets
sample_datasets="
                sample_based
                # throat_based
                "
                 
# Trim all lines, then filter out comments and empty lines
sample_datasets=$(echo "$sample_datasets" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' | grep -v '^[[:space:]]*$' | grep -v '^[[:space:]]*#')

echo "On sample datasets:"
echo "$sample_datasets"
echo ""

# Template: singularity run [SINGULARITY_OPTIONS] <sif> [COMMAND_OPTIONS]
# singularity exec /opt/images/cidacs/biome.sif  $script_for_datasets "$num_processes" "$script" "$sample_datasets"
# singularity exec /opt/images/cidacs/cidacs-jupyter-datascience-v1-r2.sif $script_for_datasets "$num_processes" "$script" "$sample_datasets"
$script_for_datasets "$num_processes" "$script" "$sample_datasets"

############################################################
#       END Job                                            #
############################################################
