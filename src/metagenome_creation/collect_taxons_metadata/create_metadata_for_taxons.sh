#!/bin/bash

# input_file=$1
# db_file="/home/pablo/nucl_gb.accession2taxid"
taxonkit_script="/home/work/aesop/taxonkit"
# output_file=$2
# delimiters='\t|,'

# echo "Executing script for file: $input_file"
# echo "Output file will be: $output_file"

# Writes the header to the output file
# echo "accession_id,accession_taxid,superkingdom,phylum,class,order,family,genus,species,"\
# "superkingdom_taxid,phylum_taxid,class_taxid,order_taxid,family_taxid,genus_taxid,"\
# "species_taxid" > "$output_file"


# INPUT DONT HAVE HEADER
# accession_id=0
unset accession_id
accession_id="AP023461.1"
taxid="9606"

# Read each line from the input file and process
# while IFS= read -r line; do
    #echo $line
    # Skip the first line (when acession_id variable is not set yet
    # if [[ -z "$accession_id" ]]; then
    #     accession_id=1
    #     continue
    # fi

    # Get the first column from the input file
    # accession_id=$(echo "$line" | awk -F"$delimiters" '{print $1}')
    echo $accession_id

    # Grep the first column in the grep file and extract the third column
    #grep_result=$(grep -m 1 -F "$accession_id" "$db_file")
    #echo $grep_result

    # taxid=$(echo "$line" | awk -F"$delimiters" '{print $1}')
    #taxid=$(echo "$grep_result" | awk -F"$delimiters" '{print $3}')
    echo $taxid
    if [[ -z "$taxid" ]]; then
        echo "taxid not found for line: ${line}"
        continue
    fi

    lineage=$(echo "$taxid" | $taxonkit_script reformat -I 1 -F -t -f ",{k},{p},{c},{o},{f},{g},{s}")
    # echo $lineage

    # Remove trailing spaces from lineage
    trimmed_lineage=$(echo "$lineage" | awk -F',' '{for(i=1;i<=NF;i++) {gsub(/^[[:space:]]+|[[:space:]]+$/,"",$i); printf "%s%s",$i,(i==NF?"":",");}}')
    echo $trimmed_lineage
    #break
    # Append the result to the output file
    echo $accession_id,$trimmed_lineage 
    #>> $output_file
# done < "$input_file"
