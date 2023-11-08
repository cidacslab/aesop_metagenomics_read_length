import os, sys, csv
from Bio import Entrez


def download_genome_from_accession_id(output_path, accession_id):
    # Email address required by Entrez
    Entrez.email = 'pablo.alessandro@gmail.com'
    Entrez.api_key = '86cf88e5ee6087442f57c78ed90336b99408'    
    handle = Entrez.efetch(db="nucleotide", id=accession_id, retmode="test", rettype="fasta")
    
    filename = os.path.join(output_path, accession_id + ".fasta")
    print(f"Writing file for: {accession_id}")
    with open(filename, 'w') as file:
        records = handle.read()
        file.write(records)
        handle.close()


def get_files_in_folder(input_path, input_extension):
    files_fullpath = []
    for root, dirs, files in os.walk(input_path):
        for file_name in files:
            if file_name.endswith(input_extension):
                file_path = os.path.join(root, file_name)
                files_fullpath.append(file_path)
    return files_fullpath



def get_all_accessions(input_path, accessions_file):
    accession_list = {}

    for file in get_files_in_folder(input_path, ".csv"):
        header = True        
        with open(file, 'r') as in_file:
            for line in in_file:
                if header:
                    header = False
                else:
                    splits = line.split(",")
                    accession_id = splits[0].strip()
                    accession_list[accession_id] = 1

    with open(accessions_file, 'w') as out_file:
        for accession in accession_list:
            out_file.write(accession + "\n")



def main():
    # input_path = r"C:\Users\pablo\Documents\github\taxonomy_analysis\data\taxons_metadata"
    # output_path = r"C:\Users\pablo\Documents\github\taxonomy_analysis\data\ncbi_genomes"
    output_path = "../../results/genomes"
    # accessions_file = os.path.join(output_path, "all_accessions.txt")
    accessions_file = "../../data/AESOP_AMVB_MOCKS_with_accession_final.csv"

    #get_all_accessions(input_path, accessions_file)
    #return
    with open(accessions_file, 'r') as file:
        csv_reader = csv.reader(file, delimiter = ",")
        next(csv_reader)
        for row in csv_reader:
            accession_id = row[0].strip()
            download_genome_from_accession_id(output_path, accession_id)
    # download_genome_from_accession_id(output_path, "GCF_002871975.1")
    
if __name__ == '__main__':
    main()
