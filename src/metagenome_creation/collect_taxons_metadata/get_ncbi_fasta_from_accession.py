import os, sys, csv
from Bio import Entrez


def download_genome_from_accession_id(output_file, accession_id):
    # Email address required by Entrez
    Entrez.email = 'pablo.alessandro@gmail.com'
    Entrez.api_key = '86cf88e5ee6087442f57c78ed90336b99408'    
    handle = Entrez.efetch(db="nucleotide", id=accession_id, retmode="test", rettype="fasta")
    
    print(f"Writing file for: {accession_id}")
    with open(output_file, 'w') as file:
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
    accession_list = set()

    for file in get_files_in_folder(input_path, ".csv"):
        with open(file, 'r') as in_file:
            csv_reader = csv.reader(in_file)
            next(csv_reader)
            for row in csv_reader:
                accession_id = row[0].strip()
                accession_list.add(accession_id)

    with open(accessions_file, 'w') as out_file:
        for accession in accession_list:
            out_file.write(accession + "\n")



def main():
    # input_path = r"C:\Users\pablo\Documents\github\taxonomy_analysis\data\taxons_metadata"
    # output_path = r"C:\Users\pablo\Documents\github\taxonomy_analysis\data\ncbi_genomes"
    output_path = "results/pipeline_mock/genomes"
    # accessions_file = os.path.join(output_path, "all_accessions.txt")
    # accessions_file = "../../data/AESOP_AMVB_MOCKS_with_accession_final.csv"
    input_path = "results/pipeline_mock/mock_composition"
    accessions_file = "results/pipeline_mock/all_accessions.txt"

    # get_all_accessions(input_path, accessions_file)
    # return
    with open(accessions_file, 'r') as file:
        csv_reader = csv.reader(file, delimiter = ",")
        next(csv_reader)
        for row in csv_reader:
            accession_id = row[0].strip()
            
            output_file = os.path.join(output_path, accession_id + ".fasta")
            if not os.path.exists(output_file):
                print(f"Trying to redownload file {output_file}")
            else: 
                print("")#f"Already downloaded file {output_file}")
                # download_genome_from_accession_id(output_file, accession_id)
    # download_genome_from_accession_id(output_path, "GCF_002871975.1")
    
if __name__ == '__main__':
    main()
