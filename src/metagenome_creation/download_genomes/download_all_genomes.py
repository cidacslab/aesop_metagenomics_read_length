from utils.utility_functions import get_all_accessions
from Bio import Entrez
import os


def download_genome_from_accession_id(output_file, accession_id):
  # Email address required by Entrez
  Entrez.email = "pablo.alessandro@gmail.com"
  Entrez.api_key = "86cf88e5ee6087442f57c78ed90336b99408"  
  handle = Entrez.efetch(db="nucleotide", id=accession_id, retmode="test", rettype="fasta")

  print(f"Writing file for: {accession_id}")
  with open(output_file, "w") as file:
    records = handle.read()
    file.write(records)
    handle.close()


def main():
  input_path = "results/mocks_sample_based/metadata"
  input_path2 = "results/mocks_throat_based/metadata"
  output_path = "data/genomes"

  accession_list = get_all_accessions(input_path, ".csv", 0)
  accession_list2 = get_all_accessions(input_path2, ".csv", 0)
  all_accessions = accession_list.union(accession_list2)
  
  os.makedirs(output_path, exist_ok=True)

  for accession_id in all_accessions:
    output_file = os.path.join(output_path, accession_id + ".fasta")
    if not os.path.exists(output_file):
      # print(f"Trying to download file {output_file}")
      download_genome_from_accession_id(output_file, accession_id)
    else: 
      print(f"Already downloaded file {output_file}")


if __name__ == '__main__':
  main()
