from utils.get_fastq_read_info import count_reads_by_sequence_id
from utils.utility_functions import get_files_in_folder
# from utils.utility_functions import get_accessions
import os


def main():
  # Example usage
  output_path = "/home/pedro/aesop/github/aesop_metagenomics_read_length/results/mocks_throat_based/composition"
  input_path = "/home/pedro/aesop/results_read_length_review_old/dataset_throat_based/0-raw_samples"
  input_extension = "_75_reads_R1.fastq"

  for file in get_files_in_folder(input_path, input_extension):
    read_abundance = count_reads_by_sequence_id(file)

    output_contents = ""
    for sequence_id, read in read_abundance.items():
      output_contents += f"{sequence_id}\t{read.abundance:.18f}\n"
      
    filename = os.path.basename(file).replace(input_extension, ".tsv")
    output_file = os.path.join(output_path, filename)
    with open(output_file, "w") as out_file:
      out_file.write(output_contents)

  # # Test file
  # file = os.path.join(input_path, "throat_with_pathogen_01" + input_extension)

  # read_abundance = count_reads_by_sequence_id(file)
  # read_accessions = set(read_abundance.keys())

  # meta_path = "/home/pedro/aesop/github/aesop_metagenomics_read_length/results/mocks_throat_based/metadata"
  # file = os.path.join(meta_path, "throat_with_pathogen_01.csv")
  # accessions = get_accessions(file)
  # taxids = get_accessions(file, -1)
  # print(f"Tax ids: {taxids}")
  # print(f"Len read_accessions {len(read_accessions)}")
  # print(f"Len tax ids {len(taxids)}")
  # print(f"Len all accessions {len(accessions)}")
  # print(f"Len diff {len(accessions.difference(read_accessions))}")

  # files = get_files_in_folder("/home/pedro/aesop/github/aesop_metagenomics_read_length/data/genomes", ".fasta")
  # files_accessions = set()
  # for file in files:
  #   filename = os.path.basename(file).replace(".fasta", "")
  #   files_accessions.add(filename)

  # print(f"Len files_accessions {len(files_accessions)}")
  # print(f"Len diff {len(accessions.difference(files_accessions))}")


if __name__ == '__main__':
  main()