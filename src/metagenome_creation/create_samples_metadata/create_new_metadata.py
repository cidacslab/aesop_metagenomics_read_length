from utils.get_ncbi_taxid_from_accession import get_taxid_from_accession
from utils.utility_functions import get_files_in_folder
import shutil, csv, os


def write_new_metadata_for_composition(composition_file, output_file):
  output_contents = "accession,tax_id\n"
  accession_tax_ids = {}
  with open(composition_file, "r") as file:
    csv_reader = csv.reader(file, delimiter="\t")
    for row in csv_reader:
      accession = row[0].strip()
      if len(accession) > 0:
        taxid,error = get_taxid_from_accession(accession)
        if error:
          print(f"Error getting taxid for accession {accession}: {error}")
        else:
          accession_tax_ids[accession] = taxid
          output_contents += f"{accession},{taxid}\n"

  with open(output_file, "w") as file:
    file.write(output_contents)


def main():
  base_path = "results/mocks_throat_based"
  # base_path = "results/mocks_sample_based"

  input_extension = ".tsv"
  input_path = f"{base_path}/composition"  
  output_path = f"{base_path}/metadata_new"
  output_extension = ".csv"    

  # shutil.rmtree(output_path)
  os.makedirs(output_path, exist_ok=True)

  for file in get_files_in_folder(input_path, input_extension):
    filename = os.path.basename(file).replace(input_extension, output_extension)
    output_file = os.path.join(output_path, filename)

    write_new_metadata_for_composition(file, output_file)

    
if __name__ == '__main__':
    main()