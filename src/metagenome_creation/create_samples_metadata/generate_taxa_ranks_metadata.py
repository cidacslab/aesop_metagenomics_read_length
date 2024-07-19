# Install the ete3 library
# !pip install ete3
from utils.get_ncbi_taxid_from_accession import get_taxid_from_accession
from utils.utility_functions import get_files_in_folder
from ete3 import NCBITaxa
import csv, os

# Initialize NCBI Taxa object
ncbi = NCBITaxa()


# Function to get taxonomic hierarchy for a given taxid
def get_taxonomic_hierarchy(taxid):

  try:
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    ranks = ncbi.get_rank(lineage)
    hierarchy = {ranks[taxid]: (names[taxid], taxid) for taxid in lineage}
    return hierarchy
  except:
    return {}
  
  
# function to get taxid of given accession
def get_taxids_from_composition_accessions(composition_file):
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
          print(f"Found taxid for accession {accession}: {taxid}")
  return accession_tax_ids


def main():
  base_path = "results/mocks_sample_based"
  # base_path = "results/mocks_sample_based"
  input_extension = ".tsv"
  input_path = f"{base_path}/composition"
  output_path = f"{base_path}/metadata"
  output_extension = ".csv"

  # shutil.rmtree(output_path)
  os.makedirs(output_path, exist_ok=True)

  for file in get_files_in_folder(input_path, input_extension):
    filename = os.path.basename(file).replace(input_extension, output_extension)
    output_file = os.path.join(output_path, filename)

    accession_tax_ids = get_taxids_from_composition_accessions(file)

    # Read the ground truth data and process each row
    with open(output_file, 'w', newline='') as outfile:
      ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
      fieldnames = ['accession', 'tax_id'] + ranks + [ rank+'_taxid' for rank in ranks]
      writer = csv.DictWriter(outfile, fieldnames=fieldnames)  
      # Write the header
      writer.writeheader()
      
      for accession, taxid in accession_tax_ids.items():
        row = {'accession': accession, 'tax_id': taxid }
        hierarchy = get_taxonomic_hierarchy(taxid)
        print(hierarchy)
        for rank in ranks:
          if rank in hierarchy:
            name,taxid = hierarchy[rank]
            row[rank] = name
            row[rank+'_taxid'] = taxid
            last_rank = rank
          else:        
            row[rank] = f"Undefined {rank} in {row[last_rank]}"
            row[rank+'_taxid'] = row[last_rank+'_taxid']
        writer.writerow(row)


if __name__ == '__main__':
    main()