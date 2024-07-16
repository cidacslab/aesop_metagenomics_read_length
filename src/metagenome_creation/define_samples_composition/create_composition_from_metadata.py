import os, datetime, random
from utils.utility_functions import get_files_in_folder
from utils.utility_functions import get_accessions
from utils.utility_functions import check_composition_abundance


def create_taxa_composition(metadata_file, include_human_contaminant=False):
  # get all accessions from metadata file
  accession_ids = get_accessions(metadata_file)

  # remove human from accessions
  if include_human_contaminant:
    human_accession = "AP023461.1"
    accession_ids.remove(human_accession)

  # Generate random values for each taxon
  random_values = [random.lognormvariate(0, 1) for _ in accession_ids]
      
  # sum all abundances to do a normalization
  sum_randoms = sum(random_values)
  total = sum_randoms
  if include_human_contaminant:
    human_perct = random.uniform(0.7, 0.9)
    human_count = (human_perct/(1-human_perct)) * sum_randoms
    total = sum_randoms + human_count
  
  # normalize random values to sum 1.0
  taxon_abundance = {taxon: round(abundance/total, 18) for taxon, abundance in zip(accession_ids, random_values)}
  
  # adjust the last abundance so that the sum is exactly 1.0  
  if include_human_contaminant:
    sum_values = sum(taxon_abundance.values())
    taxon_abundance[human_accession] = 1.0 - sum_values
  else:
    last_taxon,_ = taxon_abundance.popitem()
    sum_values = sum(taxon_abundance.values())
    taxon_abundance[last_taxon] = 1.0 - sum_values
  return taxon_abundance


def main():
  """
  Main function.
  """
  # input_path = "results/pipeline_mock/metadata"
  # output_path = "results/pipeline_mock/composition"
  random.seed(1234)

  output_path = "/home/pedro/aesop/github/aesop_metagenomics_read_length/results/mocks_throat_based/composition2"
  meta_path = "/home/pedro/aesop/github/aesop_metagenomics_read_length/results/mocks_throat_based/metadata"
  file = os.path.join(meta_path, "throat_with_pathogen_01.csv")

  # for file in get_files_in_folder(input_path, ".csv"):
  #random generate the taxa abundance
  taxon_abundance = create_taxa_composition(file)

  # Get the base name (file name with extension)
  base_name = os.path.basename(file)

  # Get the file name without the extension
  file_name = os.path.splitext(base_name)[0]

  # write output file
  output_file = os.path.join(output_path, file_name + ".tsv")

  output_contents = ""
  for k, v in taxon_abundance.items():
    output_contents += f"{k}\t{v:.18f}\n"

  with open(output_file, "w") as file:
    file.write(output_contents)
  
  file = os.path.join(output_path, "throat_with_pathogen_01.tsv")
  check_composition_abundance(file)

    
if __name__ == '__main__':
    main()