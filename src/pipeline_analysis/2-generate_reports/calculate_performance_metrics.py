from utils.utility_functions import get_files_in_folder
import os, sys, csv


def load_accession_taxids(input_file):
  accession_taxids = {}
  with open(input_file, 'r') as file:
    next(file)
    for line in file:      
      line_splits = line.strip().split(',')
      accession_id = line_splits[0].strip()
      genus_taxid = line_splits[-2].strip()
      species_taxid = line_splits[-1].strip()
      if accession_id == "":
        continue
      accession_taxids[accession_id] = (genus_taxid, species_taxid)
  return accession_taxids


def load_reported_tax_abundance(input_file):
  print(f"Load total tax id reads file: {input_file}")
  reported_taxid_abundance = {}
  with open(input_file, 'r') as file:
    csv_reader = csv.reader(file)
    next(csv_reader)
    for row in csv_reader:
      taxid = row[2].strip()
      # name = row[3]
      reads = int(row[4].strip())
      reported_taxid_abundance[taxid] = reads
  print(reported_taxid_abundance)
  return reported_taxid_abundance


def load_accession_lineage_classified(input_file):
  print(f"Load accession level abundance file: {input_file}")
  accession_lineage_taxa = {}
  lineage_taxa_counts = [{} for _ in range(10)]
  lineage_taxa_counts[0] = {"total_reads": 0}
  with open(input_file, 'r') as file:
    next(file)
    for line in file:
      line_splits = line.strip().split(',')
      accession_id = line_splits[0]
      accession_total_reads = int(line_splits[1])
      accession_lineage_taxa[accession_id] = []
      lineage_taxa_counts[0][accession_id] = accession_total_reads
      lineage_taxa_counts[0]["total_reads"] += accession_total_reads
      index = 1
      for i in range(3, 30, 3):
        taxname = line_splits[i]
        taxid = int(line_splits[i+1])
        mapped_reads = int(line_splits[i+2])
        accession_lineage_taxa[accession_id].append(taxid)
        if taxid not in lineage_taxa_counts[index]:
          lineage_taxa_counts[index][taxid] = {"tax_name": taxname, "total_reads": 0, "correct_mapped": 0}
        lineage_taxa_counts[index][taxid]["total_reads"] += accession_total_reads
        lineage_taxa_counts[index][taxid]["correct_mapped"] += mapped_reads
        index += 1
  return (accession_lineage_taxa, lineage_taxa_counts)


def get_confusion_matrix_values(sample_total_reads, total_tax_reads, total_mapped_to_tax, correct_tax_reads):
  # print(f"{sample_total_reads}, {total_tax_reads}, {total_mapped_to_tax}, {correct_tax_reads}")
  # reads from tax_id mapped to correct tax_id
  true_positive = correct_tax_reads
  # reads from tax_id not mapped to this tax_id
  false_negative = total_tax_reads - correct_tax_reads
  # reads not from tax_id mapped to this tax_id
  false_positive = total_mapped_to_tax - correct_tax_reads
  # reads not from tax_id not mapped to this tax_id
  true_negative = sample_total_reads - total_tax_reads - false_positive

  # accuracy = (true_positive + true_negative) / float(true_positive + false_negative + false_positive + true_negative)
  # sensitivity = (true_positive) / float(true_positive + false_negative)
  # specificity = (true_negative) / float(false_positive + true_negative)
  # precision = (true_positive) / float(true_positive + false_positive)
  return (true_positive, true_negative, false_positive, false_negative)


def calculate_confusion_matrix(accession_lineage_taxa, lineage_taxa_classified, reported_tax_abundance, output_file):
  print(f"Calculating confusion matrix: {output_file}")
      
  output_content = "accession_id,total_reads,genus,genus_taxid,genus_true_positive,genus_true_negative,"
  output_content += "genus_false_positive,genus_false_negative,species,species_taxid,species_true_positive,"
  output_content += "species_true_negative,species_false_positive,species_false_negative\n"
  
  sample_total_reads = lineage_taxa_classified[0]["total_reads"]

  for accession_id, lineage_taxa in accession_lineage_taxa.items():
    accession_total_reads = lineage_taxa_classified[0][accession_id]

    species_taxid= lineage_taxa[-1]
    species_name = lineage_taxa_classified[-1][species_taxid]["tax_name"]
    species_total_reads = lineage_taxa_classified[-1][species_taxid]["total_reads"]
    species_correct_reads = lineage_taxa_classified[-1][species_taxid]["correct_mapped"]
    species_total_classified = reported_tax_abundance.get(species_taxid, 0)
    
    genus_taxid = lineage_taxa[-2]
    genus_name = lineage_taxa_classified[-2][genus_taxid]["tax_name"]
    genus_total_reads = lineage_taxa_classified[-2][genus_taxid]["total_reads"]
    genus_correct_reads = lineage_taxa_classified[-2][genus_taxid]["correct_mapped"]
    genus_total_classified = reported_tax_abundance.get(genus_taxid, 0)

    print(f"{species_taxid}, {species_name}, {species_total_reads}, {species_correct_reads}, {species_total_classified}")
    print(f"{genus_taxid}, {genus_name}, {genus_total_reads}, {genus_correct_reads}, {genus_total_classified}")
    
    genus_metrics = get_confusion_matrix_values(sample_total_reads, genus_total_reads, genus_total_classified, genus_correct_reads)
    species_metrics = get_confusion_matrix_values(sample_total_reads, species_total_reads, species_total_classified, species_correct_reads)

    output_content += f"{accession_id},{accession_total_reads},{genus_name},{genus_taxid},{genus_metrics[0]},"
    output_content += f"{genus_metrics[1]},{genus_metrics[2]},{genus_metrics[3]},{species_name},{species_taxid},"
    output_content += f"{species_metrics[0]},{species_metrics[1]},{species_metrics[2]},{species_metrics[3]}\n"
      
  with open(output_file, "w") as out_file:
    out_file.write(output_content)



def main():
  base_path = "/home/work/aesop/github/aesop_metagenomics_read_length/results/mocks_throat_based"
  input_extension = '_level_abundance.csv'
  input_metadata_path = f"{base_path}/metadata"
  input_metrics_path = f"{base_path}/performance_metrics"
  
  output_path = f"{base_path}/performance_metrics"
  output_extension = "_metrics.csv"    

  os.makedirs(output_path, exist_ok=True)

  print("Starting process...")
  all_files = get_files_in_folder(input_metrics_path, input_extension)
  print(all_files)

  for file in all_files:
    print("")
    print(f"Analyzing file: {file}")

    filename = os.path.basename(file).replace(input_extension, "")
    
    # splits = filename.split("_")
    # meta_filename = "_".join(splits[0:-2])
    # metadata_file = os.path.join(input_metadata_path, meta_filename + ".csv")
    # accession_taxids = load_accession_taxids(metadata_file)

    accession_lineage_classified_file = os.path.join(input_metrics_path, filename + "_level_abundance.csv")
    accession_lineage_taxa, lineage_taxa_classified = load_accession_lineage_classified(accession_lineage_classified_file)
    
    reported_tax_abundance_file =  os.path.join(input_metrics_path, filename + "_species_abundance.csv")
    reported_tax_abundance = load_reported_tax_abundance(reported_tax_abundance_file)
    
    output_file = os.path.join(output_path, filename + output_extension)
    calculate_confusion_matrix(accession_lineage_taxa, lineage_taxa_classified, reported_tax_abundance, output_file)


if __name__ == '__main__':
    main()