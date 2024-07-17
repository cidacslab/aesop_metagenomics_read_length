import os, sys, shutil
import utils.kraken_report_parser as KrakenParser
from utils.utility_functions import get_files_in_folder
from metagenome_creation.define_samples_composition.create_composition_from_fastq import count_reads_by_sequence_id


def count_abundance_by_species(accession_taxid_counts, accession_abundance, tree_by_taxid, output_file):
  print(f"Count abundance by species, output file: {output_file}")
  output_content = "parent_taxid,level,taxid,name,kraken_classifed_reads,NT_rPM\n"

  KrakenParser.clear_abundance_from_tree(tree_by_taxid)
  for taxid_counts in accession_taxid_counts.values():
    for taxid, count in taxid_counts.items():
      if taxid not in tree_by_taxid:
        print(f"Invalid taxid: {taxid}")
        continue
      tree_by_taxid[taxid].set_abundance(count)

  total_reads = float(sum([accession_abundance[accession].count for accession in accession_abundance]))
  for taxid, node in tree_by_taxid.items():
    if node.level == "S" or node.level_enum == KrakenParser.Level.G:
      level = node.level_enum - KrakenParser.Level.G + 1
      abundance = node.acumulated_abundance
      output_content += f"{node.parent.taxid},{level},{taxid},{node.name},{abundance}\n"

  with open(output_file, "w") as file:
    file.write(output_content)


def count_abundance_by_accession_reads(accession_taxid_counts, accession_abundance, accession_taxid, tree_by_taxid, output_file):
  print(f"Count abundance by level, output file: {output_file}")
      
  output_content = "read_accession_id,read_accession_total_count,root,abundance_root"
  output_content += ",superkingdom,abundance_superkingdom,kingdom,abundance_kingdom"
  output_content += ",phylum,abundance_phylum,class,abundance_class,order,abundance_order"
  output_content += ",family,abundance_family,genus,abundance_genus,species,abundance_species\n"

  for accession_id, taxid_counts in accession_taxid_counts.items():
    KrakenParser.clear_abundance_from_tree(tree_by_taxid)
    for taxid, count in taxid_counts.items():
      if taxid not in tree_by_taxid:
        print(f"Invalid taxid: {taxid}")
        continue
      tree_by_taxid[taxid].set_abundance(count)

    abundance_by_level = [[None, 0] for _ in range(10)]
    # fill abundance level 0 as total accession reads abundance
    abundance_by_level[0][0] = accession_id
    abundance_by_level[0][1] = accession_abundance[accession_id].count
    # start at accession species tax id and fill each level of taxonomic rank classification
    tax_id = accession_taxid[accession_id]
    node = tree_by_taxid[tax_id]
    while node is not None:
      if len(node.level) == 1:
        level = node.level_enum
        abundance_by_level[level][0] = node.name
        abundance_by_level[level][1] = node.acumulated_abundance
      node = node.parent
    # get abundance from species and go up the taxonomic hierarchy filling missing ranks
    last_abundance = abundance_by_level[1][1]
    for level in range(2, 10):
      if abundance_by_level[level][0] is None:
        abundance_by_level[level][0] = "None"
        abundance_by_level[level][1] = last_abundance
      last_abundance = abundance_by_level[level][1]    

    output_content += ",".join([f"{a[0]},{a[1]}" for a in abundance_by_level])
    output_content += "\n"
                    
  with open(output_file, "w") as file:
    file.write(output_content)  


def get_accession_taxid_abundance(kout_file, output_file):
  print(f"Get accession taxid abundance from: {kout_file}")
  accession_taxid_counts = {}

  with open(kout_file, "r") as kraken_file:
    for line in kraken_file:
      line = line.strip().split()
      if len(line) >= 3 and line[0] == "C":
        accession_id = line[1].rsplit('_', 2)[0].strip()
        taxid = line[2].strip()
        if len(accession_id) == 0 or len(taxid) == 0:
          continue 
        if accession_id not in accession_taxid_counts:
          accession_taxid_counts[accession_id] = {}
        if taxid not in accession_taxid_counts[accession_id]:
          accession_taxid_counts[accession_id][taxid] = 0
        accession_taxid_counts[accession_id][taxid] += 1

  output_content = "read_accession_id,taxid,count\n"
  for accession_id, taxid_counts in accession_taxid_counts.items():
    for taxid, count in taxid_counts.items():
      output_content += f"{accession_id},{taxid},{count}\n"

  with open(output_file, "w") as out_file:
    out_file.write(output_content)

  return accession_taxid_counts


def get_accession_taxid_abundance(out_file):
  print(f"Get accession taxid abundance from: {out_file}")
  accession_taxid_counts = {}
  with open(out_file, "r") as file:
    csv_reader = csv.reader(file)
    next(csv_reader)
    for row in csv_reader:
      accession_id = row[0].strip()
      taxid = row[1].strip()
      count = int(row[2].strip())        
      if accession_id not in accession_taxid_counts:
        accession_taxid_counts[accession_id] = {}
      accession_taxid_counts[accession_id][taxid] = count
  return accession_taxid_counts


def get_accession_taxid(input_file):
  print(f"Get accession taxid from: {input_file}")
  accession_taxid = {}

  with open(input_file, "r") as file:
    csv_reader = csv.reader(file)
    next(csv_reader)
    for row in csv_reader:
      accession_id = row[0].strip()
      tax_id = row[1].strip()
      if len(accession_id) > 0:
        accession_taxid[accession_id] = tax_id

  return accession_taxid


def main():
  base_path = "results/mocks_throat_based"
  input_extension = '_1.fastq'
  input_fastq_path = f"{base_path}/mock_metagenomes"
  # input_fastq_path = f"{base_path}/pipeline_outputs/1-fastp_output"
  input_kraken_path = f"{base_path}/pipeline_outputs/2-kraken_results"
  input_metadata_path = f"{base_path}/metadata_new"
  output_path = f"{base_path}/performance_metrics"
  
  shutil.rmtree(output_path)
  os.makedirs(output_path, exist_ok=True)

  all_files = get_files_in_folder(input_fastq_path, input_extension)
  print(all_files)

  # for file in all_files:
  file = all_files[0]
  print(f"Analyzing file: {file}")
  filename = os.path.basename(file).replace(input_extension, "")

  accession_abundance = count_reads_by_sequence_id(file)

  report_file = os.path.join(input_kraken_path, filename + ".kreport")
  _, tree_by_taxid = KrakenParser.load_kraken_report_tree(report_file)

  splits = filename.split("_")
  meta_filename = "_".join(splits[0:-2])
  metadata_file = os.path.join(input_metadata_path, meta_filename + ".csv")
  accession_taxid = get_accession_taxid(metadata_file)

  # accession_taxid_abundance_file = os.path.join(input_kraken_path, filename + ".kout")
  output_class_file = os.path.join(output_path, filename + "_out.csv")
  # accession_taxid_counts = get_accession_taxid_abundance(accession_taxid_abundance_file, output_class_file)
  accession_taxid_counts = get_accession_taxid_abundance(output_class_file)

  abundance_by_level_file = os.path.join(output_path, filename + "_level_abundance.csv")
  count_abundance_by_accession_reads(accession_taxid_counts, accession_abundance, accession_taxid, tree_by_taxid, abundance_by_level_file)

  abundance_by_species_file = os.path.join(output_path, filename + "_species_abundance.csv")
  count_abundance_by_species(accession_taxid_counts, accession_abundance, tree_by_taxid, abundance_by_species_file)

  #output_file = os.path.join(output_path, filename + "_tree_out.csv")
  #write_output_file(read_counts, output_file)


if __name__ == '__main__':
    main()