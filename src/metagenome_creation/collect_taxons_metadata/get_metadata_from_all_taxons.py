import os, sys, csv

def get_files_in_folder(input_path, input_extension):
    print("Start process")
    files_fullpath = []
    for root, dirs, files in os.walk(input_path):
        for file_name in files:
            if file_name.endswith(input_extension):
                file_path = os.path.join(root, file_name)
                files_fullpath.append(file_path)
    return files_fullpath


def main():
    composition_extension = ".tsv"
    composition_path = "results/new_mocks2_composition"
    lineage_file = "results/all_taxids_lineage.csv"
    taxid_file = "data/AESOP_AMVB_MOCKS_with_accession_final.csv"
    output_path = "results/new_mocks2_composition/metadata"

    count = 0
    taxid_lineage = {}
    with open(lineage_file, "r") as file:
        for line in file:
            line_splits = line.strip().split(",", 2)
            taxid = line_splits[1]
            taxid_lineage[taxid] = line_splits[2]
            # print(f"Found for {taxid} : {taxid_lineage[taxid]}")
            count += 1
    print(f"\nNot found {count} taxids!")

    taxid_by_accession = {}
    with open(taxid_file, "r") as file:
        csv_reader = csv.reader(file, delimiter=",")
        row = next(csv_reader)
        for row in csv_reader:
            accession = row[0].strip()
            taxid = row[2].strip()
            taxid_by_accession[accession] = taxid

    for input_file in get_files_in_folder(composition_path, composition_extension):
        lineage_header = taxid_lineage['accession_taxid']
        output_content = f"accession_id,accession_taxid,{lineage_header}\n"
        with open(input_file, "r") as file:
            for line in file:
                accession = line.split("\t")[0].strip()
                taxid = taxid_by_accession[accession]
                output_content += f"{accession},{taxid},{taxid_lineage[taxid]}\n"
        
        filename = os.path.basename(input_file).split(".")
        output_file = os.path.join(output_path, filename[0] + "_metadata.csv")
        with open(output_file, "w") as file:
            file.write(output_content)



if __name__ == '__main__':
    main()