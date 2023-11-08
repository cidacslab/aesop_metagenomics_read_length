import os, sys, csv
import kraken_report_parser as KrakenParser


def get_files_in_folder(input_path, input_extension):
    print("Start process")
    files_fullpath = []
    for root, dirs, files in os.walk(input_path):
        for file_name in files:
            if file_name.endswith(input_extension):
                file_path = os.path.join(root, file_name)
                files_fullpath.append(file_path)
    return files_fullpath


def get_all_tax_ids(composition_file, taxids_set, index = 2):
    print(f"Get taxid from: {composition_file}")
    with open(composition_file, "r") as file:
        csv_reader = csv.reader(file, delimiter = ",")
        next(csv_reader) # ignore header
        for row in csv_reader:
            taxid = row[index].strip()
            taxids_set.add(taxid)



def main():

    all_taxid_set = set()

    # get_all_tax_ids("data/czid_rpip_vsp_pathogens.csv", all_taxid_set)

    get_all_tax_ids("data/AESOP_AMVB_MOCKS_with_accession_final.csv", all_taxid_set)    
    
    # for file in get_files_in_folder("data/throat_mocks_done", ".csv"):
    #     get_all_tax_ids(file, all_taxid_set)

    # for file in get_files_in_folder("data/throat_taxons_metadata", "_metadata.csv"):
    #     get_all_tax_ids(file, all_taxid_set)

    # for file in get_files_in_folder("results/mock_random", "_metadata.csv"):
    #     get_all_tax_ids(file, all_taxid_set, 1)

    with open("results/new_compositions_taxids.txt", "w") as file:
        file.write("\n".join(sorted(all_taxid_set)))
        file.write("\n")



if __name__ == '__main__':
    main()