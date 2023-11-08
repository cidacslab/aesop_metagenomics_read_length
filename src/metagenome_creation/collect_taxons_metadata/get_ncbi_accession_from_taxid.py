import os, sys, csv
from Bio import Entrez


def get_accession_from_taxid(taxid):
    # Email address required by Entrez
    Entrez.email = 'pablo.alessandro@gmail.com'
    Entrez.api_key = '86cf88e5ee6087442f57c78ed90336b99408'    
    #search_term = f"txid{taxid}[Organism:exp] AND refseq[filter] NOT mitochondrion[Title]"
    # search_term = f"txid{taxid}[Organism:exp] AND (refseq[Filter]) AND (complete genome[Title])"
    search_term = f"txid{taxid}[Organism:exp] AND (complete genome[Title])"
    # search_term = f"txid{taxid}[Organism:exp] AND (chromosome[Title])"
    print(search_term)
    handle = Entrez.esearch(db="nucleotide", retmax=10, term=search_term, idtype="acc")
    records = Entrez.read(handle)
    print(records)

    accession_id = -1
    if int(records["Count"]) > 0:
        accession_id = records["IdList"][0]
    else:
        print(f"Error taxid not found: {records}")
    return str(accession_id)


def get_files_in_folder(input_path, input_extension):
    files_fullpath = []
    for root, dirs, files in os.walk(input_path):
        for file_name in files:
            if file_name.endswith(input_extension):
                file_path = os.path.join(root, file_name)
                files_fullpath.append(file_path)
    return files_fullpath


def main():

    # input_path = r"C:\Users\pablo\Documents\github\sim_results2"
    # output_path = r"C:\Users\pablo\Documents\github\sim_results3"
    # input_file = "../../data/AESOP_AMVB_MOCKS.csv"
    input_file = "../../data/AESOP_AMVB_MOCKS_with_accession.csv"
    output_file = "../../data/AESOP_AMVB_MOCKS_with_accession2.csv"

    output_content = ""
    # output_content = "accession_id"
    with open(input_file, 'r') as file:
        csv_reader = csv.reader(file, delimiter = ",")
        # next(csv_reader)
        header = next(csv_reader)
        output_content += f"{','.join(header)}\n"
        for row in csv_reader:
            taxid = row[2].strip()
            accession_id = row[0].strip()
            print(f"Acessionid: {accession_id}")
            if accession_id == "-1":
                accession_id = get_accession_from_taxid(taxid)
            output_content += f"{accession_id},{','.join(row[1:])}\n"
            print(f"Accession for taxid: {taxid} is: {accession_id}")

    with open(output_file, 'w') as file:
        file.write(output_content)

    # for file in get_files_in_folder(input_path, ".csv"):        
    #     print(f"Analysing file: {file}")
    #     output_lines = []
    #     header = True
        
        # with open(file, 'r') as in_file:
        #     for line in in_file:
        #         if header:
        #             output_lines.append(line)
        #             header = False
        #         else:
        #             splits = line.split(",")
        #             # accession_in_file = splits[0].strip()
        #             taxid = splits[1].strip()
        #             # if accession_in_file == "-1":
        #             #     continue
        #             #accession_id = get_accession_from_taxid(taxid)
        #             #line = accession_id + "," + ",".join(splits[1:])
        #             #print(f"New accession for taxid: {taxid} is: {accession_id}")
        #             output_lines.append(line)
        #             #print(output_lines[-1])

        # filename = os.path.basename(file)
        # output_file = os.path.join(output_path, filename)
        # with open(output_file, 'w') as out_file:
        #     for line in output_lines:
        #         out_file.write(line)
        # break

    
if __name__ == '__main__':
    main()
