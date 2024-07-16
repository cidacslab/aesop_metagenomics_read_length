import os
import sys

def count_reads_by_taxid(file_path):
    accession_taxid_counts = {}

    with open(file_path, 'r') as classified_file:
        for line in classified_file:
            line = line.strip().split()
            if len(line) >= 2:
                splits = line[0][1:].split('.')
                accession_id = splits[0] + "." + splits[1].split('_')[0]
                taxid = line[1].split("|")[1]
                #print(f"Accession {accession_id}  taxid {taxid}")
                
                if accession_id not in accession_taxid_counts:
                    accession_taxid_counts[accession_id] = {}
                if taxid not in accession_taxid_counts[accession_id]:
                    accession_taxid_counts[accession_id][taxid] = 0
                accession_taxid_counts[accession_id][taxid] += 1

    return accession_taxid_counts

        
def process_files_in_folder(input_path, output_path):
    print("Start process")
    files_fullpath = []
    for root, dirs, files in os.walk(input_path):
        for file_name in files:
            if file_name.endswith('.class'):
                file_path = os.path.join(root, file_name)
                files_fullpath.append(file_path)
    print(files_fullpath)

    for file_path in files_fullpath:
        read_counts = count_reads_by_taxid(file_path)

        filename = os.path.basename(file_path).split(".")[0]
        output_file = os.path.join(output_path, filename + "_out.csv")

        with open(output_file, 'w') as out_file:
            out_file.write("read_accession_id,taxid,count\n")
            for accession_id, taxid_counts in read_counts.items():
                for taxid, count in taxid_counts.items():
                    out_file.write(f"{accession_id},{taxid},{count}\n")


def main():
    #folder_path = sys.argv[1] if len(sys.argv) > 1 else ''
    input_path = r"C:\Users\pablo\Documents\github\taxonomy_analysis\data"
    output_folder = r"C:\Users\pablo\Documents\github\taxonomy_analysis\results"
    if not os.path.isdir(input_path):
        print("Invalid input path")
        return

    process_files_in_folder(input_path, output_folder)

    
if __name__ == '__main__':
    main()