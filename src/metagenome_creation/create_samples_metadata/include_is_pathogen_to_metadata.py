import os

def get_files_in_folder(input_path, input_extension):
    files_fullpath = []
    for root, dirs, files in os.walk(input_path):
        for file_name in files:
            if file_name.endswith(input_extension):
                file_path = os.path.join(root, file_name)
                files_fullpath.append(file_path)
    return files_fullpath


def main():

    input_path = "results/new_mocks2/metadata"
    output_path = "results/new_mocks2/metadata_complete"
    pathogen_file = "data/czid_rpip_vsp_pathogens.csv"

    header = True
    pathogens = {}
    with open(pathogen_file, 'r') as file:
        for line in file:
            if header:
                header = False
            else:                
                splits = line.split(",")
                taxid = splits[2].strip()
                if taxid not in pathogens:
                    pathogens[taxid] = 1

    for file in get_files_in_folder(input_path, ".csv"):        
        print(f"\nAnalysing file: {file}")
        output_lines = []
        header = True
        
        with open(file, 'r') as in_file:
            for line in in_file:
                splits = line.split(",")                
                if header:
                    header = False
                    is_pathogen = "is_pathogen"
                else:
                    taxid = splits[1].strip()
                    is_pathogen = str(pathogens.get(taxid, 0))
                line = splits[0] + "," + is_pathogen + "," + ",".join(splits[1:])
                output_lines.append(line)

        filename = os.path.basename(file)
        output_file = os.path.join(output_path, filename)
        with open(output_file, 'w') as out_file:
            for line in output_lines:
                out_file.write(line)

    
if __name__ == '__main__':
    main()
