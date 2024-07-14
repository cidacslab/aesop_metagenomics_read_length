

def get_accession_abundance(file_path):
    accession_abundance = {}    
    with open(file_path, 'r') as fasta_file:       
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                accession_id = None
                splits = line[1:].split('.')
                accession_id = splits[0] + "." + splits[1].split('_')[0]
                
                if accession_id not in accession_abundance:
                    accession_abundance[accession_id] = 0
                accession_abundance[accession_id] += 1
    
    return accession_abundance


def main():
    # Example usage
    input_file = r"C:\Users\pablo\Documents\github\taxonomy_analysis\data\mock02.fasta"
    output_file = r"C:\Users\pablo\Documents\github\taxonomy_analysis\results\mock02_accession_abundance.csv"  

    abundance = get_accession_abundance(input_file)
    with open(output_file, 'w') as out_file:
        out_file.write("accession_id,abundance\n")
        for accession_id, counts in abundance.items():
            out_file.write(f"{accession_id},{counts}\n")

    
if __name__ == '__main__':
    main()
