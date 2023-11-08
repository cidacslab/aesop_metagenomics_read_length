from Bio import Entrez

def get_taxid_from_accession(accession_id):
    # Email address required by Entrez
    Entrez.email = 'pablo.alessandro@gmail.com'
    Entrez.api_key = '86cf88e5ee6087442f57c78ed90336b99408' 
    handle = Entrez.efetch(db='nuccore', id=accession_id, retmode='xml')
    records = Entrez.read(handle)
    
    taxid = ""
    organism = ""
    taxonomy = ""
    for subrec in records:
        if 'GBSeq_feature-table' in subrec:
            for subrec2 in subrec['GBSeq_feature-table']:
                if 'GBFeature_quals' in subrec2:
                    for subrec3 in subrec2['GBFeature_quals']:
                        if subrec3.get('GBQualifier_name', None) == 'db_xref':
                            taxid = subrec3['GBQualifier_value'].split(':')[1]
        if 'GBSeq_organism' in records[0]:
            organism = records[0]['GBSeq_organism']
        if 'GBSeq_taxonomy' in records[0]:
            taxonomy = records[0]['GBSeq_taxonomy']
    return f"{taxid},{taxonomy}; {organism}".strip()



def main():
    input_file = r"C:\Users\pablo\Documents\github\taxonomy_analysis\data\mock_random25_reads_abundance.tsv"
    output_file = r"C:\Users\pablo\Documents\github\taxonomy_analysis\results\mock_random25_accession_metadata.csv"  

    header = False
    accession_dict = {}
    with open(input_file, 'r') as file:
        for line in file:
            if header:
                header = False
                continue
            accession_id = line.strip().split()[0]
            #print(f"Accession ID: {accession_id}")
            accession_dict[accession_id] = get_taxid_from_accession(accession_id)
            print(f"{accession_id},{accession_dict[accession_id]}")

    with open(output_file, 'w') as out_file:
        out_file.write("accession_id,taxid,taxonomy\n")
        for accession_id, taxonomy in accession_dict.items():
            out_file.write(f"{accession_id},{taxonomy}\n")

    
if __name__ == '__main__':
    main()
