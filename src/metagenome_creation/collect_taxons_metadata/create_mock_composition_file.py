import os, csv, datetime, random


def get_files_in_folder(input_path, input_extension):
    files_fullpath = []
    for root, dirs, files in os.walk(input_path):
        for file_name in files:
            if file_name.endswith(input_extension):
                file_path = os.path.join(root, file_name)
                files_fullpath.append(file_path)
    return files_fullpath


def get_accessions(input_file):
    accession_list = set()
    with open(input_file, 'r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)
        for row in csv_reader:
            accession_id = row[0].strip()
            accession_list.add(accession_id)
    return accession_list


def main():
    """
    Main function.
    """
    start_time = datetime.datetime.now()
    print(f"Script start: {start_time}")

    input_path = "results/pipeline_mock/metadata"
    output_path = "results/pipeline_mock/composition"
    human_accession = "AP023461.1"

    for file in get_files_in_folder(input_path, ".csv"):
        accession_ids = get_accessions(file)

        # remove human from accessions
        accession_ids.remove(human_accession)

        # Generate random values for each taxon
        random_values = [random.lognormvariate(0, 1) for _ in accession_ids]
            
        # Normalize the values so their sum is 1
        sum_randoms = sum(random_values)
        human_perct = random.uniform(0.7, 0.9)
        human_count = (human_perct/(1-human_perct)) * sum_randoms
        total = sum_randoms + human_count

        taxon_abundance = {taxon: round(abundance/total, 18) for taxon, abundance in zip(accession_ids, random_values)}
        
        sum_values = sum(taxon_abundance.values())
        taxon_abundance[human_accession] = 1 - sum_values
    
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

    
    # Print out some timing info
    end_time = datetime.datetime.now()
    print(f"Script end: {end_time}")
    print(f"Elapsed time: {(end_time - start_time).total_seconds()}s")
    print(f"Execution finished!")
    
if __name__ == '__main__':
    main()