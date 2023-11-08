import os, csv, datetime, random


def main():
    """
    Main function.
    """
    start_time = datetime.datetime.now()
    print(f"Script start: {start_time}")

    composition_file = "data/AESOP_AMVB_MOCKS_with_accession_final.csv"
    output_path = "results/new_mocks_composition2"

    samples = []
    sample_accessions = {}
    accession_info = {}
    accession_by_genus = {}

    with open(composition_file, "r") as file:
        csv_reader = csv.reader(file, delimiter=",")
        row = next(csv_reader)
        samples = row[5:]
        for sample_name in samples:
            sample_accessions[sample_name] = {}
        for row in csv_reader:
            accession = row[0].strip()
            species = row[1].strip()
            taxid = row[2].strip()
            genus = row[3].strip()
            pathogenic = row[4].strip()
            is_pathogenic = "0" if pathogenic == "No" else "1"
            accession_info[accession] = f"{taxid},{is_pathogenic},{genus},{species}"
            for index in range(5, len(row)):
                taxon_abundance = row[index].strip()
                if taxon_abundance != "":                    
                    sample_name = samples[index-5]
                    sample_accessions[sample_name][accession] = float(taxon_abundance)
            if genus not in accession_by_genus:
                accession_by_genus[genus] = []
            accession_by_genus[genus].append(accession)


    print(samples)
    # random.seed(27092023)
    for sample_name in samples[-5:]:
        min_abundance = min(sample_accessions[sample_name].values())
        max_abundance = max(sample_accessions[sample_name].values())       
        print(f"In sample {sample_name}:  Max is {max_abundance}  Min is {min_abundance}")
        for genus in accession_by_genus:
            include_genus = random.choice([True, False, False])
            if include_genus:
                print(f"  Include genus {genus}")
                for accession in accession_by_genus[genus]:
                    abundance = random.triangular(min_abundance, max_abundance, min_abundance*10)
                    sample_accessions[sample_name][accession] = abundance
                    print(f"    For accession {accession} use {abundance}")

    last_sample = samples[-1]
    new_sample1 = f"{last_sample}_1"
    new_sample2 = f"{last_sample}_2"
    sample_accessions[new_sample1] = {}
    sample_accessions[new_sample2] = {}
    samples.append(new_sample1)
    samples.append(new_sample2)

    for sample_name in samples[-2:]:
        min_abundance = min(sample_accessions[last_sample].values())
        max_abundance = max(sample_accessions[last_sample].values())       
        print(f"In sample {sample_name}:  Max is {max_abundance}  Min is {min_abundance}")
        for genus in accession_by_genus:
            print(f"  Include genus {genus}")
            for accession in accession_by_genus[genus]:
                abundance = random.triangular(min_abundance, max_abundance, min_abundance*10)
                sample_accessions[sample_name][accession] = abundance
                print(f"    For accession {accession} use {abundance}")
    

    for sample_name in samples[20:]:
        abundance_sum = sum(sample_accessions[sample_name].values())
        final_abundance_sum = 0.0
        last_accession = ""
        for accession, value in sample_accessions[sample_name].items():
            new_value = round(value / abundance_sum, 18)
            sample_accessions[sample_name][accession] = new_value
            final_abundance_sum += new_value
            last_accession = accession
        # round final value to sum 1.0
        error = 1.0 - final_abundance_sum
        sample_accessions[sample_name][last_accession] += error

        # write output file
        output_file = os.path.join(output_path, sample_name + ".tsv")
        # output_contents = "accession_id,abundance,taxid,is_pathogenic,genus,species\n"
        output_contents = ""
        for k, v in sample_accessions[sample_name].items():
            #output_contents += f"{k},{v:.18f},{accession_info[k]}\n" 
            output_contents += f"{k}\t{v:.18f}\n"
        with open(output_file, "w") as file:
            file.write(output_contents)       


    # output_contents = {}
    # if sample_name not in output_contents:
    #     output_contents[sample_name] = "taxid,abundance,is_pathogenic,species,genus\n"
    # output_contents[sample_name] += f"{taxid},{taxon_abundance},{is_pathogenic},{species},{genus}\n"
    # for sample_name in samples:
    #     output_file = os.path.join(output_path, sample_name + ".csv")
    #     with open(output_file, "w") as file:
    #         file.write(output_contents[sample_name])       


    # taxon_list = read_report_index(report_file)
    # write_taxon_list(patric_genomes_file, taxon_list, output_file)

    # Print out some timing info
    end_time = datetime.datetime.now()
    print(f"Script end: {end_time}")
    print(f"Elapsed time: {(end_time - start_time).total_seconds()}s")
    print(f"Execution finished!")



if __name__ == '__main__':
    main()