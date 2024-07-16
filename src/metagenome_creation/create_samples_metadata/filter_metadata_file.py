## Author: Pablo Viana
## Version: 1.0 (13/07/2023)
## Created: 13/07/2023 (dd/mm/yyyy)

import datetime

# Main function.
def main():
    start_time = datetime.datetime.now()
    print("Starting script:")
    
    samples_file = r"data/final_biome_classification.csv"
    input_metadata_file = r"data/merged_metadata.csv"
    output_metadata_filename = r"results/merged_metadata.csv"

    valid_samples = set()
    with open(samples_file, "r") as file:
        # Iterate over each line in the file
        for line in file:
            sample = line.strip().split(",")[0].replace("\"", "")
            valid_samples.add(sample)


    output_content = ""
    with open(input_metadata_file, "r") as file:
        # Iterate over each line in the file
        for line in file:
            sample = line.strip().split(";")[0]
            if sample in valid_samples:
                output_content += line
                
    with open(output_metadata_filename, "w") as output_file:
        output_file.write(output_content)

    end_time = datetime.datetime.now()
    elapsed_time = (end_time - start_time).total_seconds()
    print("Script finished, processing time " + str(elapsed_time) + "s.")


# Call main function.
if __name__ == "__main__":
    main()