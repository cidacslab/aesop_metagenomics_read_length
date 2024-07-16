## Author: Eduardo HC Galvao - eduardohcgalvao@gmail.com
## Version: 4.0 (28/06/2023)
## Created: 03/01/2023 (dd/mm/yyyy)

## USAGE: python3 get_samples_metadata.py INPUT_FILE

## From a given directory, collect metadata from files with 
## the extension ".fasta" or ".fa" including 
## file_size, file number of reads, GC%, average read size, number of bases.

## INPUT_FILE ==> Main directory containing samples files.

import os, sys, datetime, statistics

# Assign input file as a variable.
input_file = sys.argv[1]

# Get actual time for log file.
def update_time():
    now = datetime.datetime.now()
    current_time = now.strftime("%H:%M:%S")
    return current_time

# Function to get the size of a file.
def get_file_size(file):
    try:
        file_size = os.path.getsize(file)
    except:
        file_size = "error_with_file_size"
    return file_size

def analyze_file(input_file):
    # Initialize variables
    GC_count = 0
    complete_sequence_length = 0
    sample_number_of_bases = 0
    total_base_pairs = 0

    read_count = 0
    read_size = 0
    total_reads = 0
    read_sizes = {}

    # Open the input file
    with open(input_file, "r") as file:
        # Iterate over each line in the file
        for line in file:
            if line.startswith(">"):
                # Count read when encountering a line starting with ">"
                read_count += 1                
                # Track the frequency of each read size
                if read_size > 0:
                    if read_size not in read_sizes:
                        read_sizes[read_size] = 1
                    else:
                        read_sizes[read_size] += 1
                    # Accumulate the total base pairs and total read count
                    total_base_pairs += read_size
                    total_reads += 1
                    read_size = 0           
            else:
                line = line.strip().upper()
                # Remove any spaces from the sequence
                complete_sequence_length += len(line.replace(" ", ""))
                # Count the number of Gs and Cs (GC content)
                GC_count += line.count("G") + line.count("C")
                # Count the total number of bases (GC, AT, N)
                sample_number_of_bases += line.count("G") + line.count("C") +\
                    line.count("A") + line.count("T") + line.count("N")
                # Determine the size of the current read
                read_size += len(line)
                
        # Track the frequency of each read size
        if read_size > 0:
            if read_size not in read_sizes:
                read_sizes[read_size] = 1
            else:
                read_sizes[read_size] += 1
            # Accumulate the total base pairs and total read count
            total_base_pairs += read_size
            total_reads += 1
            read_size = 0   

    # Calculate the GC content of the sequence
    sample_GC_content = 0.0
    if complete_sequence_length > 0:
        sample_GC_content = GC_count / float(complete_sequence_length)

    # Calculate the average read size
    average_read_size = 0.0    
    if total_reads > 0:
        average_read_size = total_base_pairs / float(total_reads)

    # Determine the mode, median, lowest, and highest read sizes
    mode = 0
    mode_frequency = 0
    read_list = []
    for number, frequency in read_sizes.items():
        read_list.extend([number] * frequency)
        if frequency > mode_frequency:
            mode_frequency = frequency
            mode = number
    sorted_read_list = sorted(read_list)
    lowest = sorted_read_list[0]
    highest = sorted_read_list[-1]
    length = len(sorted_read_list)
    middle_index = length // 2    
    median = sorted_read_list[middle_index]
    if length % 2 == 0:
        median = (median + sorted_read_list[middle_index - 1]) / 2

    # Return all the analysis results
    return [
        read_count,
        sample_GC_content,
        sample_number_of_bases,
        complete_sequence_length,
        average_read_size,
        mode,
        median,
        lowest,
        highest
    ]
    

# Main function.
def main():
    start_time = datetime.datetime.now()
    print("Starting script at " + update_time() + ".")
    #with open("samples_metadata.csv", "w") as output_file:
        #header = "sample;file_size;total_number_of_reads;GC_content;" + \
        #"total_number_of_bases;average_read_size;read_size_mode;" + \
        #"read_size_median;lowest_read_size;highest_read_size\n"
        #output_file.write(header)
        #for file in os.listdir(input_path):
    if input_file.endswith(".fasta") or input_file.endswith(".fa"):
        print("Processing file " + str(input_file) + " at " + update_time() + ".")
        #input_file = os.path.join(input_path, file)
        file_size = get_file_size(input_file)
        fasta_metrics = analyze_file(input_file)
        base_filename = os.path.basename(input_file)
        filename = os.path.splitext(base_filename)[0]
        output_filename = "samples_metadata_ " + filename + ".csv"
        with open(output_filename, "w") as output_file:
            output_file.write(str(filename) + ";")
            output_file.write(str(file_size) + ";")
            output_file.write(";".join([str(s) for s in fasta_metrics]))
            output_file.write("\n")
    end_time = datetime.datetime.now()
    elapsed_time = (end_time - start_time).total_seconds()
    print("Script finished file " + str(input_file) + " processing time " + str(elapsed_time) + "s.")


# Call main function.
if __name__ == "__main__":
    main()