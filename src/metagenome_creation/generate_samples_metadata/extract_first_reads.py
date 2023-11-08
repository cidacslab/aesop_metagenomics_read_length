## Author: Pablo Viana
## Version: 1.0 (13/07/2023)
## Created: 13/07/2023 (dd/mm/yyyy)

## USAGE: python3 extract_first_reads.py INPUT_FILE MAX_READS

## From a given directory, collect metadata from files with 
## the extension ".fasta" or ".fa" including 
## file_size, file number of reads, GC%, average read size, number of bases.

## INPUT_FILE ==> Complete path of the input file.
## MAX_READS ==> Maximum number of reads to extract.

import os, sys, datetime

# Assign input file as a variable.
input_file = sys.argv[1]
max_reads = int(sys.argv[2])


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
    complete_sequence = ""

    # Open the input file
    with open(input_file, "r") as file:
        # Iterate over each line in the file
        for line in file:
            if line.startswith(">"):
                if read_count >= max_reads:
                    break

                complete_sequence += line
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
                complete_sequence += line
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
        sample_GC_content = round(sample_GC_content, 4)

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
        highest,
        complete_sequence
    ]
    

# Main function.
def main():
    start_time = datetime.datetime.now()
    print("Starting script:")
    if input_file.endswith(".fasta"):
        print("Processing file:" + str(input_file))
        base_filename = os.path.basename(input_file)
        filename = os.path.splitext(base_filename)[0]

        fasta_metrics = analyze_file(input_file)
        print("sample;total_number_of_reads;GC_content;total_number_of_bases;"+\
              "average_read_size;read_size_mode;read_size_median;lowest_read_size;"+\
              "highest_read_size")
        print(str(filename) + ";" + ";".join([str(s) for s in fasta_metrics[0:-1]]))

        output_filename = "sample_" + filename + "_" + str(max_reads) + "_reads.fasta"
        with open(output_filename, "w") as output_file:
            output_file.write(fasta_metrics[-1])

    end_time = datetime.datetime.now()
    elapsed_time = (end_time - start_time).total_seconds()
    print("Script finished file " + str(input_file) + " processing time " + str(elapsed_time) + "s.")


# Call main function.
if __name__ == "__main__":
    main()