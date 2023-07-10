'''
Author: RICCARDO PIANEZZA

This script takes as input the output from the script "sync2matrix.py", and creates a matrix that can be directly
analyzed in R through a PCA or UMAP. Its main functions are:

- Find the SNPs in a dataset
- Invert the matrix to have every sample on a row with all its SNPs values as columns
- Process the matrix to separe the base counts and convert them into frequencies

Usage: in the command line, provide the input file, the output path, the number of requested SNPs, the minimum average coverage to take the
position into account for SNP calling, the type of sequence for which we wnat the SNPs, the file containing metadata on the dataset.

Example call: "python /Volumes/Temp1/rpianezza/PCA-SNPs-all-analysis/frequency_matrix_v15.py /Volumes/Temp1/rpianezza/PCA-SNPs-all-analysis/allte.matrix /Volumes/Temp1/rpianezza/PCA-SNPs-all-analysis/matrixes/NA/scg-cov15-1000SNPs.matrix.tsv --min_cov 15 --snps 1000 --type scg

How it works in details:

PHASE 1: FILTER

- First, the script will check if the suffix of the first column of each line matches the requested sequence type (es. _te) and if the prefix matches the requested class (es. LINE).
- For every selected line in the file, the script will first find the major allele by summing all the counts for each of the 4 bases in all the samples (columns).
The base with the higher count is considered the major allele.
- In every sample, the script is then calculating the major allele frequency (MAF) by dividing the count of the major allele by the sum of all the bases counts in the sample.
- After calculatuing all the MAFs, the script will calculate the variance of the MAF.
- After calculating the MAF variance for all the positions, it will order the positions in descending order of MAF variance, to have the most informative (variant) SNPs
at the top.
- Based on the chosen number of SNPs, the script is then writing to the output file only the requested number of lines (SNPs), starting from the top (most variant),
creating the file "output_filtered".
- "NA" assignment: the script will by default ignore positions with an averaged coverage across the samples below a certain threshold (--min_cov).
Furthermore, it will assign to a sample 4 "NA" instead of the 4 allele frequencies if the coverage for that position is below --min_cov/2. THe idea is that such a low
coverage is not realiable to call the correct allele frequencies.

PHASE 2: INVERT

- In this phase, the script will produce the output "_inverted". It simply invert the filtered file to have one line per sample with all the positions as columns.
- It's also merging the first two columns of the input into one (familyname+position) to create the new column names.

PHASE 3: PROCESS

- This is the phase in which the script produce the file "_processed". Basically, it substitute the base counts with allele frequencies and assigns "NA" with the
same criteria established in phase 1.
- It also split every column (1 position) into 4 (the 4 alleles). For example, the column "TE_1", containing base counts of the first position of a TE (es. 6:3:2:1) will be
splitted into 4 columns (TE_1_A, TE_1_T, TE_1_C, TE_1_G) with the respective allele frequencies (0.6, 0.3, 0.2, 0.1).

PHASE 4: MISSING DATA HANDLING

- Since both UMAP and PCA do not deal well with missing data, the script here impute a value for each missing data to remove the "NA". We don't want to remove the whole
SNP information cause it can be very informative.
- It will calculate the mean allele frequency for all the columns (familyname/position/base), separately for each "county" present in the metadata file. In the HGDP,
for example, there are 7 "country" in the metadata file.
- When the script finds a sample for which a particular column has a missing data, the script will assign the corresponding mean value for that specific familyname/position
in the country from which the sample is from. This should reduce the possible bias from missing data or very low coverage data.
-The final output is produced after this phase.
'''

# Import the argparse module to create a command line interface
import argparse
import csv
import os
import statistics

'''
The arguments that are included in this script and could be provided into the command line are:
- 'input': required, the input file path.
- 'output': required, the output file path.
- '--snps': the number of requested SNPs to be called.
- '--min_cov': the minimum coverage for a position of a TE to be included in the output file.
- '--type': the type of sequence we are interested in (scg, krab, te).
- '--teclass': the class of repetitive sequence we are interested in (DNA, LINE, SINE, LTR, NA, satellite)
- '--teMRCA': the most recent common ancestor of repetitive sequence we are interested in.
- '--metadata': a file containing geographic information of the samples in the input file (es. region). The default file is for the HGDP.
- '--classification': a file containing classification for each sequence in the library. The default file is for our customized human repseq library from RepBase.
'''

# Create an argparse object to define user arguments
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('input', metavar='input', type=str,
                    help='path to input file')
parser.add_argument('output', metavar='output', type=str,
                    help='path to output file')
parser.add_argument('--snps', type=int, help='Number of SNPs requested.')
parser.add_argument('--min_cov', type=int, default='15', help='Minimum coverage threshold.')
parser.add_argument('--type', type=str, default='te', help='Type of sequence ("te" for TEs, "scg" for single copy genes on autosomes and x, "krab" for KRAB-ZNF).')
parser.add_argument('--teclass', type=str, default='all', help='Class of repetitive sequences where the SNPs will be called. NA, DNA, LINE, SINE, LTR, satellite')
parser.add_argument('--teMRCA', type=str, default='all', help='Select only the TEs with the same common ancestor (es. "Eukaryota", "Homo sapiens", "Primates")')
parser.add_argument('--metadata', type=str, default='/Volumes/Temp1/rpianezza/0.old/summary-HGDP/ID-pop-country.tsv', help='Metadata file containing the sample names, population and region for each sample')
parser.add_argument('--classification', type=str, default='/Volumes/Temp1/rpianezza/GC-content/repbase_full_classification_mod.txt', help='File containing classification info of the repetitive sequence')

args = parser.parse_args()

# Read in the matrix file and determine the number of individuals in the file (the number of cell in the first line)
matrix = args.input
individuals = len(open(matrix).readline().strip().split("\t"))            

# Store the sequence type argument in the correct format
sequence2 = "_"
if args.type == "te":
    sequence = "_te"
elif args.type == "scg":
    sequence = "_scg"
    sequence2 = "_scgx"
else:
    sequence = "_krab"

# Create a dictionary containing the repetitive sequence name (key) and its class (value)
def create_dictionary(csv_file, key_col, value_col):
    with open(csv_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        result_dict = {row[key_col]: row[value_col] for row in reader}
    return result_dict

class_dict = create_dictionary(args.classification, 0, 3)
mrca_dict = create_dictionary(args.classification, 0, 2)

# Initialize a counter variable k
k = 0

# Initialize a dictionary that will contain index-maf couples (maf = major allele frequency)
maf_and_idx = {}

# Open the output file for writing
with open(args.output+"_filtered", 'w') as outfile:
    # Iterate through each line in the input file
    for line in open(args.input):
        # Initialize a list of lists to store the frequencies of each base in each sample
        freqs = [[], [], [], []]
        # Initialize a list to store the total number of bases in each sample
        bases = [0, 0, 0, 0]
        # Initialize a variable to store average coverage of the position across samples
        avg_cov = 0
        # If this is not the first line in the input file (header line)
        if k != 0:
            # Split the line into individual cells
            cells = line.strip().split('\t')
            split_familyname = cells[0].split("_")
            fn = "_".join(split_familyname[:-1])
            # If the first cell of the line ends with '_te', then it represents a TE position (skips all the krabs and scgs)
            # If the first line of the cell starts with a familyname which class is the requested one (or "all"), then we consider the line.
            if (cells[0].endswith(sequence) or cells[0].endswith(sequence2)) and ((class_dict[fn] == args.teclass) or (args.teclass == 'all')) and ((mrca_dict[fn] == args.teMRCA) or (args.teMRCA == 'all')):
                # Initialize a counter variable for counting the cells in the line
                cell_num = 0
                # Iterate through each cell in the line
                for cell in cells:
                    # Initialize a variable to store the total number of bases in the current sample
                    cell_sum = 0
                    # If this is not the first or second cell of the line (which contains "familyname" and "position")
                    if (cell_num != 0) and (cell_num != 1):
                        # Initialize a counter variable to keep track of the base number (1=A, 2=T, 3=C, 4=G, 5=NA, 6=DEL)
                        b = 0
                        # Iterate through each component in the cell
                        for component in cell.split(':'):
                            # If the current component is a count of a base (this skips NA and DEL)
                            if b < 4:
                                # Add the count to the total count of the current base
                                bases[b] += int(component)
                                # Add the count of the current base to the count of bases in the cell (specific position and specific sample)
                                cell_sum = cell_sum + int(component)
                                b += 1
                        # Re-initialize a counter variable to keep track of the base number (1=A, 2=T, 3=C, 4=G, 5=NA, 6=DEL)        
                        b = 0
                        # Iterate through each component in the cell
                        for component in cell.split(':'):
                            # If the current component is a count of a base and the total count of bases in the current sample is greater than half of the argument "min_cov"
                            # The idea is that a position with less then half of the average coverage is not reliable.
                            if b < 4 and cell_sum > (args.min_cov/2):
                                # Calculate the frequency of the current base in the current sample
                                f = round((int(component) / cell_sum), 3)
                                # Add the frequency to the list of frequencies for the current base
                                freqs[b].append(f)
                            # If the current component is a count of a base and the total count of bases in the current sample is zero
                            elif b < 4:
                                # Add "NaN" to the list of frequencies for the current base
                                freqs[b].append("NaN")
                            b += 1
                    cell_num += 1
                # Find the index of the base with the highest count (the major allele) across all the samples summed
                idx = sorted(range(len(bases)), key=lambda i: bases[i])[-1:]
                # Select the major allele frequency list based on the index previously selected
                major_allele_freqs = freqs[idx[0]]
                # Calculate the mean of the major allele frequencies
                numeric_maf = [x for x in major_allele_freqs if isinstance(x, (int, float))]
                if len(numeric_maf)>0:
                    mean_maf = sum(numeric_maf)/len(numeric_maf)
                else:
                    mean_maf = 0
                # Identify the NaN in the major allele frequencies (not the zeros, which are informative) and substitute them with the maf mean
                for i in range(len(major_allele_freqs)):
                    if major_allele_freqs[i] == "NaN":
                        major_allele_freqs[i] = mean_maf
                # Calculate the variance of the major allele frequency across all the samples
                maf_variance = statistics.variance(major_allele_freqs)
                # Calculate the average coverage for this position across all the samples
                avg_cov = sum(bases) / individuals
                # If the average coverage is greater or equal to the argument "min_cov"
                if avg_cov >= args.min_cov:
                    # Add the line index as key and the MAF as value to the dictionary "maf_and_inx"
                    maf_and_idx[k] = maf_variance
        else:
            # Write the header row to the output file
            outfile.write(line)
        # Increment the line counter
        k += 1

    # Sort the dictionary in descending order, putting the most variant MAFs at the top
    sorted_maf = dict(sorted(maf_and_idx.items(), key=lambda item: item[1], reverse=True))
    # Select the n most variable SNPs indexes in the matrix
    selected_idx = list(sorted(maf_and_idx, key=maf_and_idx.get, reverse=True))[:(args.snps)]

    index = 0
    # Iterate through each line in the input file
    for line in open(args.input):
        # If the line index is in the selected indexes, write the line to the output
        if index in selected_idx:
            outfile.write(line)
        index+=1


def invert_tsv(input_path, output_path):
    # Open the input file and read its lines
    with open(input_path, 'r') as input_file:
        lines = input_file.readlines()
        # Split each line by the tab character to get a list of columns
        rows = [line.strip().split('\t') for line in lines]
        # Add "ID" to the beginning of the header row
        rows[0].insert(0, "ID")
        # Combine the first two columns (familyname and position) into a single column called "ID" and remove the second column
        for i in range(1, len(rows)):
            rows[i][0] = f"{rows[i][0]}_{rows[i][1]}"
            rows[i].pop(1)
        # Remove the "-population" suffix from each sample in the matrix
        for i in range(1, len(rows)):
            for j in range(1, len(rows[i])):
                rows[i][j] = rows[i][j].split("-")[0]
        # Transpose the matrix by swapping rows and columns
        inverted = [[row[i] for row in rows] for i in range(len(rows[0]))]
    # Write the transposed matrix to the output file
    with open(output_path, 'w') as output_file:
        for row in inverted:
            output_file.write('\t'.join(row) + '\n')

# Call the function and invert the first output file (the one with the called SNPs)
invert_tsv(args.output+"_filtered", args.output+"_inverted")


def process_tsv(input_path, output_path):
    # Open the input file for reading
    with open(input_path, 'r') as input_file:
        # Read all lines in the file
        lines = input_file.readlines()
        # Get the headers from the first line and split them by tab
        headers = lines[0].strip().split("\t")
        # Create new headers for each possible nucleotide (ATCG)
        new_headers = [headers[0]]
        for i in range(1,len(headers)):
            for letter in ["A", "T", "C", "G"]:
                new_headers.append(headers[i]+letter)
        # Open the output file for writing
        with open(output_path, 'w') as output_file:
            # Write the new headers to the output file
            output_file.write(",".join(new_headers) + "\n")
            # Process each line (starting from the second) in the input file
            for line in lines[1:]:
                # Split the line by tab to get the columns
                columns = line.strip().split('\t')
                # Extract the first column (ID) and remove any suffixes separated by a hyphen
                first_column = columns[0].split("-")[0]
                # Create an empty list for the new columns
                new_columns = []
                # Iterate over the remaining columns (starting from the second)
                for i in range(1, len(columns)):
                    # Split the column value by column to get the individual nucleotide counts
                    values = columns[i].split(":")[:-2]
                    # Compute the total number of reads (sum of all nucleotide counts)
                    total = sum(map(int, values))
                    # If the total is not zero, compute the frequencies of each nucleotide
                    if total > (args.min_cov/2):
                        frequencies = [round(int(val)/total,3) for val in values]
                    # Otherwise, set the frequencies to "NA" for each nucleotide
                    else:
                        frequencies = ["NA" for val in values]
                    # Append the frequencies to the new columns list
                    new_columns.extend(frequencies)
                # Write the new row (ID and nucleotide frequencies) to the output file
                output_file.write(f"{first_column},{','.join(map(str,new_columns))}\n")

# Call the function to create the output "_processed"
process_tsv(args.output+"_inverted", args.output+"_processed")

def missing_data_handling(input_path, output_path, metadata_path):
    
    # Read metadata file and store it in a dictionary
    metadata = {}
    with open(metadata_path) as f:
        next(f) # skip header
        for line in f:
            fields = line.strip().split("\t")
            metadata[fields[0]] = (fields[1], fields[2])

    # Read input file and calculate mean values for each country
    with open(input_path, "r") as f:
        header = next(f).strip().split(",")
        data = {col: {} for col in header[1:]}
        for line in f:
            fields = line.strip().split(",")
            for i, value in enumerate(fields[1:]):
                if value == "NA":
                    continue
                col = header[i+1]
                sample = fields[0]
                country = metadata[sample][1]
                if country not in data[col]:
                    data[col][country] = []
                data[col][country].append(float(value))
        # Calculate the mean of each region for every position and every base, to use as imputation to replace NA
        means = {col: [] for col in header[1:]}
        for col, values in data.items():
            for country, country_values in values.items():
                mean = round((sum(country_values) / len(country_values)),3)
                means[col].append({"country": country, "mean": mean})
        row_names = []
    # Open the output file for writing
    with open(output_path, 'w') as output_file:
        # Write the headers to the output file
        output_file.write(",".join(header) + "\n")
        # Loop through each line in the file starting from the second line
        with open(input_path, "r") as f:
            next(f) # skip header
            for line in f:
                # Split the line by comma to get the values
                sample = line.strip().split(",")
                # Loop through each value and replace "NA" with the column mean
                for i, val in enumerate(sample):
                    if val == "NA":
                        colname = header[i]
                        sample_country = metadata[sample[0]][1]
                        mean_value = None
                        for country_mean in means[colname]:
                            if country_mean['country'] == sample_country:
                                mean_value = country_mean['mean']
                                break
                        if mean_value is not None:
                            sample[i] = str(mean_value)
                # Write the modified line to the output file
                output_file.write(','.join(sample) + '\n')

# Call the function to create the final output.
missing_data_handling(args.output+"_processed", args.output, args.metadata)

# Remove the other outputs from the previous code chunks
os.remove(args.output+"_filtered")
os.remove(args.output+"_processed")
os.remove(args.output+"_inverted")