'''
Author: RICCARDO PIANEZZA
'''

import argparse
import os
import statistics
import pandas as pd

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('matrix', metavar='input', type=str,
                    help='path to input SNP matrix')
parser.add_argument('sync', metavar='sync', type=str,
                    help='path to input sync file')
parser.add_argument('output', metavar='output', type=str,
                    help='path to output file')
parser.add_argument('--type', type=str, default='te', help='Type of sequence ("te" for TEs, "scg" for single copy genes on autosomes and x, "krab" for KRAB-ZNF).')
parser.add_argument('--IDs', type=str, default='/Volumes/Temp1/rpianezza/archaic-introgression/ID-afr-females.txt', help='Path to the file with the IDs of the population to compare (es. Africans)')
parser.add_argument('--SNPs', type=int, default=10000, help='Number of variants to select')
parser.add_argument('--copynumber', type=str, default='/Volumes/Temp1/rpianezza/archaic-introgression/female-copynumber.tsv', help='Path to the file with copynumber for each TE family')

args = parser.parse_args()

# Takes as input a list of IDs in a txt file, one ID per line. Return a list.
def ID_list (IDs):
    with open(IDs, "r") as infile:
        sample_list = [line.strip() for line in infile.readlines()]
        return sample_list

# Takes as input the matrix containing all the sync files. Output a matrix with selected columns of the IDs we are interested in.   
def filter_columns (input, output):
    african_indexes = [0, 1]
    with open(input) as infile:
        header = next(infile).rstrip().split('\t')
        for i, cell in enumerate(header):
            cell_part = cell.split('-')[0]
            if cell_part in africans:
                african_indexes.append(i+1)
        with open(output, 'w') as outfile:
            filtered_header = [header[i-1] for i in african_indexes[1:]]
            outfile.write('familyname\tposition' + '\t'.join(filtered_header) + '\n')
            for line in infile:
                cells = line.rstrip().split('\t')
                filtered_cells = [cells[i] for i in african_indexes]
                outfile.write('\t'.join(filtered_cells) + '\n')

# Takes as input the matrix containing all the sync files. Output a filtered matrix with only the requested sequence type (scg, te)
def filter_type (input, output):
    sequence2 = "_"
    if args.type == "te":
        sequence = "_te"
    elif args.type == "scg":
        sequence = "_scg"
        sequence2 = "_scgx"
    else:
        sequence = "_krab"
    with open(input) as infile, open(output, 'w') as outfile:
        header = next(infile)
        outfile.write(header)
        for line in infile:
            cells = line.strip().split('\t')
            if (cells[0].endswith(sequence) or cells[0].endswith(sequence2)):
                outfile.write(line)

# Takes as input the matrix containing all the sync files and a matrix containing a subset of "familyname" and "positions".
# Output: a filtered matrix that keeps only the duo familyname+positions present in the input file
def filter_rows (matrix, input, output):
    with open(input) as infile, open(matrix) as matrix, open(output, 'w') as outfile:
        header = next(matrix)
        columns = header.strip().split("\t")
        new_header = "\t".join(columns[0:])
        outfile.write("familyname\tposition\t"+new_header+"\n")
        familynames = []
        positions = []
        for line in infile:
            cells = line.strip().split('\t')
            familynames.append(cells[0])
            positions.append(cells[1])
        for line in matrix:
            cells = line.strip().split('\t')
            if cells[0] in familynames:
                for i, fam in enumerate(familynames):
                    if fam == cells[0] and positions[i] == cells[1]:
                        outfile.write(line)

# Takes as input a line of the sync matrix, calculate the average frequency of each base, return a list of 4 elements with the avg frequencies.
def avg_af(a_line):
    freqs = [[], [], [], []]
    cells = a_line.strip().split('\t')
    for cell in cells[2:]:
        components = cell.split(':')[:4]
        sam = sum(map(int, components))
        for i, component in enumerate(components):
            if sam>4:
                f = int(component) / sam
                freqs[i].append(f)
            else:
                f = "NA"
                freqs[i].append(f)
    avg_freqs = [statistics.mean([freq for freq in freq_list if freq != 'NA']) if 'NA' not in freq_list else 'NA' for freq_list in freqs]
    return avg_freqs

# Takes as input a line in the sync matrix format, calculates the allele frequencies for a speficic cell ("column") and return list of 4 allele frequencies.
def af(a_line, column):
    freqs = []
    cell = a_line.strip().split('\t')[column]
    components = cell.split(':')[:4]
    sam = sum(map(int, components))
    for i, component in enumerate(components):
        if sam>4:
            f = round((int(component) / sam),2)
            freqs.append(f)
        else:
            f = "NA"
            freqs.append(f)
    return freqs

# Takes as input a sync matrix structured file, calculate the average frequencies for all the samples and write an output file
def africa_avg(input, output):
    with open(input) as infile, open(output, 'w') as outfile:
        next(infile)
        outfile.write("familyname\tposition\tA\tT\tC\tG\n")
        for line in infile:
            cells = line.strip().split('\t')
            average_frequencies = avg_af(line)
            outfile.write(f"{cells[0]}\t{cells[1]}\t{average_frequencies[0]}\t{average_frequencies[1]}\t{average_frequencies[2]}\t{average_frequencies[3]}\n")

# Takes as input a single sync file. Output a file with allele frequencies (splitted in 4 columns for each position)
def archaic_sync(input_sync, output):
    with open(input_sync) as infile, open(output, 'w') as outfile:
        outfile.write("familyname\tposition\tA\tT\tC\tG\n")
        for line in infile:
            cells = line.strip().split('\t')
            average_frequencies = af(line, 3)
            outfile.write(f"{cells[0]}\t{cells[1]}\t{average_frequencies[0]}\t{average_frequencies[1]}\t{average_frequencies[2]}\t{average_frequencies[3]}\n")

# Merge the outputs of africa_avg(), archaic_sync() and args.copynumber.
def merge(afr_af, archaic_af, copynumber, output):
    africa = pd.read_csv(afr_af, sep="\t")
    archaic = pd.read_csv(archaic_af, sep="\t")
    copies = pd.read_csv(copynumber, sep="\t")
    merged_df1 = pd.merge(africa, archaic, on=["familyname", "position"])
    merged_df1.dropna(inplace=True)
    merged_df = pd.merge(merged_df1, copies, on=["familyname"])
    merged_df.to_csv(output, sep="\t", index=False)

# Extract diagnostic SNPs, defined as variants not found in any African genome but found in archaic in a frequency of at least 1/mean_copynumber
def extract_diagnostic(arch_var_matrix, output):
    bases = ["A", "T", "C", "G"]
    with open(arch_var_matrix) as infile, open(output, 'w') as outfile:
        next(infile)
        outfile.write("familyname\tposition\tdiagnostic\tmean_copynumber\tafrican_freq\tarchaic_freq\n")
        for line in infile:
            cells = line.strip().split('\t')
            for i,base in enumerate(bases):
                if float(cells[i+2])==0 and float(cells[i+6]) >= 1/float(cells[-1]):
                    outfile.write(f"{cells[0]}\t{cells[1]}\t{base}\t{cells[-1]}\t{cells[i+2]}\t{cells[i+6]}\n")

# Extract the diagnostic bases frequencies in each HGDP genome, write the final output file
def find_introgressed(HGDP_arch_var, diagnostic, output):
    bases = ["A", "T", "C", "G"]
    with open(HGDP_arch_var) as infile, open(diagnostic) as diag, open(output, 'w') as outfile:
        header = next(infile)
        header_pieces = header.strip().split('\t')[2:]
        outfile.write("familyname\tposition\tdiagnostic\t" + '\t'.join(header_pieces) + '\n')
        for line in infile:
            frequencies = []
            cells = line.strip().split('\t')
            for l in diag:
                c = l.strip().split('\t')
                if (c[0] == cells[0]) and (c[1] == cells[1]):
                    d = c[2]
                    diagnostic_base = bases.index(d)
            for i in range(2, len(cells)):
                freqs = af(line, i)
                freq = freqs[diagnostic_base]
                frequencies.append(str(freq))
            diag.seek(0)
            to_write = '\t'.join(frequencies)
            outfile.write(cells[0] + '\t' + cells[1] + '\t' + d + '\t' + to_write + '\n')


africans = ID_list(args.IDs)
filter_columns(args.matrix, args.output+"_columns")
filter_type(args.output+"_columns", args.output+"_type")
africa_avg(args.output+"_type", args.output+"_africa")
archaic_sync(args.sync, args.output+"_archaic")
merge(args.output+"_africa", args.output+"_archaic", args.copynumber, args.output+"_afr-arch")
extract_diagnostic(args.output+"_afr-arch", args.output+"_diagnostic")
filter_rows(args.matrix, args.output+"_diagnostic", args.output+"_HGDP-arch-var")
find_introgressed(args.output+"_HGDP-arch-var", args.output+"_diagnostic", args.output)

#os.remove(args.output+"_columns") 
#os.remove(args.output+"_type")    
#os.remove(args.output+"_africa")
#os.remove(args.output+"_archaic")
#os.remove(args.output+"_merged")
#os.remove(args.output+"_HGDP-arch-var")