'''
Author: RICCARDO PIANEZZA
'''

import argparse
import csv
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
parser.add_argument('--IDs', type=str, default='/Volumes/Temp1/rpianezza/ancient_humans/archaic-humans/analysis/archaic-variants/ID-afr', help='Path to the file with the IDs of the population to compare (es. Africans)')
parser.add_argument('--SNPs', type=int, default=10000, help='Number of variants to select')

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
    avg_freqs = [round(statistics.mean([freq for freq in freq_list if freq != 'NA']),2) if 'NA' not in freq_list else 'NA' for freq_list in freqs]
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
        header = next(infile)
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

# Merge the outputs of africa_avg() and archaic_sync(). Then, it calculates the sum of allele frequencies differences for each line.
def merge(afr_af, archaic_af, output):
    africa = pd.read_csv(afr_af, sep="\t")
    archaic = pd.read_csv(archaic_af, sep="\t")
    merged_df = pd.merge(africa, archaic, on=["familyname", "position"])
    merged_df.dropna(inplace=True)
    merged_df["A_diff"] = abs(merged_df["A_x"] - merged_df["A_y"])
    merged_df["T_diff"] = abs(merged_df["T_x"] - merged_df["T_y"])
    merged_df["C_diff"] = abs(merged_df["C_x"] - merged_df["C_y"])
    merged_df["G_diff"] = abs(merged_df["G_x"] - merged_df["G_y"])
    merged_df["total_diff"] = merged_df["A_diff"] + merged_df["T_diff"] + merged_df["C_diff"] + merged_df["G_diff"]
    sorted_df = merged_df.sort_values(by="total_diff", ascending=False)
    diff_positions = sorted_df.head(args.SNPs)
    diff_positions.to_csv(output, sep="\t", index=False)

# Takes as input the output of merge() and the output of filter_rows() and calculate, for each archaic variant,
# the difference in allele frequency with all the HGDP individuals.
def distance_matrix(arch_var_matrix, archaic, output):
    HGDP_af = pd.read_csv(arch_var_matrix, sep="\t")
    archaic_af = pd.read_csv(archaic, sep="\t")
    merged_df = pd.merge(HGDP_af, archaic_af, on=["familyname", "position"])
    merged_df.dropna(inplace=True)
    merged_df.to_csv(output+"_merged", sep="\t", index=False)
    with open(output+"_merged") as infile, open(output, "w") as outfile:
        header = next(infile)
        columns = header.strip().split("\t")
        new_header = "\t".join(columns[:-4])
        outfile.write(new_header + "\n")
        for line in infile:
            cells = line.strip().split("\t")
            diff_list = []
            for i, cell in enumerate(cells[2:-4]):
                freqs = af(line, i+2)
                if "NA" not in freqs:
                    diff = round((abs(freqs[0]-float(cells[-4])) + abs(freqs[1]-float(cells[-3])) + abs(freqs[2]-float(cells[-2])) + abs(freqs[3]-float(cells[-1]))),2)
                else:
                    diff = "NA"
                diff_list.append(str(diff))
                diff_line = "\t".join(diff_list)
            outfile.write(f"{cells[0]}\t{cells[1]}\t{diff_line}\n")


africans = ID_list(args.IDs)
filter_columns(args.matrix, args.output+"_columns")
filter_type(args.output+"_columns", args.output+"_type")
africa_avg(args.output+"_type", args.output+"_africa")
archaic_sync(args.sync, args.output+"_archaic")
merge(args.output+"_africa", args.output+"_archaic", args.output+"_afr-arch")
filter_rows(args.matrix, args.output+"_afr-arch", args.output+"_HGDP-arch-var")
distance_matrix(args.output+"_HGDP-arch-var", args.output+"_archaic", args.output)

os.remove(args.output+"_columns") 
os.remove(args.output+"_type")    
os.remove(args.output+"_africa")
os.remove(args.output+"_archaic")
os.remove(args.output+"_merged") 