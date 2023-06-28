'''
Author: RICCARDO PIANEZZA

Script that takes as input a SNP matrix created with the script frequency-matrix.py
and another sync file (in this case an archaic human). This creates a 2 rows file containing
the allele frequencies of the archaic human for all the variants selected in the SNP matrix.
'''

import argparse
import csv
import os
import statistics

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('matrix', metavar='input', type=str,
                    help='path to input SNP matrix')
parser.add_argument('sync', metavar='sync', type=str,
                    help='path to input sync file')
parser.add_argument('output', metavar='output', type=str,
                    help='path to output file')

args = parser.parse_args()


def extract_familyname_pos(input_cell):
    last_underscore = input_cell.rfind('_')
    familyname = input_cell[:last_underscore]
    position = input_cell[last_underscore+1:]
    position = position[:-1]
    return [familyname, position]

def SNP_dict(input_sync):
    with open(input_sync) as infile:
        header = infile.readline()
        cells = header.strip().split(',')
        result = []
        for i in range(1, len(cells), 4):
            fam_pos = extract_familyname_pos(cells[i])
            result.append(fam_pos)
    return(result)

def filter_sync(input, output):
    with open(output + "_filtered", "w") as outfile:
        with open(input, "r") as infile:
            for line in infile:
                cells = line.strip().split('\t')
                search_list = [cells[0],cells[1]]
                SNP_list = SNP_dict(args.matrix)
                for sublist in SNP_list:
                    if sublist == search_list:
                        outfile.write(line)

filter_sync(args.sync, args.output)

def invert_sync(input, output):
    with open(output, "w") as outfile:
        with open(input, "r") as infile:
            lines = infile.readlines()
            rows = [line.strip().split('\t') for line in lines]
            header = [f"{row[0]}_{row[1]}" for row in rows]
            values = [row[3] for row in rows]
            outfile.write('\t'.join(header) + '\n')
            outfile.write('\t'.join(values) + '\n')

invert_sync(args.output+"_filtered", args.output+"_inverted")

def process_inverted(input, output):
    # Open the input file for reading
    with open(input, 'r') as infile:
        # Read all lines in the file
        lines = infile.readlines()
        # Get the headers from the first line and split them by tab
        headers = lines[0].strip().split("\t")
        # Create new headers for each possible nucleotide (ATCG)
        new_headers = []
        for i in range(0,len(headers)):
            for letter in ["A", "T", "C", "G"]:
                new_headers.append(headers[i]+letter)
        # Open the output file for writing
        with open(output, 'w') as outfile:
            # Write the new headers to the output file
            outfile.write(",".join(new_headers) + "\n")
            columns = lines[1].strip().split('\t')
            new_columns = []
            for i in range(0, len(columns)):
                values = columns[i].split(":")[:-2]
                total = sum(map(int, values))
                if total > 5:
                    frequencies = [round(int(val)/total,3) for val in values]
                else:
                    frequencies = ["NA" for val in values]
                new_columns.extend(frequencies)
            outfile.write(','.join(map(str,new_columns))+"\n")

process_inverted(args.output+"_inverted", args.output)

os.remove(args.output+"_filtered")
os.remove(args.output+"_inverted")