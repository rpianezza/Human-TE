import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('input', metavar='input', type=str,
                    help='path to input file')
parser.add_argument('output', metavar='output', type=str,
                    help='path to output file')
args = parser.parse_args()

matrix = args.input
individuals = len(open(matrix).readline().strip().split("\t"))

k=0

for line in open(matrix):
    bases = [0,0,0,0]
    sam=0
    if (k!=0):
        cell_num = 0
        for cell in line.split('\t'):
            b=0 
            if (cell_num!=0) and (cell_num!=1):
                for component in cell.split(':'):
                    if b < 4:
                        bases[b]+=int(component)
                    b+=1
            cell_num+=1
    freq=[]
    sam=sum(bases)
    for u in bases:
        if sam!=0:
            freq.append(u/sam)
        else:
            freq.append(0)
    if ((max(freq)<0.8) and (sam/individuals)>5000) or (k==0):
        open(args.output,'a').writelines(line)   
    k+=1

def invert_tsv(input_path, output_path):
    with open(input_path, 'r') as input_file:
        lines = input_file.readlines()
        rows = [line.strip().split('\t') for line in lines]
        rows[0].insert(0, "ID")
        for i in range(1, len(rows)):
            rows[i][0] = f"{rows[i][0]}_{rows[i][1]}"
            rows[i].pop(1)
        inverted = [[row[i] for row in rows] for i in range(len(rows[0]))]

    with open(output_path, 'w') as output_file:
        for row in inverted:
            output_file.write('\t'.join(row) + '\n')

invert_tsv(args.output, args.output+"_inverted")

def process_tsv(input_path, output_path):
    with open(input_path, 'r') as input_file:
        lines = input_file.readlines()
        headers = lines[0].strip().split("\t")
        new_headers = [headers[0]]
        for i in range(1,len(headers)):
            for letter in ["A", "T", "C", "G", "N", "D"]:
                new_headers.append(headers[i]+letter)
        with open(output_path, 'w') as output_file:
            output_file.write(",".join(new_headers) + "\n")
            for line in lines[1:]:
                columns = line.strip().split('\t')
                first_column = columns[0]
                new_columns = []
                for i in range(1, len(columns)):
                    values = columns[i].split(":")
                    total = sum(map(int, values))
                    frequencies = [round(int(val)/total,2) for val in values]
                    new_columns.extend(frequencies)
                output_file.write(f"{first_column},{','.join(map(str,new_columns))}\n")

process_tsv(args.output+"_inverted", args.output+"_processed")