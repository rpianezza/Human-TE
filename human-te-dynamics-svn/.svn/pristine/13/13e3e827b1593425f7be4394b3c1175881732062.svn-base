#!/usr/bin/env python
import os
import sys
import re
import argparse
import random
import inspect



class FastaReader:
	"""
	A light-weight fasta reader;
	returns a tuple (header, sequence)
	
	"""
	def __init__(self,file):
		self.__filename=file
		self.__filehandle=open(file,"r")
		self.__prevheader=None

	def __iter__(self):
		return self
	
	def close(self):
		self.__filehandle.close()
	
	def next(self):
		line=""
		header=self.__prevheader
		seq=""
		while(1):
			line=self.__filehandle.readline()
			if line=="":					# file is empty
				if(header is not None):
					self.__prevheader=None		# return last sequence
					return (header,seq)
				else:
					raise StopIteration		# file empty and no last sequence -> STOP
			line=line.rstrip("\n")				# somethin is in the file
			if(line.startswith(">")):			# we have a header
				line=line.lstrip(">")
				if(header is None):			# if it is the first header just set the name of the header
					header=line
				else:
					self.__prevheader=line	# if it is any other header, set the previous to the current and return the sequence
					return(header,seq)
			else:
				seq+=line				# normal line, add to sequence


def mutate(seq,errorrate):
    trbase={'A':['T','C','G'], 'T':['A','C','G'],'C':['A','T','G'],'G':['A','T','C']}
    lseq=list(seq)
    for i in range(0,len(lseq)):
        # test every base if it should have an error
        if random.random()<errorrate:
            tr=trbase[lseq[i]]
            random.shuffle(tr)
            lseq[i]=tr[0]
    return "".join(lseq)


parser = argparse.ArgumentParser(description="""           
Description
-----------
    This script creates reads for TE sequences""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Version
0.1.12

Authors
-------
    Robert Kofler
    Christian Schloetterer    
""")

parser.add_argument("--te-sequences", type=str, required=True, dest="teseqs", default=None, help="A file containing the TE sequences")
parser.add_argument("--boost", type=int, required=True, dest="amp", default=None, help="Number of reads to create for every every site in the TE sequences coverage=read-length*amplification")
parser.add_argument("--read-length", type=int, required=True, dest="rl", default=None, help="The length of the reads")
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="output file, a fastq file")
parser.add_argument("--max-error-rate", type=float, required=False, dest="maxerr", default=0.05, help="the maximum 'sequencing' error rate of the generated reads")
args = parser.parse_args()

ofh=open(args.output,"w")
rl=args.rl

counter=1
for header,seq in FastaReader(args.teseqs):
    for i in range(0,len(seq)-rl):
        for k in range(0,args.amp):
            error=args.maxerr/random.random()
            read=seq[i:i+rl]
            mutread=mutate(read,error)
            ofh.write("@"+str(counter)+"\n")
            ofh.write(mutread+"\n")
            ofh.write("+"+str(counter)+"\n")
            ofh.write("I"*len(mutread)+"\n")
            counter+=1
ofh.close()
            