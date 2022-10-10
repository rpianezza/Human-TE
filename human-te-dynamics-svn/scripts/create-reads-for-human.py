#!/usr/bin/env python
import os
import sys
import re
import argparse
import random
import inspect
import gzip
import numpy as np



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
				seq+=line.upper()				# normal line, add to sequence; convert to upper case


def mutate(seq,errorrate):
	trbase={'A':['T','C','G'], 'T':['A','C','G'],'C':['A','T','G'],'G':['A','T','C']}
	# find the number of required mutations/base substitutions
	mean=errorrate*len(seq)
	errorcount=np.random.poisson(mean,1)[0]
	if(errorcount)<1: # if no base substitutions are required, just give the sequence back
		return seq
	
	# introduce mutations
	lseq=list(seq)
	posset=set([]) # set of already mutated positions; do not mutate  a position twice
	for i in range(0,errorcount):
		tomut=random.randint(0,len(seq)-1)
		while(tomut in posset):
			tomut=random.randint(0,len(seq)-1)
		base=lseq[tomut]
		if base in trbase:
			# mutate base only if the base is in the dictionary  (e.g. ATCG); otherwise, eg. N, ignore the base
			tr=trbase[base]
			random.shuffle(tr)
			lseq[tomut]=tr[0]
		posset.add(tomut) # add the mutated position to the set of mutated positions
	return "".join(lseq)


def writeread(chrname,refseq,pos,counter,rl,error,ofh):
	read=refseq[pos:pos+rl]
	mutread=mutate(read,error)
	readname="{0}:{1}/{2}".format(chrname,pos,counter)
	ofh.write("@"+readname+"\n")
	ofh.write(mutread+"\n")
	ofh.write("+"+readname+"\n")
	ofh.write("I"*len(mutread)+"\n")
	

def method_random(ip,ofh,rl,coverage,error):
	counter=1
	for header,seq in FastaReader(ip):
		rreq=int((len(seq)*coverage)/rl) # required number of reads (chr-len*coverage/read length)
		for i in range(0,rreq):
			pos=random.randint(0,len(seq)-rl-1)
			writeread(header,seq,pos,counter,rl,error,ofh)
			counter+=1


def method_uniform(ip,ofh,rl,coverage,error):
	counter=1
	for header,seq in FastaReader(ip):
		stepsize=float(rl)/float(coverage)
		pos=0.0
		while pos< len(seq)-rl:
			curpos=int(pos)
			writeread(header,seq,curpos,counter,rl,error,ofh)
			pos+=stepsize
			counter+=1


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
""")

parser.add_argument("--fasta", type=str, required=True, dest="fasta", default=None, help="A fasta-file containing the sequences")
parser.add_argument("--coverage", type=int, required=True, dest="coverage", default=None, help="Target coverage")
parser.add_argument("--read-length", type=int, required=True, dest="rl", default=None, help="The length of the reads")
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="output file, a fastq file")
parser.add_argument("--method", type=str, required=True, dest="method", default=None, help="uniform/random")
parser.add_argument("--error-rate", type=float, required=False, dest="error", default=0.0, help="the error rate of the generated reads")
args = parser.parse_args()

ofh=gzip.open(args.output,"w")
if(args.method.lower()=="uniform"):
	method_uniform(args.fasta,ofh,args.rl,args.coverage,args.error)
elif(args.method.lower()=="random"):
	method_random(args.fasta,ofh,args.rl,args.coverage,args.error)
else:
	raise Exception("Unknown method "+args.method)
ofh.close()



            