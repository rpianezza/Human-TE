import os
import sys
import re
import argparse
import random
import inspect
import collections

class Sam:
	def __init__(self,line):

		# LTR65_1	0	LTR65	1	70	50M	*	0	0	TGAGAAAGTAAAAATCTGCCCTGCCATTCATCAGGCTGGGATAACAGACA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	PG:Z:novoalign	AS:i:0	UQ:i:0	NM:i:0	MD:Z:50
		# LTR65_2	0	LTR65	2	70	50M	*	0	0	GAGAAAGTAAAAATCTGCCCTGCCATTCATCAGGCTGGGATAACAGACAC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	PG:Z:novoalign	AS:i:0	UQ:i:0	NM:i:0	MD:Z:50
		a=line.rstrip("\n").split("\t")
		readname=a[0]
		t=a[0].split("_")
		del t[-1]
		origin="_".join(t)
		flag=int(a[1])
		mapref=a[2]
		mappos=int(a[3])
		mapqual=int(a[4])
		cigar=a[5]
		self.readname=readname
		self.origin=origin
		self.mapref=mapref
		self.flag=flag
		self.mappos=mappos
		self.mapqual=mapqual
		self.cigar=cigar



def updatelh(lh, line):
	# @SQ     SN:AluSg7       LN:309  AS:humdfam.fa.nix       M5:bf8f536f0315a904020fd1bd7ce9c355
	a=line.split("\t")
	seq=a[1][3:]
	lens=int(a[2][3:])
	lh[seq]=lens
	return lh


parser = argparse.ArgumentParser(description="""           
Description
-----------
    Self mapping statistics""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
-------
    Robert Kofler
""")

#parser.add_argument("--fai", type=str, required=True, dest="fai", default=None, help="Fasta index of samtools")
parser.add_argument("--sam", type=str, required=True, dest="sam", default=None, help="mapping file in sam")
args = parser.parse_args()

lh={}
originh=collections.defaultdict(lambda:0)
maprefh=collections.defaultdict(lambda:0)

# map matrix
mm=collections.defaultdict(lambda:collections.defaultdict(lambda:0))

for l in open(args.sam):
	if l.startswith("@"):
		if l.startswith("@SQ"):
			lh=updatelh(lh,l)

	else:
		s=Sam(l)
		if s.flag & 0x04 >0: # get rid of unmapped
			continue
		originh[s.origin]+=1
		maprefh[s.mapref]+=1
		mm[s.origin][s.mapref]+=1

orkeys=sorted(originh.keys())
mapkeys=sorted(maprefh.keys())
for originte in orkeys:
	for mapte in mapkeys:
		originreads=originh[originte]
		crossmap = mm[originte][mapte]
		mapreads=maprefh[mapte]
		fracmap=float(crossmap)/float(originreads)
		fracofmap=float(crossmap)/float(mapreads)
		#sizeo=lh[originte]
		sizem=lh[mapte]

		topr=[originte,mapte, sizem, fracmap,fracofmap, crossmap,originreads]
		topr=map(str,topr)
		print "\t".join(topr)






