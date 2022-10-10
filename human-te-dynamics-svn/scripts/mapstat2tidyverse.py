#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random
import math
import collections
import fileinput

def read_tsv(fi):
	toret={}
	for l in open(fi):
		"""
		Sample name	Sex	Biosample ID	Population code	Population name	Superpopulation code	Superpopulation name	Population elastic ID	Data collections
		HGDP00003	male			Brahui		Central South Asia (HGDP)	BrahuiHGDP	Human Genome Diversity Project
		HGDP00015	male			Brahui		Central South Asia (HGDP)	BrahuiHGDP	Human Genome Diversity Project
		HGDP00096	male			Balochi		Central South Asia (HGDP)	BalochiHGDP	Human Genome Diversity Project

		"""
		if l.startswith("Sample"):
			continue
		a=l.rstrip("\n").split("\t")
		idd=a[0]
		sex=a[1]
		region=a[6][:-7]
		region=re.sub(" ","_",region)
		toret[a[0]]=(sex,region)
	return toret

parser = argparse.ArgumentParser(description="""           
Description
-----------
Summary statistics
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""

Authors
-------
    Robert Kofler
""")
parser.add_argument('--mapstat', type=str, default=None,dest="stat", required=True, help="Mapstat file")
parser.add_argument("--tsv", type=str, required=True, dest="tsv", default=None, help="a TSV file with metainfo (sex, region)")
args = parser.parse_args()

tsv=read_tsv(args.tsv)
p1=args.stat
if("/" in p1):
	p1=p1.split("/")[-1]
p2=p1.split(".")[0]
sid,popid=p2.split("-")
sex=tsv[sid][0]
region=tsv[sid][1]




for l in open(args.stat):
	"""
	summary	mapping_to_scg_weighted	396076.602730239
	summary	mapping_to_scgx_weighted	8358.3443708609
	0	1					2		3					4
	scg	chr1:916864-921016	4152	1132.1059602649	1.1058663467441125
	scg	chr1:960583-965719	5136	1232.0993377483421	0.9729569374650899                                                                 
	"""
	a=l.rstrip("\n").split("\t")
	eid=a[0]  # entry id
	if eid=="summary":
		continue
	
	topr=[sid,popid,sex,region,eid,a[1],a[2],a[3],a[4]]
	print("\t".join(topr))
	
	