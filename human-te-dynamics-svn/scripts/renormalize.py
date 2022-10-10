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


class Feature:

	def __init__(self,feature,featid,featlen,ncount):
		self.feature=feature
		self.featid=featid
		self.featlen=featlen
		self.ncount=ncount

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
parser.add_argument('--exclude', type=str, default=None,dest="exclude", required=False, help="Comma separated list of scg to exclude")
parser.add_argument("--min-threshold", type=float, required=False, dest="minth", default=0.0, help="minimum threshold for scg")
parser.add_argument("--max-threshold", type=float, required=False, dest="maxth", default=1000000.0, help="maximum threshold  for scg")
args = parser.parse_args()


# set up list of scg to exclude
toexclude=set([])
if args.exclude is not None:
	if "," in args.exclude:
		toexclude=set(args.exclude.split(","))
	else:
		toexclude.add(args.exclude)
		

ncount_control=0.0
ncount_sum=0.0

retainedfeature=[]
excluded=set([])

oldscgcov=0.0
oldscgcount=0.0
for l in open(args.stat):
	"""
	summary	mapping_to_scg_weighted	396076.602730239
	summary	mapping_to_scgx_weighted	8358.3443708609
	0	1					2		3					4
	scg	chr1:916864-921016	4152	1132.1059602649	1.1058663467441125
	scg	chr1:960583-965719	5136	1232.0993377483421	0.9729569374650899                                                                 
	"""
	l=l.rstrip("\n")
	a=l.split("\t")
	eid=a[0]  # entry id
	if eid=="summary":
		print(l)
		continue
	featid,featlen,ncount,copy_number=a[1],float(a[2]),float(a[3]),float(a[4])
	f=Feature(eid,featid,featlen,ncount)
	ncount_control+=ncount
	
	if f.feature == "scg":
		oldscgcount+=1.0
		oldscgcov+=(f.ncount/f.featlen)
		if f.featid in toexclude:
			excluded.add(f.featid)
			continue
		elif copy_number < args.minth or copy_number > args.maxth:
			excluded.add(f.featid)
			continue
	ncount_sum+=ncount
	retainedfeature.append(f)
	


# renormalize rpm; not necessary
#rescalefactor=ncount_control/ncount_sum
#for r in retainedfeature:
#	r.ncount = r.ncount * rescalefactor 


# avcov scg
meanscgcov=0.0
scgcount=0.0

for r in retainedfeature:
	if r.feature == "scg":
		coverage=r.ncount/r.featlen
		meanscgcov+=coverage
		scgcount+=1.0
meanscgcov=meanscgcov/scgcount

# report final
# print excluded
# summary renormalized excluded ""
# summary renormalization factor sumcontrol/rpmsum
# summary oldavcovscg newavcovscg

print("{0}\t{1}\t{2}".format("summary","renormalize_lost",",".join(excluded)))
print("{0}\t{1}\t{2} : {3}".format("summary","renormalize total coverage, old : new",ncount_control,ncount_sum))
print("{0}\t{1}\t{2} : {3}".format("summary","renormalize average coverage scg, old : new",oldscgcov/oldscgcount, meanscgcov))
	
# scg	chr1:960583-965719	5136	1232.0993377483421	0.9729569374650899           	
for r in retainedfeature:
	featcov=r.ncount/r.featlen
	copies =featcov/meanscgcov
	topr=[r.feature,r.featid,str(r.featlen),str(r.ncount),str(copies)]
	print("\t".join(topr))
	
	