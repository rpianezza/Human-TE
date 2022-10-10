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
args = parser.parse_args()

if os.stat(args.stat).st_size == 0:
	print("{0}\tFAILED\t{1}".format(args.stat,"empty"))
	sys.exit()

failed=[]
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
	if not (eid=="scg" or eid=="scgx"):
		continue
	coverage=float(a[4])
	if not (coverage>0.0):
		failed.append(a[1])

if len(failed)==0:
	print("{0}\tOK".format(args.stat))
else:
	print("{0}\tFAILED\t{1}".format(args.stat,",".join(failed)))
	
	