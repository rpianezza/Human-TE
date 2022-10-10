#!/bin/python3
import argparse
import collections
import re
import statistics
import sys
#from trunclib import *
#from helplib import *


#Things I want to filter by: Sex, Type (TE,SCG,KRAB)

parser = argparse.ArgumentParser(description="""           
Description
-----------
    Scan mapstat file for outlier sequences""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Prerequisites
-------------
    python version 3+

Authors
-------
    Florian Schwarz 
""")


parser.add_argument("--input", type=str, required=True, dest="inputfile", default=None, help="The input mapstat file;")
parser.add_argument("--printzero", type=bool, required=False, dest="printzero", default=False, help="The input mapstat file;")
parser.add_argument("--cutoff", type=float, required=False, dest="cutoff", default=1.5, help="The Cutoff how outlier sequences are defined. Default 1.5;")
parser.add_argument("--gender", type=str, required=False, dest="gender", default="all", help="Analyze genders jointly or seperately. WHAT ARE THE OPTIONS?!;")
parser.add_argument("--seqtype", type=str, required=False, dest="seqtype", default="all", help="Analyze te, scg, krab or all. Default: all;")

args = parser.parse_args()

inputfile,cutoff,seqtype,gender,printzero=str(args.inputfile),float(args.cutoff),str(args.seqtype),str(args.gender),bool(args.printzero)


valuehash=collections.defaultdict(lambda:0)
mainhash=collections.defaultdict(lambda:0)
with open (inputfile) as mqfile:
    for line in mqfile:
        mqline=line.rstrip('\n').split('\t')
        linename=mqline[0]
        linepop=mqline[1]
        linegender=mqline[2]
        lineregion=mqline[3]
        lineseqtype=mqline[4]
        lineseqname=mqline[5]
        linevalue=mqline[8]
        if gender!="all":
            if linegender!=gender:
                continue
        if seqtype!="all":
            if lineseqtype!=seqtype:
                continue
        if lineseqname not in mainhash.keys(): # had to include this because .append did not work otherwise, but I feel like it should not be necessary!
            valuehash[lineseqname]=[]
            mainhash[lineseqname]=[]
        valuehash[lineseqname].append(float(linevalue))
        mainhash[lineseqname].append([linename,linepop,linegender,lineregion,lineseqtype])            

#def Extract(value,lst): #Function to return the X element from each nested list
#    return [item[value] for item in lst]



#print(len(mainhash[lineseqname])) #Should always be the number of samples (170 for first test)

for key in valuehash.keys():
    meanvalue=statistics.mean(valuehash[key])
    if meanvalue==0:
        meanvalue=1
    i=-1
    for value in valuehash[key]:
        i+=1
        result=value/meanvalue
        if result<1 and result!=0:
            result=1/result
        elif result==0 and printzero:
            print(key,mainhash[key][i][0],mainhash[key][i][1],result)
        if result > cutoff:
            print(key,mainhash[key][i][0],mainhash[key][i][1],result)