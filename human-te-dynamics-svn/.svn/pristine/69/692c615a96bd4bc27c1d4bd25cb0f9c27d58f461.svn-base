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

def get_weight(cigar,readlen):
    pe=[]
    for fi in re.finditer(r"(\d+)([HSIDMN])", cigar):
        num=int(fi.group(1))
        id=fi.group(2)
        pe.append((num,id))
    matchsum=0
    for num,id in pe:
        if id=="M":
            matchsum+=num
        elif id=="D" or id=="N" or id=="I" or id=="S" or id=="H":
            pass
        else:
            raise Exception("unknown cigar"+id)
    return float(matchsum)/float(readlen)

def readfai(fai):
     telist=[]
     scglist=[]
     krablist=[]
     scgxlist=[]
     toret={}
     for l in open(fai):
          # LTR65_te	669
          a=l.rstrip("\n").split("\t")
          leng=int(a[1])
          if a[0].endswith("_te"):
               teseq=a[0][:-3]
               telist.append(teseq)
               toret[teseq]=leng
          elif a[0].endswith("_krab"):
               krabseq=a[0][:-5]
               krablist.append(krabseq)
               toret[krabseq]=leng
          elif a[0].endswith("_scg"):
               scgseq=a[0][:-4]
               scglist.append(scgseq)
               toret[scgseq]=leng
          elif a[0].endswith("_scgx"):
               scgxseq=a[0][:-5]
               scgxlist.append(scgxseq)
               toret[scgxseq]=leng
          else:
               raise Excpetion("unknown reference (supported: _te, _scg, _scgx, _krab); found: "+ref)
          

          
     return(toret,scglist,scgxlist,krablist,telist)
          



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
parser.add_argument('--sam', type=argparse.FileType('r'), default=None,dest="sam", required=True, help="A sam file")
parser.add_argument("--min-mq", type=int, required=False, dest="minmq", default=0, help="min mapping quality")
parser.add_argument("--fai", type=str, required=True, dest="fai", default=None, help="fai fasta-index (samtools faidx) of the TE database")
args = parser.parse_args()
minmq=args.minmq
lengthash,scglist,scgxlist,krablist,telist=readfai(args.fai)


hte=collections.defaultdict(lambda:0.0)
hscg=collections.defaultdict(lambda:0.0)
hkrab=collections.defaultdict(lambda:0.0)
hscgx=collections.defaultdict(lambda:0.0)

sumall=0
summapped=0
summapq=0
sumweight=0.0

sumte=0.0
sumscg=0.0
sumkrab=0.0
sumscgx=0.0


for line in args.sam:
     """
0     1         2    3    4    5   6    7    8            9                        10                  11
r1	16	M14653_te	172	70	23M	*	0	0	ATGTCGAGTTTCGTGCCGAATAA	FFFFFFFFFFFFFFFFFFBBBBB	PG:Z:novoalign	AS:i:0	UQ:i:0	NM:i:0	MD:Z:23
r2	0	M14653_te	240	70	27M	*	0	0	AACAGCTGCGGAATCGCACCGAATGCT	BBBBBFFFFFBFFFFFFFFFFFFFFFF	PG:Z:novoalign	AS:i:0	UQ:i:0	NM:i:0	MD:Z:27
     """
     if line.startswith("@"):
          continue
     line=line.rstrip("\n")
     a=line.split("\t")
     flag=int(a[1])
     
     
     
     sumall+=1
     
     if flag & 0x004 > 0:   # remove unmapped
          continue
     summapped+=1 # all weighted reads with sufficient mapping quality and 

     mq=int(a[4])
     if mq< minmq:          # remove with minmapq
          continue
     summapq+=1 # all  reads with minmapq
     
     weight=get_weight(a[5],len(a[9])) # compute weight of the read; 
     sumweight+=weight
     
     ref=a[2]
     if ref.endswith("_te"):
          teseq=ref[:-3]
          hte[teseq]+=weight
          sumte+=weight
     elif ref.endswith("_krab"):
          krabseq=ref[:-5]
          hkrab[krabseq]+=weight
          sumkrab+=weight
     elif ref.endswith("_scg"):
          scgseq=ref[:-4]
          hscg[scgseq]+=weight
          sumscg+=weight
     elif ref.endswith("_scgx"):
          scgxseq=ref[:-5]
          hscgx[scgxseq]+=weight
          sumscgx+=weight
     else:
          raise Excpetion("unknown reference (supported: _te, _scg, _scgx, _krab); found: "+ref)




meanscgcov=0.0
covscg=[]
for scg,count in hscg.items():
     length=lengthash[scg]
     coverage=float(count)/float(length)
     covscg.append(coverage)
     meanscgcov+=coverage
meanscgcov=meanscgcov/float(len(covscg))
     
 
print("{0}\t{1}\t{2}".format("summary","all_reads",sumall))
print("{0}\t{1}\t{2}".format("summary","mapped_reads",summapped))
print("{0}\t{1}\t{2}".format("summary","reads_with_mapq",summapq))
print("{0}\t{1}\t{2}".format("summary","weighted_reads_with_mapq",sumweight))
print("{0}\t{1}\t{2}".format("summary","mapping_to_te_weighted",sumte))
print("{0}\t{1}\t{2}".format("summary","mapping_to_krab_weighted",sumkrab))
print("{0}\t{1}\t{2}".format("summary","mapping_to_scg_weighted",sumscg))
print("{0}\t{1}\t{2}".format("summary","mapping_to_scgx_weighted",sumscgx))

covscg=[]
for scg in scglist:
     count=hscg[scg]
     length=lengthash[scg]
     coverage=float(count)/float(length)
     normcov=coverage/meanscgcov
     print("{0}\t{1}\t{2}\t{3}\t{4}".format("scg",scg,length,count,normcov))
for scgx in scgxlist:
     count=hscgx[scgx]
     length=lengthash[scgx]
     coverage=float(count)/float(length)
     normcov=coverage/meanscgcov
     print("{0}\t{1}\t{2}\t{3}\t{4}".format("scgx",scgx,length,count,normcov))
for krab in krablist:
     count=hkrab[krab]
     length=lengthash[krab]
     coverage=float(count)/float(length)
     normcov=coverage/meanscgcov
     print("{0}\t{1}\t{2}\t{3}\t{4}".format("krab",krab,length,count,normcov))

for te in telist:
     length=lengthash[te]
     count=hte[te]
     coverage=float(count)/float(length)
     normcov=coverage/meanscgcov
     print("{0}\t{1}\t{2}\t{3}\t{4}".format("te",te,length,count,normcov))