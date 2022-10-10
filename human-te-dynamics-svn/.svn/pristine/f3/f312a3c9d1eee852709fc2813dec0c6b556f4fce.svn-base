#!/bin/python3
import argparse
import collections
import re
import statistics
import sys


class Entry:
                # x     0   1   2   3       4       5   6       7       8
    def __init__(self,sid,pop,sex,region,feature,featid,featlen,featcov,featcopy):
        self.sid=sid
        self.pop=pop
        self.sex=sex
        self.region=region
        self.feature=feature # scg,krab,scgx
        self.featid=featid
        self.featlen=featlen
        self.featcov=featcov
        self.featcopy=featcopy
        
def load_tidy(file):
    toret=[]
    for l in open(file):
        # HGDP00001	Brahui	male	Central_South_Asia	scg	chr1:916864-921016	4152	1051.8410596026488	1.01525454783029
        # 0         1       2       3                   4   5                   6       7                   8
        a=l.rstrip("\n").split("\t")
        e=Entry(a[0],a[1],a[2],a[3],a[4],a[5],float(a[6]),float(a[7]),float(a[8]))
        toret.append(e)
    return toret

class StatEntry:
    def __init__(self,count_male,mean_male,pstdev_male,quant_male,count_fem,mean_fem,pstdev_fem,quant_fem):
        self.count_male=count_male
        self.mean_male=mean_male
        self.pstdev_male=pstdev_male
        self.quant_male=quant_male
        
        self.count_fem=count_fem
        self.mean_fem=mean_fem
        self.pstdev_fem=pstdev_fem
        self.quant_fem=quant_fem
        
    def toString(self):
        cm,cf=str(self.count_male),str(self.count_fem)
        mm,mf="{:.2f}".format(self.mean_male), "{:.2f}".format(self.mean_fem)
        sm,sf="{:.2f}".format(self.pstdev_male), "{:.2f}".format(self.pstdev_fem)
        qm,qf=",".join(["{:.2f}".format(i) for i in self.quant_male]),",".join(["{:.2f}".format(i) for i in self.quant_fem])
        return "\t".join([cm,mm,sm,qm,cf,mf,sf,qf])
            
def get_stats(featlist):
    # featlist = hash[male/female]->[entry1,entry2]
    males=[d.featcopy for d in featlist["male"]]
    females=[d.featcopy for d in featlist["female"]]
    cm,cf=len(males),len(females)
    mm,mf=statistics.fmean(males),statistics.fmean(females)
    sm,sf=statistics.pstdev(males),statistics.pstdev(females)
    qm,qf=statistics.quantiles(males,n=4),statistics.quantiles(females,n=4)
    se=StatEntry(cm,mm,sm,qm,cf,mf,sf,qf)
    return se


# todo
def get_rel_stats(featlist):
    males=[d.featcopy for d in featlist["male"]]
    females=[d.featcopy for d in featlist["female"]]
    cm,cf=len(males),len(females)
    mm,mf=statistics.fmean(males),statistics.fmean(females)
    if mm<0.00000001 or mf<0.00000001:
        return "na"
    males=[d.featcopy/mm for d in featlist["male"]]
    females=[d.featcopy/mf for d in featlist["female"]]
    sm,sf=statistics.pstdev(males),statistics.pstdev(females)
    qm,qf=statistics.quantiles(males,n=4),statistics.quantiles(females,n=4)
    se=StatEntry(cm,mm,sm,qm,cf,mf,sf,qf)
    return se


def compute_absolute(feat,minte):
    # dict[Alu/L1/...] -> dict[male/female] ->  list[entry1,entry2,entry3]

    entries=[]
    for featid,li in feat.items():
        e=get_stats(li) # requires dict[male/female] ->  list[entry1,entry2,entry3]
        if e=="na":
            continue
        if e.mean_male<minte and e.mean_fem<minte:
            continue
        e.featid=featid
        entries.append(e)
    
    # sort descending by highest variance, either male or female
    soe=sorted(entries,key=lambda i: -max(i.pstdev_fem,i.pstdev_male)) 
    return soe

def compute_relative(feat,minte):
    # dict[Alu/L1/...] -> dict[male/female] ->  list[entry1,entry2,entry3]

    entries=[]
    for featid,li in feat.items():
        e=get_rel_stats(li) # requires dict[male/female] ->  list[entry1,entry2,entry3]
        if e=="na":
            continue
        if e.mean_male<minte and e.mean_fem<minte:
            continue
        e.featid=featid
        entries.append(e)
    
    # sort descending by highest variance, either male or female
    soe=sorted(entries,key=lambda i: -max(i.pstdev_fem,i.pstdev_male)) 
    return soe
        


parser = argparse.ArgumentParser(description="""           
Description
-----------
    Summarize tidyverse""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Prerequisites
-------------
    python version 3+

Authors
-------
    Robert Kofler
""")


parser.add_argument("--tidyverse", type=str, required=True, dest="tidy", default=None, help="The input mapstat file;")
parser.add_argument("--min-te-count", type=float, required=False, dest="minte", default=1.0, help="Minimum count of TEs in both males and females")
args = parser.parse_args()

data=load_tidy(args.tidy) # load data into a long least of "Entry"s


tesh=collections.defaultdict(lambda:collections.defaultdict(lambda:[])) # te sex hash
tesrh=collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:[]))) # te sex region hash
tesph=collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:[]))) # te sex population hash

scgh=collections.defaultdict(lambda:[])
scgsh=collections.defaultdict(lambda:collections.defaultdict(lambda:[])) # scgsexhash

scgxh=collections.defaultdict(lambda:[])
scgxsh=collections.defaultdict(lambda:collections.defaultdict(lambda:[])) # scgxsexhash

krabh=collections.defaultdict(lambda:[]) 
krabsh=collections.defaultdict(lambda:collections.defaultdict(lambda:[])) # krabsexhash

# sort data
for d in data:
    if(d.feature=="te"):
        tesh[d.featid][d.sex].append(d)
    elif(d.feature=="scg"):
        scgh[d.sex].append(d)
        scgsh[d.featid][d.sex].append(d)
    elif(d.feature=="scgx"):
        scgxh[d.sex].append(d)
        scgxsh[d.featid][d.sex].append(d)
    elif(d.feature=="krab"):
        krabh[d.sex].append(d)
        krabsh[d.featid][d.sex].append(d)
    else:
        raise Exception("Invalid feature")


print("scg\toverview\t{0}".format(get_stats(scgh).toString()))
print("scgx\toverview\t{0}".format(get_stats(scgxh).toString()))
print("krab\toverview\t{0}".format(get_stats(krabh).toString()))
for i in compute_absolute(tesh,args.minte):
    print("te\tabsolute\t{0}\t{1}".format(i.featid,i.toString()))
for i in compute_relative(tesh,args.minte):
    print("te\trelative\t{0}\t{1}".format(i.featid,i.toString()))
    

for i in compute_absolute(scgsh,0.0):
    print("scg\tabsolute\t{0}\t{1}".format(i.featid,i.toString()))
for i in compute_absolute(scgxsh,0.0):
    print("scgx\tabsolute\t{0}\t{1}".format(i.featid,i.toString()))
for i in compute_absolute(krabsh,0.0):
    print("krab\tabsolute\t{0}\t{1}".format(i.featid,i.toString()))
