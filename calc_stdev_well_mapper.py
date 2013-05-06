#!/usr/bin/python
import math
import sys

mapped_file=sys.argv[1] #The file for mapped files from prev script
output=sys.argv[2]

def calc_mean(vals):
    sum=0
    count=0
    for val in vals:
        sum=sum+val
        count=count+1
    m=sum/count
    return(m)


def calc_stdev(vals):
    m=calc_mean(vals)
    sum=0
    count=0
    for val in vals:
        sum=sum+math.pow((m-val),2)
        count=count+1
    stdev=math.sqrt(sum/(count-1))
    return(stdev)

f=open(mapped_file,'r')
o=open(output,'w')
for line in f:
    if("Species" in line):
        header=line.rstrip()
        o.write(header+"\tSTDEV\n")
        next
    else:
        species,c1,c2,c3,c4,prevmean=line.rstrip().split("\t")
        Ct=[c1,c2,c3,c4]
        while("OMIT" in Ct):
            Ct.remove("OMIT")
        Ct_num=[]
        for val in Ct:
            Ct_num.append(float(val))
        sd=calc_stdev(Ct_num)
        if((sd <= 1) & (sd>0)):
            sd=str(sd)+"*"
        o.write(species+"\t"+c1+"\t"+c2+"\t"+c3+"\t"+c4+"\t"+prevmean+"\t"+str(sd)+"\n")
o.close()
f.close()