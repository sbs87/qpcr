#!/usr/bin/python
import math
import sys

def calc_mean(vals):
    sum=0
    count=0
    for val in vals:
        sum=sum+float(val)
        count=count+1
    m=sum/count
    return(m)


def calc_stdev(vals):
    m=calc_mean(vals)
    sum=0
    count=0
    for val in vals:
        sum=sum+math.pow((m-float(val)),2)
        count=count+1
    stdev=math.sqrt(sum/(count-1))
    return(stdev)

def unique(vals):
    seen=[]
    for val in vals:
        if val not in seen:
            seen.append(val)
    return(seen)

#Def for reading files, general. Reads filename into vector. Still need to figure out header issue
def readfile(file_name):
    f=open(file_name)
    lines = f.readlines()
    headers = lines[0].rstrip().split("\t")
    vectors = [line.rstrip().split("\t") for line in lines]
    f.close()
    return(vectors) #Returns entire file contents as series of vectors