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
            