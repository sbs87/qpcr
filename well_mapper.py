#!/usr/bin/python
import sys

#Define variables
#mapping files
well_ct_384={} #well # and ct value mapping
species_well={} #species and well mapping (hash of arrays)
species_order={} #species numerical order mapping
results=sys.argv[1] #The file for results
output=sys.argv[2] #name of the output file
well_mapping="/Users/stevensmith/bin/Species_to_384_96_well_mapping_file" #This SHOULD NOT CHANGE
missing_val_indicator='999'

#Read in results file, assuming results file is correctly formatted
f=open(results,'r')
for line in f:
	current={}
	try:
	    w, ct=line.rstrip().split("\t") #notice that this assumes the well # is
	except ValueError:
	    w= line.rstrip()
	    ct=missing_val_indicator
	#if(ct is "false"):
	#    print("here")
	current[w]=ct
	well_ct_384.update(current) #add to hash of arrays for well to ct mapping
f.close()

#Read in well mapping file, which is formatted as the number (1-96), corresponding species, then the four of the 384 well mapping. 
m=open(well_mapping,'r')
for line in m:
	current={}
	order={}
	num,species,w1,w2,w3,w4=line.rstrip().split("\t")
	current[species]=[w1,w2,w3,w4]
	order[num]=species
	species_well.update(current)
	species_order.update(order)
m.close()

#Output mapped Ct values
o=open(output,'w')
#Data is stored/read in the following way: a number (original sort) is mapped to a species so that the output can be properly sorted (i.e., not alphabetically)
#The current species is called with mapping the number to the species in the list. Then, each Ct value is output for the species (max 4) which is dependent on the well mapping
o.write("Species\tReplicate1\tReplicate2\tReplicate3\tReplicate4\tMean\n")
for num in sorted(species_order,key=int):
    species=species_order[num]
    o.write(species)
    mean=0
    count=0
    for well in species_well[species]:
        if(well_ct_384[well] is not missing_val_indicator):#I.e., is not missing, as defined from before
            o.write("\t" + well_ct_384[well])
            mean=mean+float(well_ct_384[well])
            count=count+1
        else:
            o.write("\t" + "OMIT")
    o.write("\t"+str(mean/count)+"\n")#Write the mean.... still figuring out a way to 
o.close()


