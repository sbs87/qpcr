import bmath
import sys

def read_ct(names,list):
    map={}
    for name in names:
        species=name[0]
        ct=float(name[1])
        if(species in list):
            map.update({species:ct})
    return map
    
def generate_control_means(Ct_controls):
    control_means={}
    control_means["bacterial"]=bmath.calc_mean([Ct_controls["Pan Bacteria 1"],Ct_controls["Pan Bacteria 2"]])
    control_means["human"]=bmath.calc_mean([Ct_controls["Hs/MmGAPDH"],Ct_controls["Hs/Mm HBB"]])
    return control_means

data_filename=sys.argv[1]
thres_filename=sys.argv[2]
bacteria_outfile=sys.argv[3]
human_outfile=sys.argv[4]

data_file=bmath.readfile(data_filename)
thres_file=bmath.readfile(thres_filename)

control_list=["Pan Bacteria 1","Pan Bacteria 2","Hs/MmGAPDH","Hs/Mm HBB"]
other_control_list=["PPC","Pan Aspergillus/Candida"]
sample_list=[]

for val in data_file:
    sample=val[0]
    ct=val[1]
    if ((sample not in control_list)&(sample not in other_control_list)):
        sample_list.append(sample)

#Read in Ct mean vals for each species
Ct_species=read_ct(data_file,sample_list) #Sample:ct_val
Ct_controls=read_ct(data_file,control_list)
Ct_controls=generate_control_means(Ct_controls)
Ct_thresholds=read_ct(thres_file,sample_list)

bacterial_out=open(bacteria_outfile,'w')
human_out=open(human_outfile,'w')
#Calc Ctthres-Ct controk for given species, taking into account constraints/threshold cutoffs
for control in ["bacterial","human"]:
    ddCt={}
    Ct_control=Ct_controls[control]
    for species in Ct_species:
        dCt=Ct_species[species]-Ct_control
        dCt_cutoff=Ct_thresholds[species]-Ct_control
        if(dCt<dCt_cutoff):
            ddCt[species]=dCt
        else:
            ddCt[species]="NP"
    #Output ranking, relative abundance, and whether species is "there" (binary)
    if(control=="bacterial"):
        out=bacterial_out
    elif(control=="human"):
        out=human_out
    out.write("Species\tdCt\tPresent\n")
    for species in sorted(ddCt):#,key=lambda dCtval: dCtval[1]):
        out.write(species+"\t"+str(ddCt[species])+"\t"+str(ddCt[species]!="NP")+"\n")
human_out.close()
bacterial_out.close()
#Calculate correlation with 16S/metadata.Plot 

