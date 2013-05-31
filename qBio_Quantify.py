import bmath
import sys
from pylab import *

def read_ct(names,list):
    map={}
    for name in names:
        species=name[0]
        ct=float(name[1])
        stdev=float(name[2])
        if(species in list):
            map.update({species:[ct,stdev]})
    return map
    
def generate_control_means(Ct_controls):
    control_means={}
    control_means["bacterial"]=[bmath.calc_mean([Ct_controls["Pan Bacteria 1"][0],Ct_controls["Pan Bacteria 2"][0]]),"NA"]
    control_means["human"]=[bmath.calc_mean([Ct_controls["Hs/MmGAPDH"][0],Ct_controls["Hs/Mm HBB"][0]]),"NA"]
    
    return control_means
def create_xy(cts,rel_abundance, intersect):
    xy={}
    for species in intersect:
        if ((cts[species]=="")|(cts[species]=="NP")):
            next
        else:
            xy[species]=[cts[species],rel_abundance[species]]
    return xy
def plot_experiments(xy):
	x=[]
	y=[]
	for val in xy:
		x.append(xy[val][0])
		y.append(xy[val][1])
	scatter(x,y)
	xlabel('qpcr')
	ylabel('16S')
	title('Method Comparison- qPCR vs 16S')
	grid(True)
	#savefig("test.png")
	show()
def read_HT_method(filename):
    file=bmath.readfile(filename)
    mapping={}
    for line in file:
        species=line[0]
        abundance=line[1]
        mapping.update({species:abundance})
    return mapping

def agreement_table(dCts,ra, intersect,control):
	agreement={}
	cat1=0
	cat2=0
	cat3=0
	cat0=0
	for species in intersect:
		val_qpcr=dCts[species]
		val_ht=ra[species]
		val_qpcr_bool=val_qpcr!="NP"
		val_ht_bool=float(val_ht)>1
		category=1*val_qpcr_bool+2*val_ht_bool
		agreement.update({species:{"qpcr":val_qpcr,"ht":val_ht,"agree_cat":category}})
		if(category==0):
		    cat0=cat0+1
		elif(category==1):
		    cat1=cat1+1
		elif(category==2):
		    cat2=cat2+1
		elif(category==3):
		    cat3=cat3+1
	summary_table=str(control+"\t\tqpcr\tqpcr\n\t\t+\t-\n16s\t+\t"+str(cat3)+"\t"+str(cat2)+"\n16s\t-\t"+str(cat1)+"\t"+str(cat0))
	return agreement,summary_table #Species val16 valq category. Seprate, formatted tablle

data_filename=sys.argv[1]
thres_filename=sys.argv[2]
outfile_prefix=sys.argv[3]
rRNA16s_filename=sys.argv[4]

data_file=bmath.readfile(data_filename)
thres_file=bmath.readfile(thres_filename)

control_list=["Pan Bacteria 1","Pan Bacteria 2","Hs/MmGAPDH","Hs/Mm HBB"]
other_control_list=["PPC","Pan Aspergillus/Candida"]
sample_list=[]

for val in data_file:
    sample=val[0]
    ct=val[1]
    stdev=val[2]
    if ((sample not in control_list)&(sample not in other_control_list)):
        sample_list.append(sample)

#Read in Ct mean vals for each species
Ct_species=read_ct(data_file,sample_list) #Sample:ct_val
Ct_controls=read_ct(data_file,control_list)
Ct_controls_combined=generate_control_means(Ct_controls)
Ct_controls.update(Ct_controls_combined)
Ct_thresholds=read_ct(thres_file,sample_list)

control_list.append("bacterial")
control_list.append("human")

outstreams={}
for control in control_list:
    if(control=="Hs/MmGAPDH"):
        control_fn="Hs_MmGAPDH"
    elif(control=="Hs/Mm HBB"):
        control_fn="Hs_MmHBB"
    else:
        control_fn=control
    outstreams.update({control:{"quant":open(str(outfile_prefix+"_quantified_"+control_fn),'w'),"agree":open(str(outfile_prefix+"_HTagreement_"+control_fn),'w'),"agree_sum":open(str(outfile_prefix+"_HTagreementSum_"+control_fn),'w')}})

ra=read_HT_method(rRNA16s_filename)  
intersect=set(Ct_species.keys())&set(ra.keys())
inter=open(str(outfile_prefix+"_agreementSpeciesSet"),'w')
inter.write(str(intersect))
inter.close()

#Calc Ctthres-Ct controk for given species, taking into account constraints/threshold cutoffs
for control in control_list:
    ddCt={}
    Ct_control=Ct_controls[control][0]
    for species in Ct_species:
        dCt=Ct_species[species][0]-Ct_control
        dCt_cutoff=Ct_thresholds[species][0]-Ct_control
        if Ct_species[species][1] in ["bacterial","human"]:
        	stdev_pass=True
        else:
        	stdev_pass=Ct_species[species][1]<=Ct_thresholds[species][1]
        if((dCt<dCt_cutoff)&(stdev_pass)):
            ddCt[species]=dCt
        else:
            ddCt[species]="NP"
    #Output ranking, relative abundance, and whether species is "there" (binary)
    outstreams[control]["quant"].write(control+"\n")
    for species in sorted(ddCt):#,key=lambda dCtval: dCtval[1]):
        outstreams[control]["quant"].write(species+"\t"+str(ddCt[species])+"\t"+str(ddCt[species]!="NP")+"\n")
    agreement,agreement_sum=agreement_table(ddCt,ra,intersect,control)
    for species in sorted(agreement):
        outstreams[control]["agree"].write(species+"\t"+str(agreement[species]["qpcr"])+"\t"+str(agreement[species]["ht"])+"\t"+str(agreement[species]["agree_cat"])+"\n")
    outstreams[control]["agree_sum"].write(agreement_sum)
    outstreams[control]["agree"].close()
    outstreams[control]["agree_sum"].close()
    #plot_experiments(xy)
    #xy=create_xy(ddCt,ra,set(ddCt.keys())&set(ra.keys()))

#Calculate correlation with 16S/metadata.Plot 

