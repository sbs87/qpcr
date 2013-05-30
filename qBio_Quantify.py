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

def agreement_table(dCts,ra, intersect):
	agreement={}
	for species in intersect:
		val_qpcr=dCts[species]
		val_ht=ra[species]
		val_qpcr_bool=val_qpcr!="NP"
		val_ht_bool=float(val_ht)>0
		category=1*val_qpcr_bool+2*val_ht_bool
		agreement.update({species:{"qpcr":val_qpcr,"ht":val_ht,"agree_cat":category}})
		#print (species+"\t"+str(val_qpcr)+"\t"+str(val_ht)+"\t"+str(category))
	return agreement#,summary_table #Species val16 valq category. Seprate, formatted tablle

data_filename=sys.argv[1]
thres_filename=sys.argv[2]
quant_outfile=sys.argv[3]
bacteria1_outfile=str(quant_outfile+"bacteria1")
bacteria2_outfile=str(quant_outfile+"bacteria2")
humanHBB_outfile=str(quant_outfile+"humanHBB")
humanGAPDH_outfile=str(quant_outfile+"humanGAPDH")
bacterial_outfile=str(quant_outfile+"bacterial")
human_outfile=str(quant_outfile+"human")
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

bacterial1_out=open(bacteria1_outfile,'w')
bacterial2_out=open(bacteria2_outfile,'w')
humanHBB_out=open(humanHBB_outfile,'w')
humanGAPDH_out=open(humanGAPDH_outfile,'w')
bacterial_out=open(bacterial_outfile,'w')
human_out=open(human_outfile,'w')
control_list.append("bacterial")
control_list.append("human")
#agreement_out=open("agreement_out",'w')
outstreams={}
for control in control_list:
	outstreams.update({control:{"agree":open("agreement_out",'w')}})
#outstream={"Pan Bacteria 1":{"quant":bacterial1_out,"agree":open("agreement_out",'w')}}    

#Calc Ctthres-Ct controk for given species, taking into account constraints/threshold cutoffs
for control in control_list:
    ddCt={}
    Ct_control=Ct_controls[control][0]#FIX.. i.e., include a bacterial and human control like before... add this to the list of controls and possible write streams
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
    if(control=="Pan Bacteria 1"):
        out=bacterial1_out
    elif(control=="Pan Bacteria 2"):
        out=bacterial2_out
    elif(control=="Hs/MmGAPDH"):
        out=humanGAPDH_out
    elif(control=="Hs/Mm HBB"):
    	out=humanHBB_out
    elif(control=="bacterial"):
    	out=bacterial_out
    elif(control=="human"):
    	out=human_out
    out.write(control+"\nSpecies\tdCt\tPresent\n")
    for species in sorted(ddCt):#,key=lambda dCtval: dCtval[1]):
        out.write(species+"\t"+str(ddCt[species])+"\t"+str(ddCt[species]!="NP")+"\n")
    ra=read_HT_method(rRNA16s_filename)
    xy=create_xy(ddCt,ra,set(ddCt.keys())&set(ra.keys()))
    agreement=agreement_table(ddCt,ra,set(ddCt.keys())&set(ra.keys()))
    blah=outstream["Pan Bacteria 1"]["agree"]
    for species in agreement:
    	blah.write(species+"\t"+str(agreement[species]["qpcr"])+"\t"+str(agreement[species]["ht"])+"\t"+str(agreement[species]["agree_cat"])+"\n")
    #plot_experiments(xy)
    
bacterial1_out.close()
bacterial2_out.close()
humanHBB_out.close()
humanGAPDH_out.close()
human_out.close()
bacterial_out.close()

#Calculate correlation with 16S/metadata.Plot 

