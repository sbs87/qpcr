import bmath
import sys
from pylab import *

#------------------
#-------Define defs
#------------------

#----Maps species:[Mean ct, Stdev] in a names list. Does this only for species in the 'list' list. Assumes names is formatted species\tmeanCt\tstdev
def read_ct(names,list):
    map={}
    for name in names: #Usually either a list of case species or control species
        species=name[0]
        ct=float(name[1])
        stdev=float(name[2])
        if(species in list): #only map species if it is in the list, i.e., the control list or regular species list
            map.update({species:[ct,stdev]})
    return map
   
#-----Calculates control Ct val means for each major control by simply taking the mean of either each bacteria control or each numan control 
def generate_control_means(Ct_controls):
    control_means={}
    control_means["bacterial"]=[bmath.calc_mean([Ct_controls["Pan Bacteria 1"][0],Ct_controls["Pan Bacteria 2"][0]]),"NA"] #control mean for "bacterial" is the mean of pan bact 1 and pan bact 2. Because the returned mapping file is expecting two values, also have to store an "NA" val 
    control_means["human"]=[bmath.calc_mean([Ct_controls["Hs/MmGAPDH"][0],Ct_controls["Hs/Mm HBB"][0]]),"NA"] #control mean for "human" is the mean of HBB and GAPDH. Because the returned mapping file is expecting two values, also have to store an "NA" val
    return control_means #The result is bacterial:[mean,"NA"], human:[mean,"NA"]
 
#-----Outdated def, used to create the x and y vectors for the eventual plot of HT_method vs qPCR.     
def create_xy(cts,rel_abundance, intersect):
    xy={}
    for species in intersect:
        if ((cts[species]=="")|(cts[species]=="NP")):
            next
        else:
            xy[species]=[cts[species],rel_abundance[species]]
    return xy
#-----Outdated def, used to plot xy vectors of HT data versus qpcr. You'd expect to see a correlation if the methods match, even though you can't directl;y compare the numbers
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

#------Def used to read in the "high throughput" (HT) method data for subsequent comparison to qpcr. Expects a filename to the correctly formatted HT data file
def read_HT_method(filename): 
    file=bmath.readfile(filename) # Read file contents using the bmath (should rename) custom python module. This simply reads all of a file's contents, and puts it into a vector of vectors, line by line
    mapping={}
    for line in file: # Read in each line one by one from the read in file
    	#Note the expected format of the HT file... the first col should be the species and the second should be the abundance or counts of metagenomic or 16s data, respectively. In other words, each time point should be separated out so that only these two fields exist. 
        #Also, the species names should match exactly to those in the qBio species list, or else they can't be compared.
        species=line[0]
        abundance=line[1] #Can also be counts, abundance is just a variable name
        mapping.update({species:abundance})
    return mapping #Returns species:number Note that the species list will be partially non overlapping with the qBio species list.

#Def for generating the agreement tables and summary agreement tables between the qPCR method and selected HT method. 

def agreement_table(dCts,ra, intersect,sample,ht_count_thres): #Expects the qPCR dCts hash (species:number), the "ra" or relative abundance hash (species:number) (could be counts), intersecting species set (species==species), sample type=[16s|meta\ used to generate output formatting, and the thresholds for selecting HT method cutoffs (file)
	agreement={}
	#intialize cats for "switch" statment below. these are agreement categories
	cat1=0
	cat2=0
	cat3=0
	cat0=0
	#For each of the overlapping species between qPCR and the HT method,
	for species in intersect:
		val_qpcr=dCts[species] #dCt value for qpcr, which can by any number or the value "NP"
		val_ht=ra[species] #abundance or count value for HT method, which can be any postiive or zero value
		val_qpcr_bool=val_qpcr!="NP" #Check to see if the qpcr val is "NP"
		val_ht_bool=float(val_ht)>ht_count_thres #Check to see if HT method value passes threshold, set to the same val for all species (currently hardcoded as 2)
		category=1*val_qpcr_bool+2*val_ht_bool #Category code from 0-3, 3 being both adetect species, 0 being both do not detect species, 1 being qpcr detects speces, 2 being ht method detecs species 
		agreement.update({species:{"qpcr":val_qpcr,"ht":val_ht,"agree_cat":category}}) #add the qpcr, ht method vals to hash as well as the catgeory code
		#The following is used to track the summary table at the end, keeping track of how many of each category there was for the entire set
		if(category==0):
		    cat0=cat0+1
		elif(category==1):
		    cat1=cat1+1
		elif(category==2):
		    cat2=cat2+1
		elif(category==3):
		    cat3=cat3+1
	#Format output for sumamry table, a 2x2 table giving the counts for each of the 4 category cases
	summary_table=str("Agreement Summary Table for "+sample+"\t\tqpcr\tqpcr\n\t\t+\t-\nmeta\t+\t"+str(cat3)+"\t"+str(cat2)+"\nmeta\t-\t"+str(cat1)+"\t"+str(cat0))
	#Output botht eh agreement table (species\tqpcr val\tht method vbal\t cat code and the agreement summary table, above
	return agreement,summary_table #Species val16 valq category. Seprate, formatted tablle

#------------------
#-------Define variables
#------------------

#---The following are variables defined by input stream. this will be inside the wrapper script, so these values have either been formatted by the .sh wrapper or are fed directly from command line
data_filename=sys.argv[1] #Pointing to source of input data, formatted by .sh script
thres_filename=sys.argv[2] #Pointing to threshold file, contained in .sh script, hardcoded  
outfile_prefix=sys.argv[3] #formatting to outfile filename prefix, formatted by .sh script as qBio_subID_WXDY_MMDDYY_
rRNA16s_filename=sys.argv[4] #Should be changed to HTMethod_filename, but this points to the HT method data filename, ie, reads or abundance values for 16s or metagenome
HT_type=sys.argv[5] #Passed from command line, this is either 16s or meta and will be used to write output files depending on the HT method 

#The folliing read in the entire filecontnes using a custom bmath package filereader. 
data_file=bmath.readfile(data_filename) 
thres_file=bmath.readfile(thres_filename)

#Define the control and other lists to be read in later. This is hardcoded because it is assumed these controls won't change. 
control_list=["Pan Bacteria 1","Pan Bacteria 2","Hs/MmGAPDH","Hs/Mm HBB"]
other_control_list=["PPC","Pan Aspergillus/Candida"]
sample_list=[]

#------------------
#-------Start main program
#------------------

#----Read in files, values 

#Read in the sample list and it's correspond mean ct val and st dev (not used) (created by wrapper script by cat the vals generated from QC script)
#This is the sample list, not the control list...
for val in data_file:
    sample=val[0]
    ct=val[1] #not used
    stdev=val[2] #not used
    #Check to see if current val is a sample or a control. If it's not in the "special" lists, add it as a sample
    if ((sample not in control_list)&(sample not in other_control_list)):
        sample_list.append(sample)

#Read in Ct mean vals for each species
Ct_species=read_ct(data_file,sample_list) #Assigns only samples (in sample list) their respective mean Ct val based on the data in data_file. Mapping is Sample:[mean ct_val, stdev]
Ct_controls=read_ct(data_file,control_list) #Assings only controls (in control_list) their respective mean Ct val based on data in data file. Mapping is Control[mean ct val, stdev]
Ct_controls_combined=generate_control_means(Ct_controls) #This takes control vals and cals mean for "bacterial" and "human" and updates the control hash
Ct_controls.update(Ct_controls_combined)
Ct_thresholds=read_ct(thres_file,sample_list) #Reads in control Ct and stdev thresholds for later
ht_count_thres=2 #Hardcoded assignment to ht_method cutoff.. CHANGE so that it is not hardcoded and changes depending on absolute counts vs relative abundance

control_list.append("bacterial") #Updates the control list accordlingly
control_list.append("human")

#Read in relative abundance file 
ra=read_HT_method(rRNA16s_filename)

#----Defines the output streams using hash. Multiple, dynamic output streams are defined based on iterative variables

outstreams={}
for control in control_list:
    #The following if statements format the GAPDH and HBB output stream so they don'thave a slash in them (interpreted as a directory path in output stream)
    if(control=="Hs/MmGAPDH"):
        control_fn="Hs_MmGAPDH"
    elif(control=="Hs/Mm HBB"):
        control_fn="Hs_MmHBB"
    else:
        control_fn=control
    #Hash for different output streams. These are defined by the outfilr prefix and the current control in forloop. os=quant:<quantification os>, agree:<agreement os>, "agree_sum" <empty>. The agree_sum key is updated after the for loop
    outstreams.update({control:{"quant":open(str("Quantification/"+outfile_prefix+"_quantified_"+control_fn),'w'),"agree":open(str("HT_Compare_"+HT_type+"/"+outfile_prefix+"_"+HT_type+"_HTagreement_"+control_fn),'w'),"agree_sum":""}})
#Update only the agree_sum os to be formatted as containing the HT_method type (incorportates all controls)
outstreams.update({"all":{"quant":"","agree":"","agree_sum":open(str("HT_Compare_"+HT_type+"/"+outfile_prefix+"_"+HT_type+"_HTagreementSum"),'w')}})
  
#Define the set of species shared by the qPCR and HT_method for use in the "agreement" def. Also output this intersecting set to output. NOTE: THIS SHOULD CHANGE DEPENDING ON THE HT_method USED!!!
intersect=set(Ct_species.keys())&set(ra.keys())
inter=open(str(outfile_prefix+"_agreementSpeciesSet"),'w')
inter.write(str(intersect))
inter.close()

#-----Calulate dCt and output results

#Calc Ctspec-Ct controk for given species, taking into account constraints/threshold cutoffs

for control in control_list:
    ddCt={}
    Ct_control=Ct_controls[control][0]
    for species in Ct_species:
        #Define dCt for the current species/control 
        dCt=Ct_species[species][0]-Ct_control
        #Define cutoff val (note technically that the dCt doesn't have to be calculated, you can define cutoffs before calc dCt
        dCt_cutoff=Ct_thresholds[species][0]-Ct_control
        #There is not stdev cutoff for bacterial or human.... so set to true. Otherwise, figue out if Ct passes stdev cutoff (in addition to Ct cutoff as well)
        if Ct_species[species][1] in ["bacterial","human"]:
        	stdev_pass=True
        else:
        	stdev_pass=Ct_species[species][1]<=Ct_thresholds[species][1]
        #If the species' Ct and stdev pass, then assign the "ddCt" (CHANGE THIS IS NOT ddCt!) to current species. Otherwise, assign a "NP" (not present)
        if((dCt<dCt_cutoff)&(stdev_pass)):
            ddCt[species]=dCt
        else:
            ddCt[species]="NP"
    #Output results
    
    #--Quantification output stream
    outstreams[control]["quant"].write(outfile_prefix+"\t"+outfile_prefix+"\t"+outfile_prefix+"\n"+control+"\t"+control+"\t"+control+"\n")
    for species in sorted(ddCt):#,key=lambda dCtval: dCtval[1]):
        outstreams[control]["quant"].write(species+"\t"+str(ddCt[species])+"\t"+str(ddCt[species]!="NP")+"\n")
    
    #---Agreement output stream
    #Generate agreement formatted outstreams using agreement_table function
    agreement,agreement_sum=agreement_table(ddCt,ra,intersect,outfile_prefix,ht_count_thres)
    #Begin writing to output stream agree (header)
    outstreams[control]["agree"].write(outfile_prefix+"\t"+outfile_prefix+"\t"+outfile_prefix+"\t"+outfile_prefix+"\n"+control+"\t"+control+"\t"+rRNA16s_filename+"\t"+control+"\n")
    #Formatted as species \t qpcr val \t ht val \t category
    for species in sorted(agreement):
        outstreams[control]["agree"].write(species+"\t"+str(agreement[species]["qpcr"])+"\t"+str(agreement[species]["ht"])+"\t"+str(agreement[species]["agree_cat"])+"\n")
    outstreams[control]["agree"].close()
#Output the agreement summary table
outstreams["all"]["agree_sum"].write(agreement_sum)
outstreams["all"]["agree_sum"].close()

#----END OF SCRIPT------