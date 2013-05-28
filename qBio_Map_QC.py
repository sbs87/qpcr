#! /Library/Frameworks/Python.framework/Versions/2.7/bin/python
import sys
import bmath
from pylab import *

#Set up defs

#Def for reading files, general. Reads filename into vector. Still need to figure out header issue
def readfile(file_name):
    f=open(file_name)
    lines = f.readlines()
    headers = lines[0].rstrip().split("\t")
    vectors = [line.rstrip().split("\t") for line in lines]
    f.close()
    return(vectors) #Returns entire file contents as series of vectors

#Def for reading in each 384 well, ct and flag value. This assumes input file is correctly formatted, which means col1=Well, col2=ct and col3=t/f regarding ct removal from plot
def read_384(well_species_map,file_name,man_fill_wells,other_wells,controls,control_rep):
    mappings={}
    ommitted=[]
    file=readfile(file_name)
    for line_384 in file:
        well=line_384[0]
        ct_val=line_384[1]
        flag_omit=line_384[2]
        if flag_omit == "false":
            flag_omit=""
        elif flag_omit=="true":
            flag_omit="X"
            ommitted.append(well)
        flag_man=""
        if well in man_fill_wells:
            flag_man="*"
        flag_other=""
        if well in other_wells:
            flag_other="#"
        species=well_species_map[well][0]
        qc_control_val="NA"
        qc_control_rep="NA"
        if species in controls:
            qc_control_val=controls[species]
            qc_control_rep=control_rep[well]
        mapping={well:{"ct":ct_val,"flags":{"man_fill":flag_man,"ommitted":flag_omit,"other":flag_other},"species":species,"Order":well_species_map[well][1],"Rep#":well_species_map[well][2],"QC_control":qc_control_val,"Control_Rep#":qc_control_rep}}
        mappings.update(mapping)
    return(mappings,ommitted)
    
def read_control_mapping(file_name):
    file=readfile(file_name)
    mapping={}
    order={}
    for line in file:
        control_sub,control_main,number=line
        mapping.update({control_sub:control_main})
        if number !="dummy":
            number=int(number)
            order.update({number:control_main})
    return mapping,order
    
def read_well_species_map(file_name):
    file=readfile(file_name)
    mapping={}
    species_order_map={}
    for well_species in file:
        order=int(well_species[0])
        species=well_species[1]
        w=[well_species[2],well_species[3],well_species[4],well_species[5]]
        if species not in species_order_map:
            species_order_map[order]=species
        rep=1
        for well in w:
            mapping.update({well:[species,order,str(rep)]})
            rep=rep+1
    return(mapping,species_order_map)
    
def unique_species_list(well_species_map):
    seen=[]
    for well in well_species_map:
        if well_species_map[well][0] not in seen:
            seen.append(well_species_map[well][0])
    return seen

def format_flags(flags):
    flag_man=flags["man_fill"]
    flag_omit=flags["ommitted"]
    flag_other=flags["other"]
    formated_flags=str(flag_man+flag_omit+flag_other)
    return(formated_flags)
    
def calculate_stats(cts):
    mean=bmath.calc_mean(cts)
    stdev=bmath.calc_stdev(cts)
    return mean,stdev

def map_and_stat(mapped_wells,well_species_map,order,controls):
    species_list=unique_species_list(well_species_map)
    well_96={}
    control={}
    for spec_key in species_list:
        well_96.update({spec_key:{"cts":[None,None,None,None],"pass":[],"stats":[]}})
        if spec_key in controls:
            control.update({controls[spec_key]:{"cts":["","","","","","","",""],"pass":[],"stats":[]}})
    for well in mapped_wells:
        species=mapped_wells[well]["species"]
        ct=mapped_wells[well]["ct"]
        rep=int(mapped_wells[well]["Rep#"])
        formated_flags=format_flags(mapped_wells[well]["flags"])
        if ((mapped_wells[well]["flags"]["ommitted"]=="")and(mapped_wells[well]["flags"]["other"]=="")):
            well_96[species]["pass"].append(ct)
            if(mapped_wells[well]["QC_control"] != "NA"):
                control[controls[species]]["pass"].append(ct)
        well_96[species]["cts"].pop(rep-1)
        well_96[species]["cts"].insert(rep-1,str(ct+formated_flags))
        if species in controls:
            rep_control=int(mapped_wells[well]["Control_Rep#"])
            control[controls[species]]["cts"].pop(rep_control-1)
            control[controls[species]]["cts"].insert(rep_control-1,str(ct+formated_flags))
    for species in well_96:
        cts=well_96[species]["cts"]
        well_96[species]["stats"]=calculate_stats(well_96[species]["pass"])
    for acontrol in control:
        cts=control[acontrol]["cts"]
        control[acontrol]["stats"]=calculate_stats(control[acontrol]["pass"])
    return well_96, control

def output_96(well_96,order, outfile):
    o=open(outfile,'w')
    o.write("Species\tRep1\tRep2\tRep3\tRep4\tMean\tSTDev\n")
    for num in sorted(order):
        species=order[num]
        o.write(species)
        for ct in well_96[species]["cts"]:
            o.write("\t"+ct)
        o.write("\t"+str(well_96[species]["stats"][0])+"\t"+str(well_96[species]["stats"][1])+"\n")
    o.write("Ct vals marked with an X (ommited due to non ideal dRn plot) or # (removed as visual outlier)\nare NOT included in mean or standard deviation calculation. Ct vals with * were manually filled.")
    o.close()  

def read_ct_thres(file_name):
    file=readfile(file_name)
    mapping={}
    for line in file:
        control, low, high=line
        mapping.update({control:{"low":float(low),"high":float(high),"std":1}})
    return mapping

def QC_call(control,mean, stdev,threshold,qc_score_code):
    mean_cal=1*(mean>threshold[control]["low"]) + 2*(mean<threshold[control]["high"])
    stdev_call=4*(stdev<threshold[control]["std"])
    score_to_text,dummy=read_control_mapping(qc_score_code)
    call=score_to_text[str(mean_cal+stdev_call)]
    return call

def output_QC(well_96,order, outfile,flagged,threshold,qc_score_code):
    o=open(outfile,'w')
    o.write("Control\tRep1\tRep2\tRep3\tRep4\tRep5\tRep6\tRep7\tRep8\tMean\tSTDev\tQC_Call\n")
    for num in sorted(order):
        control=order[num]
        o.write(control)
        mean=well_96[control]["stats"][0]
        stdev=well_96[control]["stats"][1]
        call=QC_call(control,mean,stdev,threshold,qc_score_code)
        for ct in well_96[control]["cts"]:
            o.write("\t"+ct)
        o.write("\t"+str(mean)+"\t"+str(stdev)+"\t"+call+"\n")
    o.write("Ct vals marked with an X (ommited due to non ideal dRn plot) or # (removed as visual outlier) are NOT included in mean or standard\ndeviation calculation. Ct vals with * were manually filled.\n")
    o.write("Removed wells:\n")
    for flag in flagged:
        o.write(flag+":\n")
        for well in flagged[flag]:
            o.write("\t"+well)
        o.write("\n")
    o.close() 

def generate_plot(well_96,filename):
    data=[]
    for species in sorted(well_96):
        cts_str=well_96[species]["pass"]
        cts=[]
        for ct in cts_str:
            cts.append(float(ct))
        data.append(cts)
    figure()
    matplotlib.rcParams.update({'font.size': 8})
    box=boxplot(data)
    title(filename)
    xlabel("QC_Controls")
    ylabel("Ct values")
    xticks(range(1,len(sorted(well_96))+1),sorted(well_96))
    savefig(filename)

#Variables
source_vars,order=read_control_mapping(sys.argv[5])
mapping_file=source_vars["mapping_file"] #Location of 384-96 well mapping file
results_file=sys.argv[1] #Results file name
outfile=sys.argv[2] #Output file name
outfile_control=sys.argv[3] #Output QC report file name
plotfile=sys.argv[4] #QC plot file name
control_map,control_order=read_control_mapping(source_vars["control_map"]) #This reads the control mapping file and outputs a hash that maps the controls as well as their original order
control_rep,dummy=read_control_mapping(source_vars["control_rep"]) #This file reads the replicate number for controls, since in the controls case, there are 8 instead of 4 replicates

well_species_map,order=read_well_species_map(mapping_file) #This maps the 384 wells to each species as well as retains thier order
man_fill_wells=[] #A list of manually filled wells
man_fill_wells=raw_input("Enter manually filled wells:") #As for manually filled wells

other_wells=raw_input("Enter additional wells to be removed BEFORE QC:") #Ask for any other filled wells
other_wells=other_wells.split(",")

threshold=read_ct_thres(source_vars["ct_thres_mapping"])
qc_score_code=source_vars["qc_score_code"]

flags=True #Initialize the flags var
while(flags): #Read in flags until no more are defined
    mapped_wells,omit_wells=read_384(well_species_map,results_file,man_fill_wells,other_wells,control_map,control_rep) #Read in the 384 wells, their Ct vals and ommitted flags
    well_96, control=map_and_stat(mapped_wells,well_species_map,order,control_map) #map these wells to their 96 well counterpart and calcualte mean and stdev
    output_96(well_96,order,outfile) #output the resulting mapping and statistics to perform visual QC
    flagged_wells={"manual":man_fill_wells,"other":other_wells,"omit":omit_wells} #Keep a list of flagged wells for later
    output_QC(control,control_order,outfile_control,flagged_wells,threshold,qc_score_code)  #Also output the QC for the control wells (a different function)
    generate_plot(control,plotfile) #Generate a box plot of the control wells to visualize outliers
    wells=raw_input("Enter additional wells to be removed:") #After viewing outliers and things that don't seem right, enter additional wells to be removed
    if (wells=="" or wells=="done"): #If no wells were entered or done was entered, exit loop and output final fils
        print "Removed wells:"
        print [other_wells,omit_wells]
        print "Manually filled wells:"
        print man_fill_wells
        flags=False
    else:
    	wells=wells.split(",")
    	for well in wells:
    		other_wells.append(well)
        #Otherwise, keep adding flags/wells to be removed

