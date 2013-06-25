#--------------------------------------------------
#----qBio_QQunatification_Workflow.sh-------------
#-------Created by Steven Smith -------------------
#-------Modified June 6, 2013----------------------
#--------------------------------------------------
#-------This serves as a wrapper to the QC/mapping and quantification scripts 
#--------(implemented in python) for the qBiomarker study. It takes as input 3 arguments: 
#---------File to "Results" file for qBio study (form thermocyler), pointer to high throughput method 
#--------(for comparison) and labrel for ht method (ie 16s ot meta).. These last two argyments SHOULD be 
#---------optional, but currently aren't.-----------
#--------------------------------------------------
#--------------------------------------------------

#Note: modify sourcefile pointer (SOURCEFILE) and contents, as wella s ROOT_LOCAL as needed to fit local file system. 

#Define the input and output files for result data (preformatted) 
ROOT_LOCAL="stsmith"
INFILE=$1 #should be the "result_data" file
OUTFILE=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_mappedCtVals/g' `
OUTFILE_CONTROL=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_mappedCtVals_QC/g' `
OUTFILE_PRE=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt//g' ` #File name to give to general files use in quant script
PLOTFILE=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_QCplot/g' ` #File name to give to QC plot
SOURCEFILE="/Users/$ROOT_LOCAL/bin/qpcr/source_file" #Pointer to local variable files to use used in scripts later
HT_FILE=$2 #File to 16s or metagenomic data. inthe form species \t value (counts or rel abund)
HT_TYPE=$3 #Either 16s or meta. Defines output stream variables (directory and file names)
HT_CUTOFF=$4 #Cutoff for count or metagenomic calls (present/absent_

#cut appropriate fields from results file, corresponding to well #, Ct val, and t/f hidden from Rn plot
#awk '{print $3 "\t" $7" \t" $10}' $INFILE | sed -n 10,394p > temp_formatted_results

#-----Run QC and mapping script, inputting the infile-formaated name for control and samples, as well as a pointer to the plot and source files

#python /Users/$ROOT_LOCAL/bin/qpcr/qBio_Map_QC.py temp_formatted_results $OUTFILE $OUTFILE_CONTROL $PLOTFILE $SOURCEFILE

#Do some file manipulation so that they can be opened with Excel
#cp $OUTFILE "$OUTFILE.xls"
#cp $OUTFILE_CONTROL "$OUTFILE_CONTROL.xls"

#Open results files so that user can format, save, print and paste into lab notebook
#open "$OUTFILE.xls"
#open "$OUTFILE_CONTROL.xls"
#open "$PLOTFILE.png"

#Format the QC'd, ampped files for input into the quant script/ The cut, cat takes only the mean ct vals from the samples and control files generated above
cut -f1,6,7 $OUTFILE | sed -n 2,85p > temp_sampleCt #only species:mean ct, stdev from samples
cut -f1,10,11 $OUTFILE_CONTROL | sed -n 2,7p > temp_controlCt #only species: mean ct, stdev from controls
cat temp_sampleCt temp_controlCt > Ct_summary_input #cat file of both to be read into quant script

#Make directories for use in quant script. Note that these are hardcoded in the quant script. Also assumed to be in study directory
mkdir HT_Compare_$HT_TYPE #THIS IS NOT OMPTIMAL. CHANGE. 
mkdir Quantification

#----Run quant script. Takes the Ct mean/stdev formatted file as input, along with hardcoded threshold vals, the outut file formatted, and HT file info. HT_FILE is pointer to 16s or metagenomic data. HT_TYPE is either 16s or meta. Read from command line
python /Users/$ROOT_LOCAL/bin/qpcr/qBio_Quantify.py Ct_summary_input "/Users/$ROOT_LOCAL/bin/qpcr/dCt_threshold" $OUTFILE_PRE $HT_FILE $HT_TYPE $HT_CUTOFF

#Clean up, remove temporary files. 
#rm temp_formatted_results
rm temp_sampleCt
rm temp_controlCt

#----END OF SCRIPT----


