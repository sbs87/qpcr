#Get data off of flash drive, copy to local and ssh

#Define the input and output files for result data (preformatted) 
#INFILE=`echo /Users/stevensmith/Documents/IGS/qPCR/$2/qBiomarker_Sub19_$1_$2.sdm-Result_Data.txt`
INFILE=$1 #may want to redefine so that all you have to input WXDY and MMDDYY
OUTFILE=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_mappedCtVals.xls/g' `
OUTFILE_CONTROL=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_mappedCtVals_QC.xls/g' `
PLOTFILE=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_QCplot/g' `
SOURCEFILE="/Users/stevensmith/bin/qpcr/source_file"
#cut appropriate fields from results file
awk '{print $3 "\t" $7" \t" $10}' $INFILE | sed -n 10,394p > temp_formatted_results

python /Users/stevensmith/bin/qpcr/qBio_Map_QC_Quant.py temp_formatted_results $OUTFILE $OUTFILE_CONTROL $PLOTFILE $SOURCEFILE

#Open results files
open $OUTFILE
open $OUTFILE_CONTROL
open "$PLOTFILE.png"

#Clean up
rm temp_formatted_results
