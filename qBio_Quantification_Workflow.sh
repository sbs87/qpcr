#Get data off of flash drive, copy to local and ssh

#Define the input and output files for result data (preformatted) 
ROOT_LOCAL="stsmith"
INFILE=$1 #may want to redefine so that all you have to input WXDY and MMDDYY
OUTFILE=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_mappedCtVals/g' `
OUTFILE_CONTROL=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_mappedCtVals_QC/g' `
OUTFILE_PRE=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt//g' `
#OUTFILE_QUANT_HUM=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_HumanQuantified.xls/g' `
PLOTFILE=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_QCplot/g' `
SOURCEFILE="/Users/$ROOT_LOCAL/bin/qpcr/source_file"
HT_FILE=$2
HT_TYPE=$3

#cut appropriate fields from results file
#awk '{print $3 "\t" $7" \t" $10}' $INFILE | sed -n 10,394p > temp_formatted_results

#python /Users/$ROOT_LOCAL/bin/qpcr/qBio_Map_QC.py temp_formatted_results $OUTFILE $OUTFILE_CONTROL $PLOTFILE $SOURCEFILE

#cp $OUTFILE "$OUTFILE.xls"
#cp $OUTFILE_CONTROL "$OUTFILE_CONTROL.xls"
#Open results files
#open "$OUTFILE.xls"
#open "$OUTFILE_CONTROL.xls"
#open "$PLOTFILE.png"

cut -f1,6,7 $OUTFILE | sed -n 2,85p > temp_sampleCt
cut -f1,10,11 $OUTFILE_CONTROL | sed -n 2,7p > temp_controlCt
cat temp_sampleCt temp_controlCt > Ct_summary_input
mkdir HT_Compare_$HT_TYPE #THIS IS NOT OMPTIMAL. CHANGE. 
mkdir Quantification
python /Users/$ROOT_LOCAL/bin/qpcr/qBio_Quantify.py Ct_summary_input "/Users/$ROOT_LOCAL/bin/qpcr/dCt_threshold" $OUTFILE_PRE $HT_FILE $HT_TYPE

#Clean up
#rm temp_formatted_results
rm temp_sampleCt
rm temp_controlCt


