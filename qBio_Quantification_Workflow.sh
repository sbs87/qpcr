#Get data off of flash drive, copy to local and ssh

#Define the input and output files for result data (preformatted) 
ROOT_LOCAL="stsmith"
INFILE=$1 #may want to redefine so that all you have to input WXDY and MMDDYY
OUTFILE=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_mappedCtVals.xls/g' `
OUTFILE_CONTROL=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_mappedCtVals_QC.xls/g' `
OUTFILE_QUANT_BAC=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_BacterialLoadQuantified.xls/g' `
OUTFILE_QUANT_HUM=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_HumanQuantified.xls/g' `
PLOTFILE=`echo $INFILE | sed 's/\.sdm-Result_Data\.txt/_QCplot/g' `
SOURCEFILE="/Users/$ROOT_LOCAL/bin/qpcr/source_file"

#cut appropriate fields from results file
awk '{print $3 "\t" $7" \t" $10}' $INFILE | sed -n 10,394p > temp_formatted_results

python /Users/$ROOT_LOCAL/bin/qpcr/qBio_Map_QC.py temp_formatted_results $OUTFILE $OUTFILE_CONTROL $PLOTFILE $SOURCEFILE

#Open results files
open $OUTFILE
open $OUTFILE_CONTROL
open "$PLOTFILE.png"

#cut -f1,6 $OUTFILE | sed -n 2,85p > temp_sampleCt
#cut -f1,10 $OUTFILE_CONTROL | sed -n 2,7p > temp_controlCt
#cat temp_sampleCt temp_controlCt > temp_quant_input

#python /Users/$ROOT_LOCAL/bin/qpcr/qBio_Quantify.py temp_quant_input "/Users/$ROOT_LOCAL/bin/qpcr/dCt_threshold" $OUTFILE_QUANT_BAC $OUTFILE_QUANT_HUM
#open $OUTFILE_QUANT_BAC
#open $OUTFILE_QUANT_HUM
#Clean up
rm temp_formatted_results
#rm temp_sampleCt
#rm temp_controlCt
#rm temp_quant_input

#Things to add:
#-Only incorporating acceptable Ct vals based on stdev
#-Only incorporating certain controls based on avg vals

