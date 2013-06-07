#--------------------------------------------------------
#-------copy_qBio_files.sh-------------------------------
#--------recreated June 7, 2013. (orignal was lost for some readon... checked on gitbub, could not fiind
#--------steve smith
#---------------------------------------------------------
#-------The following shell script copies the study folder originally copied from the qPCR room PC to the flash drive into the local directory. T
#--------The script takes as input the study directory, i.e., MMDDYY (060313, for example): the location of the folder wihtin the flash drive. 
#--------If multiple runs were done on the same study day, this will automatically be copied since MMDDYY->Study 1 and the same MMDDYY-> study 2
#--------------------------------------------------------
#--------------------------------------------------------

#Define variables (note they are hardcoded... flash drive and local drives assumed not to change, but may need to update time to time
STUDY_DIR=$1 #The directory folder, MMDDYY
FLASH_DIR="/Volumes/IGS_STEVE/ABS_7900HT/$STUDY_DIR"
LOCAL_DIR="/Users/stsmith/Documents/Results/qPCR/"

#Copy files
echo "Copying $FLASH_DIR to $LOCAL_DIR ..."
cp -r $FLASH_DIR $LOCAL_DIR
echo "complete!"

#-----END OF SCRIPT-----