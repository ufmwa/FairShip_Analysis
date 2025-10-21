#!/bin/bash
#######################################################################################

#ClusterId=	$1
#ProcId=$2
JOB1=$3
JOB2=$4
JOB3=$5
JOB4=$6
JOB5=$7

#######################################################################################
source /cvmfs/ship.cern.ch/24.10/setUp.sh 

export PYTHONPATH=/eos/experiment/ship/user/anupamar/NN_data/ext_pkgs:$PYTHONPATH
echo "Extra packages ready from eos"

source /afs/cern.ch/user/a/anupamar/HTCondor/configfiles/config_ECN3_2024.sh #alienv load FairShip/latest-master-release > config_<version>.sh
echo 'config sourced'
#######################################################################################


for JOB in "$JOB1" "$JOB2" "$JOB3" "$JOB4" "$JOB5"; do
    [ -z "$JOB" ] && continue   # skip empty slot (e.g. last line had fewer than 5)

    echo "===== START $JOB $(date) on $(hostname) ====="
        # Run the analysis
    python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/run_neuDIS.py -i "$JOB" --leptonrho 
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Job failed. Exiting script for safe rerun later"
        exit 1
    fi    
    
    mkdir -p "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/neuDIS_leptonrho/$JOB"
    OUTPUTDIR="/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/neuDIS_leptonrho/$JOB"
    
    xrdcp -f selectionparameters_*.root root://eospublic.cern.ch/"$OUTPUTDIR"/
	xrdcp -f selection_summary_*.csv root://eospublic.cern.ch/"$OUTPUTDIR"/
	#xrdcp -f analysis_output.txt root://eospublic.cern.ch/"$OUTPUTDIR"/

	# Clean pattern outputs so next run starts clean
    rm -f selectionparameters_*.root selection_summary_*.csv #analysis_output.txt

    echo "===== END   $JOB $(date) ====="

done

