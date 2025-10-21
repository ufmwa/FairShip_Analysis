#!/bin/bash
#######################################################################################

#ClusterId=	$1
#ProcId=$2
JOB=$3

#######################################################################################
source /cvmfs/ship.cern.ch/24.10/setUp.sh 

export PYTHONPATH=/eos/experiment/ship/user/anupamar/NN_data/ext_pkgs:$PYTHONPATH
echo "Extra packages ready from eos"

source /afs/cern.ch/user/a/anupamar/HTCondor/configfiles/config_ECN3_2024.sh #alienv load FairShip/latest-master-release > config_<version>.sh
echo 'config sourced'
#######################################################################################

python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/run_neuDIS.py -i "$JOB" --leptonrho
if [ $? -ne 0 ]; then
    echo "ERROR: Job failed. Exiting script for safe rerun later"
    exit 1
fi    

mkdir -p "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/neuDIS_leptonrho/$JOB"

OUTPUTDIR="/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/neuDIS_leptonrho/$JOB"

xrdcp -f selectionparameters_*.root root://eospublic.cern.ch/"$OUTPUTDIR"/
xrdcp -f selection_summary_*.csv root://eospublic.cern.ch/"$OUTPUTDIR"/
#xrdcp analysis_output.txt root://eospublic.cern.ch/"$OUTPUTDIR"/

rm -f selectionparameters_*.root selection_summary_*.csv #analysis_output.txt
#######################################################################################


python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/run_neuDIS.py -i "$JOB" --fullreco
if [ $? -ne 0 ]; then
    echo "ERROR: Job failed. Exiting script for safe rerun later"
    exit 1
fi    
mkdir -p "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/neuDIS_fullreco/$JOB"

OUTPUTDIR="/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/neuDIS_fullreco/$JOB"

xrdcp -f selectionparameters_*.root root://eospublic.cern.ch/"$OUTPUTDIR"/
xrdcp -f selection_summary_*.csv root://eospublic.cern.ch/"$OUTPUTDIR"/

rm -f selectionparameters_*.root selection_summary_*.csv #analysis_output.txt

#######################################################################################

python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/run_neuDIS.py -i "$JOB" --partialreco
if [ $? -ne 0 ]; then
    echo "ERROR: Job failed. Exiting script for safe rerun later"
    exit 1
fi    
mkdir -p "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/neuDIS/$JOB"

OUTPUTDIR="/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/neuDIS/$JOB"

xrdcp -f selectionparameters_*.root root://eospublic.cern.ch/"$OUTPUTDIR"/
xrdcp -f selection_summary_*.csv root://eospublic.cern.ch/"$OUTPUTDIR"/

rm -f selectionparameters_*.root selection_summary_*.csv #analysis_output.txt
#######################################################################################
