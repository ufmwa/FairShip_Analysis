#!/bin/bash
#######################################################################################

#ClusterId=	$1
#ProcId=$2
JOB=$3
OUTDIR=$4
SCRIPTDIR=$5

#--------------------------------------------------------------------------------------

# source /cvmfs/ship.cern.ch/24.10/setUp.sh 

# export PYTHONPATH=/eos/experiment/ship/user/anupamar/NN_data/ext_pkgs:$PYTHONPATH
# echo "Extra packages ready from eos"

# source /afs/cern.ch/user/a/anupamar/HTCondor/configfiles/config_ECN3_2024.sh #alienv load FairShip/latest-master-release > config_<version>.sh
# echo 'config sourced'

#######################################################################################

python "$SCRIPTDIR/BackgroundRejection_Studies/run_neuDIS.py" -i "$JOB" --leptonrho
if [ $? -ne 0 ]; then
    echo "ERROR: Job failed. Exiting script for safe rerun later"
    exit 1
fi    

mkdir -p "$OUTDIR/neuDIS/leptonrho/$JOB"

OUTPUTDIR="$OUTDIR/neuDIS/leptonrho/$JOB"

xrdcp selectionparameters_*.root root://eospublic.cern.ch/"$OUTPUTDIR"/
xrdcp selection_summary_*.csv root://eospublic.cern.ch/"$OUTPUTDIR"/

rm selectionparameters_*.root selection_summary_*.csv 

#--------------------------------------------------------------------------------------

python "$SCRIPTDIR/BackgroundRejection_Studies/run_neuDIS.py" -i "$JOB" --fullreco
if [ $? -ne 0 ]; then
    echo "ERROR: Job failed. Exiting script for safe rerun later"
    exit 1
fi    
mkdir -p "$EOSDIR/neuDIS_fullreco/$JOB"

OUTPUTDIR="$EOSDIR/neuDIS_fullreco/$JOB"

xrdcp selectionparameters_*.root root://eospublic.cern.ch/"$OUTPUTDIR"/
xrdcp selection_summary_*.csv root://eospublic.cern.ch/"$OUTPUTDIR"/

rm selectionparameters_*.root selection_summary_*.csv 

#--------------------------------------------------------------------------------------

python "$SCRIPTDIR/BackgroundRejection_Studies/run_neuDIS.py" -i "$JOB" --partialreco
if [ $? -ne 0 ]; then
    echo "ERROR: Job failed. Exiting script for safe rerun later"
    exit 1
fi    
mkdir -p "$EOSDIR/neuDIS/$JOB"

OUTPUTDIR="$EOSDIR/neuDIS/$JOB"

xrdcp selectionparameters_*.root root://eospublic.cern.ch/"$OUTPUTDIR"/
xrdcp selection_summary_*.csv root://eospublic.cern.ch/"$OUTPUTDIR"/

rm selectionparameters_*.root selection_summary_*.csv 

#######################################################################################
