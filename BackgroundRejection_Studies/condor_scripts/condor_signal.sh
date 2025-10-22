#!/bin/bash
#######################################################################################

#ClusterId=	$1
#ProcId=$2
JOB=$3
EOSDIR=$4
SCRIPTDIR=$5
KEYWORD=$6 #mumuv or mupi or erho or murho 

#######################################################################################
source /cvmfs/ship.cern.ch/24.10/setUp.sh 

export PYTHONPATH=/eos/experiment/ship/user/anupamar/NN_data/ext_pkgs:$PYTHONPATH
echo "Extra packages ready from eos"

source /afs/cern.ch/user/a/anupamar/HTCondor/configfiles/config_ECN3_2024.sh #alienv load FairShip/latest-master-release > config_<version>.sh
echo 'config sourced'
#######################################################################################

python "$SCRIPTDIR/BackgroundRejection_Studies/run_${KEYWORD}.py" -i "$JOB"
if [ $? -ne 0 ]; then
    echo "ERROR: Step failed. Exiting script."
    exit 1
fi

OUTPUTDIR="$EOSDIR/signalEventCalc/${KEYWORD}/$JOB"
mkdir -p "$OUTPUTDIR"

xrdcp selectionparameters_*.root root://eospublic.cern.ch/"$OUTPUTDIR"/
xrdcp selection_summary_*.csv root://eospublic.cern.ch/"$OUTPUTDIR"/

rm selectionparameters_*.root selection_summary_*.csv 
#######################################################################################