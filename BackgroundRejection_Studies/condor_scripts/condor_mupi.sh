#!/bin/bash
#######################################################################################

#ClusterId=	$1
#ProcId=$2
JOB=$3

#######################################################################################
source /cvmfs/ship.cern.ch/24.10/setUp.sh 
################################################################################
# (a) point Python to the extra packages
export PYTHONPATH=/eos/experiment/ship/user/anupamar/NN_data/ext_pkgs:$PYTHONPATH
echo "Extra packages ready from eos"
################################################################################
source /afs/cern.ch/user/a/anupamar/HTCondor/configfiles/config_ECN3_2024.sh #alienv load FairShip/latest-master-release > config_<version>.sh
echo 'config sourced'
#######################################################################################

python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/run_fullreco.py -i "$JOB"
if [ $? -ne 0 ]; then
    echo "ERROR: Step failed. Exiting script."
    exit 1
fi

mkdir -p "/eos/experiment/ship/user/anupamar/BackgroundStudies/mupi_EventCalc/$JOB"

OUTPUTDIR="/eos/experiment/ship/user/anupamar/BackgroundStudies/mupi_EventCalc/$JOB"

xrdcp selectionparameters_*.root root://eospublic.cern.ch/"$OUTPUTDIR"/
xrdcp selection_summary_*.csv root://eospublic.cern.ch/"$OUTPUTDIR"/
