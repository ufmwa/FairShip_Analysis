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


python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/run_erho.py -i "$JOB"

mkdir -p "/eos/experiment/ship/user/anupamar/BackgroundStudies/erho_EventCalc/$JOB"

OUTPUTDIR="/eos/experiment/ship/user/anupamar/BackgroundStudies/erho_EventCalc/$JOB"

xrdcp selectionparameters_*.root root://eospublic.cern.ch/"$OUTPUTDIR"/
xrdcp selection_summary_*.csv root://eospublic.cern.ch/"$OUTPUTDIR"/
#only 1GeV representative here