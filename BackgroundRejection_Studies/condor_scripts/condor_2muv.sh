#!/bin/bash
#######################################################################################

#ClusterId=	$1
#ProcId=$2
JOB=$3

#######################################################################################
source /cvmfs/ship.cern.ch/24.10/setUp.sh 
source /afs/cern.ch/user/a/anupamar/HTCondor/configfiles/config_ECN3_2024.sh #alienv load FairShip/latest-master-release > config_<version>.sh
echo 'config sourced'
#######################################################################################

python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/run_partialreco.py -i "$JOB"

mkdir -p "/eos/experiment/ship/user/anupamar/BackgroundStudies/2muv_EventCalc/$JOB"

OUTPUTDIR="/eos/experiment/ship/user/anupamar/BackgroundStudies/2muv_EventCalc/$JOB"

xrdcp selectionparameters_*.root root://eospublic.cern.ch/"$OUTPUTDIR"/
xrdcp selection_summary_*.csv root://eospublic.cern.ch/"$OUTPUTDIR"/
