#!/bin/bash
#######################################################################################

#ClusterId	=	$1
#ProcId		=	$2
#jobid	=	$3

#######################################################################################
source /cvmfs/ship.cern.ch/24.10/setUp.sh 
source /afs/cern.ch/user/a/anupamar/HTCondor/configfiles/config_ECN3_2024.sh #alienv load FairShip/latest-master-release > config_<version>.sh
echo 'config sourced'
#######################################################################################
python /afs/cern.ch/user/a/anupamar/Analysis/FalseVetoProbability_BasicSBTveto/adapted_code/falsevetoprobability_v3.py --nEvents 100 --jobID "$2"  2>&1 | tee terminal_output.txt

#should always have the whole path for condor to work

mkdir -p /afs/cern.ch/work/a/anupamar/FalseVetoProbability/job_"$2"

mv *.root /afs/cern.ch/work/a/anupamar/FalseVetoProbability/job_"$2"/
mv *.csv /afs/cern.ch/work/a/anupamar/FalseVetoProbability/job_"$2"/
mv *.txt /afs/cern.ch/work/a/anupamar/FalseVetoProbability/job_"$2"/


