#!/bin/bash
#######################################################################################

#ClusterId= $1
#ProcId=    $2
#Comment=   $3
#runnumber= $4

#######################################################################################
source /cvmfs/ship.cern.ch/24.10/setUp.sh 
#alienv load FairShip/latest-master-release > config_<version>.sh
source /afs/cern.ch/user/a/anupamar/HTCondor/configfiles/config_ECN3_2024.sh
echo 'config sourced'
#######################################################################################

python /afs/cern.ch/user/a/anupamar/Analysis/SBTGeometry_optimisation/adapted_code/MuonBack_hitrates.py 
