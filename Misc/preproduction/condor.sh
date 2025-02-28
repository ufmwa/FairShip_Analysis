#!/bin/bash

#######################################################################################
source /cvmfs/ship.cern.ch/24.10/setUp.sh 

#alienv load FairShip/latest-master-release > config_<version>.sh
source /afs/cern.ch/user/a/anupamar/HTCondor/configfiles/config_ECN3_2024.sh
echo 'config sourced'

#######################################################################################

python /afs/cern.ch/user/a/anupamar/Analysis/preproduction/muon_weights.py > preproduction_stat.txt
