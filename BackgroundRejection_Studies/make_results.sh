#!/bin/bash
#######################################################################################
source /cvmfs/ship.cern.ch/24.10/setUp.sh 
source /afs/cern.ch/user/a/anupamar/HTCondor/configfiles/config_ECN3_2024.sh #alienv load FairShip/latest-master-release > config_<version>.sh
echo 'config sourced'
#######################################################################################

echo "Merging muonDIS now.."
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS --all > muonDIS_all.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS --vesselCase > muonDIS_vesselCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS --heliumCase > muonDIS_heliumCase.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS --all
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS --vesselCase
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS --heliumCase
#condor_scripts/condor_hadd.sh "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/muonDIS" "muonDIS" 
echo "completed."

echo "Merging muonDIS_fullreco now.."
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS_fullreco --all > muonDIS_fullreco_all.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS_fullreco --vesselCase > muonDIS_fullreco_vesselCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS_fullreco --heliumCase > muonDIS_fullreco_heliumCase.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS_fullreco --all
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS_fullreco --vesselCase
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS_fullreco --heliumCase
#condor_scripts/condor_hadd.sh "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/muonDIS_fullreco" "muonDIS_fullreco"
wait
echo "completed."

echo "Merging muonDIS lepton-rho now.."
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS_leptonrho --all > muonDIS_leptonrho_all.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS_leptonrho --vesselCase > muonDIS_leptonrho_vesselCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS_leptonrho --heliumCase > muonDIS_leptonrho_heliumCase.txt 

#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS_leptonrho --all
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS_leptonrho --vesselCase
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS_leptonrho --heliumCase
#condor_scripts/condor_hadd.sh "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/muonDIS_leptonrho" "muonDIS_leptonrho" 
wait
echo "completed."

exit 0

#---------------------scaling neuDIS-------------------------------------------------------------------------

python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/scale_neuDIS.py --neuDIS --all > neuDIS_all.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/scale_neuDIS.py --neuDIS --vesselCase > neuDIS_vesselCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/scale_neuDIS.py --neuDIS --heliumCase > neuDIS_heliumCase.txt 

python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/scale_neuDIS.py --neuDIS_fullreco --all > neuDIS_fullreco_all.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/scale_neuDIS.py --neuDIS_fullreco --vesselCase > neuDIS_fullreco_vesselCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/scale_neuDIS.py --neuDIS_fullreco --heliumCase > neuDIS_fullreco_heliumCase.txt 

python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/scale_neuDIS.py --neuDIS_leptonrho --all > neuDIS_leptonrho_all.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/scale_neuDIS.py --neuDIS_leptonrho --vesselCase > neuDIS_leptonrho_vesselCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/scale_neuDIS.py --neuDIS_leptonrho --heliumCase > neuDIS_leptonrho_heliumCase.txt 



#------------------------------------------------------------------------------------------------------------
echo "Merging mupi now.."
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --mupi --all > mupi_all.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --mupi --vesselCase > mupi_vesselCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --mupi --heliumCase > mupi_heliumCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --mupi --all
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --mupi --vesselCase
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --mupi --heliumCase
#condor_scripts/condor_hadd.sh "/eos/experiment/ship/user/anupamar/BackgroundStudies/mupi_EventCalc" "mupi_EventCalc"
echo "completed."

echo "Merging 2muv now.."
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --mumuv --all > mumuv_all.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --mumuv --vesselCase > mumuv_vesselCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --mumuv --heliumCase > mumuv_heliumCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --mumuv --all
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --mumuv --vesselCase
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --mumuv --heliumCase
#condor_scripts/condor_hadd.sh "/eos/experiment/ship/user/anupamar/BackgroundStudies/2muv_EventCalc" "2muv_EventCalc"
echo "completed."

echo "Merging erho now.."
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --erho --all > erho_all.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --erho --vesselCase > erho_vesselCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --erho --heliumCase > erho_heliumCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --erho --all
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --erho --vesselCase
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --erho --heliumCase
#condor_scripts/condor_hadd.sh "/eos/experiment/ship/user/anupamar/BackgroundStudies/erho_EventCalc" "erho_EventCalc"
echo "completed."

echo "Merging murho now.."
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --murho --all > murho_all.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --murho --vesselCase > murho_vesselCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --murho --heliumCase > murho_heliumCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --murho --all
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --murho --vesselCase
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --murho --heliumCase
#condor_scripts/condor_hadd.sh "/eos/experiment/ship/user/anupamar/BackgroundStudies/murho_EventCalc" "murho_EventCalc"
echo "completed."


echo "Merging neuDIS now .."
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --neuDIS --all > neuDIS_all.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --neuDIS --vesselCase > neuDIS_vesselCase.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --neuDIS --heliumCase > neuDIS_heliumCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --neuDIS --all
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --neuDIS --vesselCase
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --neuDIS --heliumCase
#condor_scripts/condor_hadd.sh "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/neuDIS" "neuDIS"
echo "completed."


echo "Merging neuDIS_fullreco now .."
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --neuDIS_fullreco --all > neuDIS_fullreco_all.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --neuDIS_fullreco --vesselCase > neuDIS_fullreco_vesselCase.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --neuDIS_fullreco --heliumCase > neuDIS_fullreco_heliumCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --neuDIS_fullreco --all
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --neuDIS_fullreco --vesselCase
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --neuDIS_fullreco --heliumCase
#condor_scripts/condor_hadd.sh "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/neuDIS_fullreco" "neuDIS_fullreco"
echo "completed."


echo "Merging neuDIS leptonrho now.."
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --neuDIS_leptonrho --all > neuDIS_leptonrho_all.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --neuDIS_leptonrho --vesselCase > neuDIS_leptonrho_vesselCase.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --neuDIS_leptonrho --heliumCase > neuDIS_leptonrho_heliumCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --neuDIS_leptonrho --all
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --neuDIS_leptonrho --vesselCase
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --neuDIS_leptonrho --heliumCase
#condor_scripts/condor_hadd.sh "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/neuDIS_leptonrho" "neuDIS_leptonrho" 
echo "completed."



echo "Merging muonDIS now.."
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS --all > muonDIS_all.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS --vesselCase > muonDIS_vesselCase.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS --heliumCase > muonDIS_heliumCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS --all
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS --vesselCase
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS --heliumCase
#condor_scripts/condor_hadd.sh "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/muonDIS" "muonDIS" 
echo "completed."

echo "Merging muonDIS_fullreco now.."
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS_fullreco --all > muonDIS_fullreco_all.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS_fullreco --vesselCase > muonDIS_fullreco_vesselCase.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS_fullreco --heliumCase > muonDIS_fullreco_heliumCase.txt 
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS_fullreco --all
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS_fullreco --vesselCase
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS_fullreco --heliumCase
#condor_scripts/condor_hadd.sh "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/muonDIS_fullreco" "muonDIS_fullreco"
wait
echo "completed."

echo "Merging muonDIS lepton-rho now.."
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS_leptonrho --all > muonDIS_leptonrho_all.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS_leptonrho --vesselCase > muonDIS_leptonrho_vesselCase.txt 
#python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/merge.py --muonDIS_leptonrho --heliumCase > muonDIS_leptonrho_heliumCase.txt 

python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS_leptonrho --all
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS_leptonrho --vesselCase
python /afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/surviving_xyz.py --muonDIS_leptonrho --heliumCase
#condor_scripts/condor_hadd.sh "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/muonDIS_leptonrho" "muonDIS_leptonrho" 
wait
echo "completed."

echo "All jobs completed"