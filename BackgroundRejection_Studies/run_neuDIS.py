#!/usr/bin/env python3
"""Wrapper Script to run selection checks on neuDIS samples."""

#-----------------------------------------------------------------------------------------------------------

from selectioncheck import main
import sys, argparse


def calcweight_neuDIS(event,SHiP_running=15,N_gen=100000*98): #Each file has 100k events each change N_gen according to files(1) used for analysis, and 98 successful jobs
    
    w_DIS    =  event.MCTrack[0].GetWeight()
    nPOTinteraction     =(2.e+20)*(SHiP_running/5)
    nPOTinteraction_perspill =5.e+13

    n_Spill  = nPOTinteraction/nPOTinteraction_perspill #number of spill in SHiP_running(default=15) years
    
    nNu_perspill=4.51e+11       #number of neutrinos in a spill.
    
    N_nu=nNu_perspill*n_Spill   #Expected number of neutrinos in 15 years

    w_nu=nNu_perspill/N_gen     #weight of each neutrino considered scaled to a spill such that sum(w_nu)=(nNu_perspill/N_gen)*N_gen= nNu_perspill = number of neutrinos in a spill.
    
    N_A=6.022*10**23
    E_avg=2.57 #GeV
    sigma_DIS=7*(10**-39)*E_avg*N_A  #cross section cm^2 per mole
    
    return w_DIS*sigma_DIS*w_nu*n_Spill  #(rho_L*N_nu*N_A*neu_crosssection*E_avg)/N_gen     #returns the number of the DIS interaction events of that type in SHiP running(default=5) years.   #DIS_multiplicity=1 here


p = argparse.ArgumentParser(description=__doc__)
p.add_argument("-p", "--path", default="/eos/experiment/ship/user/Iaroslava/train_sample_N2024_big/")

known, rest = p.parse_known_args(sys.argv[1:])

# Pass the parsed path plus any *remaining* CLI args to the core.
sys.argv = [sys.argv[0], *rest, "-p", known.path]
main(weight_function=calcweight_neuDIS)