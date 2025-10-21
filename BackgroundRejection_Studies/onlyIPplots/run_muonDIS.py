#!/usr/bin/env python3
"""Wrapper Script to run selection checks on muonDIS samples."""

#-----------------------------------------------------------------------------------------------------------

import sys, argparse

def calcweight_muonDIS(event,SHiP_running=15,w_DIS=None):
    """Calculate event weight in 15 years."""    
    
    w_mu=event.MCTrack[0].GetWeight()  #weight of the incoming muon*DIS multiplicity normalised to a full spill   sum(w_mu) = nMuons_perspill = number of muons in a spill. w_mu is not the same as N_muperspill/N_gen, where N_gen = nEvents*DISmultiplicity ( events enhanced in Pythia to increase statistics) .

    cross=event.CrossSection
    
    if w_DIS==None:
        rho_l=event.MCTrack[2].GetWeight()
    else:
        rho_l=w_DIS

    N_a=6.022e+23 

    sigma_DIS=cross*1e-27*N_a #cross section cm^2 per mole
    
    nPOTinteraction     =(2.e+20)*(SHiP_running/5) #in years
    nPOTinteraction_perspill =5.e+13
    
    n_Spill  = nPOTinteraction/nPOTinteraction_perspill  #Number of Spills in SHiP running( default=5) years  
        
    weight_i = rho_l*sigma_DIS*w_mu*n_Spill 

    return weight_i    


p = argparse.ArgumentParser(description=__doc__)
p.add_argument("-p", "--path", default="/eos/experiment/ship/simulation/bkg/MuonDIS_2024helium/8070735")

g = p.add_mutually_exclusive_group(required=True)

g.add_argument("--leptonrho",   dest="channel", action="store_const", const="leptonrho",
               help="Run studies for the partially reconstructed semi-leptonic (l ρ) channel")
g.add_argument("--partialreco", dest="channel", action="store_const", const="partialreco",
               help="Run studies for the partially reconstructed (l l ν) channel")
g.add_argument("--fullreco", dest="channel", action="store_const", const="fullreco",
               help="Run studies for the fully reconstructed (l π) channel")


known, rest = p.parse_known_args(sys.argv[1:])

# ---------- choose the analysis module -----------------------------------
if known.channel == "leptonrho":
    from selectioncheck_leptonrho_muonDIS import main
    ipcut=250
elif known.channel == "partialreco":
    from selectioncheck_muonDIS import main
    ipcut=250
elif known.channel == "fullreco":
    from selectioncheck_muonDIS import main
    ipcut=10
else:
    raise RuntimeError("Bug: unhandled channel flag")

# Pass the parsed path plus any *remaining* CLI args to the core.
sys.argv = [sys.argv[0], *rest, "-p", known.path]
main(IP_CUT=ipcut,weight_function=calcweight_muonDIS)