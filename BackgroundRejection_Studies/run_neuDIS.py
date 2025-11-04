#!/usr/bin/env python3
"""Wrapper Script to run selection checks on neuDIS samples."""

#-----------------------------------------------------------------------------------------------------------
import numpy as np
import ROOT
import shipunit as u
import sys, argparse

from helperfunctions import torch_available


def calcweight_neuDIS(event,SHiP_running=15,N_gen=6000*19969,w_DIS=None):#6k events per job, 19.993k jobs #For Iaroslava productions 2024,N_gen=100000*98): #Each file has 100k events each change N_gen according to files(1) used for analysis, and 98 successful jobs
    
    if w_DIS==None:
        w_DIS    =  event.MCTrack[0].GetWeight()
    
    nPOTinteraction     =(2.e+20)*(SHiP_running/5)
    nPOTinteraction_perspill =5.e+13

    n_Spill  = nPOTinteraction/nPOTinteraction_perspill #number of spill in SHiP_running(default=15) years
    
    nNu_perspill=4.51e+11       #number of neutrinos in a spill.
    
    #N_nu=nNu_perspill*n_Spill   #Expected number of neutrinos in 15 years

    w_nu=nNu_perspill/N_gen     #weight of each neutrino considered scaled to a spill such that sum(w_nu)=(nNu_perspill/N_gen)*N_gen= nNu_perspill = number of neutrinos in a spill.
    
    N_A=6.022*10**23
    E_avg=2.57 #GeV
    sigma_DIS=7*(10**-39)*E_avg*N_A  #cross section cm^2 per mole
    
    return w_DIS*sigma_DIS*w_nu*n_Spill  #(rho_L*N_nu*N_A*neu_crosssection*E_avg)/N_gen     #returns the number of the DIS interaction events of that type in SHiP running(default=5) years.   #DIS_multiplicity=1 here

def TDC_correction(event,candidate):#resolve time bug in neuDIS production. to be removed for new productions post 2024
    #print("TDC_correction called")

    def define_t_vtx():

        t0=event.ShipEventHeader.GetEventTime()

        candidatePos = ROOT.TLorentzVector()
        candidate.ProductionVertex(candidatePos)

        d1, d2 = candidate.GetDaughter(0), candidate.GetDaughter(1)
        d1_mc, d2_mc = event.fitTrack2MC[d1], event.fitTrack2MC[d2]

        time_vtx_from_strawhits=[]

        for hit in event.strawtubesPoint:

            if not (int( str( hit.GetDetectorID() )[:1]) ==1 or int( str( hit.GetDetectorID() )[:1]) ==2) : continue #if hit.GetZ() > ( ShipGeo.TrackStation2.z + 0.5*(ShipGeo.TrackStation3.z - ShipGeo.TrackStation2.z) ): continue #starwhits only from T1 and T2 before the SHiP magnet .

            if not (hit.GetTrackID()==d1_mc or hit.GetTrackID()==d2_mc) : continue

            t_straw    = event.MCTrack[0].GetStartT()/1e4+(hit.GetTime()-event.MCTrack[0].GetStartT())
            
            d_strawhit  = [hit.GetX(),hit.GetY(),hit.GetZ()]

            dist     = np.sqrt( (candidatePos.X()-hit.GetX() )**2+( candidatePos.Y() -hit.GetY())**2+ ( candidatePos.Z()-hit.GetZ() )**2) #distance to the vertex #in cm

            Mom          = event.MCTrack[hit.GetTrackID()].GetP()/u.GeV
            mass         = event.MCTrack[hit.GetTrackID()].GetMass()
            v            = u.c_light*Mom/np.sqrt(Mom**2+(mass)**2)

            t_vertex   = t_straw-(dist/v)

            time_vtx_from_strawhits.append(t_vertex)

        t_vtx=np.average(time_vtx_from_strawhits)+t0

        return t_vtx

    ElossPerDetId   = {}
    tOfFlight       = {}
    
    for aMCPoint in event.vetoPoint:

        detID=aMCPoint.GetDetectorID()
        Eloss=aMCPoint.GetEnergyLoss()

        if detID not in ElossPerDetId:
            
            ElossPerDetId[detID]=0
            tOfFlight[detID]=[]

        ElossPerDetId[detID] += Eloss

        hittime = event.MCTrack[0].GetStartT()/1e4+(aMCPoint.GetTime()-event.MCTrack[0].GetStartT()) 

        tOfFlight[detID].append(hittime)

    digiSBT=[]
    
    t0=event.ShipEventHeader.GetEventTime()
    
    for detID in ElossPerDetId:

        aHit = ROOT.vetoHit(detID,ElossPerDetId[detID])
        aHit.SetTDC(min( tOfFlight[detID] ) + t0 )    
        if ElossPerDetId[detID]<0.045:    aHit.setInvalid()  
        digiSBT.append(aHit)
        #index=index+1

    candidate_time = define_t_vtx()
    
    return digiSBT,candidate_time

p = argparse.ArgumentParser(description=__doc__)
#p.add_argument("-p", "--path", default="/eos/experiment/ship/user/Iaroslava/train_sample_N2024_big/")
p.add_argument("-p", "--path", default="/eos/experiment/ship/simulation/bkg/NeutrinoDIS_2024helium/10864335/")
p.add_argument("--no-gnn", action="store_true",
               help="Disable torch-based SBT GNN even if torch is available.")

g = p.add_mutually_exclusive_group(required=True)

g.add_argument("--leptonrho",   dest="channel", action="store_const", const="leptonrho",
               help="Run studies for the partially reconstructed semi-leptonic (l ρ) channel")
g.add_argument("--partialreco", dest="channel", action="store_const", const="partialreco",
               help="Run studies for the partially reconstructed (l l ν) channel")
g.add_argument("--fullreco", dest="channel", action="store_const", const="fullreco",
               help="Run studies for the fully reconstructed (l π) channel")


known, rest = p.parse_known_args(sys.argv[1:])

USE_GNN = (not known.no_gnn) and torch_available()
if not USE_GNN:
    print("[SBT-GNN] torch not available or --no-gnn set → using basic SBT veto (Edep>45 MeV)")


# ---------- choose the analysis module -----------------------------------
from selectioncheck import main

print("Warning: final numbers may need rescaling, use merge.py for final results")

if known.channel == "leptonrho":
    
    print(f"Partial Reco. (l ρ) channel Analysis starts now ")
    ipcut=(10,250)
    finalstate='semileptonic'

elif known.channel == "partialreco":

    print(f"Partial Reco. (l l ν) channel Analysis starts now ")
    ipcut=250
    finalstate='dileptonic'

elif known.channel == "fullreco":

    print(f"Fully Reco. (l π) channel Analysis starts now ")
    ipcut=10
    finalstate='semileptonic'

else:
    raise RuntimeError("Unknown channel flag")


sys.argv = [sys.argv[0], *rest, "-p", known.path]# Pass the parsed path plus any remaining args
main(IP_CUT=ipcut,weight_function=calcweight_neuDIS,finalstate=finalstate,use_gnn=USE_GNN)
