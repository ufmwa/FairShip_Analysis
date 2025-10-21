#!/usr/bin/env python
"""Example script for usage of the analysis_toolkit for signal selection."""

from argparse import ArgumentParser
import os
import ROOT
import rootUtils as ut
#from experimental import analysis_toolkit
import math
from collections import defaultdict
import pandas as pd
#from pathlib import Path

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

def helium_rhoL(event, track_index=0,
                rho_he=1.78e-4, tol=6e-5, eps=1e-4, max_hops=256):
    """
    Return helium-only rho*L (g/cm^2) for a DIS vertex in any DV shape.
    Uses TGeoNavigator to follow the track in both directions from the vertex.
    """

    nav = ROOT.gGeoManager.GetCurrentNavigator()
    if not nav:
        raise RuntimeError("No TGeoNavigator; load geometry first.")

    # --- get vertex & direction ---
    v = ROOT.TVector3()
    event.MCTrack[track_index].GetStartVertex(v)
    xv, yv, zv = v.X(), v.Y(), v.Z()

    px = event.MCTrack[track_index].GetPx()
    py = event.MCTrack[track_index].GetPy()
    pz = event.MCTrack[track_index].GetPz()
    pm = math.sqrt(px*px + py*py + pz*pz) or 1.0
    dx, dy, dz = px/pm, py/pm, pz/pm

    # --- helper: is this point in helium? ---
    def in_helium(x, y, z):
        nav.SetCurrentPoint(x, y, z)
        nav.FindNode()
        node = nav.GetCurrentNode()
        if not node: return False
        try:
            rho = node.GetMedium().GetMaterial().GetDensity()
        except Exception:
            return False
        return abs(rho - rho_he) <= tol

    # quick check that vertex is in helium (nudge if needed)
    if not (in_helium(xv, yv, zv) or
            in_helium(xv + dx*eps, yv + dy*eps, zv + dz*eps) or
            in_helium(xv - dx*eps, yv - dy*eps, zv - dz*eps)):
        return 0.0

    # --- helper: integrate helium length from point along direction ---
    def helium_len_from(x0, y0, z0, dx, dy, dz):
        nav.SetCurrentPoint(x0 + dx*eps, y0 + dy*eps, z0 + dz*eps)
        nav.SetCurrentDirection(dx, dy, dz)
        nav.FindNode()
        total, seen_he = 0.0, False
        for _ in range(max_hops):
            node = nav.GetCurrentNode()
            if not node: break
            try:
                rho = node.GetMedium().GetMaterial().GetDensity()
            except Exception:
                rho = -1
            in_he = abs(rho - rho_he) <= tol
            nav.FindNextBoundaryAndStep()
            step = nav.GetStep()
            if in_he:
                total += step
                seen_he = True
            elif seen_he:
                break
            if nav.IsOutside(): break
            cp = nav.GetCurrentPoint()
            nav.SetCurrentPoint(cp[0] + dx*eps, cp[1] + dy*eps, cp[2] + dz*eps)
        return total

    # --- integrate forward + backward from vertex ---
    L_fwd = helium_len_from(xv, yv, zv,  dx,  dy,  dz)
    L_bwd = helium_len_from(xv, yv, zv, -dx, -dy, -dz)
    L_he  = L_fwd + L_bwd

    return rho_he * L_he  # g/cm^2

def ip_category(ip_elem):
    """Return which sub-sample the event belongs to."""
    
    if ip_elem.startswith(("LiSc", "VetoInnerWall", "VetoOuterWall","VetoVerticalRib","VetoLongitRib")):
        return "vesselCase"
    if ip_elem.startswith("DecayVacuum"):
        return "heliumCase"
    
    return "other"       



def main():
    
    """Sample function to analyse the pre-selection parameters."""
    parser = ArgumentParser(description=__doc__)
    
    parser.add_argument("--path", help="Path to simulation file", default="/eos/experiment/ship/simulation/bkg/MuonDIS_2024helium/8070735")
    parser.add_argument(     "--test"    ,dest="testing_code",help="Run Test on 100 events of the input file"              ,  action="store_true")
    options = parser.parse_args()

    geo_file=None
    
    numbermuons,weightedcross,weightedrho_l,numberdis,numberdis_vtx={},{},{},{},{}
    globaleventnr=-1
    
    #numberDIS_permuon={'all':[],'heliumCase':[],'vesselCase':[],'other':[]}
    #numbermuons_rescaled={'all':[],'heliumCase':[],'vesselCase':[],'other':[]}
    
    NUMBERMUONS={'all':[]}
    
    EVENT_TO_MUONMAPPING={}
    
    NUMBERDISPERMUON={'all':{},'heliumCase':{},'vesselCase':{},'other':{}}
    muon_index=-1
    sum_w={'all':0,'heliumCase':0,'vesselCase':0,'other':0}
    sum_w2={'all':0,'heliumCase':0,'vesselCase':0,'other':0}

    for muon_folder in os.listdir(options.path):
        
        if not os.path.isdir(f"{options.path}/{muon_folder}"):
            continue
        
        for job_nr,job_folder in enumerate(os.listdir(f"{options.path}/{muon_folder}")):
            


            print(muon_folder,"\t",job_nr,job_folder)

            if options.testing_code and job_nr>2:
                    break
            
            try:
            
                f = ROOT.TFile.Open(f"{options.path}/{muon_folder}/{job_folder}/ship.conical.muonDIS-TGeant4_rec.root", "read")
                tree = f.Get("cbmsim")
            
            except:
                continue

            if geo_file==None:
                geo_file = ROOT.TFile.Open(
                    f"{options.path}/{muon_folder}/{job_folder}/geofile_full.conical.muonDIS-TGeant4.root", "read"
                )
                
                ROOT.geometry_manager = geo_file.Get("FAIRGeom")
                #unpickler = Unpickler(geo_file)
                #ship_geo = unpickler.load("ShipGeo")
                #import helperfunctions as analysis_toolkit
                #ctx       = analysis_toolkit.AnalysisContext(tree, geo_file)
                #selection = analysis_toolkit.selection_check(geo_file)
                #inspector = analysis_toolkit.event_inspector()

            mu_Pz,mu_Px,mu_Py=0,0,0

            for event_nr, event in enumerate(tree):
                
                globaleventnr+=1

                if options.testing_code and event_nr>200:
                    break
        

                # -----------------------------------------------------------------
                # interaction-point category
                interaction_point = ROOT.TVector3()
                event.MCTrack[0].GetStartVertex(interaction_point)
                try:
                
                    node  = ROOT.gGeoManager.FindNode(interaction_point.X(),
                                                      interaction_point.Y(),
                                                      interaction_point.Z())
                    ip_elem = node.GetVolume().GetName()
                    
                except Exception:
                    ip_elem = ""                 # falls back to global-only

                category = ip_category(ip_elem)       # "all", "vesselCase", "heliumCase"
                # -----------------------------------------------------------------
                
                if not (mu_Pz==tree.MCTrack[0].GetPz() and mu_Px==tree.MCTrack[0].GetPx() and mu_Py==tree.MCTrack[0].GetPy()):
                    muon_index+=1
                    #this is a new muon case.

                    mu_Pz=tree.MCTrack[0].GetPz()
                    mu_Py=tree.MCTrack[0].GetPy()
                    mu_Px=tree.MCTrack[0].GetPx()

                    NUMBERMUONS['all'].append(muon_index)
                    NUMBERDISPERMUON['all'][muon_index]=0
                    NUMBERDISPERMUON['vesselCase'][muon_index]=0
                    NUMBERDISPERMUON['heliumCase'][muon_index]=0
                    NUMBERDISPERMUON['other'][muon_index]=0
                
                mu_weight   = event.MCTrack[0].GetWeight()
                cross       = event.CrossSection
                
                EVENT_TO_MUONMAPPING[(f"{muon_folder}/{job_folder}",event_nr)] = muon_index # each entry added to this contains the f"{muon_folder}/{job_folder}", event_nr and the muonid in a dictionary. 
                
                for cat in {"all",category}:
                    
                    NUMBERDISPERMUON[cat][muon_index]+=1    

                    if cat not in numbermuons:
                        numbermuons[cat],weightedcross[cat],weightedrho_l[cat],numberdis[cat]=0,0,0,0
                        numberdis_vtx[cat]={}
                    
                    if cat=='heliumCase':
                        rho_l       = helium_rhoL(event)
                        print(f"before:{event.MCTrack[2].GetWeight()},after:{rho_l}")
                    else:
                        rho_l       = event.MCTrack[2].GetWeight()
                    
                    numberdis[cat]       +=calcweight_muonDIS(event,w_DIS=rho_l)
                    numbermuons[cat]     +=mu_weight
                    weightedcross[cat]   +=mu_weight*cross
                    weightedrho_l[cat]   +=mu_weight*rho_l

                    sum_w[cat]           +=calcweight_muonDIS(event,w_DIS=rho_l)
                    sum_w2[cat]          +=(calcweight_muonDIS(event,w_DIS=rho_l))**2
                    
                    if len(event.Particles)>0:
                        numberdis_vtx[cat][globaleventnr]=calcweight_muonDIS(event,w_DIS=rho_l)


    for cat in ("all","other", "vesselCase", "heliumCase"):
        print(f"\n\n--------------------{cat}---------------------------")
        avg_rhol=weightedrho_l[cat]/numbermuons[cat]
        avg_cross=weightedcross[cat]/numbermuons[cat]
        print(f"\n\n\taverage rho_l={weightedrho_l[cat]}/{numbermuons[cat]}\t={weightedrho_l[cat]/numbermuons[cat]}")
        print(f"\taverage cross={weightedcross[cat]}/{numbermuons[cat]}\t={weightedcross[cat]/numbermuons[cat]}")
        print(f"\tnumbermuons={numbermuons[cat]}")
        print(f"\tnumberdis={numberdis[cat]}")

        print(f"\tnumberdis_vtx={sum(numberdis_vtx[cat].values())}[generated={len(numberdis_vtx[cat])}]")
        
        N_a=6.022e+23 
        print(f"\n\tN_A*avgrho_l*avg_cross=\t{N_a*avg_rhol*avg_cross*1e-27}\t={N_a*avg_rhol*avg_cross*1e-27:.2e}")
        
        n_eff=(sum_w[cat]**2)/sum_w2[cat]
        
        print("\n\nn_eff=",n_eff)

        continue

        UL_CL = (2.3/n_eff) #* N_a*avg_rhol*avg_cross*1e-27 * sum(numberdis_vtx[cat].values())/len(numberdis_vtx[cat])#UpperLimit of (90%CL)

        print(f"\n90%Confidence limit :{UL_CL}")

        Inefficiency=UL_CL/sum(numberdis_vtx[cat].values())

        print("Inefficiency=",Inefficiency)

        print("Inefficiencywith LoI case=",166/sum(numberdis_vtx[cat].values()))
        print("\n\n")

    #------------------------------------------------------------------
    print("SUMMARYY")

    rows = []
    for (job_id, eventnr), muon_id in EVENT_TO_MUONMAPPING.items():
        rows.append({
            "job_folder": job_id,
            "eventNr": eventnr,
            "muon_id": muon_id,
            "nDIS_all":    NUMBERDISPERMUON["all"].get(muon_id, 0),
            "nDIS_helium": NUMBERDISPERMUON["heliumCase"].get(muon_id, 0),
            "nDIS_vessel": NUMBERDISPERMUON["vesselCase"].get(muon_id, 0),
            "nDIS_other":  NUMBERDISPERMUON["other"].get(muon_id, 0),
        })


    df = pd.DataFrame(rows, columns=[
        "job_folder","eventNr","muon_id","nDIS_all","nDIS_helium","nDIS_vessel","nDIS_other"
    ]).sort_values(["job_folder","eventNr","muon_id"])

    #df.to_csv("muon_dis_summary.csv", index=False)
    #print("saved: muon_dis_summary.csv")


if __name__ == "__main__":
    main()
