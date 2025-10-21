#!/usr/bin/env python
"""Example script for usage of the analysis_toolkit for signal selection."""

from argparse import ArgumentParser
import os
import ROOT
import rootUtils as ut
from experimental import analysis_toolkit
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


dis = (pd.read_csv("/afs/cern.ch/user/a/anupamar/FairShip_Analysis/BackgroundRejection_Studies/muon_dis_summary.csv")#contains the number of DIS in helium /vessel for each muon tagged to the eventNr and job_id; easy lookup!
         .rename(columns={"muon_folder/job_folder": "job_folder"})   # harmless if already named job_folder
         .astype({"eventNr": int})
         .set_index(["job_folder", "eventNr"]))

def ndis_rescale(job_folder,event_nr, ip_cat):
    
    row = dis.loc[(job_folder, int(event_nr))]
    num = int(row["nDIS_all"])
    
    if ip_cat=='vesselCase':
        den = int(row["nDIS_vessel"])
    if ip_cat=='heliumCase':
        den = int(row["nDIS_helium"])
    if ip_cat=='all':
        den = int(row["nDIS_all"])

    return num / den


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
    
    return "all"



def main():
    
    """Sample function to analyse the pre-selection parameters."""
    parser = ArgumentParser(description=__doc__)
    
    parser.add_argument("--path", help="Path to simulation file", default="/eos/experiment/ship/simulation/bkg/MuonDIS_2024helium/8070735")
    parser.add_argument(     "--test"    ,dest="testing_code",help="Run Test on 100 events of the input file"              ,  action="store_true")
    options = parser.parse_args()

    
    cats = ("all", "vesselCase", "heliumCase")

    hist_dict = {}

    for c1 in cats:
        
        hist_dict[c1]={}
        pre=f"{c1}_"

        ut.bookHist(hist_dict[c1], pre+"impact_parameter", "Impact parameter; cm", 300, 0, 300)
        ut.bookHist(hist_dict[c1], pre+"impact_parameter_vs_zvtx", "Impact parameter vs z-vtx; cm", 6000,-3000,3000,300, 0, 300)
        ut.bookHist(hist_dict[c1], pre+"impact_parameter_vs_z_interac_point", "Impact parameter vs z-interaction point; cm", 6000,-3000,3000,300, 0, 300)

        ut.bookHist(hist_dict[c1], pre+"wsel_impact_parameter_vs_zvtx", "Impact parameter vs z-vtx; cm", 6000,-3000,3000,300, 0, 300)
        ut.bookHist(hist_dict[c1], pre+"wsel_impact_parameter_vs_z_interac_point", "Impact parameter vs z-interaction point; cm", 6000,-3000,3000,300, 0, 300)
        ut.bookHist(hist_dict[c1], pre+ "dist_to_innerwall", "Distance to inner wall; cm", 200, 0, 100)
        ut.bookHist(hist_dict[c1], pre+ "dist_to_innerwall_unweighted", "Distance to inner wall; cm", 200, 0, 100)
        ut.bookHist(hist_dict[c1], pre+ "wsel_dist_to_innerwall", "Distance to inner wall; cm", 200, 0, 100)

    globaleventnr=-1

    ncandidates=0
    
    for muon_folder in os.listdir(options.path):
        
        if not os.path.isdir(f"{options.path}/{muon_folder}"):
            continue
        
        for job_nr,job_folder in enumerate(os.listdir(f"{options.path}/{muon_folder}")):
            
            print(muon_folder,"\t",job_nr,job_folder)

            #if job_nr>20:
            #        break
            if ncandidates>5000:
                break            
            
            try:
            
                f = ROOT.TFile.Open(f"{options.path}/{muon_folder}/{job_folder}/ship.conical.muonDIS-TGeant4_rec.root", "read")
                tree = f.Get("cbmsim")
                geo_file=None

            except:
                continue

            if geo_file==None:
                geo_file = ROOT.TFile.Open(
                    f"{options.path}/{muon_folder}/{job_folder}/geofile_full.conical.muonDIS-TGeant4.root", "read"
                )
                
                ROOT.geometry_manager = geo_file.Get("FAIRGeom")
                selection = analysis_toolkit.selection_check(geo_file)

            for event_nr, event in enumerate(tree):
                
                globaleventnr+=1

                if options.testing_code and event_nr>15:
                    break
        

                selection.access_event(event)
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

                cross       = event.CrossSection

                scalefactor={}

                for c2 in {"all",category}:

                    if c2 not in scalefactor:            

                        scalefactor[c2] = ndis_rescale(f"{muon_folder}/{job_folder}",event_nr,c2)

                for candidate_id_in_event, signal in enumerate(event.Particles):
                    
                    ncandidates+=1

                    for cat in {"all",category}:
                        
                        if cat=='heliumCase':
                            rho_l       = helium_rhoL(event)
                        else:
                            rho_l       = event.MCTrack[2].GetWeight()
                        
                        event_weight = calcweight_muonDIS(event,w_DIS=rho_l)*scalefactor[cat]

                        signal_pos = ROOT.TLorentzVector()
                        signal.ProductionVertex(signal_pos)
                        
                        hist_dict[cat][cat+"_impact_parameter"].Fill(selection.impact_parameter(signal),event_weight)
                        hist_dict[cat][cat+"_impact_parameter_vs_z_interac_point"].Fill(interaction_point.Z(),selection.impact_parameter(signal),event_weight)
                        hist_dict[cat][cat+"_impact_parameter_vs_zvtx"].Fill(signal_pos.Z(),selection.impact_parameter(signal),event_weight)
                        hist_dict[cat][cat+"_dist_to_innerwall"       ].Fill(selection.dist_to_innerwall(signal),event_weight)
                        hist_dict[cat][cat+"_dist_to_innerwall_unweighted"       ].Fill(selection.dist_to_innerwall(signal))

                        pre_ok = selection.preselection_cut(signal, IP_cut=250, show_table=False)
                        
                        if selection.is_in_fiducial(signal):
                            hist_dict[cat][cat+"_wsel_dist_to_innerwall"       ].Fill(selection.dist_to_innerwall(signal),event_weight)

                        if pre_ok:

                            hist_dict[cat][cat+"_wsel_impact_parameter_vs_z_interac_point"].Fill(interaction_point.Z(),selection.impact_parameter(signal),event_weight)
                            hist_dict[cat][cat+"_wsel_impact_parameter_vs_zvtx"].Fill(signal_pos.Z(),selection.impact_parameter(signal),event_weight)
                            

    for c in cats:
        ut.writeHists(hist_dict[c], f"ip_plots_muDIS_{c}_wdist2innerwall.root")


if __name__ == "__main__":
    main()
