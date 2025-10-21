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
    
    parser.add_argument("--path", help="Path to simulation file", default="/eos/experiment/ship/simulation/bkg/NeutrinoDIS_2024helium/10864335/")
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
    
    for job_nr,job_folder in enumerate(os.listdir(f"{options.path}")):
        
        print("\t",job_nr,job_folder)

        if job_nr>100:
                break
        
        #if ncandidates>10000:
        #    break

        try:
        
            f = ROOT.TFile.Open(f"{options.path}/{job_folder}/ship.conical.Genie-TGeant4_rec.root", "read")
            tree = f.Get("cbmsim")
            geo_file=None

        except:
            continue

        if geo_file==None:
            geo_file = ROOT.TFile.Open(
                f"{options.path}/{job_folder}/geofile_full.conical.Genie-TGeant4.root", "read"
            )
            
            ROOT.geometry_manager = geo_file.Get("FAIRGeom")
            #unpickler = Unpickler(geo_file)
            #ship_geo = unpickler.load("ShipGeo")
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

            for candidate_id_in_event, signal in enumerate(event.Particles):
                

                for cat in {"all",category}:
                    
                    if cat=='heliumCase':
                        rho_l       = helium_rhoL(event)
                    else:
                        rho_l       = event.MCTrack[2].GetWeight()
                    
                    event_weight = rho_l #calcweight_neuDIS(event,w_DIS=rho_l)

                    signal_pos = ROOT.TLorentzVector()
                    signal.ProductionVertex(signal_pos)
                    
                    hist_dict[cat][cat+"_impact_parameter"].Fill(selection.impact_parameter(signal),event_weight)
                    hist_dict[cat][cat+"_impact_parameter_vs_z_interac_point"].Fill(interaction_point.Z(),selection.impact_parameter(signal),event_weight)
                    hist_dict[cat][cat+"_impact_parameter_vs_zvtx"].Fill(signal_pos.Z(),selection.impact_parameter(signal),event_weight)
                    hist_dict[cat][cat+"_dist_to_innerwall"       ].Fill(selection.dist_to_innerwall(signal),event_weight)
                    hist_dict[cat][cat+"_dist_to_innerwall_unweighted"       ].Fill(selection.dist_to_innerwall(signal))

                    if selection.is_in_fiducial(signal):
                        hist_dict[cat][cat+"_wsel_dist_to_innerwall"       ].Fill(selection.dist_to_innerwall(signal),event_weight)

                    pre_ok = selection.preselection_cut(signal, IP_cut=250, show_table=False)
                    
                    if pre_ok:
                        
                        ncandidates+=1

                        hist_dict[cat][cat+"_wsel_impact_parameter_vs_z_interac_point"].Fill(interaction_point.Z(),selection.impact_parameter(signal),event_weight)
                        hist_dict[cat][cat+"_wsel_impact_parameter_vs_zvtx"].Fill(signal_pos.Z(),selection.impact_parameter(signal),event_weight)

    for c in cats:
        #ut.writeHists(hist_dict[c], f"ip_plots_neuDIS_{c}.root")
        ut.writeHists(hist_dict[c], f"ip_plots_neuDIS_{c}_wdist2innerwall.root")


if __name__ == "__main__":
    main()
