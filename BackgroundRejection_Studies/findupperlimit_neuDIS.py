#!/usr/bin/env python
"""Example script for usage of the analysis_toolkit for signal selection."""

from argparse import ArgumentParser
import os
import ROOT
import rootUtils as ut
import math
#from experimental import analysis_toolkit

def calcweight_neuDIS(event,SHiP_running=15,N_gen=1.2e8,w_DIS=None) #6000*50,w_DIS=None):#6k events per job, 19.993k jobs #For Iaroslava productions 2024,N_gen=100000*98): #Each file has 100k events each change N_gen according to files(1) used for analysis, and 98 successful jobs
    
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
    
    return "other"       

def main():
    
    """Sample function to analyse the pre-selection parameters."""
    parser = ArgumentParser(description=__doc__)

    parser.add_argument("--path", help="Path to simulation file", default="/eos/experiment/ship/user/Iaroslava/NeutrinoSample2023/")#default="/eos/experiment/ship/simulation/bkg/NeutrinoDIS_2024helium/10864335")
    parser.add_argument(     "--test"    ,dest="testing_code",help="Run Test on 100 events of the input file"              ,  action="store_true")
    options = parser.parse_args()

    geo_file=None
    
    numberneutrinos,weightedcross,weightedrho_l,numberdis,numberdis_vtx={},{},{},{},{}
    globaleventnr=-1
    successjobnr=0
    for job_nr,job_folder in enumerate(os.listdir(f"{options.path}")):

        print("\t",job_nr,job_folder)

        #if options.testing_code and job_nr>500:
        #        break
    
        try:
        
            f = ROOT.TFile.Open(f"{options.path}/{job_folder}/ship.conical.Genie-TGeant4_rec.root", "read")
            tree = f.Get("cbmsim")
            successjobnr+=1
        except:
            continue

        if successjobnr>5: break

        if geo_file==None:
            geo_file = ROOT.TFile.Open(
                f"{options.path}/{job_folder}/geofile_full.conical.Genie-TGeant4.root", "read"
            )
            
            ROOT.geometry_manager = geo_file.Get("FAIRGeom")
            
            #unpickler = Unpickler(geo_file)
            #ship_geo = unpickler.load("ShipGeo")
            #import helperfunctions as analysis_toolkit
            #ctx       = analysis_toolkit.AnalysisContext(tree, geo_file)
            #selection = analysis_toolkit.selection_check(geo_file)
            #inspector = analysis_toolkit.event_inspector()

        
        for event_nr, event in enumerate(tree):
            
            globaleventnr+=1

            #if options.testing_code and event_nr>50:
            #    break

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

            cat = ip_category(ip_elem)       # "all", "vesselCase", "heliumCase"
            # -----------------------------------------------------------------
            
            #print(cat)
            
            if cat not in numberneutrinos:
                numberneutrinos[cat],weightedcross[cat],weightedrho_l[cat],numberdis[cat]=0,0,0,0
                numberdis_vtx[cat]={}
            
            nu_weight   = 4.51e+11/(1.2e8)

            
            if cat=='heliumCase':
                rho_l       = helium_rhoL(event)
            else:
                rho_l       = event.MCTrack[0].GetWeight()
            
            numberdis[cat]       +=calcweight_neuDIS(event,w_DIS=rho_l)
            numberneutrinos[cat] +=nu_weight
            weightedrho_l[cat]   +=nu_weight*rho_l

            #print(nu_weight,rho_l,calcweight_neuDIS(event))
            
            if len(event.Particles)>0:
                numberdis_vtx[cat][globaleventnr]=calcweight_neuDIS(event,w_DIS=rho_l)


    for cat in ("other", "vesselCase", "heliumCase"):
        
        if  cat not in numberdis_vtx : continue

        if len(numberdis_vtx[cat])==0: continue

        print(f"\n\n--------------------{cat}---------------------------")
        avg_rhol=weightedrho_l[cat]/numberneutrinos[cat]

        print(f"\n\n\taverage rho_l={weightedrho_l[cat]}/{numberneutrinos[cat]}\t={weightedrho_l[cat]/numberneutrinos[cat]}")
            
        N_A=6.022*10**23
        E_avg=2.57 #GeV
        sigma_DIS=7*(10**-39)*E_avg*N_A  #cross section cm^2 per mole
        
        print(f"\tcross={sigma_DIS}")
        print(f"\tnumberneutrinos={numberneutrinos[cat]}")
        print(f"\tnumberdis={numberdis[cat]}")

        print(f"\tnumberdis_vtx={sum(numberdis_vtx[cat].values())}[generated={len(numberdis_vtx[cat])}]")
        
        print(f"\n\tN_A*avgrho_l*avg_cross=\t{avg_rhol*sigma_DIS}\t={avg_rhol*sigma_DIS:.2e}")
        

        UL_CL = 2.3 * avg_rhol*sigma_DIS * sum(numberdis_vtx[cat].values())/len(numberdis_vtx[cat])#UpperLimit of (90%CL)

        print(f"\nConfidence limit :{UL_CL}")

        Inefficiency=UL_CL/sum(numberdis_vtx[cat].values())

        print("Inefficiency=",Inefficiency)

        print("\n\n")

if __name__ == "__main__":
    main()
