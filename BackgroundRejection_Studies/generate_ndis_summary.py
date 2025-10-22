#!/usr/bin/env python
from argparse import ArgumentParser
import os
import ROOT
import pandas as pd

def ip_category(ip_elem):
    """Return which sub-sample the event belongs to."""
    if ip_elem.startswith(("LiSc", "VetoInnerWall", "VetoOuterWall","VetoVerticalRib","VetoLongitRib")):
        return "vesselCase"
    if ip_elem.startswith("DecayVacuum"):
        return "heliumCase"
    return "other"       

def main():
    
    """Generate metadata for n_dis per muon in different medium (He/SBT) for rescaling."""
    parser = ArgumentParser(description=__doc__)
    
    parser.add_argument("--path", help="Path to simulation file", default="/eos/experiment/ship/simulation/bkg/MuonDIS_2024helium/8070735")
    parser.add_argument("--test"    ,dest="testing_code",help="Run Test on few events of the input file"              ,  action="store_true")
    options = parser.parse_args()

    geo_file=None
    
    event_to_muonmapping={}
    
    ndis_per_muon={'all':{},'heliumCase':{},'vesselCase':{},'other':{}}

    muon_index=-1

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

            mu_Pz,mu_Px,mu_Py=0,0,0

            for event_nr, event in enumerate(tree):
                
                if options.testing_code and event_nr>200:
                    break
        
                # -----------------------------------------------------------------
                
                interaction_point = ROOT.TVector3()
                event.MCTrack[0].GetStartVertex(interaction_point)
                
                try:
                
                    node  = ROOT.gGeoManager.FindNode(interaction_point.X(),
                                                      interaction_point.Y(),
                                                      interaction_point.Z())
                    ip_elem = node.GetVolume().GetName()
                    
                except Exception:
                    ip_elem = ""                 # falls back to global-only

                category = ip_category(ip_elem)       # "other", "vesselCase", "heliumCase"
                
                # -----------------------------------------------------------------
                
                if not (mu_Pz==tree.MCTrack[0].GetPz() and mu_Px==tree.MCTrack[0].GetPx() and mu_Py==tree.MCTrack[0].GetPy()):
                    muon_index+=1
                    #this is a new muon case.

                    mu_Pz=tree.MCTrack[0].GetPz()
                    mu_Py=tree.MCTrack[0].GetPy()
                    mu_Px=tree.MCTrack[0].GetPx()

                    ndis_per_muon['all'][muon_index]=0
                    ndis_per_muon['vesselCase'][muon_index]=0
                    ndis_per_muon['heliumCase'][muon_index]=0
                    ndis_per_muon['other'][muon_index]=0
                
                
                event_to_muonmapping[(f"{muon_folder}/{job_folder}",event_nr)] = muon_index 
                
                for cat in {"all",category}:
                    
                    ndis_per_muon[cat][muon_index]+=1    



    #------------------------------------------------------------------

    rows = []
    for (job_id, eventnr), muon_id in event_to_muonmapping.items():
        rows.append({
            "job_folder": job_id,
            "eventNr": eventnr,
            "muon_id": muon_id,
            "nDIS_all":    ndis_per_muon["all"].get(muon_id, 0),
            "nDIS_helium": ndis_per_muon["heliumCase"].get(muon_id, 0),
            "nDIS_vessel": ndis_per_muon["vesselCase"].get(muon_id, 0),
            "nDIS_other":  ndis_per_muon["other"].get(muon_id, 0),
        })


    df = pd.DataFrame(rows, columns=[
        "job_folder","eventNr","muon_id","nDIS_all","nDIS_helium","nDIS_vessel","nDIS_other"
    ]).sort_values(["job_folder","eventNr","muon_id"])

    df.to_csv("ndis_summary.csv", index=False)
    print("saved: ndis_summary.csv")


if __name__ == "__main__":
    main()
