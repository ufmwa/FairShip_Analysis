#!/usr/bin/env python3
"""Script to evaluate signal selection on different samples (Studies separated into : all,vesselCase,heliumCase) """

#-----------------------------------------------------------------------------------------------------------

from argparse import ArgumentParser
import ROOT
import rootUtils as ut
import shipunit as u
from collections import defaultdict
import helperfunctions as analysis_toolkit
from tabulate import tabulate
import numpy as np
import pandas as pd
import pathlib, datetime 
import glob

other_ip=set()

def dump_summary_csv(job_id,event_stats,pass_stats,out_dir= ".",**meta):
    """Dump Analysis Summary onto csv file."""  

    baseline = dict(job=job_id,time=datetime.datetime.utcnow().isoformat(timespec="seconds"),**meta)

    def pack(tag, ev_map):
        return {**baseline,"tag": tag,"nCandidates": len(ev_map),"nEvents15y": round(sum(ev_map.values()), 3)}

    rows = [pack(t, event_stats[t]) for t in ("simulated", "reconstructed")] \
         + [pack(t, m) for t, m in pass_stats.items()]

    path = pathlib.Path(out_dir) / f"selection_summary_{job_id.replace('/', '_')}.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    print("summary saved in", path)

def ip_category(ip_elem):
    """Return which sub-sample the event belongs to."""
    
    if ip_elem.startswith(("LiSc", "VetoInnerWall", "VetoOuterWall","VetoVerticalRib","VetoLongitRib")):
        return "vesselCase"
    if ip_elem.startswith("DecayVacuum"):
        return "heliumCase"
    
    return "all"       

def fixwidth_tabulate(rows, headers, *, width=50, **kw):
    """Pad only the first column to `first_width` spaces."""
    pad = lambda s: f"{s:<{width}}"
    rows2    = [[pad(str(r[0])), *r[1:]] for r in rows]
    headers2 = [pad(str(headers[0])), *headers[1:]]
    return tabulate(rows2,
                    headers=headers2,
                    **kw)

def main(weight_function,IP_CUT = 250,fixTDC=None):
    
    """Main function to analyse the selection efficiency of different cuts."""
    
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("-p", "--path"  ,dest="path"         ,help="Path to simulation file",required=True)
    parser.add_argument("-i","--jobDir"  ,dest="jobDir"      ,help="job name of input file",  type=str,required=True)
    parser.add_argument(     "--test"    ,dest="testing_code",help="Run Test on 100 events of the input file"              ,  action="store_true")

    options = parser.parse_args()

    file_name = glob.glob(f"{options.path}/{options.jobDir}/ship.conical*_rec.root")[0]
    
    f = ROOT.TFile.Open(file_name,"read")
    tree = f.Get("cbmsim")

    geofile_name = glob.glob(f"{options.path}/{options.jobDir}/geofile_full.conical*.root")[0]
    geo_file = ROOT.TFile.Open(geofile_name,"read")#f"{options.path}/{options.jobDir}/geofile_full.conical.muonDIS-TGeant4.root", "read")
    
    ctx       = analysis_toolkit.AnalysisContext(tree, geo_file)

    selection = analysis_toolkit.selection_check(ctx)
    inspector = analysis_toolkit.event_inspector(ctx)
    veto_ship = analysis_toolkit.veto_tasks(ctx)


    pre_tags = [
                "n_particles",
                "fiducial",
                "dist2innerwall",
                "dist2vesselentrance",
                "impact_par",
                "doca",
                "n_dof",
                "reduced_chi2",
                "d_mom",
                "preselection",
                ]
    veto_tags = [
                "BasicSBT@45MeV",
                "BasicSBT@90MeV",
                "BasicSBT@0MeV",
                "AdvSBT@45MeV",
                "AdvSBT@90MeV",
                "TOFSBT@45MeV", 
                "TOFSBT@90MeV", 
                "UBT",
                ]

    other_tags= ["inv_mass"
                ,"PID"]


    combinedveto_tags = ["UBT+BasicSBT@45MeV",
                         "UBT+BasicSBT@90MeV",
                         "UBT+AdvSBT@45MeV",
                         "UBT+AdvSBT@90MeV",
                         "UBT+TOFSBT@45MeV",
                         "UBT+TOFSBT@90MeV",
                            ]


    combined_Basic45       = ["preselection+UBT",
                        "preselection+UBT+BasicSBT@45MeV",
                        "preselection+UBT+BasicSBT@45MeV+PID",
                        "preselection+UBT+BasicSBT@45MeV+PID+inv_mass"]

    combined_Basic90       = ["preselection+UBT",
                        "preselection+UBT+BasicSBT@90MeV",
                        "preselection+UBT+BasicSBT@90MeV+PID",
                        "preselection+UBT+BasicSBT@90MeV+PID+inv_mass"]
    
    combined_Adv45       = ["preselection+UBT",
                        "preselection+UBT+AdvSBT@45MeV",
                        "preselection+UBT+AdvSBT@45MeV+PID",
                        "preselection+UBT+AdvSBT@45MeV+PID+inv_mass"]

    combined_Adv90       = ["preselection+UBT",
                        "preselection+UBT+AdvSBT@90MeV",
                        "preselection+UBT+AdvSBT@90MeV+PID",
                        "preselection+UBT+AdvSBT@90MeV+PID+inv_mass"]


    combined_TOF45       = ["preselection+UBT",
                        "preselection+UBT+TOFSBT@45MeV",
                        "preselection+UBT+TOFSBT@45MeV+PID",
                        "preselection+UBT+TOFSBT@45MeV+PID+inv_mass"]

    combined_TOF90       = ["preselection+UBT",
                        "preselection+UBT+TOFSBT@90MeV",
                        "preselection+UBT+TOFSBT@90MeV+PID",
                        "preselection+UBT+TOFSBT@90MeV+PID+inv_mass"]

    combined_other    = [
                        "preselection+AdvSBT@45MeV",
                        "preselection+AdvSBT@90MeV",
                        "preselection+TOFSBT@45MeV",
                        "preselection+TOFSBT@90MeV",

                        "preselection+BasicSBT@45MeV",
                        "preselection+BasicSBT@90MeV",
                        "preselection+inv_mass",
                        "preselection+UBT+inv_mass",
                        "preselection+UBT+PID+inv_mass",
                        "preselection+UBT+BasicSBT@45MeV+inv_mass",
                        "preselection+UBT+BasicSBT@90MeV+inv_mass",
                        "preselection+UBT+TOFSBT@45MeV+inv_mass",
                        "preselection+UBT+TOFSBT@90MeV+inv_mass",

                        "preselection+UBT+AdvSBT@45MeV+inv_mass",
                        "preselection+UBT+AdvSBT@90MeV+inv_mass",
                        ]
    
    from itertools import chain

    ALL_TAGS = list(chain(
        pre_tags,
        veto_tags,
        other_tags,
        combinedveto_tags,
        combined_Basic45, combined_Basic90,
        combined_Adv45,   combined_Adv90,
        combined_TOF45,   combined_TOF90,
        combined_other
    ))

    cats = ("all", "vesselCase", "heliumCase")

    pass_stats = {
        c: defaultdict(lambda: {}, {t: {} for t in ALL_TAGS})
        for c in cats
    }

    event_stats = {c: {"simulated": {}, "reconstructed": {}} for c in cats}

    hist_dict = {}
    
    for c in cats:
        
        hist_dict[c]={}
        
        pre = f"{c}_"
        
        ut.bookHist(hist_dict[c], pre +"rho_l", "rho_l" , 1000, 0, 1000)
        ut.bookHist(hist_dict[c], pre +"candidate_time","candidate time @ decay vertex; ns",300,0, 300)
        ut.bookHist(hist_dict[c], pre +"impact_parameter", "Impact parameter; cm", 500, 0, 500)
        ut.bookHist(hist_dict[c], pre +"dist_to_innerwall", "Distance to inner wall; cm", 200, 0, 100)
        ut.bookHist(hist_dict[c], pre +"dist_to_vesselentrance","Distance to Decay Vessel Entrance; cm", 500, 0,5000)
        ut.bookHist(hist_dict[c], pre +"inv_mass", "Invariant mass; GeV", 500, 0, 5)
        ut.bookHist(hist_dict[c], pre +"DOCA", "Distance of closest approach; cm", 1000, 0, 10)
        ut.bookHist(hist_dict[c], pre +"len_Particles", "len(tree.Particles); Number of candidates per event", 5,0,5)
        ut.bookHist(hist_dict[c], pre +"d_mom", "momentum of daughters; d1 (GeV); d2 (GeV)",400,0,400,400,0,400)
        ut.bookHist(hist_dict[c], pre +"nDOF", "nDOF; d1; d2", 30, 0, 30, 30, 0, 30)
        ut.bookHist(hist_dict[c], pre +"chi2nDOF", "chi2nDOF; d1; d2", 10, 0, 10, 10, 0, 10)


    for event_nr, event in enumerate(tree):
                
        if options.testing_code and event_nr>99: continue
        
        print(f"Event {event_nr}")

        #inspector.dump_event(track_p_threshold=0.5)  # in GeV


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

        event_weight = weight_function(event)
    
        for c in {"all", cat}:
            pre = f"{c}_"

            hist_dict[c][pre+"rho_l"].Fill(event.MCTrack[2].GetWeight())
    
            event_stats[c]["simulated"][event_nr] = event_weight
    
            if len(event.Particles):
                event_stats[c]["reconstructed"][event_nr] = event_weight

        for candidate_id_in_event, signal in enumerate(event.Particles):
            
            selection_list =  defaultdict(dict) 
            
            if cat=="all":
                other_ip.add(ip_elem)

            for c in {"all", cat}:
                pre = f"{c}_"

                hist_dict[c][pre+"impact_parameter"        ].Fill(selection.impact_parameter(signal),event_weight)   
                hist_dict[c][pre+"dist_to_innerwall"       ].Fill(selection.dist_to_innerwall(signal),event_weight)
                hist_dict[c][pre+"dist_to_vesselentrance"  ].Fill(selection.dist_to_vesselentrance(signal),event_weight)
                hist_dict[c][pre+"DOCA"                    ].Fill(selection.DOCA(signal),event_weight)
                hist_dict[c][pre+"len_Particles"           ].Fill(len(tree.Particles),event_weight)
                hist_dict[c][pre+"d_mom"                   ].Fill(*selection.daughtermomentum(signal),event_weight)
                hist_dict[c][pre+"nDOF"                    ].Fill(*selection.nDOF(signal),event_weight)
                hist_dict[c][pre+"chi2nDOF"                ].Fill(*selection.chi2nDOF(signal),event_weight)
                hist_dict[c][pre+"inv_mass"                ].Fill(selection.invariant_mass(signal),event_weight)
                #hist_dict[c][pre+"candidate_time"          ].Fill(selection.define_candidate_time(signal),event_weight)
                
                
            #--------------------------Preselection---------------------------------

            if len(event.Particles) ==1:
                selection_list["n_particles"] = True

            if selection.is_in_fiducial(signal):
                selection_list["fiducial"] = True
                
            if selection.dist_to_innerwall(signal) > 5 * u.cm:
                selection_list["dist2innerwall"] = True
                        
            if selection.dist_to_vesselentrance(signal) > 100 * u.cm:
                selection_list["dist2vesselentrance"] = True

            if selection.impact_parameter(signal) < IP_CUT * u.cm:
                selection_list["impact_par"] = True

            if selection.DOCA(signal) < 1 * u.cm:
                selection_list["doca"] = True

            if np.all(selection.nDOF(signal) > 25):
                selection_list["n_dof"] = True

            if np.all(selection.chi2nDOF(signal) < 5):
                selection_list["reduced_chi2"] = True

            if np.all(selection.daughtermomentum(signal) > 1 * u.GeV):
                selection_list["d_mom"] = True


            pre_ok = selection.preselection_cut(signal, IP_cut=IP_CUT, show_table=False)

            #if pre_ok:
            #    print(f"Event:{event_nr} Candidate_index: {candidate_id_in_event} <--passes the pre-selection\n\n")

            selection_list['preselection'] = pre_ok

            #-------------------------------------------------------------------------
            #-------------------------Veto Decisions----------------------------------
            
            BasicSBT45_veto ,wBasicSBT45,HitsSBT45  =   veto_ship.SBT_decision(threshold=45)
            selection_list['BasicSBT@45MeV'        ]   = not(BasicSBT45_veto)

            
            BasicSBT90_veto ,wBasicSBT90,HitsSBT90  =   veto_ship.SBT_decision(threshold=90)
            selection_list['BasicSBT@90MeV'        ]   = not(BasicSBT90_veto)
            
            
            BasicSBT0_veto      =   bool(len(event.Digi_SBTHits)) #any sbt activity
            selection_list['BasicSBT@0MeV'         ]   = not(BasicSBT0_veto)


            UBT_veto,ubthits    =   veto_ship.UBT_decision()
            selection_list['UBT'                   ]   = not(UBT_veto)

            
            xs, ys, zs,bestHits=[],[],[],[]
            
            AdvSBT90_veto,AdvSBT45_veto,AdvSBT0_veto=False,False,False

            track_index_first,track_index_last  =  signal.GetDaughter(0),signal.GetDaughter(1)
            
            for tr in [track_index_first,track_index_last]:

                bestHit,xs_, ys_, zs_=veto_ship.extrapolateTrackToSBT(tr)
                
                xs.append(xs_)
                ys.append(ys_)
                zs.append(zs_)
                if len(bestHit):
                    bestHits.extend(bestHit)
            
            for hit in bestHits:
                AdvSBT0_veto=True
                ELoss    = hit.GetEloss()
                if ELoss>=90*0.001:
                    AdvSBT90_veto=True
                if ELoss>=45*0.001:
                    AdvSBT45_veto=True
            
            selection_list['AdvSBT@45MeV'          ]   = not(AdvSBT45_veto)
            selection_list['AdvSBT@90MeV'          ]   = not(AdvSBT90_veto)
            
            if fixTDC: #fix timing bug for neuDIS
            
                Digi_SBTHits,candidate_time=fixTDC(event,signal)

                hits_matched,TOFSBT45_veto=veto_ship.pointing_to_vertex(candidate=signal,threshold=45,Digi_SBTHits=Digi_SBTHits,t_vtx=candidate_time)
                hits_matched,TOFSBT90_veto=veto_ship.pointing_to_vertex(candidate=signal,threshold=90,Digi_SBTHits=Digi_SBTHits,t_vtx=candidate_time)
            
            else:
                hits_matched,TOFSBT45_veto=veto_ship.pointing_to_vertex(candidate=signal,threshold=45)
                hits_matched,TOFSBT90_veto=veto_ship.pointing_to_vertex(candidate=signal,threshold=90)
            
            selection_list['TOFSBT@45MeV'          ]   = not(TOFSBT45_veto)
            selection_list['TOFSBT@90MeV'          ]   = not(TOFSBT90_veto)
            #-------------------------------------------------------------------------
            #-----------------------------Other Cuts----------------------------------
            inv_mass_pass=selection.invariant_mass(signal)  > 0.15*u.GeV

            selection_list['inv_mass']   = inv_mass_pass

            pid= selection.pid_decision(candidate=signal)

            if IP_CUT > 10*u.cm:   pid_pass  = (pid==1 or pid==3) #dileptonic final state
            elif IP_CUT <= 10*u.cm:  pid_pass  = (pid==2 or pid==3) #semileptonic final state
            
            selection_list['PID']   = pid_pass
            
            #--------------------------------Combined Cuts-----------------------------------------            
            
            #combined veto efficiency
            selection_list['UBT+'+ 'BasicSBT@45MeV' ]   =   selection_list['UBT'] and selection_list['BasicSBT@45MeV']
            selection_list['UBT+'+ 'BasicSBT@90MeV' ]   =   selection_list['UBT'] and selection_list['BasicSBT@90MeV']
            selection_list['UBT+'+ 'AdvSBT@45MeV'   ]   =   selection_list['UBT'] and selection_list['AdvSBT@45MeV']
            selection_list['UBT+'+ 'AdvSBT@90MeV'   ]   =   selection_list['UBT'] and selection_list['AdvSBT@90MeV']            
            selection_list['UBT+'+ 'TOFSBT@45MeV'   ]   =   selection_list['UBT'] and selection_list['TOFSBT@45MeV']
            selection_list['UBT+'+ 'TOFSBT@90MeV'   ]   =   selection_list['UBT'] and selection_list['TOFSBT@90MeV']      

            selection_list['preselection+'+ 'UBT']   = selection_list['preselection'] and selection_list['UBT']
            
            #--Basic@45--
            selection_list['preselection+'+ 'UBT+'+ 'BasicSBT@45MeV'                    ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['BasicSBT@45MeV']
            selection_list['preselection+'+ 'UBT+'+ 'BasicSBT@45MeV+'+'PID'             ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['BasicSBT@45MeV']  and selection_list['PID']
            selection_list['preselection+'+ 'UBT+'+ 'BasicSBT@45MeV+'+'PID+'+'inv_mass' ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['BasicSBT@45MeV']  and selection_list['PID']  and selection_list['inv_mass']

            #--Basic@45--
            selection_list['preselection+'+ 'UBT+'+ 'BasicSBT@90MeV'                    ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['BasicSBT@90MeV']
            selection_list['preselection+'+ 'UBT+'+ 'BasicSBT@90MeV+'+'PID'             ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['BasicSBT@90MeV']  and selection_list['PID']
            selection_list['preselection+'+ 'UBT+'+ 'BasicSBT@90MeV+'+'PID+'+'inv_mass' ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['BasicSBT@90MeV']  and selection_list['PID']  and selection_list['inv_mass']

            #--Adv@45MeV--
            selection_list['preselection+'+ 'UBT+'+ 'AdvSBT@45MeV'                      ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['AdvSBT@45MeV']
            selection_list['preselection+'+ 'UBT+'+ 'AdvSBT@45MeV+'   +'PID'            ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['AdvSBT@45MeV']    and selection_list['PID']
            selection_list['preselection+'+ 'UBT+'+ 'AdvSBT@45MeV+'   +'PID+'+'inv_mass']   = selection_list['preselection'] and selection_list['UBT'] and selection_list['AdvSBT@45MeV']    and selection_list['PID']  and selection_list['inv_mass']

            #--Adv@90MeV--
            selection_list['preselection+'+ 'UBT+'+ 'AdvSBT@90MeV'                      ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['AdvSBT@90MeV']
            selection_list['preselection+'+ 'UBT+'+ 'AdvSBT@90MeV+'   +'PID'            ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['AdvSBT@90MeV']    and selection_list['PID']
            selection_list['preselection+'+ 'UBT+'+ 'AdvSBT@90MeV+'   +'PID+'+'inv_mass']   = selection_list['preselection'] and selection_list['UBT'] and selection_list['AdvSBT@90MeV']    and selection_list['PID']  and selection_list['inv_mass']

            #--TOF@45MeV--
            selection_list['preselection+'+ 'UBT+'+ 'TOFSBT@45MeV'                      ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['TOFSBT@45MeV']
            selection_list['preselection+'+ 'UBT+'+ 'TOFSBT@45MeV+'   +'PID'            ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['TOFSBT@45MeV']    and selection_list['PID']
            selection_list['preselection+'+ 'UBT+'+ 'TOFSBT@45MeV+'   +'PID+'+'inv_mass']   = selection_list['preselection'] and selection_list['UBT'] and selection_list['TOFSBT@45MeV']    and selection_list['PID']  and selection_list['inv_mass']

            #--TOF@90MeV--
            selection_list['preselection+'+ 'UBT+'+ 'TOFSBT@90MeV'                      ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['TOFSBT@90MeV']
            selection_list['preselection+'+ 'UBT+'+ 'TOFSBT@90MeV+'   +'PID'            ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['TOFSBT@90MeV']    and selection_list['PID']
            selection_list['preselection+'+ 'UBT+'+ 'TOFSBT@90MeV+'   +'PID+'+'inv_mass']   = selection_list['preselection'] and selection_list['UBT'] and selection_list['TOFSBT@90MeV']    and selection_list['PID']  and selection_list['inv_mass']

            
            # other cuts as backup info
            selection_list['preselection+'+ 'UBT+'+ '[AdvSBT + TOFSBT ]@45MeV '   +'PID+'+'inv_mass']   = selection_list['preselection'] and selection_list['UBT'] and (selection_list['AdvSBT@45MeV'] and selection_list['TOFSBT@45MeV'])   and selection_list['PID']  and selection_list['inv_mass']
            selection_list['preselection+'+ 'UBT+'+ '[AdvSBT + TOFSBT ]@90MeV '   +'PID+'+'inv_mass']   = selection_list['preselection'] and selection_list['UBT'] and (selection_list['AdvSBT@90MeV'] and selection_list['TOFSBT@90MeV'])    and selection_list['PID']  and selection_list['inv_mass']
            
            selection_list['preselection+'                                  +'inv_mass' ]   = selection_list['preselection']                                                                                            and selection_list['inv_mass']
            selection_list['preselection+'+ 'UBT+'                          +'inv_mass' ]   = selection_list['preselection'] and selection_list['UBT']                                                                  and selection_list['inv_mass']
            selection_list['preselection+'+ 'UBT+'+ 'PID+'                  +'inv_mass' ]   = selection_list['preselection'] and selection_list['UBT']                                       and selection_list['PID']  and selection_list['inv_mass']
            selection_list['preselection+'+ 'UBT+'+ 'BasicSBT@45MeV+'       +'inv_mass' ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['BasicSBT@45MeV']                             and selection_list['inv_mass']
            selection_list['preselection+'+ 'UBT+'+ 'BasicSBT@90MeV+'       +'inv_mass' ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['BasicSBT@90MeV']                             and selection_list['inv_mass']
            selection_list['preselection+'+ 'UBT+'+ 'AdvSBT@45MeV+'         +'inv_mass' ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['AdvSBT@45MeV']                               and selection_list['inv_mass']
            selection_list['preselection+'+ 'UBT+'+ 'AdvSBT@90MeV+'         +'inv_mass' ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['AdvSBT@90MeV']                               and selection_list['inv_mass']
            selection_list['preselection+'+ 'UBT+'+ 'TOFSBT@45MeV+'         +'inv_mass' ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['TOFSBT@45MeV']                               and selection_list['inv_mass']
            selection_list['preselection+'+ 'UBT+'+ 'TOFSBT@90MeV+'         +'inv_mass' ]   = selection_list['preselection'] and selection_list['UBT'] and selection_list['TOFSBT@90MeV']                               and selection_list['inv_mass']


            selection_list['preselection+'        + 'BasicSBT@45MeV'                    ]   = selection_list['preselection']                           and selection_list['BasicSBT@45MeV']
            selection_list['preselection+'        + 'BasicSBT@90MeV'                    ]   = selection_list['preselection']                           and selection_list['BasicSBT@90MeV']
            selection_list['preselection+'        + 'AdvSBT@45MeV'                      ]   = selection_list['preselection']                           and selection_list['AdvSBT@45MeV']
            selection_list['preselection+'        + 'AdvSBT@90MeV'                      ]   = selection_list['preselection']                           and selection_list['AdvSBT@90MeV']
            selection_list['preselection+'        + 'TOFSBT@45MeV'                      ]   = selection_list['preselection']                           and selection_list['TOFSBT@45MeV']
            selection_list['preselection+'        + 'TOFSBT@90MeV'                      ]   = selection_list['preselection']                           and selection_list['TOFSBT@90MeV']

            #-------------------------------------------------------------------------

            for selection_name, passed in selection_list.items():

                if passed:
                    for c in {"all", cat}:                       # update both dicts
                        pass_stats[c][selection_name][event_nr] = event_weight

    for cat in cats:

        print(f"----------------------------------------------------------{cat}---------------------------------------------------------")
        
        printed_tags =set()

        #--Table 1------------------------
        rows_event_stats = [
                                [tag,
                                len(ev_dict),                 # nCandidates
                                sum(ev_dict.values())]        # nEvents in 15 y
                            for tag, ev_dict in event_stats[cat].items()
                            ]

        
        recon_table = fixwidth_tabulate(
                        rows_event_stats,     
                        headers=[" ", "  nEvents  ", "nEvents in 15 y"],
                        floatfmt=".3f",
                        tablefmt="rounded_grid",
                    )

        print(recon_table)
        printed_tags.update(event_stats[cat])

        #--Table 2------------------------
        rows_pass_stats_presel  = [
                                        [tag,
                                        len(pass_stats[cat].get(tag, {-99:-99})),                 # nCandidates
                                        sum((pass_stats[cat].get(tag, {-99:-99})).values())]      # nEvents in 15 y
                                    for tag in pre_tags
                                    ]

        presel_table = fixwidth_tabulate(
                        rows_pass_stats_presel,
                        headers=["Pre-Selection Cut", "nCandidates", "nEvents in 15 y"],
                        floatfmt=".3f", tablefmt="rounded_grid"
                    )

        print(presel_table)
        printed_tags.update(pre_tags)

        #--Table 3--------------------------
        rows_pass_stats_vetosystems  = [
                                            [tag,
                                            len(pass_stats[cat].get(tag, {-99:-99})),                 # nCandidates
                                            sum((pass_stats[cat].get(tag, {-99:-99})).values())]        # nEvents in 15 y
                                        for tag in veto_tags
                                        ]

        veto_table = fixwidth_tabulate(
            rows_pass_stats_vetosystems,
            headers=["Veto System",   "nCandidates", "nEvents in 15 y"],
            floatfmt=".3f", tablefmt="rounded_grid"
        )
        
        print(veto_table)
        printed_tags.update(veto_tags)


        #--Tables 4-7--------------------------
        
        table_specs = [
                        ("Combined vetosystem efficiency (UBT x SBT)"   , combinedveto_tags),
                        ("Ordered cuts (BasicSBT@45 MeV threshold)"     , combined_Basic45),
                        ("Ordered cuts (BasicSBT@90 MeV threshold)"     , combined_Basic90),
                        ("Ordered cuts (AdvSBT@45 MeV threshold)"       , combined_Adv45),
                        ("Ordered cuts (AdvSBT@90 MeV threshold)"       , combined_Adv90),
                        ("Ordered cuts (TOFSBT@45 MeV threshold)"       , combined_TOF45),
                        ("Ordered cuts (TOFSBT@90 MeV threshold)"       , combined_TOF90),
                        ]
        

        all_tables = []

        for title, tag_order in table_specs:
            rows = []
            printed_tags.update(tag_order)
            for tag in tag_order:
                ev_map = pass_stats[cat].get(tag, {-99:-99})       # empty dict if tag never filled
                n_ev   = len(ev_map)                          # how many distinct events
                w_sum  = sum(ev_map.values())                 # total weighted events
                rows.append([tag, n_ev, f"{w_sum:.3f}"])

            tbl = fixwidth_tabulate(
                rows,
                headers=[title, "nCandidates", "nEvents in 15 y"],
                floatfmt=".3f",
                tablefmt="rounded_grid",
            )
            all_tables.append(tbl)



        for tbl in all_tables:
            print(tbl, end="\n\n")


        #--Table 8------------------------
        
        fallback_tags = [
            tag for tag, ev_map in pass_stats[cat].items()
            if tag not in printed_tags                      #and ev_map if  at least one event passed
        ]


        if fallback_tags:
            rows = [
                [tag, len(pass_stats[cat][tag]), f"{sum(pass_stats[cat][tag].values()):.3f}"]
                for tag in sorted(fallback_tags)
            ]
            fallback_table=fixwidth_tabulate(rows,
                                  headers=["Other cuts ", "nCandidates", "nEvents in 15 y"],
                                  floatfmt=".3f",
                                  tablefmt="rounded_grid")

            print(fallback_table)
        #---------------------------------------------------------------------
        
        ut.writeHists(hist_dict[cat], f"selectionparameters_{cat}.root")

        #dump_summary_csv(options.jobDir, event_stats, pass_stats,out_dir=".", ship_version="2024_helium")
    
        dump_summary_csv(
            job_id=f"{cat}_{options.jobDir}",     # file name gets the suffix
            event_stats=event_stats[cat],
            pass_stats =pass_stats[cat],
            out_dir    =".",
            sample     =cat,                      # extra metadata column if you like
            ship_version="2024_helium"
        )
    print(other_ip)

if __name__ == "__main__":
    main()