"""
Calculate the invariant mass for HNL->rho lepton channel using a rho constraint using files generated using Fairship.

USAGE: python HNLinvmass_Fairship.py

"""
import HNLinvmass_EventCalc as functions
import ROOT
from tabulate import tabulate
pdg = ROOT.TDatabasePDG.Instance()
import pythia8_conf
pythia8_conf.addHNLtoROOT()
import shipunit as u
import os


def dump(event,mom_threshold=0):

    headers=['#','particle','pdgcode','mother_id','Momentum [Px,Py,Pz] (GeV/c)','StartVertex[x,y,z] (m)','Process', 'GetWeight()', ]
    
    event_table=[]
    for trackNr,track in enumerate(event.MCTrack): 
        
        if track.GetP()/u.GeV < mom_threshold :  continue
        
        try: particlename=pdg.GetParticle(track.GetPdgCode()).GetName()
        except: particlename='----'

        event_table.append([trackNr,
                        particlename,
                        track.GetPdgCode(),
                        track.GetMotherId(),
                        f"[{track.GetPx()/u.GeV:7.3f},{track.GetPy()/u.GeV:7.3f},{track.GetPz()/u.GeV:7.3f}]",
                        f"[{track.GetStartX()/u.m:7.3f},{track.GetStartY()/u.m:7.3f},{track.GetStartZ()/u.m:7.3f}]",
                        track.GetProcName().Data(),
                        track.GetWeight()
                        ])
    
    print(tabulate(event_table,headers=headers,floatfmt=".3f",tablefmt='simple_outline'))

def process_hnl_fairshipevents_with_rho_constraint(path,nametag):

    h= functions.book_histograms(nametag)
    
    nSuccessful_candidates = 0
    nReconstructed_candidates = 0
    nNegative_radicands =0

    mass_list={ '11':0.511*0.001,
                '211':139.57039*0.001,
                '13':105.66*0.001}

    truemassdist, calcmassdist1,calcmassdist2, invmassdist = [], [], [], []
    truth_pion0_pt,calc_pion0_pt=[],[]
    truth_pion0_pp,sol1_calc_pion0_pp,sol2_calc_pion0_pp=[],[],[]
    
    for inputFolder in os.listdir(path):
        try:
            file = ROOT.TFile.Open(f"{path}/{inputFolder}/ship.conical.Pythia8-TGeant4_rec.root")
            sTree = file.cbmsim
        except:
            continue

        for eventNr,event in enumerate(sTree):

            if not len(event.Particles): continue #only look at events with candidates

            nReconstructed_candidates += 1

            #dump(event,mom_threshold=0)#GeV

            for candidate_id_in_event,signal in enumerate(event.Particles):


                HNLtrack = next((track for track in event.MCTrack if pdg.GetParticle(track.GetPdgCode()).GetName() == 'N2'), None)

                if HNLtrack is None:
                    continue  # Skip if no HNL track found

                # Retrieve daughter tracks and masses

                d1_rec = event.FitTracks[signal.GetDaughter(0)]
                d1_mc  = event.MCTrack[event.fitTrack2MC[signal.GetDaughter(0)]]

                d2_rec = event.FitTracks[signal.GetDaughter(1)]
                d2_mc  = event.MCTrack[event.fitTrack2MC[signal.GetDaughter(1)]]

                lepton, pion = None, None

                #print(f"d1: MCTrack #{event.fitTrack2MC[signal.GetDaughter(0)]} pdg:{d1_rec.getFittedState().getPDG()}")
                #print(f"d2: MCTrack #{event.fitTrack2MC[signal.GetDaughter(1)]} pdg:{d2_rec.getFittedState().getPDG()}")

                for d_ in [d1_rec, d2_rec]:
                    d_pdg = d_.getFittedState().getPDG()
                    if abs(d_pdg) in [13,11]: lepton = d_
                    if abs(d_pdg) ==211: pion = d_

                if not (lepton and pion):
                    continue

                production_vertex = ROOT.TVector3(HNLtrack.GetStartX(), HNLtrack.GetStartY(), HNLtrack.GetStartZ())
                h[f"{nametag}_HNL_prodvertex_z"].Fill(HNLtrack.GetStartZ())

                candidatePos = ROOT.TLorentzVector()
                signal.ProductionVertex(candidatePos)
                decay_vertex = ROOT.TVector3(candidatePos.X(), candidatePos.Y(), candidatePos.Z())

                h[f"{nametag}_HNL_decayvertex_z"].Fill(candidatePos.Z())
                h[f"{nametag}_HNL_decayvertex_xy"].Fill(candidatePos.X(),candidatePos.Y())


                # Prepare momentum 4-vectors
                lepton_pdg=lepton.getFittedState().getPDG()
                lepton_vec = ROOT.TLorentzVector()
                lepton_vec.SetPtEtaPhiM(lepton.getFittedState().getMom().Pt(),
                                        lepton.getFittedState().getMom().Eta(),
                                        lepton.getFittedState().getMom().Phi(),
                                        mass_list[f'{abs(lepton_pdg)}']
                                        )
                lepton_rotated  =functions.rotate_momentum_to_hnl_frame(decay_vertex, production_vertex, lepton_vec)


                pion_pdg = pion.getFittedState().getPDG()
                pion_vec = ROOT.TLorentzVector()
                pion_vec.SetPtEtaPhiM(  pion.getFittedState().getMom().Pt(),
                                        pion.getFittedState().getMom().Eta(),
                                        pion.getFittedState().getMom().Phi(),
                                        mass_list[f'{abs(pion_pdg)}']
                                        )
                pion_rotated    =functions.rotate_momentum_to_hnl_frame(decay_vertex, production_vertex, pion_vec)
                
                try:
                    truth_pi0track = next((track for track in event.MCTrack if pdg.GetParticle(track.GetPdgCode()).GetName() == 'pi0'), None)

                    truth_pion0_vec=ROOT.TLorentzVector()
                    truth_pion0_vec.SetPtEtaPhiM(truth_pi0track.GetPt(),
                                            truth_pi0track.GetRapidity(),
                                            ROOT.TMath.ATan2(truth_pi0track.GetPy(), truth_pi0track.GetPx()),
                                            truth_pi0track.GetMass()
                                            )
                    truth_pion0_rotated=functions.rotate_momentum_to_hnl_frame(decay_vertex, production_vertex, truth_pion0_vec)
                    h[f"{nametag}_truth_pion0_pt"].Fill(truth_pion0_rotated.Pt())
                    h[f"{nametag}_truth_pion0_pp"].Fill(truth_pion0_rotated.Pz())

                    rho_track=event.MCTrack[truth_pi0track.GetMotherId()]
                    h[f"{nametag}_truth_rho_M"].Fill(rho_track.GetMass())

                except:
                    pass

                calc_pion0_rotated = -(pion_rotated + lepton_rotated)

                pion0_pt_rotated = calc_pion0_rotated.Pt()

                h[f"{nametag}_calc_pion0_pt"].Fill(calc_pion0_rotated.Pt())

                calc_pion0_pt.append(calc_pion0_rotated.Pt())

                sol1_pion0_pp,sol2_pion0_pp=functions.calculate_pion0_parallel(lepton_rotated,pion_rotated,calc_pion0_rotated)

                if ROOT.TMath.IsNaN(sol1_pion0_pp) or ROOT.TMath.IsNaN(sol2_pion0_pp):
                    rho_mass_min=functions.calculate_min_rho_mass(pion_rotated,lepton_rotated, calc_pion0_rotated)
                    rho_median=functions.rho_median_mass_from_min(rho_mass_min)
                    sol1_pion0_pp,sol2_pion0_pp=functions.calculate_pion0_parallel(lepton_rotated,pion_rotated,calc_pion0_rotated,rho_mass=rho_median)

                h[f"{nametag}_sol1_calc_pion0_pp"].Fill(sol1_pion0_pp)
                h[f"{nametag}_sol2_calc_pion0_pp"].Fill(sol2_pion0_pp)

                HNLmass1=functions.reconstruct_hnl_mass(sol1_pion0_pp,lepton_rotated, pion_rotated,calc_pion0_rotated)
                HNLmass2=functions.reconstruct_hnl_mass(sol2_pion0_pp,lepton_rotated, pion_rotated,calc_pion0_rotated)

                if ROOT.TMath.IsNaN(HNLmass1) or ROOT.TMath.IsNaN(HNLmass2):
                    nNegative_radicands+=1
                    continue

                truemass = HNLtrack.GetMass()
                h[f"{nametag}_truemassdist"].Fill(truemass)

                h[f"{nametag}_calcmassdist1"].Fill(HNLmass1)

                h[f"{nametag}_calcmassdist2"].Fill(HNLmass2)

                inv_mass = signal.GetMass()
                h[f"{nametag}_invmassdist"].Fill(inv_mass)

                nSuccessful_candidates += 1

    print("nEvents available:",nReconstructed_candidates,"\tEvents Successful=",nSuccessful_candidates)
    
    if nNegative_radicands:
        print("nNegative_radicands!!", nNegative_radicands)

    return h


h_lrho =process_hnl_fairshipevents_with_rho_constraint('/eos/experiment/ship/user/anupamar/signal/erho','e+rho')#1GeV
h_lpi  =process_hnl_fairshipevents_with_rho_constraint('/eos/experiment/ship/user/anupamar/signal/epi','e+pi')#1GeV


functions.plot_distributions(h_lrho,h_lpi,filename="plots_HNLMass_Fairship.root")

