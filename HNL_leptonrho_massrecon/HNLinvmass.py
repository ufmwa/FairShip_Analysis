"""
Calculate the invariant mass for HNL->rholepton channel using Maral's method.

USAGE: python HNLinvmass.py --signalpath <locationtosignalfile>

"""
import os,ROOT
from argparse import ArgumentParser
import numpy as np
from matplotlib.colors import LogNorm
from tabulate import tabulate
import rootUtils as ut
import matplotlib.pyplot as plt

pdg = ROOT.TDatabasePDG.Instance()
import pythia8_conf
pythia8_conf.addHNLtoROOT()
import shipunit as u

parser = ArgumentParser()

parser.add_argument("-sg", "--signalpath", dest="signal_path", help="Path to Signal Files", required=False,type=str, default= '/afs/cern.ch/user/a/anupamar/Analysis/HNLMass/10000_rhoe')#HNLrhoe2GeV/'#
options = parser.parse_args()

if "10000_rhoe" in options.signal_path: nametag="e+rho_1GeV"
if "rhoe2GeV" in options.signal_path:   nametag="e+rho_2GeV"

if "rhomu1GeV" in options.signal_path:  nametag="mu+rho_1GeV"
if "rhomu2GeV" in options.signal_path:  nametag="mu+rho_2GeV"


pion0_mass = 134.9768*0.001     # Mass of the neutral pion.
rho_mass = 775.45*0.001         # Mass of the rho meson.

def writeplots(h1,h2):
    # Create an output ROOT file to store all the histograms

    output_file = ROOT.TFile("comparison_plots.root", "RECREATE")
    
    for h in [h1,h2]:
        for akey in h:
            if not hasattr(h[akey],'Class'): continue
            cln = h[akey].Class().GetName()
            if not cln.find('TH')<0 or not cln.find('TP')<0:   h[akey].Write()
    
    canvas = ROOT.TCanvas("HNLMass", "", 800, 600)
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.85)
    ROOT.gPad.SetLogy()
    histograms = []

    nametag = "e+rho_1GeV"
    histograms.append(h1[f"{nametag}_calcmassdist1"])
    histograms.append(h1[f"{nametag}_calcmassdist2"])
    histograms.append(h1[f"{nametag}_invmassdist"])

    nametag = "e+pi_1GeV"
    histograms.append(h2[f"{nametag}_calcmassdist1"])
    histograms.append(h2[f"{nametag}_calcmassdist2"])
    histograms.append(h2[f"{nametag}_invmassdist"])

    max_value = max(h.GetMaximum() for h in histograms)
    scale_factor = 1.2  # Scale y-axis by 20% to avoid cutoff

    for h in histograms:
        h.SetStats(0)
        h.SetMaximum(max_value * scale_factor)

    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kCyan, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kYellow+2]

    for i, h in enumerate(histograms):
        h.SetLineColor(colors[i])
        h.SetLineWidth(3)
        h.SetFillColor(colors[i])
        #h.SetFillStyle(3205+i*10)
        h.SetFillStyle(3005)
        legend.AddEntry(h, h.GetTitle(), "l")
        h.SetTitle("")  
        if i == 0:
            h.GetXaxis().SetTitle("Reconstructed HNL Mass (GeV)")
            h.Draw("HIST")  # Draw first histogram normally
        else:
            h.Draw("HIST SAME")  # Overlay the rest

    legend.Draw()
    canvas.Update()
    canvas.Write()

    # Define nametags and corresponding histograms
    nametags = ['e+rho_1GeV', 'e+pi_1GeV']

    # Define a list of colors for multiple histograms
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kOrange+7]

    for nametag in nametags:
        h = h1 if nametag == 'e+rho_1GeV' else h2  # Select corresponding histogram container
            # Define histogram groups for comparison
        histogram_groups = [
            (["calc_pion0_pt","truth_pion0_pt"], f"{nametag}: Pion0 Transverse Momentum", "Momentum (GeV/c)", "Entries"),
            (["pos_calc_pion0_pz", "neg_calc_pion0_pz","truth_pion0_pz"], f"{nametag}: Pion0 Parallel Momentum", "Momentum (GeV/c)", "Entries"),
            (["calcmassdist1", "calcmassdist2"], f"{nametag}: HNL Mass Calculation",  "Mass (GeV/c^2)", "Entries")
        ]   

        for hist_keys, plot_title, x_label, y_label in histogram_groups:
            
            canvas = ROOT.TCanvas(f"Comparison_{nametag}_{plot_title}", plot_title, 800, 600)
            
            legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.85)  # Legend in top-right corner
            first_hist = True  # Flag to handle first histogram drawing
            histograms = []  # Store histograms for later adjustments

            for i, hist_key in enumerate(hist_keys):
                hist_name = f"{nametag}_{hist_key}"
                if hist_name in h:
                    hist = h[hist_name]
                    hist.SetStats(0)
                    hist.SetLineColor(colors[i % len(colors)])  # Cycle through colors
                    hist.SetLineWidth(2)
                    
                    if first_hist:
                        hist.Draw("HIST")  # Draw first histogram
                        hist.GetXaxis().SetTitle(x_label)
                        hist.GetYaxis().SetTitle(y_label)
                        hist.SetTitle(f"Comparison of {plot_title}")
                        first_hist = False
                    else:
                        hist.Draw("HIST SAME")  # Overlay subsequent histograms

                    legend.AddEntry(hist, hist_key.replace("_", " "), "l")
                    histograms.append(hist)  # Store for further manipulation
            
            legend.Draw()
            canvas.Write()  # Save the canvas

    output_file.Close()
    print("All plots have been saved in 'comparison_plots.root'.")

def dump(event,mom_threshold=0):

    headers=['#','particle','pdgcode','mother_id','Momentum [Px,Py,Pz] (GeV/c)','StartVertex[x,y,z] (m)','Process', 'GetWeight()', ]#'px','py','pz','vx','vy','vz'
    
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

def pp(pt, theta):
    """Calculate the parallel momentum component (p_parallel) along the beam axis."""
    sin_theta = ROOT.TMath.Sin(theta)
    
    if sin_theta == 0: #particle passes parallel to the beam axis
        return float('inf')  
    return (pt / sin_theta) * ROOT.TMath.Cos(theta)

def calculate_pion0_parallel(lepton_new,pion_new,pion0_vec):

    pion_pt     = pion_new.Pt()
    pion_phi    = pion_new.Phi()
    pion_mass   = pion_new.M()
    pion_pz     = pion_new.Pz() 
    pion_energy = pion_new.E()

    lepton_pt       = lepton_new.Pt()
    lepton_phi      = lepton_new.Phi()
    
    pion0_pt    =  pion0_vec.Pt()   
    

    energy_squared_diff = pion_energy**2 - pion_pz**2  #E^2-pp^2 = pT^2 + m^2
    cos_phi_diff =  ROOT.TMath.Cos(lepton_phi - pion_phi) 
    
    # Compute the term inside the square root
    sqrt_inner_term = (
                        pion_energy**2 * (
                                          pion_mass**4 
                                          - (4*pion_energy**2 * pion0_mass**2) 
                                          + (2*pion_mass**2 * pion0_mass**2) 
                                          + pion0_mass**4 - (2*pion_mass**2*rho_mass**2) 
                                          - 2 * pion0_mass**2 * rho_mass**2
                                          + rho_mass**4
                                          + 4 *pion0_mass**2 * pion_pz**2
                                          + 4 * pion_mass**2 * pion_pt**2
                                          + 4 * pion0_mass**2 * pion_pt**2
                                          - 4 * rho_mass**2 * pion_pt**2
                                          + 2 * lepton_pt**2 * pion_pt**2
                                          + 4 * pion_pt**4
                                          - 4 * pion_energy**2 * pion0_pt**2
                                          + 4 * pion_pz**2 * pion0_pt**2
                                          + 4 * pion_pt * lepton_pt * (pion_mass**2 + pion0_mass**2 - rho_mass**2 + 2 * pion_pt**2) * cos_phi_diff
                                          + 2 * pion_pt**2 * lepton_pt**2 *  ROOT.TMath.Cos(2 * (lepton_phi - pion_phi))
                                          )
                        )##CHECK THIS AGAIN
    
    sqrt_argument = max(sqrt_inner_term, 0)
    term2_part2 =  ROOT.TMath.Sqrt(sqrt_argument)

    # Compute the terms
    term1 = (
        pion_mass**2 * pion_pz + pion0_mass**2 * pion_pz - rho_mass**2 * pion_pz + 
        2 * pion_pz * pion_pt**2 + 2 * pion_pz * pion_pt * lepton_pt * cos_phi_diff
    )
    term2 = term2_part2 / energy_squared_diff

    # Sum the terms and compute P_pi_parallel
    P_pi_parallel = -((1 / (2 * energy_squared_diff)) * (term1 + term2_part2))
    P_pi_parallel1 = -((1 / (2 * energy_squared_diff)) * (term1 - term2_part2))

    return P_pi_parallel,P_pi_parallel1

def transform_to_HNL_frame(decay_vertex, production_vertex,vec):

    def calculate_momentum_in_new_coordinates(vector, x, y, z):
        return ROOT.TVector3(vector.Dot(x), vector.Dot(y), vector.Dot(z))   

    def get_unit_vector(vector):
        return vector * (1.0 / vector.Mag()) 
    
    # Define the z-axis as the direction from the production to the decay vertex of HNL
    z_axis = get_unit_vector(decay_vertex - production_vertex)

    # Choose a reference axis that is not collinear with z-axis
    orthogonal_to_z = ROOT.TVector3(1, 0, 0)
    if abs(orthogonal_to_z.Dot(z_axis)) > 0.9999:  
        orthogonal_to_z = ROOT.TVector3(0, 1, 0) 

    # Define the x and y axes for the perpendicular plane
    x_axis = get_unit_vector(z_axis.Cross(orthogonal_to_z))
    y_axis = get_unit_vector(z_axis.Cross(x_axis))

    # Calculate momentum for vector in the new coordinate system
    transformed_momentum = calculate_momentum_in_new_coordinates(vec.Vect(), x_axis, y_axis, z_axis)    
    
    recalc_E = ROOT.TMath.Sqrt(transformed_momentum.Mag2() + vec.M()**2)  # E' = sqrt(p^2 + m^2) #this is never used

    # Create a new TLorentzVector with transformed momentum and Energy
    transformed_vec = ROOT.TLorentzVector()
    transformed_vec.SetPxPyPzE(transformed_momentum.X(), transformed_momentum.Y(), transformed_momentum.Z(), recalc_E)

    return transformed_vec


def calculate_HNLMass(pion0_pz,lepton_new, pion_new, pion0_vec):

    lepton_pt       = lepton_new.Pt()
    lepton_mass     = lepton_new.M()
    lepton_energy   = lepton_new.E()
    lepton_pz = lepton_new.Pz()
    
    pion_energy = pion_new.E()
    pion_pz     = pion_new.Pz()
    
    pion_lepton_dot = pion_new.Dot(lepton_new)
    
    pion0_pt =  pion0_vec.Pt()   
    pion0_p = np.sqrt(pion0_pt**2 + pion0_pz**2)
    pion0_energy = np.sqrt(pion0_p**2 + pion0_mass**2)
    
    pion0_lepton_dot = pion0_vec.Dot(lepton_new)
    
    
    term1 = lepton_energy * pion0_energy
    term2 = abs(pion_pz)*abs(lepton_pz) + abs(pion0_pz) * abs(lepton_pz)
    
    term3 = 2*lepton_energy*pion_energy
    
    term4 = 2*(lepton_pt**2)

    radicand = (
                lepton_mass**2 
                + rho_mass**2
                + 2*term1
                + term3
                - (2 * pion0_lepton_dot)
                - (2 * pion_lepton_dot)
                )
    
    sqrt_arg = max(radicand, 0)

    HNL_reconmass = np.sqrt( 
                            lepton_mass**2 
                            + rho_mass**2 
                            + 2 * term1 
                            - 2 * term2 
                            + term3 
                            + term4
                            )

    return HNL_reconmass
def bookhists(nametag):
    h={}

    ut.bookHist(h,  f"{nametag}_HNL_decayvertex_xy"," ; x (cm); y (cm)", 1200,-600,600,1200,-600,600)
    ut.bookHist(h,  f"{nametag}_HNL_decayvertex_z"," ; z (cm); ", 1000,-3000,3000)
    ut.bookHist(h,  f"{nametag}_HNL_prodvertex_z"," ; z (cm); ", 146,-5888,-5742) #length of the target
    

    ut.bookHist(h,f"{nametag}_truth_pion0_pt"  ," ; p_{T} (GeV/c) ; ", 100,0,2)
    ut.bookHist(h,f"{nametag}_calc_pion0_pt"   ," ; p_{T} (GeV/c); ", 100,0,2)

    ut.bookHist(h,f"{nametag}_truth_pion0_pz"      ," ;p_{\\parallel} (GeV/c) ; ", 300,-1500,1500)
    ut.bookHist(h,f"{nametag}_pos_calc_pion0_pz"   ," ;p_{\\parallel} (GeV/c) ; ", 300,-1500,1500)
    ut.bookHist(h,f"{nametag}_neg_calc_pion0_pz"   ," ;P_{\\parallel}(GeV/c) ; ", 300,-1500,1500)

    ut.bookHist(h,f"{nametag}_calcmassdist1"   ,f"{nametag} mass (sol 1);GeV/c^{2} ; ", 75,0,3)
    ut.bookHist(h,f"{nametag}_calcmassdist2"   ,f"{nametag} mass (sol 2) ;GeV/c^{2} ; ", 75,0,3)
    ut.bookHist(h,f"{nametag}_invmassdist"     ,f"{nametag} invmass ;GeV/c^{2} ; ", 75,0,3)
    ut.bookHist(h,f"{nametag}_truemassdist"    ,f"{nametag}_truemass ;GeV/c^{2} ; ", 75,0,3)
    
    return h

def Main_function(path,nametag):
    """Maral's method. Check her mathematica script for details."""

    h= bookhists(nametag)
    
    nSuccessful_candidates = 0
    nReconstructed_candidates = 0
    nNegative_mass =0

    file = ROOT.TFile.Open(f"{path}/ship.conical.Pythia8-TGeant4_rec.root")
    sTree = file.cbmsim
    
    mass_list={ '11':0.511*0.001,
                '211':139.57039*0.001,
                '13':105.66*0.001}

    truemassdist, calcmassdist1,calcmassdist2, invmassdist = [], [], [], []
    truth_pion0_pt,calc_pion0_pt=[],[]
    truth_pion0_pz,pos_calc_pion0_pz,neg_calc_pion0_pz=[],[],[]
    
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
            lepton_new  =transform_to_HNL_frame(decay_vertex, production_vertex, lepton_vec)
            

            pion_pdg = pion.getFittedState().getPDG()
            pion_vec = ROOT.TLorentzVector()
            pion_vec.SetPtEtaPhiM(  pion.getFittedState().getMom().Pt(),
                                    pion.getFittedState().getMom().Eta(),
                                    pion.getFittedState().getMom().Phi(),
                                    mass_list[f'{abs(pion_pdg)}']
                                    )
            pion_new    =transform_to_HNL_frame(decay_vertex, production_vertex, pion_vec)
            
            try:
                truth_pi0track = next((track for track in event.MCTrack if pdg.GetParticle(track.GetPdgCode()).GetName() == 'pi0'), None)

                truth_pion0_vec=ROOT.TLorentzVector()
                truth_pion0_vec.SetPtEtaPhiM(truth_pi0track.GetPt(),
                                        truth_pi0track.GetRapidity(),
                                        ROOT.TMath.ATan2(truth_pi0track.GetPy(), truth_pi0track.GetPx()),
                                        truth_pi0track.GetMass()
                                        )
                truth_pion0_new=transform_to_HNL_frame(decay_vertex, production_vertex, truth_pion0_vec)
                h[f"{nametag}_truth_pion0_pt"].Fill(truth_pion0_new.Pt())
                h[f"{nametag}_truth_pion0_pz"].Fill(truth_pion0_new.Pz())
                truth_pion0_pt.append(truth_pion0_new.Pt()) # has to be in the rotated frame
                truth_pion0_pz.append(truth_pion0_new.Pz())
            
            except:
                pass        
                
            
            # Create the HNL momentum vector
            #HNL_vec = ROOT.TLorentzVector()
            #HNL_vec.SetPtEtaPhiM(   HNLtrack.GetPt(),
            #                        HNLtrack.GetRapidity(),
            #                        ROOT.TMath.ATan2(HNLtrack.GetPy(),HNLtrack.GetPx()),
            #                        HNLtrack.GetMass()
            #                        )
            #HNL_new     =transform_to_HNL_frame(decay_vertex, production_vertex, HNL_vec)

            pion0_vec = -(pion_new + lepton_new)
            pion0_pt_rotated = pion0_vec.Pt()

            h[f"{nametag}_calc_pion0_pt"].Fill(pion0_vec.Pt())            
            
            #print("pi0momentum in the rotated frame PT",pion0_vec.Pt(),"\ttruth",truth_pion0_pp.Pt())
            
            calc_pion0_pt.append(pion0_vec.Pt())
            
            sol1_pion0_Parallel,sol2_pion0_Parallel=calculate_pion0_parallel(lepton_new,pion_new,pion0_vec)

            
            h[f"{nametag}_pos_calc_pion0_pz"].Fill(sol1_pion0_Parallel)
            h[f"{nametag}_neg_calc_pion0_pz"].Fill(sol2_pion0_Parallel)
            
            pos_calc_pion0_pz.append(sol1_pion0_Parallel)
            neg_calc_pion0_pz.append(sol2_pion0_Parallel)

            HNLmass1=calculate_HNLMass(sol1_pion0_Parallel,lepton_new, pion_new,pion0_vec)
            HNLmass2=calculate_HNLMass(sol2_pion0_Parallel,lepton_new, pion_new,pion0_vec)
            
            
            if HNLmass1<0 or HNLmass2 <0: 
                nNegative_mass+=1

            truemass = HNLtrack.GetMass()
            truemassdist.append(truemass)
            h[f"{nametag}_truemassdist"].Fill(truemass)

            calcmassdist1.append(HNLmass1)
            h[f"{nametag}_calcmassdist1"].Fill(HNLmass1)

            calcmassdist2.append(HNLmass2)
            h[f"{nametag}_calcmassdist2"].Fill(HNLmass2)

            inv_mass = signal.GetMass()
            invmassdist.append(inv_mass)
            h[f"{nametag}_invmassdist"].Fill(inv_mass)

            nSuccessful_candidates += 1

    print("Events with succesful recon=",nReconstructed_candidates)
    print("Events Analysis successful=",nSuccessful_candidates)
    print("Negative masses", nNegative_mass)

    return h


h_lrho=Main_function(options.signal_path,'e+rho_1GeV')
h_lpi=Main_function("/afs/cern.ch/user/a/anupamar/Analysis/HNLMass/10000_pie/",'e+pi_1GeV')

writeplots(h_lrho,h_lpi)

