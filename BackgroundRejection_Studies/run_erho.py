#!/usr/bin/env python3
"""Wrapper Script to run selection checks on Fully recon. samples (EventCalc)."""

#-----------------------------------------------------------------------------------------------------------
from selectioncheck import main
import sys, argparse
import os,ROOT
import pandas as pd
import shipunit as u
import numpy as np

from helperfunctions import torch_available
p = argparse.ArgumentParser(description=__doc__)
p.add_argument("-p", "--path", default="/eos/experiment/ship/user/anupamar/Signal_EventCalc/2pie/12968088/HNL_1.000e+00_7.133e+04_4.808e-02_7.692e-01_1.827e-01",help="Path to simulation file") #mumunu 1 GeV
p.add_argument("--no-gnn", action="store_true",
               help="Disable torch-based SBT GNN even if torch is available.")

known, rest = p.parse_known_args(sys.argv[1:])

USE_GNN = (not known.no_gnn) and torch_available()
if not USE_GNN:
    print("[SBT-GNN] torch not available or --no-gnn set → using basic SBT veto (Edep>45 MeV)")

def define_time_till_vtx(event):
	
	Mom = event.MCTrack[0].GetP()/u.GeV
	mass = event.MCTrack[0].GetMass()

	v= u.c_light*Mom/np.sqrt(Mom**2+(mass)**2)
	
	Target_Point = ROOT.TVector3(0., 0., -5814.25)  # production point
	Decay_Vertex = ROOT.TVector3(event.MCTrack[0].GetStartX(),  event.MCTrack[0].GetStartY(),  event.MCTrack[0].GetStartZ())  
	
	r_vec = Decay_Vertex - Target_Point          
	dist_from_target     = r_vec.Mag()              # cm   
	
	t_ns  = dist_from_target / v
	# return the time taken for the particle to reach the X,YZ from the start point
	return t_ns

def define_weight_EventCalc(event,path=known.path,w_DIS=None): #w_DIS is kept to maintain the similiarity in code
	
	w_event=event.MCTrack[0].GetWeight()

	eventcalcdata_path='/eos/experiment/ship/user/anupamar/EventCalc_data/FairShip_benchmarkcoupling/'

	channel=path.split('/')[-3] 
	foldername=os.path.basename(path)
	
	
	masstag = float(foldername.split('_')[1])

	#print(f"channel:\t{channel},\tmass:\t{masstag}")

	datfile=f'{eventcalcdata_path}/{channel}_sample/HNL/total/HNL_4.808e-02_7.692e-01_1.827e-01_total.txt'
	inputfile=f'{eventcalcdata_path}/rootfiles/{channel}/{foldername}_data.root'
	
	#print(f"{path}\n\t should match \n {inputfile}\n\t does it?")

	file = ROOT.TFile.Open(inputfile)    
	tree=file.Get("LLP_tree")

	#print("nEntries",tree.GetEntries())
	n_events=tree.GetEntries()

	df = pd.read_csv(datfile, delim_whitespace=True)

	row = df.loc[df['mass'] == masstag]

	if row.empty:
	    raise ValueError(f"No entry for HNL mass = {masstag} in {datfile} EVENT WEIGHTS ARE AN ISSUE")

	row = row.iloc[0]

	# Extract factors
	N_LLP_tot      = row['N_LLP_tot']
	eps_pol        = row['epsilon_polar']
	eps_azi        = row['epsilon_azimuthal']
	BR_vis         = row['Br_visible']

	base = N_LLP_tot * eps_pol * eps_azi * BR_vis

	# 5) apply flavour‐specific scaling
	if 'e' in channel:
		expected   = base * (w_event / n_events)
	if 'mu' in channel:
		expected   = base * (w_event / n_events)

	return expected


sys.argv = [sys.argv[0], *rest, "-p", known.path] # Pass the parsed path plus any remaining CLI args

print(f"Partial Reco. (l ρ) channel Analysis starts now ")
ipcut=(10,250)
finalstate='semileptonic'

main(IP_CUT=ipcut,weight_function=define_weight_EventCalc,fix_candidatetime=define_time_till_vtx,finalstate=finalstate,use_gnn=USE_GNN)
