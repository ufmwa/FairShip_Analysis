"""
Check the false veto probability  for a signal candidate + MuonBackground scenario for different SBT vetoes

INUPUT:

"""
import sys
sys.path.insert(0, "/eos/experiment/ship/user/anupamar/NN_data/ext_pkgs")
import typing_extensions

from argparse import ArgumentParser
import numpy as np
import os,ROOT
from rootpyPickler import Unpickler
from tabulate import tabulate
import rootUtils as ut
import shipunit as u
import csv
import time
import pythia8_conf


pdg = ROOT.TDatabasePDG.Instance()
pythia8_conf.addHNLtoROOT()

parser = ArgumentParser()
parser.add_argument("-jobID", "--jobID",dest="jobID",help="job ID ( Process number) unique index to save output", type=int, default=0)
parser.add_argument("-n", "--nEvents",dest="nEvents",help="nEvents per root file", type=int, default=1000)
parser.add_argument("-sg", "--signalpath", dest="signal_path", help="Path to Signal Files", required=False, default= '/eos/experiment/ship/user/anupamar/signal/mumunu',type=str)
parser.add_argument("-bg", "--Backgroundpath", dest="bg_path", help="Path to MuonBack Files", required=False, default='/eos/experiment/ship/simulation/bkg/MuonBack_2024helium/8070735',type=str)
parser.add_argument("--test" , dest="testing_code" , help="Run Test"              , required=False, action="store_true",default=False)
options = parser.parse_args()


from tabulate import tabulate

class FalseVetoProbability:
	def __init__(self,signal_path,bg_path):

		self.bg_path=bg_path
		self.signal_path=signal_path
		
		self.SBTefficiency=1 #0.99
		self.random=ROOT.TRandom()
		
		if options.testing_code:
			seed_value = int(123456)
		else:
			seed_value = int(time.time())
		
		print(f"Setting Seed: {seed_value}")
		self.random.SetSeed(seed_value)

		self.threshold_list=[0,10,20,30,45,60,90]

		
		self.vetoed_cases={}
		self.n_cases={}
		
		#self.geo_file=None
		self.build_embgchain()

		self.f_signal = ROOT.TFile.Open(f"{self.signal_path}/ship.conical.Pythia8-TGeant4_rec.root","read")
		self.signal_tree = self.f_signal.cbmsim
		self.signal_entries =self.signal_tree.GetEntries()

		geofile_name=os.path.join(f"{self.signal_path}/geofile_full.conical.Pythia8-TGeant4.root")
		self.geo_file = ROOT.TFile.Open(geofile_name,"read")
					


	def load_geofile(self):
		"""
		Load the geometry file and set the global geometry manager.
		"""
		try:
			fgeo = ROOT.TFile(self.geo_file)  	
			self.fGeo = fgeo.FAIRGeom  
			ROOT.gGeoManager = self.fGeo  

			upkl    = Unpickler(fgeo)
			self.ShipGeo = upkl.load('ShipGeo')
			print(f"Loaded geometry file: {self.geo_file}")

			import geomGeant4
			
			if hasattr(self.ShipGeo.Bfield,"fieldMap"):
			  self.fieldMaker = geomGeant4.addVMCFields(self.ShipGeo, '', True, withVirtualMC = False)
			else:
			  print("no fieldmap given, geofile too old, not anymore support")
			  exit(-1)

			sGeo   = fgeo.FAIRGeom
			self.fM=None

		except Exception as e:
			raise FileNotFoundError(f"Error loading geo file: {self.geo_file}. Error: {e}")

	def digitizecombinedSBT(self,candidate_t0):
	    
	    index=0 
	    digiSBT={}
	    
	    for detID in self.ElossPerDetId:
	        aHit = ROOT.vetoHit(detID,self.ElossPerDetId[detID])
	        aHit.SetTDC(min( self.tOfFlight[detID] )+ candidate_t0 )    
	        if self.ElossPerDetId[detID]<0.045:    aHit.setInvalid()  
	        digiSBT[index] = aHit
	        index=index+1
	    return digiSBT		

	def generate_event_time_in_spill(self,eventweight,starttime=0,endtime=10**9):
		return np.array([self.random.Uniform(starttime, endtime) for _ in range(int(eventweight))])  
	
	def assign_event_time_signal(self,nevents):

		print("Assigning t0 time for candidate events now...")

		self.signal_eventtimes=self.generate_event_time_in_spill(eventweight=nevents,starttime=self.timeframe[0]+75,endtime=self.timeframe[1]-75)

		print(f"Signal times assigned,{len(self.signal_eventtimes)} entries available")		
	
	def build_embgchain(self):

		self.embg_chain = ROOT.TChain("cbmsim")

		for jobNr,job_folder in enumerate(os.listdir(self.bg_path)):
		
			if not job_folder.startswith('job'): continue # to ignore the README

			if "geofile_full.conical.MuonBack-TGeant4.root" not in os.listdir(f"{self.bg_path}/{job_folder}"): continue
			
			if options.testing_code and jobNr>2: break 
			
			try:

				file_path=f"{self.bg_path}/{job_folder}/ship.conical.MuonBack-TGeant4_rec.root"
				
				with ROOT.TFile.Open(file_path,"read") as f:
					
					test=f.Get("cbmsim")
					
					if test:
						self.embg_chain.Add(file_path)
						print(f"{file_path} added to TChain")
					
					
			except Exception as e:
				print(f"build_embgchain error:{e}")

		print(f"Number of events in the MuBack sample {self.embg_chain.GetEntries()}")
		# Disable all branches then enable only those used in the analysis
		self.embg_chain.SetBranchStatus("*", 0)
		self.embg_chain.SetBranchStatus("vetoPoint*", 1)
		self.embg_chain.SetBranchStatus("Digi_SBTHits*", 1)
		self.embg_chain.SetBranchStatus("MCTrack*", 1)
		self.embg_chain.SetCacheSize(10000000) 
		
	def append_vetoPoints(self,vetoPoints,time_offset):

		for aMCPoint in vetoPoints:

			detID=aMCPoint.GetDetectorID()
			Eloss=aMCPoint.GetEnergyLoss()

			if detID not in self.ElossPerDetId:
			    self.ElossPerDetId[detID]=0
			    self.tOfFlight[detID]=[]

			self.ElossPerDetId[detID] += Eloss

			hittime = aMCPoint.GetTime()

			self.tOfFlight[detID].append(hittime+time_offset)

	def assign_event_time_bg(self):
		
		print("Assigning t0 time for MuBack events now...")
		
		self.bg_eventtimes=[]

		print(f"time frame used:{self.timeframe}")
		
		for embgNr,self.embg_event in enumerate(self.embg_chain):
			
			for track in self.embg_event.MCTrack: 
				if track.GetPdgCode() in [-13,13]:
					self.weight_i=track.GetWeight()
					break
			
			if options.testing_code:
				eventtimes=self.generate_event_time_in_spill(self.weight_i,starttime=self.timeframe[0],endtime=self.timeframe[1]) 
			else:
				eventtimes=self.generate_event_time_in_spill(self.weight_i) 
			

			valid_mask = (eventtimes > self.timeframe[0]) & (eventtimes < self.timeframe[1])
			valid_times = eventtimes[valid_mask]
			
			if valid_times.size == 0:
			    continue # if no event time falls within the timeframe
			
			#print(f"MuBack Event {embgNr} falls within the timeframe")			  
			
			# Fill histogram and record the valid times
			for t in valid_times:
			    self.h['bg_t0'].Fill(t)
			    self.bg_eventtimes.append({"entry": embgNr, "t0": t})
			    self.h["nSBThits"].Fill(len(self.embg_event.Digi_SBTHits)) #weighted

		print(f"BG times assigned,{len(self.bg_eventtimes)} entries available")

	def check_veto(self):

		combine_table_header = ['Signal Event Index','Candidate t0 (ns)','Candidate_nVetoPoints(=0)','MuBack Event Index','MuBack t0 (ns)', 'MuBack nVetoPoints','Total_cumulate_nVetoPoints ','Total_Digi_SBTHits@0MeV']
		veto_table_header = ['Signal Event Index', "[Extrapolate+GNN]Veto@0MeV","[Extrapolate+GNN]Veto@45MeV","[Extrapolate+GNN]Veto@90MeV",  f"nSBTHits>45MeV",f"nSBTHits>0 MeV",f"nSBTHits>10 MeV",f"nSBTHits>20 MeV",f"nSBTHits>30 MeV",f"nSBTHits>60 MeV",f"nSBTHits>90 MeV"]

		# Once before the loop
		bg_event_t0s = np.array([e["t0"] for e in self.bg_eventtimes])
		bg_event_entries = np.array([e["entry"] for e in self.bg_eventtimes])
		
		veto_table=[]
		
		import helperfunctions as analysis_toolkit #torch ROOT 6.32 crash workaround, import torch AFTER initialising ROOT

		ctx       = analysis_toolkit.AnalysisContext(self.signal_tree, self.geo_file)

		selection = analysis_toolkit.selection_check(ctx)
		inspector = analysis_toolkit.event_inspector(ctx)
		veto_ship = analysis_toolkit.veto_tasks(ctx)


		for signalNr,signal_t0 in enumerate(self.signal_eventtimes):
			
			self.signal_tree.GetEntry( self.random.Integer(self.signal_entries) )
			#self.signal_tree.GetEntry(signalNr)
			
			self.signal_event=self.signal_tree
			
			if not len(self.signal_event.Particles): continue #only look if a reconstructed particle exists

			combine_table=[]	
				
			self.h['signal_t0'].Fill(signal_t0) 
			
			self.ElossPerDetId    = {}
			self.tOfFlight        = {}
				
			mask = np.abs(bg_event_t0s - signal_t0) <= 75
			matching_entries = bg_event_entries[mask]
			matching_t0s = bg_event_t0s[mask]

			for bg_t0, entry in zip(matching_t0s, matching_entries):

				self.embg_chain.GetEntry(entry)
				
				self.embg_event = self.embg_chain
								
				self.append_vetoPoints(self.embg_event.vetoPoint,bg_t0-signal_t0)

				combine_table.append([
							signalNr,  # Signal Event Index
							round(signal_t0, 3),  # Signal Time (ns)
							0,
							entry,  # Background Event Number
							round(bg_t0, 3),  # BG Time (ns)
							len(self.embg_event.vetoPoint), #number of vetoPoints in the EMBG event
							sum(len(vetopoints) for vetopoints in self.tOfFlight.values()),
							len(self.ElossPerDetId),#number of digihits
							])
			
			combined_Digi_SBTHits=self.digitizecombinedSBT(self.signal_event.ShipEventHeader.GetEventTime())
			
			print(tabulate(combine_table, headers=combine_table_header, tablefmt="pretty"))

			self.h["ncombinedSBThits"].Fill(len(combined_Digi_SBTHits))
			
			
			for signal in self.signal_event.Particles: 
	
				track_index_first,track_index_last = signal.GetDaughter(0),signal.GetDaughter(1)

				xs, ys, zs=[],[],[]
				bestHits=[]

				for tr in [track_index_first,track_index_last]:

					bestHit,xs_, ys_, zs_=veto_ship.extrapolateTrackToSBT(tr,Digi_SBTHits=combined_Digi_SBTHits.values())

					xs.append(xs_)
					ys.append(ys_)
					zs.append(zs_)
					
					if len(bestHit):
						bestHits.extend(bestHit)
						

				veto, w, hitSegments={},{},{}
				AdvSBT_veto={}

				for threshold in self.threshold_list:

					self.n_cases[self.iterationNr][threshold].add(signalNr) #+=1 #considering only one particle per event!

					ExtSBT_veto=False #default value
				
					veto[threshold], w[threshold], hitSegments[threshold]=veto_ship.SBT_decision(Digi_SBTHits=combined_Digi_SBTHits.values(),threshold=threshold)
					
					for hit in bestHits:
						
						ELoss    = hit.GetEloss()

						if ELoss>=threshold*0.001:
							ExtSBT_veto=True
							break #vetoed by one of the tracks, no need to repeat

					GNN_veto, pBG = veto_ship.Veto_decision_GNNbinary_wdeltaT(threshold=0.6,offset=0)
					
					AdvSBT_veto[threshold]= GNN_veto or ExtSBT_veto #(if either options tells to veto it means veto)
				
					if AdvSBT_veto[threshold]:
						self.vetoed_cases[self.iterationNr][threshold].add(signalNr)
								
				veto_table.append([
									signalNr,  # Signal Event Index
									AdvSBT_veto[0],
									AdvSBT_veto[45],
									AdvSBT_veto[90],
									len(hitSegments[45]),
									len(hitSegments[0]),
									len(hitSegments[10]),
									len(hitSegments[20]),
									len(hitSegments[30]),
									len(hitSegments[60]),
									len(hitSegments[90]),
									])

		print(tabulate(veto_table, headers=veto_table_header, tablefmt="pretty"))

	def calculate_probability(self):
		
		print("\n\nCalculating probability now:")
		
		prob_header=['Threshold','nEvents checked','nEvents vetoed','Prob(%)']
		
		Total_vetoed_cases = {thresh: 0 for thresh in self.threshold_list}
		Total_n_cases = {thresh: 0 for thresh in self.threshold_list}
		
		for iteration in range(1):
			
			prob_data=[]
			
			for threshold in self.threshold_list:

				Total_vetoed_cases[threshold]+=len(self.vetoed_cases[iteration][threshold])
				Total_n_cases[threshold]+=len(self.n_cases[iteration][threshold])

				Prob=len(self.vetoed_cases[iteration][threshold])/float(len(self.n_cases[iteration][threshold]))
		
				prob_data.append([threshold,len(self.n_cases[iteration][threshold]),len(self.vetoed_cases[iteration][threshold]),Prob*100])

			print(f"\n\nIteration {iteration}\n{tabulate(prob_data,headers=prob_header, tablefmt='pretty')}")
			print("------------------------------------------------------------------------------------------------------------------------\n")
		
		Total_prob_data=[]

		for threshold in self.threshold_list:

			Prob=Total_vetoed_cases[threshold]/float(Total_n_cases[threshold])
	
			Total_prob_data.append([threshold,Total_n_cases[threshold],Total_vetoed_cases[threshold],Prob*100])
		print("\n\n------------------------------------------------------------------------------------------------------------------------")
		print("------------------------------------------------------------------------------------------------------------------------")
		print(f"\n\nCumulative results: \n{tabulate(Total_prob_data,headers=prob_header, tablefmt='pretty')}")

		with open(f"falseveto_results_{options.jobID}.csv", mode='w', newline='') as file:
			writer = csv.writer(file)
			writer.writerow(prob_header)  
			writer.writerows(Total_prob_data)   

	def run_analysis(self):
		
		for iteration in range(1):
			
			self.iterationNr=iteration
			
			print("\n\n------------------------------------------------------------------------------------------------------------------------\n\n")

			timeblockreference=self.random.Uniform(0+500, (10**9)-(500)) #reference time in the spill to set the timeframe
		
			self.timeframe=[timeblockreference-500,timeblockreference+500] #timeframe of 1/millionth of a spill

			self.h={}

			ut.bookHist(self.h,  "signal_t0","; signal t0 (ns); nEvents", 1000,self.timeframe[0],self.timeframe[1])
			ut.bookHist(self.h,  "bg_t0","; MuonBack t0 (ns); nEvents", 1000,self.timeframe[0],self.timeframe[1])
			ut.bookHist(self.h,  "nSBThits","; nSBThits@0MeV in the BG sample event (weighted); nEvents", 1000,0,1000)
			ut.bookHist(self.h,  "ncombinedSBThits","; ncombinedSBThits per checked event (no threshold); nEvents", 1000,0,1000)

			for threshold in self.threshold_list:
				ut.bookHist(self.h,  f"Edep_{threshold}",f" E deposition per cell for events checked (@ Threshold= {threshold} MeV) ; E deposition per cell (GeV);", 100,0,1) #"nEvents checked" is after weighting and then sampling

			self.vetoed_cases[iteration] = {thresh: set() for thresh in self.threshold_list}
			self.n_cases[iteration] = {thresh: set() for thresh in self.threshold_list}
		
			print(f"Iteration {iteration}")
		
			self.assign_event_time_bg()

			nevents_to_check = min(options.nEvents,self.signal_entries)
		
			print(f"nevents_to_check: {nevents_to_check}")

			self.assign_event_time_signal(nevents_to_check)
			
			self.check_veto()
			
			print("------------------------------------------------------------------------------------------------------------------------")

		self.calculate_probability()
		ut.writeHists(self.h,"sbtplots.root")


import random

subdirs = [
    d for d in os.listdir(options.signal_path)
    if os.path.isdir(os.path.join(options.signal_path, d))
]

signal_subfolder = random.choice(subdirs)
print(f"Candidate file in {options.signal_path}/{signal_subfolder}")

# Initialize and run the analysis
falsevetoprobability = FalseVetoProbability(f"{options.signal_path}/{signal_subfolder}",options.bg_path)
falsevetoprobability.run_analysis()
