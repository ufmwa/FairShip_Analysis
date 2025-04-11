"""
Check the false veto probability of SBT for a signal candidate 

v2, for a given signal candidate what is the likelihood that it gets vetoed.

INUPUT:
USAGE: python falsevetoprobability.py -bg <path to BG file> -sg <path to signal file> 
"""
from argparse import ArgumentParser
import numpy as np
import os,ROOT
from rootpyPickler import Unpickler
from tabulate import tabulate
import rootUtils as ut
import shipunit as u
import csv
import time

parser = ArgumentParser()
parser.add_argument("-jobID", "--jobID",dest="jobID",help="job ID ( Process number) unique index to save output", type=int, default=0)
parser.add_argument("-n", "--nEvents",dest="nEvents",help="nEvents per root file", type=int, default=100)
parser.add_argument("-bg", "--Backgroundpath", dest="bg_path", help="Path to MuonBack Files", required=False, default='/eos/experiment/ship/simulation/bkg/MuonBack_2024helium/8070735',type=str)
parser.add_argument("--test"          , dest="testing_code" , help="Run Test"              , required=False, action="store_true",default=False)
options = parser.parse_args()


class FalseVetoProbability:
	def __init__(self,bg_path,nIterations=1):

		self.bg_path=bg_path
		
		self.SBTefficiency=1 #0.99
		self.random=ROOT.TRandom()
		
		if options.testing_code:
			seed_value = int(123456)
		else:
			seed_value = int(time.time())
		
		print(f"Setting Seed: {seed_value}")
		self.random.SetSeed(seed_value)
		
		print(f"nIterations: {nIterations}")
		self.nIterations=nIterations

		self.threshold_list=[0,10,20,30,45,60,90]

		
		self.vetoed_cases={}
		self.n_cases={}
		
		self.geo_file=None
		self.build_embgchain()

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

	def SBT_decision(self,Digi_SBTHits,mcParticle=None,detector='liquid',threshold=45,candidate=None):
		
		# if mcParticle >0, only count hits with this particle
		# if mcParticle <0, do not count hits with this particle
		################################

		#hitSegments = 0

		hitSegments = []
		index = -1
		fdetector = detector=='liquid'

		for aDigi in Digi_SBTHits.values():
		 
		 index+=1 
		 detID    = aDigi.GetDetectorID()

		 if fdetector and detID > 999999:continue
		 if not fdetector and not detID > 999999:continue 
		 if mcParticle:
		    found = False
		    for mcP in self.sTree.digiSBT2MC[index]: 
		     if mcParticle>0 and mcParticle != mcP : found=True
		     if mcParticle<0 and abs(mcParticle) == mcP : found=True
		    if found: continue
		 
		 position = aDigi.GetXYZ()
		 ELoss    = aDigi.GetEloss()
		 
		 if ELoss>=threshold*0.001: 
		 	self.h[f"Edep_{threshold}"].Fill(ELoss) 
		 	hitSegments.append(index)#hitSegments += 1 
		 
		 #if aDigi.isValid(): hitSegments += 1 #threshold of 45 MeV per segment

		w = (1-self.SBTefficiency)**len(hitSegments)  
		veto = self.random.Rndm() > w

		#print 'SBT :',hitSegments
		return veto, w, hitSegments#hitSegments contain the index of the Digihit that causes the veto

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
					
					if self.geo_file is None:
						self.geo_file=os.path.join(f"{self.bg_path}/{job_folder}/geofile_full.conical.MuonBack-TGeant4.root")
						self.load_geofile()
					
			except Exception as e:
				print(f"build_embgchain error:{e}")

		print(f"Number of events in the MuBack sample {self.embg_chain.GetEntries()}")
		# Disable all branches then enable only those used in the analysis
		self.embg_chain.SetBranchStatus("*", 0)
		self.embg_chain.SetBranchStatus("vetoPoint*", 1)
		self.embg_chain.SetBranchStatus("Digi_SBTHits*", 1)
		self.embg_chain.SetBranchStatus("MCTrack*", 1)
		self.embg_chain.SetCacheSize(10000000)  # 10 MB cache, adjust as needed
		

	def append_vetoPoints(self,vetoPoints):

		for aMCPoint in vetoPoints:

			detID=aMCPoint.GetDetectorID()
			Eloss=aMCPoint.GetEnergyLoss()

			if detID not in self.ElossPerDetId:
			    self.ElossPerDetId[detID]=0
			    self.tOfFlight[detID]=[]

			self.ElossPerDetId[detID] += Eloss

			hittime = aMCPoint.GetTime()

			self.tOfFlight[detID].append(hittime)

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
		veto_table_header = ['Signal Event Index', "Veto@45MeV",  f"nSBTHits>45MeV",f"nSBTHits>0 MeV",f"nSBTHits>10 MeV",f"nSBTHits>20 MeV",f"nSBTHits>30 MeV",f"nSBTHits>60 MeV",f"nSBTHits>90 MeV"]

		# Once before the loop
		bg_event_t0s = np.array([e["t0"] for e in self.bg_eventtimes])
		bg_event_entries = np.array([e["entry"] for e in self.bg_eventtimes])
		
		veto_table=[]

		for signalNr,signal_t0 in enumerate(self.signal_eventtimes):
			
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
								
				self.append_vetoPoints(self.embg_event.vetoPoint)

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
			
			
			combined_Digi_SBTHits=self.digitizecombinedSBT(self.random.Rndm()*1*u.microsecond)#candidate_event.ShipEventHeader.GetEventTime())
			
			print(tabulate(combine_table, headers=combine_table_header, tablefmt="pretty"))

			self.h["ncombinedSBThits"].Fill(len(combined_Digi_SBTHits))

			veto, w, hitSegments={},{},{}
			
			for threshold in self.threshold_list:

				veto[threshold], w[threshold], hitSegments[threshold]=self.SBT_decision(combined_Digi_SBTHits,threshold=threshold)
						
				if veto[threshold]:
					#self.vetoed_cases[threshold]+=1
					self.vetoed_cases[self.iterationNr][threshold].add(signalNr)#+=1
				#self.n_cases[threshold]+=1
				self.n_cases[self.iterationNr][threshold].add(signalNr)#+=1

			veto_table.append([
								signalNr,  # Signal Event Index
								veto[45],
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
		
		for iteration in range(self.nIterations):
			
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
		
		for iteration in range(self.nIterations):
			
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

			nevents_to_check = options.nEvents
		
			print(f"nevents_to_check: {nevents_to_check}")

			self.assign_event_time_signal(nevents_to_check)
			self.check_veto()
			print("------------------------------------------------------------------------------------------------------------------------")

		self.calculate_probability()
		ut.writeHists(self.h,"sbtplots.root")

# Initialize and run the analysis
falsevetoprobability = FalseVetoProbability(options.bg_path)
falsevetoprobability.run_analysis()

