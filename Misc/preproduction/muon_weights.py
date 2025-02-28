import os
import ROOT
from tabulate import tabulate
import shipunit as u
PDG = ROOT.TDatabasePDG.Instance()

def printMCTrack(n,MCTrack):
    mcp = MCTrack[n]
    # ANSI escape codes for colors
    RED = '\033[91m'  # Red color for highlighting
    RESET = '\033[0m'  # Reset to default color

    try:
       # Get the particle name
       particle_name = PDG.GetParticle(mcp.GetPdgCode()).GetName()

       # Check if the particle is mu+ or mu-
       if particle_name == 'mu+' or particle_name == 'mu-':
            particle_name = f"{RED}{particle_name}{RESET}       "  # Apply red color to particle name

       print(' %6s %-10s %10i %6.3F %6.3F %7.3F %7.3F %7.3F %7.3F %6s %10.3F %28s'%(n,particle_name,mcp.GetPdgCode(),mcp.GetPx()/u.GeV,mcp.GetPy()/u.GeV,mcp.GetPz()/u.GeV, \
                          mcp.GetStartX()/u.m,mcp.GetStartY()/u.m,mcp.GetStartZ()/u.m,mcp.GetMotherId(),mcp.GetWeight(),mcp.GetProcName().Data()    ))
    except: 
        print(' %6s %-10s %10i %6.3F %6.3F %7.3F %7.3F %7.3F %7.3F %6s %10.3F %28s'%(n,"----",mcp.GetPdgCode(),mcp.GetPx()/u.GeV,mcp.GetPy()/u.GeV,mcp.GetPz()/u.GeV, \
                      mcp.GetStartX()/u.m,mcp.GetStartY()/u.m,mcp.GetStartZ()/u.m,mcp.GetMotherId(),mcp.GetWeight(),mcp.GetProcName().Data()    ))

def dump(event,pcut=0,print_whole_event=True):
       #print("print_whole_event",print_whole_event)
       if print_whole_event:
           print('\n %6s %-10s %10s %6s %6s %7s %7s %7s %7s %6s %10s %18s'%( '#','particle','pid','px','py','pz','vx','vy','vz','mid', 'w', 'Process'))
           print(' %6s %10s %10s %6s %6s %7s %7s %7s %7s %6s %10s %18s\n '%(' ','--------','---','--','--','--','--','--','--','---','---','-------'))
       n=-1 
       daughters='->'
       for mcp in event.MCTrack: 
         n+=1
         if mcp.GetP()/u.GeV < pcut :  continue
         if print_whole_event: printMCTrack(n,event.MCTrack)
         
         #if mcp.GetMotherId()==-1: mother=PDG.GetParticle(mcp.GetPdgCode()).GetName() #excited proton 9902210
         #if mcp.GetMotherId()==0: daughters = daughters +' '+PDG.GetParticle(mcp.GetPdgCode()).GetName()
       
       return #mother+daughters

path='/eos/experiment/ship/data/Mbias/background-prod-2018/'
data=[]

for index in range(67):#os.listdir(path):
	
	inputFile=f"pythia8_Geant4_10.0_withCharmandBeauty{index*1000}_mu.root"
	
	#if index !=0: continue
	#print(inputFile,index)

	f = ROOT.TFile.Open(           os.path.join(path, inputFile),            "read")
		
	tree=f.cbmsim
	
	muon_weight={}
	nMuons=0
	nEvents=tree.GetEntries()

	for eventNr,event in enumerate(tree,start=1):
		#dump(event)
		#dummy=input("continue?")
		#if eventNr>10: break
		#print(eventNr)
		for track in event.MCTrack:
			pdgcode=track.GetPdgCode()
			weight=track.GetWeight()

			if not abs(pdgcode)==13: continue
			nMuons+=1
			if eventNr not in muon_weight:
				muon_weight[eventNr]=[]	
			
			muon_weight[eventNr].append(weight)
		
	
	first_index_sum = sum(arr[0] for arr in muon_weight.values() if arr)
	#first_index_sum = sum(arrays[0] for arrays in muon_weight.values() if arrays)

	total_sum = sum(value for arrays in muon_weight.values() for value in arrays)
	
	data.append([index,nEvents,nMuons,first_index_sum,total_sum])

	
	


	
            
	
# Print the table using tabulate
headers=["index","nEvents","nMuons","sum(weight)","sum(weight*nMuons)"]
print(tabulate(data, headers=headers, tablefmt="grid"))
