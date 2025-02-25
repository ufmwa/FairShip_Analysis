"""
USAGE<within FairSHiP>:

python MuonBack_hitrates.py --test

"""

#-----------------------------------------------------------------------------------------------------------

import rootUtils as ut
from rootpyPickler import Unpickler
import ROOT, os,sys
from argparse import ArgumentParser
import numpy as np 
from array import array
from tabulate import tabulate
import shipunit as u

PDGData = ROOT.TDatabasePDG.Instance()
parser = ArgumentParser()
parser.add_argument('--test'    	, dest='testing_code' 	, help='Run Test'   , required=False, action='store_true',default=False)
parser.add_argument('--path'     , dest='path' 			, help='path to the MuonBack files'		, required=False, action='store_true',default='/eos/experiment/ship/simulation/bkg/MuonBack_2024helium/8070735')

options = parser.parse_args()


if options.testing_code:    directory = './test_hitrates'    
else: 						directory = './'
	

tag='ECN3_2024_newgeo'
path=options.path


h={}

ut.bookHist(h, 'weights',' ; muon weights; '							,2000,6,6)
ut.bookHist(h, 'ubtpoint_nHits','UBT 	; nHits per event; '	,1000,0,100)
ut.bookHist(h, 'ubtpoint_hitZ', 'UBT 	; z (cm); '				,600,-2600,-2400)
ut.bookHist(h, 'ubtpoint_hitXY','UBT	; x (cm);y (cm);'			,200,-200,200,200,-200,200)
ut.bookHist(h, 'ubtpoint_momentum','UBT	; hit momentum (GeV/c); '		,1000,0,10)

ut.bookHist(h, 'strawtube_detID',   'SST			 	; detectorID; '		,500,0.5e7,5e7)
ut.bookHist(h, 'strawtube_momentum','SST			 	; momentum (GeV/c); ',10000,0,2)
ut.bookHist(h, 'strawtube_momentum_muons','SST (muons) ; momentum (GeV/c); ',10000,0,2)

for i in range(1,5):
	ut.bookHist(h, f'strawtube_station{i}_hitXY'		,f'SST station{i} 			; x (cm) ;y (cm);	'		,100,-300,300,100,-400,400)
	ut.bookHist(h, f'strawtube_station{i}_muons_hitXY'	,f'SST station{i} (muons)	; x (cm) ;y (cm);	'		,100,-300,300,100,-400,400)

ut.bookHist(h, 'strawtubepoint_pdg'				,'SST		; pdg  		;'				,20,-0.5,20.5)
ut.bookHist(h, 'strawtubepoint_pdg_vs_mom'		,'SST		; pdg 		; momentum (GeV/c)'	,20,-0.5,20.5,100000,1e-9,1)

ut.bookHist(h, 'vetopoint_topology_phi'					,'SBT (vetoPoint info) hitrate ; z(cm) ; #phi 	'			,100,-3000.,3000.,36,0,360)

ut.bookHist(h, 'vetopoint_spatial_dist'					,'SBT (vetoPoint info) position of particle hit within the SBT cell ; x(cm) ; z(cm)	; y(cm) '				,100,-200.,200.,100,-50,50,100,-300.,300.)

ut.bookHist(h, 'vetopoint_energydeposition'				,'SBT (vetoPoint info)					; Energy deposition per particle hit(GeV); 		'		,1000,0.,1)
ut.bookHist(h, 'vetopoint_energydeposition_shapewise'	,'SBT (vetoPoint info)					; Shape ID; Energy deposition per particle hit(GeV)	;'	,6,0.5,6.5,1000,0,1)#2D plot
ut.bookHist(h, 'vetopoint_pdg'							,'SBT (vetoPoint info)					; pdg  		; 	'				,20,-0.5,20.5)
ut.bookHist(h, 'vetopoint_pdg_vs_energydeposition'		,'SBT (vetoPoint info)					; pdg 		; Energy deposition (GeV) '		,20,-0.5,20.5,1000,0,1)
ut.bookHist(h, 'vetopoint_min_energydeposition_muons'			,'SBT (vetoPoint info) min(energy deposition) of muons 	; z(cm) ; #phi ;energy deposition(GeV)'	,100,-3000.,3000.,36,0,360)

h['vetopoint_min_energydeposition_muons'].GetZaxis().SetTitleOffset(-0.5);  
h['vetopoint_min_energydeposition_muons'].GetZaxis().SetTitleSize(0.03);   

threshold_list=[0,10,20,30,45,50,60,90]

for threshold in threshold_list:
	
	ut.bookHist(h, f'{threshold}_digihit_topology'				,f'SBT (Digitised hits @ {threshold}MeV threshold ) hitrate	; x(cm)	; z(cm)	; y(cm) 	',100,-500.,500.,100,-3000.,3000.,100,-500.,500.)
	ut.bookHist(h, f'{threshold}_digihit_topology_phi'			,f'SBT (Digitised hits @ {threshold}MeV threshold ) hitrate 	; z(cm) ; #phi ;			',100,-3000.,3000.,36,0,360)
	
	ut.bookHist(h, f'{threshold}_vetopoint_multiplicity'		,f'SBT (Digitised hits @ {threshold}MeV threshold )	vetoPoint multiplicity		;Number of vetoPoints hitting the SBT cell	; 	',100,0,100)
	ut.bookHist(h, f'{threshold}_z_vs_vetopoint_multiplicity'	,f'SBT (Digitised hits @ {threshold}MeV threshold )			; z(cm)				;Number of vetoPoints hitting the SBT cell 	',100,-3000.,3000.,100,0,100)

	ut.bookHist(h, f'{threshold}_maxenergydeposition'			,f'SBT (Digitised hits @ {threshold}MeV threshold )			; max(Energy deposition) (GeV)'				,1000,0,1)
	ut.bookHist(h, f'{threshold}_n_maxenergydeposition'			,f'SBT (Digitised hits @ {threshold}MeV threshold )			; Number of cells with max(Energy deposition) '	,100,0,100)
	ut.bookHist(h, f'{threshold}_cell_maxenergydeposition'		,f'SBT (Digitised hits @ {threshold}MeV threshold )			; Max(E_deposition) in a cell (GeV) 	;Number of cells with max(E_deposition) '	,1000,0.,1,50,0,50)
	
	ut.bookHist(h, f'{threshold}_digihit_max_edepval_topology_phi'		,f'SBT (Digitised hits @ {threshold}MeV threshold ) min( max(energy deposition) per event ) 	; z(cm) ; #phi ;energy deposition(GeV)',100,-3000.,3000.,36,0,360)
	
	h[f'{threshold}_digihit_max_edepval_topology_phi'].GetZaxis().SetTitleOffset(-0.5);  
	h[f'{threshold}_digihit_max_edepval_topology_phi'].GetZaxis().SetTitleSize(0.03);   

	ut.bookHist(h, f'{threshold}_digihit_energydeposition'			,f'SBT (Digitised hits @ {threshold}MeV threshold ) ; E_deposition per digihit (GeV) 				;										',1000,0.,1)
	ut.bookHist(h, f'{threshold}_digihit_multiplicity'				,f'SBT (Digitised hits @ {threshold}MeV threshold ) ; Number of triggered SBT cells in an event 	;										',2000,-0.5,2000.5)
	ut.bookHist(h, f'{threshold}_digihit_rate_shapewise'			,f'SBT (Digitised hits @ {threshold}MeV threshold ) ; Shape ID 										; Digitised hit rate					',6,0.5,6.5)
	ut.bookHist(h, f'{threshold}_digihit_rate_cellwise'				,f'SBT (Digitised hits @ {threshold}MeV threshold ) ; Cellwise digihit rate  						; 										',2000,100e3,700e3)
	ut.bookHist(h, f'{threshold}_digihit_energydeposition_shapewise',f'SBT (Digitised hits @ {threshold}MeV threshold ) ; Shape ID 										; Energy deposition per Digihit(GeV)	',6,0.5,6.5,1000,0,1)
	ut.bookHist(h, f'{threshold}_digihit_z'							,f'SBT (Digitised hits @ {threshold}MeV threshold ) ; z(cm) 										; 										',100,-3000.,3000)
	
def Phicalc(x, y):
	"""Calculate the azimuthal angle phi in degrees with a 90° offset."""

	r = ROOT.TMath.Sqrt(x*x + y*y)

	if r == 0:  
	    return np.inf #Prevent division by zero

	if(y>=0):   phi =   ROOT.TMath.ACos(x/r)
	else:       phi =-1*ROOT.TMath.ACos(x/r)+2*ROOT.TMath.Pi()

	phi = phi*180/ ROOT.TMath.Pi()# Convert radians to degrees

	phi = (phi + 90) % 360  #+90 offset to start reading from bottom center, phi range= [0,360)

	return phi  
    

def dump(event,mom_threshold=0):

    headers=['#','particle','pdgcode','mother_id','Momentum [Px,Py,Pz] (GeV/c)','StartVertex[x,y,z] (m)','Process', 'GetWeight()', ]#'px','py','pz','vx','vy','vz'
    
    event_table=[]
    for trackNr,track in enumerate(event.MCTrack): 
        
        if track.GetP()/u.GeV < mom_threshold :  continue
        
        try: particlename=PDGData.GetParticle(track.GetPdgCode()).GetName()
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

def print_result_new(tag):
    global h, Event_weight, SBT_Event_weight, digihitrate, sst_hitrate

    # Save histograms
    ut.writeHists(h, directory + tag + '.root')

    # Open readme file
    with open(directory + tag + '_readme.txt', 'w') as readme:
        
        # Muon BG Statistics
        muon_stats = [
            ["Muon BG Statistics", len(Event_weight), sum(Event_weight.values())],
            ["Muon BG Statistics with SBT activity", len(SBT_Event_weight), sum(SBT_Event_weight.values())]
        ]
        print("\n\n" + tabulate(muon_stats, headers=["Description", "Generated", "Scaled to one spill"], tablefmt="pretty"))
        readme.write(tabulate(muon_stats, headers=["Description", "Generated", "Scaled to one spill"], tablefmt="pretty") + "\n\n")

        # Digihit Multiplicity Section
        print("\n  =========================== SBT ===========================\n")
        readme.write("\n  =========================== SBT ===========================\n\n")

        digihit_table = []
        for threshold in threshold_list:
            total_digihit = round(sum(digihitrate[f'{threshold}MeV'].values()) * 1e-6, 4)
            max_cell, max_digihit = max(digihitrate[f'{threshold}MeV'].items(), key=lambda k: k[1])
            digihit_table.append([f"{threshold} MeV", f"{total_digihit} MHz", f"{max_digihit:.2f} Hz", f"(Detector ID: {max_cell})"])

        print(tabulate(digihit_table, headers=["Threshold", "Total Digihit-Rate", "Max Digihit-Rate", "Max Detector"], tablefmt="pretty"))
        readme.write(tabulate(digihit_table, headers=["Threshold", "Total Digihit-Rate", "Max Digihit-Rate", "Max Detector"], tablefmt="pretty") + "\n\n")

        # Tracking Station Hit Rates
        print("\n  ==================== TRACKING STATIONS ====================\n")
        readme.write("\n  ==================== TRACKING STATIONS ====================\n\n")

        TOTALHITRATE = 0
        muonhitrate = {}

        for station in range(1, 5):
            muonhitrate[station] = 0
            print(f"\n  PLANE {station}\n  ---------------------")
            readme.write(f"\n  PLANE {station}\n  ---------------------\n")

            if station not in sst_hitrate:
                continue

            station_table = []
            for particle_name, hitrate in sst_hitrate[station].items():
                if particle_name.startswith("mu"):
                    muonhitrate[station] += float(hitrate)
                station_table.append([particle_name, f"{hitrate:.5f} Hz"])

            print(tabulate(station_table, headers=["Particle", "Hitrate"], tablefmt="pretty"))
            readme.write(tabulate(station_table, headers=["Particle", "Hitrate"], tablefmt="pretty") + "\n\n")

            station_hitrate = sum(sst_hitrate[station].values())
            print(f"  Total hitrate in station {station}: {station_hitrate:.5f} Hz\n")
            readme.write(f"  Total hitrate in station {station}: {station_hitrate:.5f} Hz\n\n")
            TOTALHITRATE += station_hitrate

        # Summary of Total Hit Rates
        print("\n  =========================================================\n")
        total_summary = [
            ["Total hitrate in all trackers", f"{TOTALHITRATE:.5f} Hz"],
            ["Total muon hitrate in all trackers", f"{sum(muonhitrate.values()):.5f} Hz"]
        ]
        print(tabulate(total_summary, tablefmt="pretty"))
        readme.write("\n" + tabulate(total_summary, tablefmt="pretty") + "\n")
        print("\n  =========================================================\n")
        readme.write("\n  =========================================================\n")


def print_result(tag):
	
	global h,Event_weight,SBT_Event_weight,digihitrate,sst_hitrate
	ut.writeHists(h, directory +tag+'.root')

	with open(directory +tag+'_readme.txt', 'w') as readme: 
				
		print("\n\n\n")
		

		
		print(" {:46} Generated: {:10.5}\t Scaled to one spill: {:10.5}".format('Muon BG Statistics',float(len(Event_weight)),float(sum(Event_weight.values()))))
		print(" {:46} Generated: {:10.5}\t Scaled to one spill: {:10.5}".format('Muon BG Statistics with SBT activity',float(len(SBT_Event_weight)),float(sum(SBT_Event_weight.values()))))
		
		print("\n\n  ================================================================================")
		print("\n  SBT")
		print("  --------------------------------------------------------------------------------\n\n")
		

		
		readme.write("\n {:46} Generated: {:10.5}\t Scaled to one spill: {:10.5}".format('Muon BG Statistics',float(len(Event_weight)),float(sum(Event_weight.values()))))
		readme.write("\n {:46} Generated: {:10.5}\t Scaled to one spill: {:10.5}".format('Muon BG Statistics with SBT activity',float(len(SBT_Event_weight)),float(sum(SBT_Event_weight.values()))))
		
		readme.write("\n\n  ================================================================================")
		readme.write("\n  SBT")
		readme.write("\n  --------------------------------------------------------------------------------\n\n")
		
		
		print("  {:10}\t {:10}\t\t{:10}".format('THRESHOLD','TOTAL DIGIHIT-RATE',"MAXIMUM DIGIHIT-RATE IN A CELL"))
		print("  --------------------------------------------------------------------------------\n")
		
		readme.write( "\n\n  {:10}\t {:10}\t\t{:10}\n------------------------------------------------------------------------------------------------------------------------\n".format('THRESHOLD','TOTAL DIGIHIT-RATE',"MAXIMUM DIGIHIT-RATE IN A CELL"))
		for threshold in threshold_list:
			print(         " {:5} MeV\t {:10} MHz\t\t {:10.2} Hz --->{:10}".format(threshold,round( sum(digihitrate[f'{threshold}MeV'].values())*1e-6,4),max(digihitrate[f'{threshold}MeV'].items(), key=lambda k: k[1])[1]," ( Detector ID: "+str(max(digihitrate[f'{threshold}MeV'].items(), key=lambda k: k[1])[0])+")"))
			readme.write("\n {:5} MeV\t {:10} MHz\t\t {:10.2} Hz --->{:10}".format(threshold,round( sum(digihitrate[f'{threshold}MeV'].values())*1e-6,4),max(digihitrate[f'{threshold}MeV'].items(), key=lambda k: k[1])[1]," ( Detector ID: "+str(max(digihitrate[f'{threshold}MeV'].items(), key=lambda k: k[1])[0])+")"))
		
		print(       "\n  --------------------------------------------------------------------------------\n")
		readme.write("\n  --------------------------------------------------------------------------------\n")


		print("\n\n  ================================================================================")
		print("\n  TRACKING STATIONS")
		print("  ---------------------------------------------------------\n")
		print("\n  {:5} {:20}      {:10}".format("","","HITRATE"))		
		print( "  {:5}  {:20}     {:10}".format("","","--------"))		
		readme.write("\n\n  ================================================================================\n")
		readme.write("\n\n  {:5} {:20}      {:10}".format("","","HITRATE"))		
		readme.write("\n  {:5}  {:20}     {:10}".format("","","--------"))		

		TOTALHITRATE=0
		muonhitrate={}
		for station in range(1,5):
			if station not in muonhitrate:muonhitrate[station]=0
			print("   PLANE ",station)
			print( "  {:20}  {:30}     {:10}".format("-------------------","",""))		
			readme.write("\n   PLANE "+str(station))
			readme.write( "\n  {:5}  {:20}     {:10}".format("-------------------","",""))		

			if not station in sst_hitrate: continue
			for particle_name in sst_hitrate[station]:
				if particle_name.startswith("mu"):muonhitrate[station]+=float(sst_hitrate[station][particle_name])
				print( 			"  {:5} {:20}  {:10.5} Hz".format("",particle_name,float(sst_hitrate[station][particle_name])))	
				readme.write( "\n  {:5} {:20}  {:10.5} Hz".format("",particle_name,float(sst_hitrate[station][particle_name])))	
			print("  -----------------------------------------\n")
			print(" {:10}  {:10.5} Hz\n\n".format(" Total hitrate in station "+str(station)+" =",float(sum(sst_hitrate[station].values()))))			
			readme.write("\n  -----------------------------------------\n")
			readme.write("\n {:10}  {:10.5} Hz\n\n".format(" Total hitrate in station "+str(station)+" =",float(sum(sst_hitrate[station].values()))))			
			TOTALHITRATE+=sum(sst_hitrate[station].values())
		
		print("  ---------------------------------------------------------\n")
		print( " {:20}{:10}{:10.5} Hz".format(" Total hitrate in all the trackers:","",float(TOTALHITRATE)))
		print( " {:20}{:5}{:10.5} Hz".format(" Total muon hitrate in all the trackers:","",float(sum(muonhitrate.values()))))
		print("\n  ================================================================================")
		
		readme.write("\n  ================================================================================\n")
		readme.write("\n  {:20}{:10}{:10.5} Hz".format(" Total hitrate in all the trackers:","",float(TOTALHITRATE)))
		readme.write("\n  {:20}{:5}{:10.5} Hz".format(" Total muon hitrate in all the trackers:","",float(sum(muonhitrate.values()))))
		readme.write("\n  ================================================================================\n")
			
def print_SBTcell_relative_pos(vetoPoint):
    # Initialize the navigator
    nav = ROOT.gGeoManager.GetCurrentNavigator()
    if not nav:
        print("Navigator could not be initialized.")
        return

    # Define master (global) coordinates for the vetoPoint
    master_point = array('d', [vetoPoint.GetX(), vetoPoint.GetY(), vetoPoint.GetZ()])

    # Set the current point to the master coordinates in the navigator
    nav.SetCurrentPoint(master_point[0], master_point[1], master_point[2])

    # Find the node at this point, which moves the navigator to the correct volume
    current_node = nav.FindNode()
    if not current_node:
        print("No node found at the given coordinates.")
        return

    # Print the name of the node where the point is located
    #print(f"Found node: {current_node.GetName()}")

    # Prepare an array to hold the local coordinates
    local_coords = array('d', [0, 0, 0])

    # Perform the transformation from master (global) to local coordinates
    nav.MasterToLocal(master_point, local_coords)

    # Print the result
    #print(f"Master coordinates: {master_point}")
    #print(f"Local coordinates in '{nav.GetCurrentNode().GetName()}': {local_coords}")

    return local_coords[0],local_coords[2],local_coords[1]

def Main_function():	
	
	global h,Event_weight,SBT_Event_weight,digihitrate,sst_hitrate
	
	Event_weight,SBT_Event_weight ={},{}
	digihitrate ={}
	total_particlehitrate=0
	files=0
	global_event_id=-1
	sbt_pdg_list={}
	sst_pdg_list={}
	sbt_pdg_index = 0
	sst_pdg_index = 0
	sst_hitrate = {}

	min_maxEloss_array = {} 
	for threshold in threshold_list:
		min_maxEloss_array[threshold]=np.full((100, 36), np.inf)  # Create a 2D array or dictionary to store minimum eLoss values per (z, phi) bin, initialized with inf
	
	muon_min_eloss_array = np.full((100, 36), np.inf)  # Initialize with infinity
	
	sGeo=None

	exception_issues={}

	for jobDir in os.listdir(path):
		
		#if jobDir.startswith('README'): continue
		f,fgeo=None,None
		
		try:             
		
			inputFile=f'{path}/{jobDir}/ship.conical.MuonBack-TGeant4_rec.root' 
					
			f = ROOT.TFile.Open(inputFile)
			tree = f.cbmsim

			geoFile=f'{path}/{jobDir}/geofile_full.conical.MuonBack-TGeant4.root'
			fgeo = ROOT.TFile(geoFile)
			
			if not sGeo:
				upkl    = Unpickler(fgeo) #load Shipgeo dictionary written by run_simSfcript.py
				ShipGeo = upkl.load('ShipGeo')
				sGeo   = fgeo.FAIRGeom

			if options.testing_code and files>6: break

			print(files,jobDir)
			
			files+=1
						
			for eventNr,event in enumerate(tree):

				global_event_id+=1

				for track in event.MCTrack:	
					if (track.GetPdgCode() in [-13,13]):
						Event_weight[global_event_id]=track.GetWeight()		
						h['weights'].Fill(Event_weight[global_event_id])
						break		
				
				weight=Event_weight[global_event_id]
				
				#------------------------UBT------------------------------------------------
				
				h['ubtpoint_nHits'].Fill(len(event.UpstreamTaggerPoint),weight)	
				
				mom = ROOT.TVector3()
				for ubt_MCPoint in event.UpstreamTaggerPoint:
					ubt_MCPoint.Momentum(mom)
					h[ 'ubtpoint_momentum'].Fill(mom.Mag(), weight)
					h[ 'ubtpoint_hitXY'	  ].Fill(ubt_MCPoint.GetX(),ubt_MCPoint.GetY(), weight)
					h[ 'ubtpoint_hitZ' 	  ].Fill(ubt_MCPoint.GetZ(), weight)
				
				#------------------------------------------------------------------------
				#-----------------------SST----------------------------------------------
				
				for sst_MCPoint in event.strawtubesPoint:
					
					detID = sst_MCPoint.GetDetectorID()					
					Eloss = sst_MCPoint.GetEnergyLoss()
					x,y,z = sst_MCPoint.GetX(), sst_MCPoint.GetY(), sst_MCPoint.GetZ()
					P = ROOT.TMath.Sqrt(sst_MCPoint.GetPx() ** 2 + sst_MCPoint.GetPy() ** 2 + sst_MCPoint.GetPz() ** 2)
					
					station = detID // 10000000  # will need to update if strawtubes detID format is changed

					if station not in sst_hitrate: sst_hitrate[station]={}
						
					pdgCode = event.MCTrack[sst_MCPoint.GetTrackID()].GetPdgCode()
					
					try:	
						particle_name = PDGData.GetParticle(pdgCode).GetName()
					except:	
						particle_name = 'others'
						
					if particle_name not in sst_pdg_list.values():
					    sst_pdg_list[sst_pdg_index] = particle_name
					    sst_pdg_index+=1
					    
					if particle_name not in sst_hitrate[station]:
						sst_hitrate[station][particle_name]=0

					sst_pdg_key = [k for k, v in sst_pdg_list.items() if v == particle_name][0]

					h[ 'strawtube_momentum'					].Fill(P, weight)
					h[ 'strawtube_detID'					].Fill(detID, weight)
					h[ f'strawtube_station{station}_hitXY'	].Fill(x,y, weight)
					h[ 'strawtubepoint_pdg'					].Fill(sst_pdg_key, weight)
					h[ 'strawtubepoint_pdg_vs_mom'			].Fill(sst_pdg_key,P, weight)
					
					if particle_name.startswith('mu'):
						h[ f'strawtube_station{station}_muons_hitXY'].Fill(x,y, weight)
						h[ 'strawtube_momentum_muons'				 ].Fill(P, weight)
					
					sst_hitrate[station][particle_name]+=Event_weight[global_event_id]

				#------------------------------------------------------------------------
				#-----------------------SBT----------------------------------------------


				ElossPerDetId    = {}
				listOfVetoPoints = {}				
				#tOfFlight        = {}

				for key,veto_MCPoint in enumerate(event.vetoPoint):
					
					if global_event_id not in SBT_Event_weight:
						SBT_Event_weight[global_event_id]=Event_weight[global_event_id] #saving nEvents with SBT activity
					
					total_particlehitrate+=Event_weight[global_event_id] 

					pdgCode = event.MCTrack[veto_MCPoint.GetTrackID()].GetPdgCode()
					detID 	= veto_MCPoint.GetDetectorID()
					shape_nr= detID//100000

					vetopoint_z,vetopoint_x,vetopoint_y = veto_MCPoint.GetZ(),veto_MCPoint.GetX(),veto_MCPoint.GetY()
					
					Eloss = veto_MCPoint.GetEnergyLoss()

					if detID not in ElossPerDetId: 
						ElossPerDetId[detID]=0
						listOfVetoPoints[detID]=[]
						#tOfFlight[detID]=[]
						
					ElossPerDetId[detID] += Eloss
					listOfVetoPoints[detID].append(key)
					#tOfFlight[detID].append(veto_MCPoint.GetTime())
					
					try:	particle_name=PDGData.GetParticle(pdgCode).GetName()
					except:	particle_name='others'

					if particle_name not in sbt_pdg_list.values():
					    sbt_pdg_list[sbt_pdg_index]=particle_name
					    sbt_pdg_index+=1

					veto_pdg_key = [k for k, v in sbt_pdg_list.items() if v == particle_name][0]

					h[ 'vetopoint_pdg'							].Fill(veto_pdg_key,weight)
					h[ 'vetopoint_pdg_vs_energydeposition'		].Fill(veto_pdg_key,Eloss,weight)
					h[ 'vetopoint_energydeposition'			 	].Fill(Eloss,weight)
					h[ 'vetopoint_energydeposition_shapewise'	].Fill(shape_nr,Eloss,weight)
					h[ 'vetopoint_topology_phi'				 	].Fill(vetopoint_z,Phicalc(vetopoint_x,vetopoint_y),weight) 

					if particle_name.startswith('mu'):

						z_bin 	= h['vetopoint_min_energydeposition_muons'].GetXaxis().FindBin(vetopoint_z)
						phi_bin = h['vetopoint_min_energydeposition_muons'].GetYaxis().FindBin(Phicalc(vetopoint_x,vetopoint_y))
						
						if (Eloss/0.001) < muon_min_eloss_array[z_bin-1, phi_bin-1]:  # -1 to adjust for array index
							muon_min_eloss_array[z_bin-1, phi_bin-1] = Eloss/0.001
						
					relative_x,relative_z,relative_y=print_SBTcell_relative_pos(veto_MCPoint)
					h['vetopoint_spatial_dist'].Fill(relative_x,relative_z,relative_y,weight)
								
				
				#Explicit Digitisation 
				digiSBT={}
				#digihit_multiplicity={} 
				digihit_multiplicity = {threshold: 0 for threshold in threshold_list} #nSBT cells fired per event

				for index,detID in enumerate(ElossPerDetId):
					
					aHit = ROOT.vetoHit(detID,ElossPerDetId[detID])
					#aHit.SetTDC(min( tOfFlight[detID] ) + event.t0 )
					#if ElossPerDetId[detID]<0.045:    aHit.setInvalid()  
					
					digiSBT[index] = aHit
									
					for threshold in threshold_list:
						
						#if threshold not in digihit_multiplicity:
						#	digihit_multiplicity[threshold]=0
						
						if ElossPerDetId[detID]<0.001*threshold:	continue
						
						if f'{threshold}MeV' not in digihitrate: 		
							digihitrate[f'{threshold}MeV']={}
						
						if detID not in digihitrate[f'{threshold}MeV']: 	
							digihitrate[f'{threshold}MeV'][detID]=0
						
						digihit_multiplicity[threshold] 	 +=1
						digihitrate[f'{threshold}MeV'][detID]+=Event_weight[global_event_id]
						
						h[f'{threshold}_vetopoint_multiplicity'			].Fill(len(listOfVetoPoints[detID]),weight) 	#how many vetopoints per digitised hit	
						h[f'{threshold}_z_vs_vetopoint_multiplicity'	].Fill(aHit.GetZ(),len(listOfVetoPoints[detID]),weight)


				#Reading Digitised Data

				#maxeLoss={}
				#for threshold in threshold_list:
				#	maxeLoss[threshold]=-1

				maxeLoss = {threshold: -1 for threshold in threshold_list} #maximum energy deposition percell
				
				nmaxcells={}
				max_z={}
				max_phi={}

				for aHit in digiSBT.values():
						
					#if not aHit.isValid(): continue #only hit which pass the threshold of 45 MeV

					x 		=aHit.GetX()
					y 		=aHit.GetY()
					z 	 	=aHit.GetZ()
					eLoss  	=aHit.GetEloss()
					detID 	=aHit.GetDetectorID()
					shape_nr=int(ROOT.TMath.Floor(detID/100000))
					
					for threshold in threshold_list:
					
						if eLoss<0.001*threshold: continue

						if eLoss>maxeLoss[threshold]:
							nmaxcells[threshold]=0
							maxeLoss[threshold]= eLoss
							max_z[threshold]   = z
							max_phi[threshold] = Phicalc(x,y)

						#print(ElossPerDetId[detID],maxeLoss)
						
						if eLoss==maxeLoss[threshold]:
							nmaxcells[threshold]+=1
						
						h[ f'{threshold}_digihit_topology'					].Fill(x,z,y,weight)
						h[ f'{threshold}_digihit_topology_phi'				].Fill(z,Phicalc(x,y),weight)
						h[ f'{threshold}_digihit_z'							].Fill(z,weight)
						h[ f'{threshold}_digihit_rate_shapewise'			].Fill(shape_nr,weight)
						h[ f'{threshold}_digihit_rate_cellwise'				].Fill(detID,weight)
						h[ f'{threshold}_digihit_energydeposition'			].Fill(eLoss,weight)
						h[ f'{threshold}_digihit_energydeposition_shapewise'].Fill(shape_nr,eLoss,weight)

						# Update the min_maxEloss_array if the current eLoss is smaller than the stored value
						#h[ f'{threshold}_digihit_max_edepval_topology_phi'			  ].Fill(z,Phicalc(x,y),eLoss/0.001)
						#print('ElossPerDetId[detID]',maxeLoss,nmaxcells)

				for threshold in threshold_list:

					h[f'{threshold}_digihit_multiplicity'].Fill(digihit_multiplicity[threshold],weight) #how many digihits per event
					
					if maxeLoss[threshold]==-1: continue
				
					h[f'{threshold}_cell_maxenergydeposition'	].Fill(maxeLoss[threshold],nmaxcells[threshold],weight)#2Dplot
					h[f'{threshold}_n_maxenergydeposition'		].Fill(nmaxcells[threshold],weight)
					h[f'{threshold}_maxenergydeposition'		].Fill(maxeLoss[threshold],weight)

					z_bin 	= h[f'{threshold}_digihit_max_edepval_topology_phi'].GetXaxis().FindBin(max_z[threshold])
					phi_bin = h[f'{threshold}_digihit_max_edepval_topology_phi'].GetYaxis().FindBin(max_phi[threshold])
					
					if maxeLoss[threshold]/0.001 < min_maxEloss_array[threshold][z_bin-1, phi_bin-1]:  # -1 to adjust for array index
						min_maxEloss_array[threshold][z_bin-1, phi_bin-1] = maxeLoss[threshold]/0.001
			f.Close()
			fgeo.Close()
		except Exception as e:      
			if f:
				f.Close() 
			if fgeo:
				fgeo.Close() 
			#print(f"Except called for (Reason :{e})")
			exception_issues[jobDir]=e
			continue		
	for key in h:
		h[key].SetOption('HIST')
	
	for index in sbt_pdg_list:
	
	    bin_id = index +1
	    h[ 'vetopoint_pdg'			].GetXaxis().SetBinLabel(bin_id,str(sbt_pdg_list[index]))	
	    h[ 'vetopoint_pdg_vs_energydeposition'].GetXaxis().SetBinLabel(bin_id,str(sbt_pdg_list[index]))	
	
	for index in sst_pdg_list:
	    
	    bin_id = index +1
	    h[ 'strawtubepoint_pdg'			 ].GetXaxis().SetBinLabel(bin_id,str(sst_pdg_list[index]))	
	    h[ 'strawtubepoint_pdg_vs_mom'].GetXaxis().SetBinLabel(bin_id,str(sst_pdg_list[index]))	

	
	# Fill the histogram with the minimum eLoss values
	for z_bin in range(1,101):
		for phi_bin in range(1,37):
			for threshold in threshold_list:
			    min_eloss = min_maxEloss_array[threshold][z_bin-1, phi_bin-1]
			    if min_eloss != np.inf:  # Only fill if there's a valid min eLoss
			        h[f'{threshold}_digihit_max_edepval_topology_phi'].SetBinContent(z_bin, phi_bin, min_eloss)
			    
			min_eloss_veto = muon_min_eloss_array[z_bin-1, phi_bin-1]
			if min_eloss_veto != np.inf:  # Only fill if there's a valid min eLoss
			    h['vetopoint_min_energydeposition_muons'].SetBinContent(z_bin, phi_bin, min_eloss_veto)	


	print_result(tag)
	print('Exceptions:\n', exception_issues)

Main_function()


