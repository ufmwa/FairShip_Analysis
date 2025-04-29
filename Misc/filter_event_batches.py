#!/usr/bin/env python
"""
Script to only filter and save events with a reconstructed candidate within a given simulation path in batches of 100 events per filtered file.
"""
from argparse import ArgumentParser
import glob
import ROOT
import os

parser = ArgumentParser(description=__doc__);
group = parser.add_mutually_exclusive_group(required=True)

group.add_argument('--path', dest='path'         , help='parent path to simulation files'	, required=False, type=str)
group.add_argument("--test", dest="testing_code" , help="run test"              , required=False, action="store_true")
parser.add_argument(
    "-o", "--output_file",
    dest="output_file",
    help="output file name",
    type=str,
    default="events_withcandidates_filtered.root",
)

options = parser.parse_args()

if options.testing_code:
    test_path='/eos/experiment/ship/user/anupamar/signal/mumunu/'
    options.path=test_path

def filter_events_withcandidates(path,BATCH_SIZE=100):

    n_file=0
    chain_input = ROOT.TChain("cbmsim")
    
    for job_folder in os.listdir(path):
        
        input_file=glob.glob(os.path.join(path,job_folder,"ship.conical*_rec.root"))[0]
        n_file+=1
        if options.testing_code and n_file>5: break
        chain_input.Add(input_file)
        print(f"{input_file} added to TChain")
    
    chain_input.SetBranchStatus("Digi_UpstreamTaggerHits", 0)   # skip the bad branch 
    
    n_events = chain_input.GetEntries()      
    
    print("Tree has", n_events, "entries")
    
    if not n_events: 
        print("Exit Normally")
        exit(0)
    
    output_file = options.output_file
    
    n_events_withcandidate=0
    n_events_in_batch=0
    batch_id=-1
    index=0
    f_out=None
    end_of_file=False

    while not end_of_file:
        
        print(index)
        
        if n_events_in_batch==0:
            
            if f_out:
                f_out.Write()
                f_out.Close()

            batch_id+=1
            
            print(f"\tBatch:{batch_id}")
            
            output_dir  = f"job_{batch_id}"
            os.makedirs(output_dir, exist_ok=True)  #create a folder job_{} and save the file inside it
            
            output_file = os.path.join(output_dir, options.output_file)

            f_out = ROOT.TFile(output_file, "recreate")
            f_out.cd()                
            newtree = chain_input.CloneTree(0)
        
        for i in range(index,n_events):
            
            print(f"\tReading entry:{i}")
            
            chain_input.GetEntry(i)
            if chain_input.Particles.GetEntries()>0:
               n_events_withcandidate+=1 
               n_events_in_batch+=1
               print("\t\t has recon candidate")
               newtree.Fill()      
            if n_events_in_batch>=BATCH_SIZE:
                print("Batch size limit reached, proceeding to generate next batch")
                index=i+1
                n_events_in_batch=0
                break
            
            if i==n_events-1:
                end_of_file=True

    if f_out:
        f_out.Write()
        f_out.Close()

    print(f"Filtering successful.:\n\t saved {n_events_withcandidate} candidate events into {batch_id+1} batches of up to {BATCH_SIZE} events per file.")

def inspect_files(path):
    
    n_events_withcandidate=0
    
    for output_dir in os.listdir(path):
        
        if not output_dir.startswith("job_"): continue

        output_file = os.path.join(output_dir, options.output_file)
        f = ROOT.TFile(output_file, "read")
        newtree = f.cbmsim
        n_events=newtree.GetEntries()
        for i in range(n_events):
            newtree.GetEntry(i)
            if newtree.Particles.GetEntries()>0:
               n_events_withcandidate+=1 
               
    print(f'Number of events with candidates={n_events_withcandidate}')

filter_events_withcandidates(options.path)

#inspect_files("./")

