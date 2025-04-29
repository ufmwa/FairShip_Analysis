#!/usr/bin/env python
"""
Script to only filter and save events with a reconstructed candidate within a given folder
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
    default="filtered_events_withcandidates.root",      
)
options = parser.parse_args()

if options.testing_code:
    test_path='/eos/experiment/ship/user/anupamar/signal/mumunu/'
    options.path=test_path

def filter_events_withcandidates(path):

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
    
    f_out = ROOT.TFile(output_file, "recreate")
    f_out.cd()                
    newtree = chain_input.CloneTree(0)
    n_events_withcandidate=0
    for i in range(n_events):
        chain_input.GetEntry(i)
        if chain_input.Particles.GetEntries()>0:
           n_events_withcandidate+=1 
           newtree.Fill()        
    
    f_out.Write()
    print(f'Filtering successful. \n{n_events_withcandidate} events saved in {output_file}')

def inspect_file(file):

    f = ROOT.TFile(file, "read")
    newtree = f.cbmsim
    n_events_withcandidate=0
    n_events=newtree.GetEntries()
    for i in range(n_events):
        newtree.GetEntry(i)
        if newtree.Particles.GetEntries()>0:
           n_events_withcandidate+=1 
           
    print(f'Number of events with candidates={n_events_withcandidate}, totalentries in filtered file {newtree.GetEntries()}')


filter_events_withcandidates(options.path)

inspect_file(options.output_file)

