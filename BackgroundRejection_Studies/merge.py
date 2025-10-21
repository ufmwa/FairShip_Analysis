#!/usr/bin/env python3
"""Script to combine results from run_<background>.py and interactively query contributing jobs by number."""

import glob
from pathlib import Path
import pandas as pd
from tabulate import tabulate
from argparse import ArgumentParser
import sys


parser = ArgumentParser(description=__doc__)
group1 = parser.add_mutually_exclusive_group(required=True)

group1.add_argument("--muonDIS", dest="muonDIS",help="Summarise muonDIS background studies", action="store_true")
group1.add_argument("--neuDIS",  dest="neuDIS" ,help="Summarise neuDIS background studies",  action="store_true")

group1.add_argument("--muonDIS_fullreco", dest="muonDIS_fullreco",help="Summarise muonDIS background studies for fully reco.", action="store_true")
group1.add_argument("--muonDIS_leptonrho", dest="muonDIS_leptonrho",help="Summarise muonDIS background studies for leptonrho", action="store_true")

group1.add_argument("--neuDIS_fullreco",  dest="neuDIS_fullreco" ,help="Summarise neuDIS background studies  for fully reco.",  action="store_true")
group1.add_argument("--neuDIS_leptonrho", dest="neuDIS_leptonrho",help="Summarise neuDIS background studies for leptonrho", action="store_true")

group1.add_argument("--mupi",    dest="mupi"   ,help="Summarise full. reco studies",         action="store_true")

group1.add_argument("--erho",    dest="erho"   ,help="Summarise erho studies",         action="store_true")

group1.add_argument("--murho",    dest="murho"   ,help="Summarise erho studies",         action="store_true")

group1.add_argument("--mumuv",    dest="mumuv"   ,help="Summarise partial. reco studies",      action="store_true")

group2 = parser.add_mutually_exclusive_group(required=True)
group2.add_argument("--all"        , dest="all"       ,help="All interactions"                            , action="store_true")
group2.add_argument("--vesselCase" , dest="vesselCase",help="Interactions only in the SBT Vessel"         , action="store_true")
group2.add_argument("--heliumCase" , dest="heliumCase",help="Interactions only in the DecayVolume (helium)", action="store_true")

parser.add_argument("--test", dest="testing_code" , help="Run Test" , required=False, action="store_true",default=False)
parser.add_argument("--dump", dest="dump", help="Write merged and aggregated DataFrames to CSV", action="store_true")

options = parser.parse_args()

if options.all:         keyword = "all"
elif options.vesselCase:  keyword = "vesselCase"
elif options.heliumCase:  keyword = "heliumCase"


main_path='/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/'

if options.muonDIS:
    pathlist = [
        main_path+'/muonDIS/SBT',
        main_path+'/muonDIS/Tr'
    ]

if options.muonDIS_fullreco:
    pathlist = [
        main_path+'/muonDIS_fullreco/SBT',
        main_path+'/muonDIS_fullreco/Tr'
    ]

if options.neuDIS:
    pathlist = [main_path+'/neuDIS/']

if options.neuDIS_fullreco:
    pathlist = [main_path+'/neuDIS_fullreco/']


if options.neuDIS_leptonrho:
    pathlist = [main_path+'/neuDIS_leptonrho/']

if options.muonDIS_leptonrho:
    pathlist = [
        main_path+'/muonDIS_leptonrho/SBT',
        main_path+'/muonDIS_leptonrho/Tr'
    ]

if options.mupi:
    pathlist = [main_path+'/mupi_EventCalc/']

if options.mumuv:
    pathlist = [main_path+'/2muv_EventCalc/']

if options.erho:
    pathlist = [main_path+'/erho_EventCalc/']

if options.murho:
    pathlist = [main_path+'/murho_EventCalc/']


# tag definitions
pre_tags = [
    "n_particles", "fiducial", "dist2innerwall", "dist2vesselentrance",
    "impact_par", "doca", "n_dof", "reduced_chi2", "d_mom", "preselection",
]
veto_tags = [
    "BasicSBT@45MeV", "BasicSBT@90MeV", "BasicSBT@0MeV",
    "AdvSBT@45MeV", "AdvSBT@90MeV","GNNSBT@45MeV", "TOFSBT@45MeV", "TOFSBT@90MeV", "UBT",
]
combinedveto_tags = ["UBT+BasicSBT@45MeV","UBT+BasicSBT@90MeV","UBT+AdvSBT@45MeV","UBT+AdvSBT@90MeV","UBT+GNNSBT@45MeV","UBT+TOFSBT@45MeV","UBT+TOFSBT@90MeV"]
combined_Basic45 = ["preselection+UBT","preselection+UBT+BasicSBT@45MeV","preselection+UBT+BasicSBT@45MeV+PID","preselection+UBT+BasicSBT@45MeV+PID+inv_mass"]
combined_Basic90 = ["preselection+UBT","preselection+UBT+BasicSBT@90MeV","preselection+UBT+BasicSBT@90MeV+PID","preselection+UBT+BasicSBT@90MeV+PID+inv_mass"]
combined_Adv45   = ["preselection+UBT","preselection+UBT+AdvSBT@45MeV","preselection+UBT+AdvSBT@45MeV+PID","preselection+UBT+AdvSBT@45MeV+PID+inv_mass"]
combined_GNN45   = ["preselection+UBT","preselection+UBT+GNNSBT@45MeV","preselection+UBT+GNNSBT@45MeV+PID","preselection+UBT+GNNSBT@45MeV+PID+inv_mass"]
combined_Adv90   = ["preselection+UBT","preselection+UBT+AdvSBT@90MeV","preselection+UBT+AdvSBT@90MeV+PID","preselection+UBT+AdvSBT@90MeV+PID+inv_mass"]
combined_TOF45   = ["preselection+UBT","preselection+UBT+TOFSBT@45MeV","preselection+UBT+TOFSBT@45MeV+PID","preselection+UBT+TOFSBT@45MeV+PID+inv_mass"]
combined_TOF90   = ["preselection+UBT","preselection+UBT+TOFSBT@90MeV","preselection+UBT+TOFSBT@90MeV+PID","preselection+UBT+TOFSBT@90MeV+PID+inv_mass"]
combined_45   = ["preselection+UBT","preselection+UBT+[AdvSBT + GNNSBT ]@45MeV","preselection+UBT+[AdvSBT + GNNSBT ]@45MeV+ PID","preselection+UBT+[AdvSBT + GNNSBT ]@45MeV+ PID+inv_mass"]

table_specs = [
    ("Combined vetosystem efficiency (UBT x SBT)", combinedveto_tags),
    ("Ordered cuts (BasicSBT@45 MeV threshold)", combined_Basic45),
    ("Ordered cuts (BasicSBT@90 MeV threshold)", combined_Basic90),
    ("Ordered cuts (AdvSBT@45 MeV threshold)",   combined_Adv45),
    ("Ordered cuts (AdvSBT@90 MeV threshold)",   combined_Adv90),
    ("Ordered cuts (GNNSBT@45 MeV threshold)",   combined_GNN45),
    ("Ordered cuts (TOFSBT@45 MeV threshold)",  combined_TOF45),
    ("Ordered cuts (TOFSBT@90 MeV threshold)",  combined_TOF90),
    ("Ordered cuts ([AdvSBT + GNNSBT] @45MeV threshold)",  combined_45),
]

def fmt(value, denom, pct_fmt=".2f"):
    if denom:
        pct = 100 * value / denom
        return f"{value:{'.2e'}} ({pct:{pct_fmt}} %)"
    return f"{value:{'.2e'}}"


def load_csvs(pathlist, keyword):
    csvs = [p for base in pathlist
                  for p in Path(base).glob(f"job_*/selection_summary_{keyword}*.csv")]

    if not csvs:
        raise FileNotFoundError(f"No CSVs matching {keyword!r} under {pathlist}")

    frames = []
    for f in csvs:
        try:
            frames.append(pd.read_csv(f))
        except (FileNotFoundError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
            print(f"[warn] skipped {f}: {e}")

    if not frames:
        raise ValueError("All discovered CSVs were empty or unreadable.")

    return pd.concat(frames, ignore_index=True)


df = load_csvs(pathlist, keyword)

agg = df.groupby('tag')[['nCandidates','nEvents15y']].sum().sort_index()

if options.dump:
    raw_out = "merged_rows.csv"
    agg_out = "aggregated_by_tag.csv"
    df.to_csv(raw_out, index=False)
    agg.reset_index().to_csv(agg_out, index=False)
    
    print(f"[debug] Wrote raw merged rows → {raw_out}")
    print(f"[debug] Wrote aggregated summary → {agg_out}")

rec_nc, rec_n15 = agg.loc["reconstructed", ["nCandidates", "nEvents15y"]]
sim_nc, sim_n15 = agg.loc["simulated", ["nCandidates", "nEvents15y"]]

# print summary blocks
def _block(title, tags):
    rows = []
    for t in tags:
        nc = int(agg.at[t,'nCandidates']) if t in agg.index else 0
        n15 = agg.at[t,'nEvents15y'] if t in agg.index else 0.0

        if t == "simulated" : 
            rows.append([t, 
                        f"{nc:.2e}",
                        f"{n15:.2e}"])
        
        elif t == "reconstructed" :                       
            
            rows.append([t,
                         fmt(nc,  sim_nc),       # nEvents
                         fmt(n15, sim_n15)])     # nEvents15y
                         
        else:
            # all other rows: % relative to 'reconstructed'
            rows.append([t,
                         fmt(nc,  rec_nc),       # nEvents
                         fmt(n15, rec_n15)])     # nEvents15y


    if title=='Event Stats':
        print(tabulate(rows, headers=[title,'nEvents generated(/nSim in %)','nEvents in 15 y (/nSim in %)'], tablefmt='rounded_grid', floatfmt='.2e'))
    else:
        print(tabulate(rows, headers=[title,'nEvents generated(/nReco in %)','nEvents in 15 y (/nReco in %)'], tablefmt='rounded_grid', floatfmt='.2e'))
    print()

# 1) Event stats
_block('Event Stats',['simulated','reconstructed'])

# 2) Pre-selection
_block('Pre-Selection Cut', pre_tags)

# 3) Veto systems
_block('Veto System', veto_tags)

# 4) Combined / ordered
for title, tags in table_specs:
    _block(title, tags)

# 5) Other cuts
printed = set(['simulated','reconstructed']) | set(pre_tags) | set(veto_tags) | set(t for _,g in table_specs for t in g)
others = [t for t in agg.index if t not in printed]
if others: _block('Other cuts', others)
"""
#---------------------------------------------------------------------------------------------------------------------------------------
# indexed menu clustered for further interactive debugging
all_tags = []
clustered = [('Event Stats',['simulated','reconstructed']),
             ('Pre-Selection',pre_tags),('Veto',veto_tags)] + table_specs + [('Other cuts',others)]
for _, tags in clustered:
    for t in tags:
        all_tags.append(t)
menu = [[i,t,int(agg.at[t,'nCandidates']) if t in agg.index else 0, f"{agg.at[t,'nEvents15y']}" if t in agg.index else '0.000']
        for i,t in enumerate(all_tags)]
print(tabulate(menu, headers=['#','tag','nEvents generated','nEvents15y'], tablefmt='rounded_grid'))
print()

# prompt

try:
    choice = input('Enter tag number (blank to exit): ').strip()
    if not choice: sys.exit(0)
    idx = int(choice)
    tag = all_tags[idx]
except:
    print('Invalid.'); sys.exit(1)
print(f"Selected [{idx}]: {tag}\n")

# lookup
sel = df[(df.tag==tag)&(df.nCandidates>0)][['job',"nCandidates", "nEvents15y"]].reset_index(drop=True)
if sel.empty:
    print('No jobs with non-zero for',tag)
else:
    sel["nCandidates"] = sel["nCandidates"].apply(lambda x: fmt(x, rec_nc))
    sel["nEvents15y"]  = sel["nEvents15y"].apply(lambda x: fmt(x, rec_n15))
    print(tabulate(sel.values.tolist(), headers=['job','nEvents generated','nEvents15y'], tablefmt='rounded_grid', floatfmt='.2e'))

#---------------------------------------------------------------------------------------------------------------------------------------
"""
#---------------------------------------------------------------------------------------------------------------------------------------
# indexed menu clustered for further interactive debugging

def _build_all_tags():
    all_tags = []
    clustered = [('Event Stats',['simulated','reconstructed']),
                 ('Pre-Selection', pre_tags), ('Veto', veto_tags)] + table_specs + [('Other cuts', others)]
    for _, tags in clustered:
        for t in tags:
            all_tags.append(t)
    return all_tags

def _build_menu(all_tags):
    return [[i, t,
             int(agg.at[t, 'nCandidates']) if t in agg.index else 0,
             f"{agg.at[t, 'nEvents15y']}" if t in agg.index else '0.000']
            for i, t in enumerate(all_tags)]

while True:
    all_tags = _build_all_tags()
    menu = _build_menu(all_tags)
    print(tabulate(menu, headers=['#','tag','nEvents generated','nEvents15y'],
                   tablefmt='rounded_grid'))
    print()

    try:
        choice = input('Enter tag number (blank to exit): ').strip()
        if not choice:
            break
        idx = int(choice)
        if idx < 0 or idx >= len(all_tags):
            print('Invalid index.\n')
            continue
        tag = all_tags[idx]
    except (ValueError, KeyboardInterrupt, EOFError):
        print('Invalid.\n')
        continue

    print(f"Selected [{idx}]: {tag}\n")

    # lookup
    sel = df[(df.tag == tag) & (df.nCandidates > 0)][['job', 'nCandidates', 'nEvents15y']].reset_index(drop=True)
    if sel.empty:
        print('No jobs with non-zero for', tag)
    else:
        sel["nCandidates"] = sel["nCandidates"].apply(lambda x: fmt(x, rec_nc))
        sel["nEvents15y"]  = sel["nEvents15y"].apply(lambda x: fmt(x, rec_n15))
        print(tabulate(sel.values.tolist(),
                       headers=['job','nEvents generated','nEvents15y'],
                       tablefmt='rounded_grid', floatfmt='.2e'))

    try:
        again = input("Do you wish to continue? [yes/no] ").strip().lower()
    except (KeyboardInterrupt, EOFError):
        break
    if again in ('', 'y', 'yes'):
        print()
        continue
    else:
        break
#---------------------------------------------------------------------------------------------------------------------------------------
