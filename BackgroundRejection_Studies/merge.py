#!/usr/bin/env python3
"""Script to combine results from run_<background>.py."""

#-----------------------------------------------------------------------------------------------------------
import glob
from pathlib import Path
import pandas as pd
from tabulate import tabulate
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)

group1 = parser.add_mutually_exclusive_group(required=True)

group1.add_argument("--muonDIS"   ,dest="muonDIS",help="Summarise muonDIS background studies", action="store_true")
group1.add_argument("--neuDIS"    ,dest="neuDIS" ,help="Summarise neuDIS background studies",  action="store_true")

group2 = parser.add_mutually_exclusive_group(required=True)

group2.add_argument("--all"           ,dest="all"        ,help="All interactions"                               , action="store_true")
group2.add_argument("--vesselCase"    ,dest="vesselCase" ,help="Interactions only in the SBT Vessel"            , action="store_true")
group2.add_argument("--heliumCase"    ,dest="heliumCase" ,help="Interactions only in the DecayVolume(helium)"   , action="store_true")

options = parser.parse_args()


if options.all:         keyword="all"
if options.vesselCase:  keyword="vesselCase"
if options.heliumCase:  keyword="heliumCase"

if options.muonDIS:
    pathlist=['/eos/experiment/ship/user/anupamar/BackgroundStudies/muonDIS/SBT','/eos/experiment/ship/user/anupamar/BackgroundStudies/muonDIS/Tr']
if options.neuDIS:
    pathlist=['/eos/experiment/ship/user/anupamar/BackgroundStudies/neuDIS/']

pre_tags = [
    "n_particles", "fiducial", "dist2innerwall", "dist2vesselentrance",
    "impact_par", "doca", "n_dof", "reduced_chi2", "d_mom", "preselection",
]

veto_tags = [
    "BasicSBT@45MeV", "BasicSBT@90MeV", "BasicSBT@0MeV",
    "AdvSBT@45MeV", "AdvSBT@90MeV", "UBT",
]

combinedveto_tags = [
    "UBT+BasicSBT@45MeV", "UBT+BasicSBT@90MeV",
    "UBT+AdvSBT@45MeV",   "UBT+AdvSBT@90MeV",
]

combined_Basic45 = [
    "preselection+UBT",
    "preselection+UBT+BasicSBT@45MeV",
    "preselection+UBT+BasicSBT@45MeV+PID",
    "preselection+UBT+BasicSBT@45MeV+PID+inv_mass",
]
combined_Basic90 = [
    "preselection+UBT",
    "preselection+UBT+BasicSBT@90MeV",
    "preselection+UBT+BasicSBT@90MeV+PID",
    "preselection+UBT+BasicSBT@90MeV+PID+inv_mass",
]
combined_Adv45 = [
    "preselection+UBT",
    "preselection+UBT+AdvSBT@45MeV",
    "preselection+UBT+AdvSBT@45MeV+PID",
    "preselection+UBT+AdvSBT@45MeV+PID+inv_mass",
]
combined_Adv90 = [
    "preselection+UBT",
    "preselection+UBT+AdvSBT@90MeV",
    "preselection+UBT+AdvSBT@90MeV+PID",
    "preselection+UBT+AdvSBT@90MeV+PID+inv_mass",
]

table_specs = [
    ("Combined vetosystem efficiency (UBT x SBT)", combinedveto_tags),
    ("Ordered cuts (BasicSBT@45 MeV threshold)",   combined_Basic45),
    ("Ordered cuts (BasicSBT@90 MeV threshold)",   combined_Basic90),
    ("Ordered cuts (AdvSBT@45 MeV threshold)",     combined_Adv45),
    ("Ordered cuts (AdvSBT@90 MeV threshold)",     combined_Adv90),
]


csv_files = [
    str(csv_path)
    for base in pathlist
    for csv_path in Path(base).glob(
        f"job_*/selection_summary_{keyword}*.csv")   
]


df = pd.concat((pd.read_csv(f) for f in csv_files), ignore_index=True)

agg = (df.groupby("tag")[["nCandidates", "nEvents15y"]]
         .sum()
         .sort_index())

def _vals(tag):
    """return nCand, n15y (0, 0) if tag never appears"""
    if tag in agg.index:
        row = agg.loc[tag]
        return int(row.nCandidates), row.nEvents15y
    return 0, 0.0


def _block(title, tags, show_title=True):
    rows = []
    for t in tags:
        nC, n15 = _vals(t)
        rows.append([t, nC, f"{n15:.3f}"])
    print(tabulate(rows,
                   headers=[title if show_title else " ",
                            "nCandidates", "nEvents in 15 y"],
                   tablefmt="rounded_grid",
                   floatfmt=".3f"),"\n\n")

# Event statistics --------------------------------------------------------
_block("", ["simulated", "reconstructed"], show_title=False)

# Pre-selection -----------------------------------------------------------
_block("Pre-Selection Cut", pre_tags)

# Veto systems ------------------------------------------------------------
_block("Veto System", veto_tags)

# Combined / ordered tables ----------------------------------------------
for title, tags in table_specs:
    _block(title, tags)

# Everything that never made it into a named table ------------------------
printed = (set(["simulated", "reconstructed"])
           | set(pre_tags) | set(veto_tags)
           | {t for _, g in table_specs for t in g})
others = [t for t in agg.index if t not in printed]

if others:
    _block("Other cuts", others)
