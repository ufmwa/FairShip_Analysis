#!/usr/bin/env python3
"""
Merge *_pos TTrees across many jobs and generate one figure per selection (tree name).
- CLI matches your CSV merge script flags.
- Saves plots into a dedicated folder: plots_<sample>/.
- Filenames use ONLY the TTree name (sanitized), e.g.:
    plots_muonDIS_fullreco/all_preselection_UBT_BasicSBT@45MeV_PID_inv_mass_pos.png
"""

import sys, glob, re
from pathlib import Path
from argparse import ArgumentParser

import ROOT
ROOT.gROOT.SetBatch(True)

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Polygon, Rectangle

# ---------------- style & geometry ----------------
def load_style():
    # Keep your house style for cohesive look
    style_path = Path('/afs/cern.ch/user/a/anupamar/Analysis/Tools/thesis.mplstyle')
    try:
        plt.style.use(str(style_path))
    except Exception:
        plt.style.use('default')

def add_geometry(ax_zx, ax_zy, ax_xy, show_detectors=True):
    """Draw the SHiP DV views exactly like your eventDisplay_mini."""
    z_start, z_end = -2500, 2500

    # x vs z (top view)
    ax_zx.add_patch(Polygon([(z_start,-74),(z_start,74),(z_end,224),(z_end,-224)],
                            fill=True, facecolor='#e6f2fa', edgecolor='#6baed6',
                            linewidth=2, linestyle='-'))
    # y vs z (side view)
    ax_zy.add_patch(Polygon([(z_start,-159),(z_start,159),(z_end,324),(z_end,-324)],
                            fill=True, facecolor='#e6f2fa', edgecolor='#6baed6',
                            linewidth=2, linestyle='-'))
    # x vs y (back view): mirrored trapezoids + connectors
    back = [(-74,-159),(-74,159),(-224,324),(-224,-324)]
    backR = [(-x,y) for x,y in back]
    ax_xy.add_patch(Polygon(back,  fill=True, facecolor='#e6f2fa', edgecolor='#6baed6', linewidth=2))
    ax_xy.add_patch(Polygon(backR, fill=True, facecolor='#e6f2fa', edgecolor='#6baed6', linewidth=2))
    ax_xy.add_patch(Polygon([(-74,159),(74,159),(224,324),(-224,324)],
                            fill=True, facecolor='#e6f2fa', edgecolor='#6baed6', linewidth=2))
    ax_xy.add_patch(Polygon([(-74,-159),(74,-159),(224,-324),(-224,-324)],
                            fill=True, facecolor='#e6f2fa', edgecolor='#6baed6', linewidth=2))

    if not show_detectors:
        return

    detectors = {
        'Tr1_1':[2588,2608,-241.57,241.57,-333.28,333.28],
        'Tr2_2':[2788,2808,-241.57,241.57,-333.28,333.28],
        'ShipMagnet_1':[2902,3234,-274.66,274.66,-335.93,335.93],
        'Tr3_3':[3328,3348,-274.66,274.66,-335.93,335.93],
        'Tr4_4':[3528,3548,-274.66,274.66,-335.93,335.93],
        'Timing Detector_1':[3605,3610,-259,259,-331,331],
    }
    colors = {'Tr1_1':'#e6550d','Tr2_2':'#e6550d','ShipMagnet_1':'gray',
              'Tr3_3':'#e6550d','Tr4_4':'#e6550d','Timing Detector_1':'#31a354'}
    for n,(z0,z1,x0,x1,y0,y1) in detectors.items():
        ax_zx.add_patch(Rectangle((z0,x0), z1-z0, x1-x0, fill=False,
                                  edgecolor=colors.get(n,'black'), linewidth=2, alpha=0.4))
        ax_zy.add_patch(Rectangle((z0,y0), z1-z0, y1-y0, fill=False,
                                  edgecolor=colors.get(n,'black'), linewidth=2, alpha=0.4))

def safe_filename(s: str) -> str:
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', s)

# ---------------- CLI (same flags as your CSV merge) ----------------
parser = ArgumentParser(description=__doc__)
g1 = parser.add_mutually_exclusive_group(required=True)
g1.add_argument("--muonDIS",             action="store_true", help="Summarise muonDIS background studies")
g1.add_argument("--neuDIS",              action="store_true", help="Summarise neuDIS background studies")
g1.add_argument("--muonDIS_fullreco",    action="store_true", help="Summarise muonDIS background (fully reco.)")
g1.add_argument("--muonDIS_leptonrho",   action="store_true", help="Summarise muonDIS background (leptonrho)")
g1.add_argument("--neuDIS_fullreco",     action="store_true", help="Summarise neuDIS background (fully reco.)")
g1.add_argument("--neuDIS_leptonrho",    action="store_true", help="Summarise neuDIS background (leptonrho)")
g1.add_argument("--mupi",                action="store_true", help="Summarise full. reco studies")
g1.add_argument("--erho",                action="store_true", help="Summarise erho studies")
g1.add_argument("--murho",               action="store_true", help="Summarise murho studies")
g1.add_argument("--mumuv",               action="store_true", help="Summarise partial. reco studies")

g2 = parser.add_mutually_exclusive_group(required=True)
g2.add_argument("--all",         dest="all",        action="store_true", help="All interactions")
g2.add_argument("--vesselCase",  action="store_true", help="Interactions in SBT Vessel")
g2.add_argument("--heliumCase",  action="store_true", help="Interactions in DecayVolume (helium)")

parser.add_argument("--test", dest="testing_code", action="store_true", default=False,
                    help="Process a small subset of files (quick check).")
parser.add_argument("--dump", action="store_true",
                    help="Also write a CSV summary of merged counts per selection.")
opts = parser.parse_args()

# ---------------- resolve sample → paths/tag ----------------
#main_path = '/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/'
pathlist, sample_tag = [], None

if opts.muonDIS:
    main_path = '/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/'
    pathlist   = [main_path+'/muonDIS/SBT', main_path+'/muonDIS/Tr']
    sample_tag = "muonDIS"
elif opts.muonDIS_fullreco:
    main_path = '/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/'
    pathlist   = [main_path+'/muonDIS_fullreco/SBT', main_path+'/muonDIS_fullreco/Tr']
    sample_tag = "muonDIS_fullreco"
elif opts.neuDIS:
    main_path = '/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/'
    pathlist   = [main_path+'/neuDIS/']
    sample_tag = "neuDIS"
elif opts.neuDIS_fullreco:
    main_path = '/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/'
    pathlist   = [main_path+'/neuDIS_fullreco/']
    sample_tag = "neuDIS_fullreco"
elif opts.neuDIS_leptonrho:
    main_path = '/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/'
    pathlist   = [main_path+'/neuDIS_leptonrho/']
    sample_tag = "neuDIS_leptonrho"
elif opts.muonDIS_leptonrho:
    main_path = '/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/'
    pathlist   = [main_path+'/muonDIS_leptonrho/SBT', main_path+'/muonDIS_leptonrho/Tr']
    sample_tag = "muonDIS_leptonrho"
elif opts.mupi:
    main_path = '/eos/experiment/ship/user/anupamar/BackgroundStudies/'
    pathlist   = [main_path+'/mupi_EventCalc/']
    sample_tag = "mupi"
elif opts.mumuv:
    main_path = '/eos/experiment/ship/user/anupamar/BackgroundStudies/'
    pathlist   = [main_path+'/2muv_EventCalc/']
    sample_tag = "mumuv"
elif opts.erho:
    main_path = '/eos/experiment/ship/user/anupamar/BackgroundStudies/'
    pathlist   = [main_path+'/erho_EventCalc/']
    sample_tag = "erho"
elif opts.murho:
    main_path = '/eos/experiment/ship/user/anupamar/BackgroundStudies/'
    pathlist   = [main_path+'/murho_EventCalc/']
    sample_tag = "murho"

if not sample_tag:
    print("Unknown sample flag."); sys.exit(1)

if opts.all: keyword = "all"
elif opts.vesselCase: keyword = "vesselCase"
elif opts.heliumCase: keyword = "heliumCase"
else:
    print("Pick one of --all/--vesselCase/--heliumCase."); sys.exit(1)

# Output directory
outdir = Path(f"plots_{sample_tag}")
outdir.mkdir(parents=True, exist_ok=True)

# ---------------- merge reader ----------------
def collect_groups(pathlist, keyword, limit_files=None):
    """
    Return dict:
      groups[treename] = {
          'pts': {'reco': {'zx':[Z,X], 'zy':[Z,Y], 'xy':[X,Y]},
                  'ip':   {'zx':[Z,X], 'zy':[Z,Y], 'xy':[X,Y]}},
          'has_ip': bool,
          'files': set([...])
      }
    """
    groups = {}
    files = []
    for base in pathlist:
        files.extend(glob.glob(f"{base}/job_*/selectionparameters_{keyword}.root"))
    files = sorted(files)
    if limit_files:
        files = files[:limit_files]

    if not files:
        print(f"(no selectionparameters_{keyword}.root found under {pathlist})")
        return groups

    for path in files:
        f = ROOT.TFile.Open(path, "READ")
        if not f or f.IsZombie():
            print(f"! skip unreadable: {path}")
            continue
        for key in f.GetListOfKeys():
            obj = key.ReadObj()
            if not obj.InheritsFrom("TTree"):
                continue
            name = obj.GetName()
            if not name.endswith("_pos"):
                continue   # only the position TTrees

            t = obj
            brs = [b.GetName() for b in t.GetListOfBranches()]
            has_ip = all(b in brs for b in ("IP_x","IP_y","IP_z"))

            if name not in groups:
                groups[name] = {
                    'pts': {
                        'reco': {'zx':[[],[]], 'zy':[[],[]], 'xy':[[],[]]},
                        'ip':   {'zx':[[],[]], 'zy':[[],[]], 'xy':[[],[]]},
                    },
                    'has_ip': has_ip,
                    'files': set(),
                }
            g = groups[name]
            g['has_ip'] = g['has_ip'] and has_ip
            g['files'].add(path)

            # Fill points
            for ev in t:
                x, y, z = float(ev.x), float(ev.y), float(ev.z)
                g['pts']['reco']['zx'][0].append(z); g['pts']['reco']['zx'][1].append(x)
                g['pts']['reco']['zy'][0].append(z); g['pts']['reco']['zy'][1].append(y)
                g['pts']['reco']['xy'][0].append(x); g['pts']['reco']['xy'][1].append(y)
                if has_ip:
                    IPx, IPy, IPz = float(ev.IP_x), float(ev.IP_y), float(ev.IP_z)
                    g['pts']['ip']['zx'][0].append(IPz); g['pts']['ip']['zx'][1].append(IPx)
                    g['pts']['ip']['zy'][0].append(IPz); g['pts']['ip']['zy'][1].append(IPy)
                    g['pts']['ip']['xy'][0].append(IPx); g['pts']['ip']['xy'][1].append(IPy)
        f.Close()
    return groups

# ---------------- plotting ----------------
def draw_one_merged(name, group, outdir: Path):
    load_style()
    pts = group['pts']

    fig = plt.figure(figsize=(15,7), constrained_layout=True)
    gs  = gridspec.GridSpec(2,8, figure=fig)

    ax_zx = fig.add_subplot(gs[0,0:5]); ax_zx.set_title('top view')
    ax_zx.set_ylabel('x (cm)'); ax_zx.set_xlim(-3000,3800); ax_zx.set_ylim(-600,600)

    ax_zy = fig.add_subplot(gs[1,0:5]); ax_zy.set_title('side view')
    ax_zy.set_xlabel('z (cm)'); ax_zy.set_ylabel('y (cm)')
    ax_zy.set_xlim(-3000,3800); ax_zy.set_ylim(-600,600)

    ax_xy = fig.add_subplot(gs[0:2,5:7]); ax_xy.set_title('back view')
    ax_xy.set_xlabel('x (cm)'); ax_xy.set_ylabel('y (cm)')
    ax_xy.set_xlim(-250,250);  ax_xy.set_ylim(-400,400)

    add_geometry(ax_zx, ax_zy, ax_xy, show_detectors=True)

    def scat(ax, X, Y, **kw):
        if X: ax.scatter(X, Y, **kw)

    # True IP (black star)
    scat(ax_zx, pts['ip']['zx'][0], pts['ip']['zx'][1], marker='*', s=100, color='black')
    scat(ax_zy, pts['ip']['zy'][0], pts['ip']['zy'][1], marker='*', s=100, color='black')
    scat(ax_xy, pts['ip']['xy'][0], pts['ip']['xy'][1], marker='*', s=100, color='black', label='DIS Interaction Point')

    # Reco vertex (red x)
    scat(ax_zx, pts['reco']['zx'][0], pts['reco']['zx'][1], color='red', marker='x', alpha=1, s=100)
    scat(ax_zy, pts['reco']['zy'][0], pts['reco']['zy'][1], color='red', marker='x', alpha=1, s=100)
    scat(ax_xy, pts['reco']['xy'][0], pts['reco']['xy'][1], color='red', marker='x', alpha=1, s=100, label='Recon. Candidate Vertex')

    # Legend (dedup)
    handles, labels = ax_xy.get_legend_handles_labels()
    uniq = dict(zip(labels, handles))
    fig.legend(list(uniq.values()), list(uniq.keys()),
               loc='lower center', ncol=2, bbox_to_anchor=(0.45,-0.08),
               columnspacing=1.5, fontsize=14, frameon=True)

    # Save using ONLY the treename into plots_<sample>/
    safe = safe_filename(name)
    base = outdir / safe
    fig.savefig(str(base) + ".png", dpi=300, bbox_inches='tight', pad_inches=0.4)
    #fig.savefig(str(base) + ".pdf",            bbox_inches='tight', pad_inches=0.4)
    plt.close(fig)
    print(f"saved: {base}.png  (from {len(group['files'])} files)")

# ---------------- run ----------------
limit = 5 if opts.testing_code else None
groups = collect_groups(pathlist, keyword, limit_files=limit)

if not groups:
    print("(no matching *_pos trees found)")
    sys.exit(0)

# Plot per selection (tree name)
for name, grp in sorted(groups.items()):
    n_reco = len(grp['pts']['reco']['xy'][0])
    n_ip   = len(grp['pts']['ip']['xy'][0])
    if n_reco + n_ip == 0:
        continue
    draw_one_merged(name, grp, outdir)

# Optional CSV summary (saved in the same folder)
if opts.dump:
    import csv
    outcsv = outdir / f"merged_counts_{sample_tag}_{keyword}.csv"
    with open(outcsv, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["treename","n_reco","n_ip","n_files"])
        for name, g in sorted(groups.items()):
            w.writerow([name,
                        len(g['pts']['reco']['xy'][0]),
                        len(g['pts']['ip']['xy'][0]),
                        len(g['files'])])
    print("[dump] wrote", outcsv)

exit()

#!/usr/bin/env python3
# Read *_pos trees (x,y,z[,w,IP_x,IP_y,IP_z]) and draw z:x, z:y, x:y
# one FIGURE PER TREE. Filenames include the tree name.

import glob, re
import argparse
from pathlib import Path
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Polygon, Rectangle

ROOT.gROOT.SetBatch(True)

def load_style():
    tools_dir = Path('/afs/cern.ch/user/a/anupamar/Analysis/Tools')
    style_path = tools_dir / 'thesis.mplstyle'
    try:
        plt.style.use(str(style_path))
    except Exception:
        plt.style.use('default')

def add_geometry(ax_zx, ax_zy, ax_xy, show_detectors=True):
    z_start, z_end = -2500, 2500
    side_coords_x = [(z_start, -74), (z_start,  74), (z_end, 224), (z_end, -224)]
    ax_zx.add_patch(Polygon(side_coords_x, fill=True, facecolor='#e6f2fa',
                            edgecolor='#6baed6', linewidth=2, linestyle='-'))
    side_coords_y = [(z_start, -159), (z_start, 159), (z_end, 324), (z_end, -324)]
    ax_zy.add_patch(Polygon(side_coords_y, fill=True, facecolor='#e6f2fa',
                            edgecolor='#6baed6', linewidth=2, linestyle='-'))

    back_coords = [(-74, -159), (-74, 159), (-224, 324), (-224, -324)]
    back_coords_mirror = [(-x, y) for x, y in back_coords]
    ax_xy.add_patch(Polygon(back_coords, fill=True, facecolor='#e6f2fa',
                            edgecolor='#6baed6', linewidth=2, linestyle='-'))
    ax_xy.add_patch(Polygon(back_coords_mirror, fill=True, facecolor='#e6f2fa',
                            edgecolor='#6baed6', linewidth=2, linestyle='-'))
    top_patch = Polygon([(-74, 159), (74, 159), (224, 324), (-224, 324)],
                        fill=True, facecolor='#e6f2fa', edgecolor='#6baed6',
                        linewidth=2, linestyle='-')
    bottom_patch = Polygon([(-74, -159), (74, -159), (224, -324), (-224, -324)],
                           fill=True, facecolor='#e6f2fa', edgecolor='#6baed6',
                           linewidth=2, linestyle='-')
    ax_xy.add_patch(top_patch); ax_xy.add_patch(bottom_patch)

    if not show_detectors:
        return
    detectors = {
        'Tr1_1':             [2588.0, 2608.0,  -241.57, 241.57, -333.28, 333.28],
        'Tr2_2':             [2788.0, 2808.0,  -241.57, 241.57, -333.28, 333.28],
        'ShipMagnet_1':      [2902.0, 3234.0,  -274.66, 274.66, -335.93, 335.93],
        'Tr3_3':             [3328.0, 3348.0,  -274.66, 274.66, -335.93, 335.93],
        'Tr4_4':             [3528.0, 3548.0,  -274.66, 274.66, -335.93, 335.93],
        'Timing Detector_1': [3605.0, 3610.0,  -259.0,  259.0,  -331.0,  331.0],
    }
    colors = {
        'Tr1_1': '#e6550d','Tr2_2': '#e6550d','ShipMagnet_1': 'gray',
        'Tr3_3': '#e6550d','Tr4_4': '#e6550d','Timing Detector_1': '#31a354',
    }
    for name, (z0, z1, x0, x1, y0, y1) in detectors.items():
        ax_zx.add_patch(Rectangle((z0, x0), z1-z0, x1-x0, fill=False,
                                  edgecolor=colors.get(name, 'black'),
                                  linewidth=2, linestyle='-', alpha=0.4, label=name))
        ax_zy.add_patch(Rectangle((z0, y0), z1-z0, y1-y0, fill=False,
                                  edgecolor=colors.get(name, 'black'),
                                  linewidth=2, linestyle='-', alpha=0.4, label=name))

def safe_filename(s):
    """Keep it readable but filesystem-safe."""
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', s)

def read_all_trees(input_glob, cat_filter, contains):
    """
    Yield dicts: { 'file': path, 'cat': cat, 'name': treename, 'title': title, 'pts': {...} }
    pts = {'reco': {'zx': [Z,X], 'zy':[Z,Y], 'xy':[X,Y]}, 'ip': {...} }
    """
    cats_ok = {"all","vesselCase","heliumCase"}
    for path in glob.glob(input_glob):
        f = ROOT.TFile.Open(path, "READ")
        if not f or f.IsZombie():
            print("! Skipping unreadable file:", path); continue
        for key in f.GetListOfKeys():
            obj = key.ReadObj()
            if not obj.InheritsFrom("TTree"): continue
            name = obj.GetName()
            if not name.endswith("_pos"): continue
            if contains and contains not in name: continue

            # Extract category from name prefix
            tree_cat = name.split("_", 1)[0] if "_" in name else "all"
            if tree_cat not in cats_ok: tree_cat = "all"
            if cat_filter != "any" and tree_cat != cat_filter: 
                continue

            t = obj
            branches = [b.GetName() for b in t.GetListOfBranches()]
            has_ip = all(b in branches for b in ("IP_x","IP_y","IP_z"))

            pts = {'reco': {'zx': [[],[]], 'zy': [[],[]], 'xy': [[],[]]},
                   'ip':   {'zx': [[],[]], 'zy': [[],[]], 'xy': [[],[]]}}

            for ev in t:
                x, y, z = float(ev.x), float(ev.y), float(ev.z)
                pts['reco']['zx'][0].append(z); pts['reco']['zx'][1].append(x)
                pts['reco']['zy'][0].append(z); pts['reco']['zy'][1].append(y)
                pts['reco']['xy'][0].append(x); pts['reco']['xy'][1].append(y)
                if has_ip:
                    IPx, IPy, IPz = float(ev.IP_x), float(ev.IP_y), float(ev.IP_z)
                    pts['ip']['zx'][0].append(IPz); pts['ip']['zx'][1].append(IPx)
                    pts['ip']['zy'][0].append(IPz); pts['ip']['zy'][1].append(IPy)
                    pts['ip']['xy'][0].append(IPx); pts['ip']['xy'][1].append(IPy)

            yield {
                'file': path,
                'cat' : tree_cat,
                'name': name,
                'title': "",#t.GetTitle() if hasattr(t, "GetTitle") else "",
                'pts' : pts
            }
        f.Close()

def draw_tree_one_figure(treeinfo, outprefix):
    """Make the 3-panel figure for a single tree."""
    load_style()
    cat   = treeinfo['cat']
    name  = treeinfo['name']
    title = treeinfo['title']
    pts   = treeinfo['pts']

    # Figure & axes (same as your eventDisplay_mini)
    fig = plt.figure(figsize=(15, 7), constrained_layout=True)
    gs  = gridspec.GridSpec(2, 8, figure=fig)
    ax_zx = fig.add_subplot(gs[0, 0:5]); ax_zx.set_title('top view')
    ax_zx.set_ylabel('x (cm)'); ax_zx.set_xlim(-3000, 3800); ax_zx.set_ylim(-600, 600)
    ax_zy = fig.add_subplot(gs[1, 0:5]); ax_zy.set_title('side view')
    ax_zy.set_xlabel('z (cm)'); ax_zy.set_ylabel('y (cm)')
    ax_zy.set_xlim(-3000, 3800); ax_zy.set_ylim(-600, 600)
    ax_xy = fig.add_subplot(gs[0:2, 5:7]); ax_xy.set_title('back view')
    ax_xy.set_xlabel('x (cm)'); ax_xy.set_ylabel('y (cm)')
    ax_xy.set_xlim(-250, 250); ax_xy.set_ylim(-400, 400)

    # Geometry overlays
    add_geometry(ax_zx, ax_zy, ax_xy, show_detectors=True)

    # Plot points (match your marker choices)
    def scatter_pair(ax, X, Y, **kw):
        if len(X): ax.scatter(X, Y, **kw)

    # True IP as black star
    scatter_pair(ax_zx, pts["ip"]["zx"][0], pts["ip"]["zx"][1], marker='*', s=100, color='black')
    scatter_pair(ax_zy, pts["ip"]["zy"][0], pts["ip"]["zy"][1], marker='*', s=100, color='black')
    scatter_pair(ax_xy, pts["ip"]["xy"][0], pts["ip"]["xy"][1], marker='*', s=100, color='black', label='DIS Interaction Point')

    # Reco vertex as red 'x'
    scatter_pair(ax_zx, pts["reco"]["zx"][0], pts["reco"]["zx"][1], color='red', marker='x', alpha=1, s=100)
    scatter_pair(ax_zy, pts["reco"]["zy"][0], pts["reco"]["zy"][1], color='red', marker='x', alpha=1, s=100)
    scatter_pair(ax_xy, pts["reco"]["xy"][0], pts["reco"]["xy"][1], color='red', marker='x', alpha=1, s=100, label='Recon. Candidate Vertex')

    # Legend
    handles, labels = ax_xy.get_legend_handles_labels()
    uniq = dict(zip(labels, handles))
    fig.legend(list(uniq.values()), list(uniq.keys()),
               loc='lower center', ncol=2, bbox_to_anchor=(0.45, -0.08),
               columnspacing=1.5, fontsize=14, frameon=True)

    # Title shows cat + tree title/name
    shown_title = title if title else name
    #fig.suptitle("Category: {}  |  {}".format(cat, shown_title), y=0.99, fontsize=16)

    # Save: include category and tree name in the filename (filesystem-safe)
    safe = safe_filename(name)
    base = "{}_{}_{}".format(outprefix, cat, safe)
    fig.savefig(base + ".png", dpi=300, bbox_inches='tight', pad_inches=0.4)
    #fig.savefig(base + ".pdf",            bbox_inches='tight', pad_inches=0.4)
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default="selectionparameters_*.root", help="Glob of ROOT files")
    ap.add_argument("--cat", choices=["all","vesselCase","heliumCase","any"], default="any")
    ap.add_argument("--contains", default="", help="Substring filter for tree names (optional)")
    ap.add_argument("--outprefix", default="combo_positions", help="Output file prefix")
    args = ap.parse_args()

    count = 0
    for treeinfo in read_all_trees(args.inputs, args.cat, args.contains):
        # Skip empty trees (no reco + no IP)
        pts = treeinfo['pts']
        npts = sum(len(pts[k][pl][0]) for k in ('reco','ip') for pl in ('zx','zy','xy'))
        #if npts == 0:
        #    continue
        draw_tree_one_figure(treeinfo, args.outprefix)
        count += 1

    if count == 0:
        print("(no matching trees found or all empty)")

if __name__ == "__main__":
    main()
