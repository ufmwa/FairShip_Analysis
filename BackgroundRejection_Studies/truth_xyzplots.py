#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
truth_xyzplots.py
Plot truth interaction vertices (MCTrack[0] start vertex) in SHiP coordinates,
in the same 3-view layout as surviving_xyzplots.py, but BEFORE any cuts.

Typical use for your hypothesis:
- plot air-like truth vertices upstream of the DV entrance (hole region check)
  by selecting:
    --select density:air --select region:UpstreamOfDV --z-window -3300 -2450

Output:
- PNGs in --outputdir
- Optional JSON summary
"""

import os, glob, re, json, random
from pathlib import Path
from argparse import ArgumentParser
from collections import Counter, defaultdict

import ROOT
ROOT.gROOT.SetBatch(True)

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Polygon, Rectangle
import re
# ---------------- style & geometry (borrowed from surviving_xyzplots.py) ----------------
def beamline_segments(zmin, zmax, zstep, x0, y0):
    """Return list of (volname, z_start, z_end) along beamline."""
    segs = []
    last = None
    seg_start = zmin
    z = zmin
    while z <= zmax:
        node = ROOT.gGeoManager.FindNode(x0, y0, z)
        name = "OUTSIDE"
        if node and not ROOT.gGeoManager.IsOutside():
            name = node.GetVolume().GetName()
        if last is None:
            last = name
            seg_start = z
        elif name != last:
            segs.append((last, seg_start, z - zstep))
            last = name
            seg_start = z
        z += zstep
    if last is not None:
        segs.append((last, seg_start, zmax))
    return segs

def halfwidth_for_exact_volume(volname, z, axis, step=2.0, maxw=1200.0):
    """
    Scan outward from (0,0,z) until we leave *exactly* volname.
    Returns halfwidth or None.
    """
    w = 0.0
    while w <= maxw:
        x = w if axis == "x" else 0.0
        y = w if axis == "y" else 0.0
        node = ROOT.gGeoManager.FindNode(x, y, z)
        cur = "OUTSIDE"
        if node and not ROOT.gGeoManager.IsOutside():
            cur = node.GetVolume().GetName()
        if cur != volname:
            break
        w += step
    if w <= step:
        return None
    return w - step

def build_beamline_overlays(segs, include_re, exclude_re, min_dz, hw_step, hw_max):
    """
    For each segment, estimate x/y halfwidth and return drawable overlays:
      [(name,z0,z1,xw,yw), ...]
    """
    overlays = []
    for name, z0, z1 in segs:
        dz = z1 - z0
        if dz < min_dz:
            continue
        if exclude_re and re.search(exclude_re, name):
            continue
        if include_re and not re.search(include_re, name):
            continue

        zmid = 0.5*(z0+z1)
        xw = halfwidth_for_exact_volume(name, zmid, "x", step=hw_step, maxw=hw_max)
        yw = halfwidth_for_exact_volume(name, zmid, "y", step=hw_step, maxw=hw_max)
        if xw is None or yw is None:
            continue
        overlays.append((name, z0, z1, xw, yw))
    return overlays

def draw_beamline_overlays(ax_zx, ax_zy, overlays, labels=False):
    """Draw rectangles for each overlay segment in x–z and y–z."""
    for name, z0, z1, xw, yw in overlays:
        ax_zx.add_patch(Rectangle((z0, -xw), z1-z0, 2*xw, fill=False, linewidth=1, alpha=0.7))
        ax_zy.add_patch(Rectangle((z0, -yw), z1-z0, 2*yw, fill=False, linewidth=1, alpha=0.7))
        if labels:
            ax_zx.text(0.5*(z0+z1), xw, name, fontsize=7, rotation=90, va="bottom", ha="center")


def load_style():
    # keep house style if available; fallback to default
    style_path = Path('/afs/cern.ch/user/a/anupamar/Analysis/Tools/thesis.mplstyle')
    try:
        plt.style.use(str(style_path))
    except Exception:
        plt.style.use('default')

def safe_filename(s: str) -> str:
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', s)

def add_geometry_fixed(ax_zx, ax_zy, ax_xy, show_detectors=True):
    """Hardcoded DV & detector overlay exactly like surviving_xyzplots.py."""
    z_start, z_end = -2500, 2500

    # x vs z (top view)
    ax_zx.add_patch(Polygon([(z_start,-74),(z_start,74),(z_end,224),(z_end,-224)],
                            fill=True, linewidth=2, linestyle='-',
                            edgecolor='#6baed6', facecolor='#e6f2fa'))
    # y vs z (side view)
    ax_zy.add_patch(Polygon([(z_start,-159),(z_start,159),(z_end,324),(z_end,-324)],
                            fill=True, linewidth=2, linestyle='-',
                            edgecolor='#6baed6', facecolor='#e6f2fa'))

    # x vs y (back view): mirrored trapezoids + connectors
    back  = [(-74,-159),(-74,159),(-224,324),(-224,-324)]
    backR = [(-x,y) for x,y in back]
    ax_xy.add_patch(Polygon(back,  fill=True, linewidth=2, edgecolor='#6baed6', facecolor='#e6f2fa'))
    ax_xy.add_patch(Polygon(backR, fill=True, linewidth=2, edgecolor='#6baed6', facecolor='#e6f2fa'))
    ax_xy.add_patch(Polygon([(-74,159),(74,159),(224,324),(-224,324)],
                            fill=True, linewidth=2, edgecolor='#6baed6', facecolor='#e6f2fa'))
    ax_xy.add_patch(Polygon([(-74,-159),(74,-159),(224,-324),(-224,-324)],
                            fill=True, linewidth=2, edgecolor='#6baed6', facecolor='#e6f2fa'))

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

def find_z_range_for_prefix(prefix: str, zmin: float, zmax: float, zstep: float, x0: float, y0: float):
    """Scan along z at fixed (x0,y0), find first/last z where volume startswith(prefix)."""
    z_first, z_last = None, None
    z = zmin
    while z <= zmax:
        node = ROOT.gGeoManager.FindNode(x0, y0, z)
        if node and not ROOT.gGeoManager.IsOutside():
            v = node.GetVolume().GetName()
            if v.startswith(prefix):
                if z_first is None:
                    z_first = z
                z_last = z
        z += zstep
    return z_first, z_last

def find_halfwidth(prefix: str, z: float, axis: str, step: float = 1.0, maxw: float = 4000.0):
    """
    Approximate DV half-width in x or y at given z by scanning outward from 0.
    Returns float or None.
    """
    if axis not in ("x","y"):
        return None
    w = 0.0
    while w <= maxw:
        x = w if axis == "x" else 0.0
        y = w if axis == "y" else 0.0
        node = ROOT.gGeoManager.FindNode(x, y, z)
        if not node or ROOT.gGeoManager.IsOutside():
            break
        v = node.GetVolume().GetName()
        if not v.startswith(prefix):
            break
        w += step
    if w <= step:
        return None
    return w - step

def add_geometry_auto(ax_zx, ax_zy, ax_xy, dv_prefix: str, z_first: float, z_last: float, show_detectors=True):
    """
    Draw DV trapezoids from geometry by estimating half-widths at z_first and z_last.
    Falls back to fixed overlay if auto fails.
    """
    # half widths at entrance/exit
    x0 = find_halfwidth(dv_prefix, z_first, "x")
    x1 = find_halfwidth(dv_prefix, z_last,  "x")
    y0 = find_halfwidth(dv_prefix, z_first, "y")
    y1 = find_halfwidth(dv_prefix, z_last,  "y")

    if any(v is None for v in (x0,x1,y0,y1)):
        add_geometry_fixed(ax_zx, ax_zy, ax_xy, show_detectors=show_detectors)
        return False, None

    # x vs z
    ax_zx.add_patch(Polygon([(z_first,-x0),(z_first,x0),(z_last,x1),(z_last,-x1)],
                            fill=True, linewidth=2, linestyle='-',
                            edgecolor='#6baed6', facecolor='#e6f2fa'))
    # y vs z
    ax_zy.add_patch(Polygon([(z_first,-y0),(z_first,y0),(z_last,y1),(z_last,-y1)],
                            fill=True, linewidth=2, linestyle='-',
                            edgecolor='#6baed6', facecolor='#e6f2fa'))

    # x vs y back view
    back  = [(-x0,-y0),(-x0,y0),(-x1,y1),(-x1,-y1)]
    backR = [( x0,-y0),( x0,y0),( x1,y1),( x1,-y1)]
    ax_xy.add_patch(Polygon(back,  fill=True, linewidth=2, edgecolor='#6baed6', facecolor='#e6f2fa'))
    ax_xy.add_patch(Polygon(backR, fill=True, linewidth=2, edgecolor='#6baed6', facecolor='#e6f2fa'))
    ax_xy.add_patch(Polygon([(-x0,y0),(x0,y0),(x1,y1),(-x1,y1)],
                            fill=True, linewidth=2, edgecolor='#6baed6', facecolor='#e6f2fa'))
    ax_xy.add_patch(Polygon([(-x0,-y0),(x0,-y0),(x1,-y1),(-x1,-y1)],
                            fill=True, linewidth=2, edgecolor='#6baed6', facecolor='#e6f2fa'))

    if show_detectors:
        # keep detector rectangles identical to surviving_xyzplots.py
        add_geometry_fixed(ax_zx, ax_zy, ax_xy, show_detectors=True)

    return True, {"x0":x0,"x1":x1,"y0":y0,"y1":y1}

# ---------------- classification ----------------
def density_bin(rho: float, he_thr: float, air_thr: float) -> str:
    if rho != rho:
        return "nan"
    if rho < he_thr:
        return "he"
    if rho < air_thr:
        return "air"
    return "solid"

def pick_first(patterns, base_dir):
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(base_dir, pat)))
        if hits:
            return hits[0]
    return None

def load_geom(geofile):
    fg = ROOT.TFile.Open(geofile, "READ")
    if not fg or fg.IsZombie():
        raise RuntimeError(f"Cannot open geofile: {geofile}")
    _ = fg.Get("FAIRGeom")
    if not ROOT.gGeoManager:
        raise RuntimeError("gGeoManager not initialised (FAIRGeom missing?)")
    return fg

# ---------------- plotting helpers ----------------
def plot_points(points, outpath, title, dv_info, axes_cfg, show_detectors,
                mark_z=None, plot_mode="auto",
                beam_overlays=None, overlay_labels=False):
    """
    points: dict with keys 'zx','zy','xy' each as [Xlist,Ylist]
            where zx=[Z,X], zy=[Z,Y], xy=[X,Y]
    """
    load_style()

    fig = plt.figure(figsize=(15,7), constrained_layout=True)
    gs  = gridspec.GridSpec(2,8, figure=fig)

    ax_zx = fig.add_subplot(gs[0,0:5]); ax_zx.set_title('top view')
    ax_zx.set_ylabel('x (cm)')
    ax_zx.set_xlim(axes_cfg["zlim"][0], axes_cfg["zlim"][1])
    ax_zx.set_ylim(axes_cfg["xlim"][0], axes_cfg["xlim"][1])

    ax_zy = fig.add_subplot(gs[1,0:5]); ax_zy.set_title('side view')
    ax_zy.set_xlabel('z (cm)'); ax_zy.set_ylabel('y (cm)')
    ax_zy.set_xlim(axes_cfg["zlim"][0], axes_cfg["zlim"][1])
    ax_zy.set_ylim(axes_cfg["ylim"][0], axes_cfg["ylim"][1])

    ax_xy = fig.add_subplot(gs[0:2,5:7]); ax_xy.set_title('back view')
    ax_xy.set_xlabel('x (cm)'); ax_xy.set_ylabel('y (cm)')
    ax_xy.set_xlim(axes_cfg["xylim_x"][0], axes_cfg["xylim_x"][1])
    ax_xy.set_ylim(axes_cfg["xylim_y"][0], axes_cfg["xylim_y"][1])

    # Geometry overlay
    if dv_info["mode"] == "auto" and dv_info.get("z_first") is not None and dv_info.get("z_last") is not None:
        ok, hw = add_geometry_auto(ax_zx, ax_zy, ax_xy, dv_info["prefix"], dv_info["z_first"], dv_info["z_last"], show_detectors=show_detectors)
        if not ok:
            # fallback already drawn inside
            pass
    else:
        add_geometry_fixed(ax_zx, ax_zy, ax_xy, show_detectors=show_detectors)
                    
    # Mark vertical z lines if requested (DV entrance etc.)
    if mark_z:
        for zline in mark_z:
            ax_zx.axvline(zline, linewidth=1)
            ax_zy.axvline(zline, linewidth=1)

    # Decide plot mode
    def draw(ax, X, Y):
        n = len(X)
        if n == 0:
            return
        if plot_mode == "scatter" or (plot_mode == "auto" and n <= 50000):
            ax.scatter(X, Y, s=4, alpha=0.15)
        else:
            # fast & readable for large n
            ax.hexbin(X, Y, gridsize=220, bins='log', mincnt=1)

    draw(ax_zx, points["zx"][0], points["zx"][1])
    draw(ax_zy, points["zy"][0], points["zy"][1])
    draw(ax_xy, points["xy"][0], points["xy"][1])

    fig.suptitle(title, fontsize=14)
    fig.savefig(outpath, dpi=250, bbox_inches='tight', pad_inches=0.4)
    plt.close(fig)

# ---------------- main ----------------
def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument("--overlay-beamline", action="store_true", default=False,
                help="Overlay volumes intersecting the beamline (x=beam_x,y=beam_y) from geofile.")
    ap.add_argument("--overlay-zmin", type=float, default=-6000.0)
    ap.add_argument("--overlay-zmax", type=float, default= 5000.0)
    ap.add_argument("--overlay-zstep", type=float, default=   10.0)
    ap.add_argument("--overlay-include", default=None,
                    help="Regex for volumes to include (e.g. 'Magn|Upstream_Tagger|muon').")
    ap.add_argument("--overlay-exclude", default=r"^(cave|DecayVacuum.*|gas.*|OUTSIDE)$",
                    help="Regex for volumes to exclude.")
    ap.add_argument("--overlay-min-dz", type=float, default=20.0)
    ap.add_argument("--overlay-hw-step", type=float, default=2.0)
    ap.add_argument("--overlay-hw-max", type=float, default=1200.0)
    ap.add_argument("--overlay-labels", action="store_true", default=False)

    ap.add_argument("--inputdir", required=True, help="Directory containing job_* (or a single job directory).")
    ap.add_argument("--outputdir", required=True, help="Output directory for plots.")
    ap.add_argument("--jobs", default="job_*", help="Job glob under inputdir (default: job_*)")
    ap.add_argument("--max-jobs", type=int, default=None)
    ap.add_argument("--max-events", type=int, default=None, help="Per job limit (default: all).")
    ap.add_argument("--max-points", type=int, default=500000, help="Global cap for stored points (per selection).")

    ap.add_argument("--dv-prefix", default="DecayVacuum", help="DV volume name prefix (default: DecayVacuum)")
    ap.add_argument("--zscan-min", type=float, default=-5000.0)
    ap.add_argument("--zscan-max", type=float, default= 5000.0)
    ap.add_argument("--zscan-step", type=float, default=  10.0)
    ap.add_argument("--beam-x", type=float, default=0.0)
    ap.add_argument("--beam-y", type=float, default=0.0)

    ap.add_argument("--he-thr", type=float, default=3e-4, help="rho < he_thr => He-like")
    ap.add_argument("--air-thr", type=float, default=3e-3, help="he_thr <= rho < air_thr => air-like")

    ap.add_argument("--select", action="append", default=[],
                    help="Selections like: density:air|he|solid|any  and/or region:UpstreamOfDV|DV|DownstreamOrSide|any. "
                         "Can be given multiple times. Default is any/any.")
    ap.add_argument("--z-window", nargs=2, type=float, default=None, help="Only keep events with zmin<=z<=zmax")
    ap.add_argument("--subsample", type=int, default=None, help="Randomly keep at most N points (after filters).")
    ap.add_argument("--seed", type=int, default=13)

    ap.add_argument("--geometry", choices=["auto","fixed"], default="auto",
                    help="DV overlay from geofile (auto) or hardcoded (fixed, like surviving_xyzplots).")
    ap.add_argument("--show-detectors", action="store_true", default=False)

    ap.add_argument("--plot-mode", choices=["auto","scatter","hexbin"], default="auto",
                    help="auto chooses scatter for small N, hexbin for large N.")
    ap.add_argument("--write-summary", action="store_true", default=True)
    ap.add_argument("--dump-z-transitions", action="store_true", default=False,
                    help="Scan along z at (beam_x,beam_y) and print volume transitions (helps locate muon-shield end / gap).")

    # axis limits (same defaults as surviving_xyzplots)
    ap.add_argument("--zlim", nargs=2, type=float, default=[-3000, 3800])
    ap.add_argument("--xlim", nargs=2, type=float, default=[-600, 600])
    ap.add_argument("--ylim", nargs=2, type=float, default=[-600, 600])
    ap.add_argument("--xylim-x", nargs=2, type=float, default=[-250, 250])
    ap.add_argument("--xylim-y", nargs=2, type=float, default=[-400, 400])

    args = ap.parse_args()
    random.seed(args.seed)

    in_dir = Path(args.inputdir).resolve()
    out_dir = Path(args.outputdir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    beam_overlays = None  # set near top of main()

    # selection parsing
    sel_density = "any"
    sel_region  = "any"
    for s in args.select:
        if s.startswith("density:"):
            sel_density = s.split(":",1)[1]
        if s.startswith("region:"):
            sel_region = s.split(":",1)[1]

    # job dirs
    if list(in_dir.glob("*_rec.root")):
        job_dirs = [in_dir]
    else:
        job_dirs = sorted([p for p in in_dir.glob(args.jobs) if p.is_dir()])

    if args.max_jobs is not None:
        job_dirs = job_dirs[:max(0, args.max_jobs)]

    if not job_dirs:
        raise RuntimeError(f"No job dirs found in {in_dir} with pattern {args.jobs}")

    # global containers (merged across jobs)
    pts = {"zx":[[],[]], "zy":[[],[]], "xy":[[],[]]}
    counters = Counter()
    region_counts = Counter()
    dv_info_global = {"prefix": args.dv_prefix, "z_first": None, "z_last": None, "mode": args.geometry}

    for jd in job_dirs:
        # discover files
        rec = pick_first(["*_rec.root","ship.*_rec.root","ship.conical*_rec.root"], str(jd))
        geo = pick_first(["geofile*.root","*geofile*.root","geofile_full*.root"], str(jd))
        if not rec or not geo:
            print(f"! skip {jd} (missing rec/geofile)")
            continue

        frec = ROOT.TFile.Open(rec, "READ")
        if not frec or frec.IsZombie():
            print(f"! skip unreadable rec: {rec}")
            continue
        tree = frec.Get("cbmsim")
        if not tree:
            print(f"! skip {rec} (no cbmsim)")
            frec.Close()
            continue

        fgeo = load_geom(geo)
        if args.overlay_beamline and beam_overlays is None:
            segs = beamline_segments(args.overlay_zmin, args.overlay_zmax, args.overlay_zstep, args.beam_x, args.beam_y)
            beam_overlays = build_beamline_overlays(
                segs,
                include_re=args.overlay_include,
                exclude_re=args.overlay_exclude,
                min_dz=args.overlay_min_dz,
                hw_step=args.overlay_hw_step,
                hw_max=args.overlay_hw_max,
            )
        # determine DV z-range by scanning
        z_first, z_last = find_z_range_for_prefix(args.dv_prefix, args.zscan_min, args.zscan_max, args.zscan_step, args.beam_x, args.beam_y)
        if dv_info_global["z_first"] is None:
            dv_info_global["z_first"] = z_first
            dv_info_global["z_last"]  = z_last

        if args.dump_z_transitions:
            # print volume name transitions along z (beamline)
            last = None
            z = args.zscan_min
            print(f"\n--- z transitions for job {jd.name} at (x,y)=({args.beam_x},{args.beam_y}) ---")
            while z <= args.zscan_max:
                node = ROOT.gGeoManager.FindNode(args.beam_x, args.beam_y, z)
                name = "OUTSIDE"
                if node and not ROOT.gGeoManager.IsOutside():
                    name = node.GetVolume().GetName()
                if name != last:
                    print(f"z={z:8.1f}  {name}")
                    last = name
                z += args.zscan_step

        n_entries = int(tree.GetEntries())
        n_scan = n_entries if args.max_events is None else min(n_entries, int(args.max_events))

        v = ROOT.TVector3()

        for i in range(n_scan):
            tree.GetEntry(i)
            tree.MCTrack[0].GetStartVertex(v)
            x, y, z = float(v.X()), float(v.Y()), float(v.Z())

            if args.z_window is not None:
                if not (args.z_window[0] <= z <= args.z_window[1]):
                    continue

            node = ROOT.gGeoManager.FindNode(x, y, z)
            if not node or ROOT.gGeoManager.IsOutside():
                vol = "OUTSIDE"
                rho = float("nan")
            else:
                vol = node.GetVolume().GetName()
                mat = node.GetVolume().GetMaterial()
                rho = float(mat.GetDensity()) if mat else float("nan")

            dens = density_bin(rho, args.he_thr, args.air_thr)

            # region logic (needs z_first)
            in_dv = vol.startswith(args.dv_prefix)
            if z_first is None:
                region = "DV" if in_dv else "NonDV"
            else:
                if in_dv:
                    region = "DV"
                elif z < z_first:
                    region = "UpstreamOfDV"
                else:
                    region = "DownstreamOrSide"

            counters[dens] += 1
            region_counts[region] += 1

            # apply selections
            if sel_density != "any" and dens != sel_density:
                continue
            if sel_region != "any" and region != sel_region:
                continue

            # store points (cap)
            if len(pts["zx"][0]) >= args.max_points:
                continue

            pts["zx"][0].append(z); pts["zx"][1].append(x)
            pts["zy"][0].append(z); pts["zy"][1].append(y)
            pts["xy"][0].append(x); pts["xy"][1].append(y)

        frec.Close()
        fgeo.Close()

    # optional subsample (after filtering)
    if args.subsample is not None and len(pts["zx"][0]) > args.subsample:
        idx = list(range(len(pts["zx"][0])))
        random.shuffle(idx)
        idx = idx[:args.subsample]
        idx_set = set(idx)

        def filt(a):
            return [v for k,v in enumerate(a) if k in idx_set]

        pts = {
            "zx":[filt(pts["zx"][0]), filt(pts["zx"][1])],
            "zy":[filt(pts["zy"][0]), filt(pts["zy"][1])],
            "xy":[filt(pts["xy"][0]), filt(pts["xy"][1])],
        }

    # Build title + output name
    title = f"Truth vertices | density={sel_density} region={sel_region}"
    if args.z_window is not None:
        title += f" | z in [{args.z_window[0]}, {args.z_window[1]}] cm"

    outname = f"truthVtx_density-{sel_density}_region-{sel_region}"
    if args.z_window is not None:
        outname += f"_z{int(args.z_window[0])}to{int(args.z_window[1])}"
    outname = safe_filename(outname) + ".png"
    outpath = out_dir / outname

    axes_cfg = {
        "zlim": args.zlim,
        "xlim": args.xlim,
        "ylim": args.ylim,
        "xylim_x": args.xylim_x,
        "xylim_y": args.xylim_y,
    }

    mark_z = []
    if dv_info_global.get("z_first") is not None:
        mark_z.append(dv_info_global["z_first"])
    if dv_info_global.get("z_last") is not None:
        mark_z.append(dv_info_global["z_last"])

    plot_mode = args.plot_mode
    if plot_mode == "hexbin":
        plot_mode = "hexbin"  # mapped in draw() via else branch
    plot_points(
    pts, str(outpath), title, dv_info_global, axes_cfg,
    show_detectors=args.show_detectors,
    mark_z=mark_z,
    plot_mode=("scatter" if args.plot_mode=="scatter" else "auto"),
    beam_overlays=beam_overlays,
    overlay_labels=args.overlay_labels
    )



    print(f"\nSaved: {outpath}")
    print(f"Stored points (after filters): {len(pts['zx'][0])}")
    print(f"Density counts (before filters): {dict(counters)}")
    print(f"Region counts  (before filters): {dict(region_counts)}")
    #print(f"DV scan z-range: {dv_info_global.get('z_first')} .. {dv_info_global.get('z_last')}  (prefix={args.dv_prefix})")

    if args.write_summary:
        summary = {
            "inputdir": str(in_dir),
            "jobs_used": [p.name for p in job_dirs],
            "selection": {"density": sel_density, "region": sel_region, "z_window": args.z_window},
            "counters_density": dict(counters),
            "counters_region": dict(region_counts),
            "dv_prefix": args.dv_prefix,
            "dv_z_first": dv_info_global.get("z_first"),
            "dv_z_last": dv_info_global.get("z_last"),
            "stored_points": len(pts["zx"][0]),
            "plot": str(outpath),
        }
        jout = out_dir / (safe_filename(outname.replace(".png","")) + ".json")
        with open(jout, "w") as f:
            json.dump(summary, f, indent=2, sort_keys=True)
        print(f"Summary: {jout}")

if __name__ == "__main__":
    main()
