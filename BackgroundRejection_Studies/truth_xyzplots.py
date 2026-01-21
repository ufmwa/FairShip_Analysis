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

Tip:
- Provide --detector-boxes-csv (or keep subsystem_boxes.csv nearby) to auto-show Target/MuonShield
  and enable auto-limits that include upstream components.
- Supported detector boxes CSV schemas:
  * full-bbox: subsystem,xmin_cm,xmax_cm,ymin_cm,ymax_cm,zmin_cm,zmax_cm
  * z-only:   subsystem,zmin_cm,zmax_cm (x/y inferred from axis limits)
"""

import os, glob, re, json, random, csv, sys
from pathlib import Path
from argparse import ArgumentParser, BooleanOptionalAction
from collections import Counter, defaultdict

import ROOT
ROOT.gROOT.SetBatch(True)

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Polygon, Rectangle, Patch
from matplotlib.lines import Line2D

DEFAULT_LIMITS = {
    "zlim": [-3500, 3800],
    "xlim": [-2000, 2000],
    "ylim": [-2000, 2000],
    "xylim_x": [-2000, 2000],
    "xylim_y": [-2000, 2000],
}

# ---------------- load geometry from geofile via (subsystem_boxes.C) ----------------

def load_subsystem_boxes_csv(path):
    """
    Expected columns:
      subsystem,xmin_cm,xmax_cm,ymin_cm,ymax_cm,zmin_cm,zmax_cm,...
      or subsystem,zmin_cm,zmax_cm,...
    """
    boxes = {}
    with open(path, "r", newline="") as fh:
        r = csv.DictReader(fh)
        fieldnames = [f.strip() for f in (r.fieldnames or [])]
        has_full = "xmin_cm" in fieldnames
        has_z = "zmin_cm" in fieldnames and "zmax_cm" in fieldnames
        if not has_z:
            raise RuntimeError(f"Detector boxes CSV missing zmin_cm/zmax_cm: {path}")
        mode = "full-bbox" if has_full else "z-only"
        for row in r:
            name = row["subsystem"]
            entry = {
                "z0": float(row["zmin_cm"]),
                "z1": float(row["zmax_cm"]),
            }
            if has_full:
                entry.update({
                    "x0": float(row["xmin_cm"]),
                    "x1": float(row["xmax_cm"]),
                    "y0": float(row["ymin_cm"]),
                    "y1": float(row["ymax_cm"]),
                })
            boxes[name] = entry
    return boxes, mode

def filter_detector_boxes(boxes, include=None, exclude=None):
    filtered = {}
    for name, b in boxes.items():
        if exclude and re.search(exclude, name):
            continue
        if include and not re.search(include, name):
            continue
        filtered[name] = b
    return filtered

def draw_detector_boxes(ax_zx, ax_zy, ax_xy, boxes, include=None, exclude=None,
                        labels=False, zorder=4, mode="full-bbox"):
    """
    Draws rectangles for each subsystem box.
    include/exclude are regex strings on subsystem names.
    """
    if mode == "z-only":
        zlabel_idx = 0
        for name, b in boxes.items():
            if exclude and re.search(exclude, name):
                continue
            if include and not re.search(include, name):
                continue
            z0, z1 = b["z0"], b["z1"]
            ax_zx.axvspan(z0, z1, alpha=0.08, color="#444444", zorder=zorder)
            ax_zy.axvspan(z0, z1, alpha=0.08, color="#444444", zorder=zorder)
            if labels:
                offset = 10 * (zlabel_idx % 5)
                ax_zx.text(0.5*(z0+z1), ax_zx.get_ylim()[1] + offset, name, fontsize=7,
                           rotation=90, va="bottom", ha="center")
                zlabel_idx += 1
        return

    label_idx = 0
    for name, b in boxes.items():
        if exclude and re.search(exclude, name):
            continue
        if include and not re.search(include, name):
            continue

        z0, z1 = b["z0"], b["z1"]
        x0, x1 = b["x0"], b["x1"]
        y0, y1 = b["y0"], b["y1"]

        ax_zx.add_patch(Rectangle((z0, x0), z1-z0, x1-x0, fill=False, linewidth=1.6, alpha=0.55, zorder=zorder))
        ax_zy.add_patch(Rectangle((z0, y0), z1-z0, y1-y0, fill=False, linewidth=1.6, alpha=0.55, zorder=zorder))
        ax_xy.add_patch(Rectangle((x0, y0), x1-x0, y1-y0, fill=False, linewidth=1.6, alpha=0.55, zorder=zorder))

        if labels:
            offset = 10 * (label_idx % 5)
            ax_zx.text(0.5*(z0+z1), x1 + offset, name, fontsize=7, rotation=90, va="bottom", ha="center")
            label_idx += 1


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

def draw_beamline_overlays(ax_zx, ax_zy, overlays, labels=False, zorder=3):
    """Draw rectangles for each overlay segment in x–z and y–z."""
    label_idx = 0
    for name, z0, z1, xw, yw in overlays:
        ax_zx.add_patch(Rectangle((z0, -xw), z1-z0, 2*xw, fill=False, linewidth=1, alpha=0.7, zorder=zorder))
        ax_zy.add_patch(Rectangle((z0, -yw), z1-z0, 2*yw, fill=False, linewidth=1, alpha=0.7, zorder=zorder))
        if labels:
            offset = 10 * (label_idx % 5)
            ax_zx.text(0.5*(z0+z1), xw + offset, name, fontsize=7, rotation=90, va="bottom", ha="center")
            label_idx += 1


def load_style():
    # keep house style if available; fallback to default
    style_path = Path('/afs/cern.ch/user/a/anupamar/Analysis/Tools/thesis.mplstyle')
    try:
        plt.style.use(str(style_path))
    except Exception:
        plt.style.use('default')

def safe_filename(s: str) -> str:
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', s)

def add_geometry_fixed(
    ax_zx, ax_zy, ax_xy,
    show_detectors=True, show_dv=True, zorder=2,
    dv_edgecolor="#6baed6", dv_facecolor="#e6f2fa", dv_alpha=0.25):
    """Hardcoded DV & detector overlay exactly like surviving_xyzplots.py."""
    z_start, z_end = -2500, 2500
    if show_dv:

        # x vs z (top view)
        ax_zx.add_patch(Polygon([(z_start,-74),(z_start,74),(z_end,224),(z_end,-224)],
                                fill=True, linewidth=2, linestyle='-',
                                edgecolor=dv_edgecolor, facecolor=dv_facecolor, alpha=dv_alpha, zorder=zorder))
        # y vs z (side view)
        ax_zy.add_patch(Polygon([(z_start,-159),(z_start,159),(z_end,324),(z_end,-324)],
                                fill=True, linewidth=2, linestyle='-',
                                edgecolor=dv_edgecolor, facecolor=dv_facecolor, alpha=dv_alpha, zorder=zorder))

        # x vs y (back view): mirrored trapezoids + connectors
        back  = [(-74,-159),(-74,159),(-224,324),(-224,-324)]
        backR = [(-x,y) for x,y in back]
        ax_xy.add_patch(Polygon(back,  fill=True, linewidth=2, edgecolor=dv_edgecolor, facecolor=dv_facecolor, alpha=dv_alpha, zorder=zorder))
        ax_xy.add_patch(Polygon(backR, fill=True, linewidth=2, edgecolor=dv_edgecolor, facecolor=dv_facecolor, alpha=dv_alpha, zorder=zorder))
        ax_xy.add_patch(Polygon([(-74,159),(74,159),(224,324),(-224,324)],
                                fill=True, linewidth=2, edgecolor=dv_edgecolor, facecolor=dv_facecolor, alpha=dv_alpha, zorder=zorder))
        ax_xy.add_patch(Polygon([(-74,-159),(74,-159),(224,-324),(-224,-324)],
                                fill=True, linewidth=2, edgecolor=dv_edgecolor, facecolor=dv_facecolor, alpha=dv_alpha, zorder=zorder))

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
                                  edgecolor=colors.get(n,'black'), linewidth=2, alpha=0.4, zorder=zorder+1))
        ax_zy.add_patch(Rectangle((z0,y0), z1-z0, y1-y0, fill=False,
                                  edgecolor=colors.get(n,'black'), linewidth=2, alpha=0.4, zorder=zorder+1))

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

def add_geometry_auto(
    ax_zx, ax_zy, ax_xy,
    dv_prefix: str, z_first: float, z_last: float,
    show_detectors=True, zorder=2,
    dv_edgecolor="#6baed6", dv_facecolor="#e6f2fa", dv_alpha=0.25):
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
        add_geometry_fixed(
        ax_zx, ax_zy, ax_xy,
        show_detectors=show_detectors, zorder=zorder,
        dv_edgecolor=dv_edgecolor, dv_facecolor=dv_facecolor, dv_alpha=dv_alpha)
        return False, None

    # x vs z
    ax_zx.add_patch(Polygon([(z_first,-x0),(z_first,x0),(z_last,x1),(z_last,-x1)],
                            fill=True, linewidth=2, linestyle='-',
                            edgecolor='#6baed6', facecolor='#e6f2fa', alpha=0.25, zorder=zorder))
    # y vs z
    ax_zy.add_patch(Polygon([(z_first,-y0),(z_first,y0),(z_last,y1),(z_last,-y1)],
                            fill=True, linewidth=2, linestyle='-',
                            edgecolor='#6baed6', facecolor='#e6f2fa', alpha=0.25, zorder=zorder))

    # x vs y back view
    back  = [(-x0,-y0),(-x0,y0),(-x1,y1),(-x1,-y1)]
    backR = [( x0,-y0),( x0,y0),( x1,y1),( x1,-y1)]
    ax_xy.add_patch(Polygon(back,  fill=True, linewidth=2, edgecolor='#6baed6', facecolor='#e6f2fa', alpha=0.25, zorder=zorder))
    ax_xy.add_patch(Polygon(backR, fill=True, linewidth=2, edgecolor='#6baed6', facecolor='#e6f2fa', alpha=0.25, zorder=zorder))
    ax_xy.add_patch(Polygon([(-x0,y0),(x0,y0),(x1,y1),(-x1,y1)],
                            fill=True, linewidth=2, edgecolor='#6baed6', facecolor='#e6f2fa', alpha=0.25, zorder=zorder))
    ax_xy.add_patch(Polygon([(-x0,-y0),(x0,-y0),(x1,-y1),(-x1,-y1)],
                            fill=True, linewidth=2, edgecolor='#6baed6', facecolor='#e6f2fa', alpha=0.25, zorder=zorder))

    if show_detectors:
        add_geometry_fixed(ax_zx, ax_zy, ax_xy, show_detectors=True, show_dv=False, zorder=zorder)


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

def find_detector_boxes_csv(search_dirs):
    patterns = ["subsystem_boxes.csv", "subsystem_boxes*.csv"]
    for d in search_dirs:
        for pat in patterns:
            hits = sorted(Path(d).glob(pat))
            if hits:
                return hits[0]
    return None

def expand_limits(lo, hi, frac=0.08):
    if lo is None or hi is None:
        return None, None
    span = hi - lo
    if span <= 0:
        span = abs(hi) if hi != 0 else 1.0
    pad = span * frac
    return lo - pad, hi + pad

def limits_from_boxes(boxes):
    if not boxes:
        return None
    z0 = min(b["z0"] for b in boxes.values())
    z1 = max(b["z1"] for b in boxes.values())
    zlim = expand_limits(z0, z1)
    if "x0" in next(iter(boxes.values())):
        x0 = min(b["x0"] for b in boxes.values())
        x1 = max(b["x1"] for b in boxes.values())
        y0 = min(b["y0"] for b in boxes.values())
        y1 = max(b["y1"] for b in boxes.values())
        xlim = expand_limits(x0, x1)
        ylim = expand_limits(y0, y1)
        return {"zlim": zlim, "xlim": xlim, "ylim": ylim}
    return {"zlim": zlim}

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
                beam_overlays=None, overlay_labels=False, detector_boxes=None, detector_boxes_cfg=None,
                label_mode="panel", panel_info=None, hexbin_cfg=None,
                point_color=None, point_alpha=0.15, point_size=4.0):

    """
    points: dict with keys 'zx','zy','xy' each as [Xlist,Ylist]
            where zx=[Z,X], zy=[Z,Y], xy=[X,Y]
    """
    load_style()
    if hexbin_cfg is None:
        hexbin_cfg = {"gridsize": 100, "mincnt": 1, "bins": "log"}

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

    # Decide plot mode
    def draw(ax, X, Y):
        n = len(X)
        if n == 0:
            return
        if plot_mode == "hexbin":
            ax.hexbin(
                X, Y,
                gridsize=hexbin_cfg["gridsize"],
                bins=hexbin_cfg["bins"],
                mincnt=hexbin_cfg["mincnt"],
                zorder=1,
            )
        elif plot_mode == "scatter":
            ax.scatter(X, Y, s=point_size, alpha=point_alpha,
                    color=point_color if point_color else None,
                    zorder=1)
        else:  # auto
            if n <= 50000:
                ax.scatter(X, Y, s=4, alpha=0.15, zorder=1)
            else:
                ax.hexbin(
                    X, Y,
                    gridsize=hexbin_cfg["gridsize"],
                    bins=hexbin_cfg["bins"],
                    mincnt=hexbin_cfg["mincnt"],
                    zorder=1,
                )


    draw(ax_zx, points["zx"][0], points["zx"][1])
    draw(ax_zy, points["zy"][0], points["zy"][1])
    draw(ax_xy, points["xy"][0], points["xy"][1])

    # Geometry overlay
    if dv_info["mode"] == "auto" and dv_info.get("z_first") is not None and dv_info.get("z_last") is not None:
        if dv_info["mode"] == "auto" and dv_info.get("z_first") is not None and dv_info.get("z_last") is not None:
            ok, hw = add_geometry_auto(
                ax_zx, ax_zy, ax_xy,
                dv_info["prefix"], dv_info["z_first"], dv_info["z_last"],
                show_detectors=show_detectors, zorder=2,
                dv_edgecolor=dv_info.get("edgecolor", "#6baed6"),
                dv_facecolor=dv_info.get("facecolor", "#e6f2fa"),
                dv_alpha=dv_info.get("alpha", 0.25),
            )
            if not ok:
                pass
        else:
            add_geometry_fixed(
                ax_zx, ax_zy, ax_xy,
                show_detectors=show_detectors, zorder=2,
                dv_edgecolor=dv_info.get("edgecolor", "#6baed6"),
                dv_facecolor=dv_info.get("facecolor", "#e6f2fa"),
                dv_alpha=dv_info.get("alpha", 0.25),
            )


    # optional beamline overlays
    if beam_overlays:
        draw_beamline_overlays(ax_zx, ax_zy, beam_overlays, labels=overlay_labels, zorder=3)

    # optional detector subsystem boxes
    if detector_boxes:
        inc = detector_boxes_cfg.get("include") if detector_boxes_cfg else None
        exc = detector_boxes_cfg.get("exclude") if detector_boxes_cfg else None
        lab = detector_boxes_cfg.get("labels") if detector_boxes_cfg else False
        mode = detector_boxes_cfg.get("mode") if detector_boxes_cfg else "full-bbox"
        draw_detector_boxes(ax_zx, ax_zy, ax_xy, detector_boxes, include=inc, exclude=exc,
                            labels=lab, zorder=4, mode=mode)

    # Mark vertical z lines if requested (DV entrance etc.)
    zline_handles = []
    if mark_z:
        for zline, label in mark_z:
            ax_zx.axvline(zline, linewidth=1, linestyle="--", color="#444444", zorder=6)
            ax_zy.axvline(zline, linewidth=1, linestyle="--", color="#444444", zorder=6)
            zline_handles.append(Line2D([0], [0], color="#444444", linestyle="--", linewidth=1, label=label))

    if label_mode == "panel":
        ax_panel = fig.add_subplot(gs[:,7])
        ax_panel.set_axis_off()
        handles = []
        if show_detectors or dv_info.get("mode"):
            handles.append(Patch(facecolor="#e6f2fa", edgecolor="#6baed6", alpha=0.25, label="DV polygon"))
        if detector_boxes:
            handles.append(Line2D([0], [0], color="black", linewidth=1.6, label="Detector boxes"))
        if beam_overlays:
            handles.append(Line2D([0], [0], color="black", linewidth=1.0, linestyle="-", alpha=0.7, label="Beamline overlays"))
        handles.extend(zline_handles)
        if handles:
            ax_panel.legend(handles=handles, loc="upper left", frameon=False, fontsize=9)
        if panel_info:
            ax_panel.text(0.02, 0.55, panel_info, ha="left", va="top", fontsize=9)
    elif zline_handles:
        ax_zx.legend(handles=zline_handles, loc="upper right", frameon=False, fontsize=9)

    fig.suptitle(title, fontsize=14)
    fig.savefig(outpath, dpi=250, bbox_inches='tight', pad_inches=0.4)
    plt.close(fig)

# ---------------- main ----------------
def main():
    ap = ArgumentParser(description=__doc__)

    ap.add_argument("--detector-boxes-csv", default=None,
                help="CSV with subsystem,zmin_cm,zmax_cm (+ optional xmin/xmax/ymin/ymax for full bbox).")

    ap.add_argument("--point-color", default=None, help="Scatter point color (e.g. 'tab:orange' or '#ff8800').")
    ap.add_argument("--point-alpha", type=float, default=0.15, help="Scatter alpha.")
    ap.add_argument("--point-size", type=float, default=4.0, help="Scatter marker size.")

    ap.add_argument("--dv-edgecolor", default="#6baed6", help="DV polygon edge color.")
    ap.add_argument("--dv-facecolor", default="#e6f2fa", help="DV polygon fill color.")
    ap.add_argument("--dv-alpha", type=float, default=0.25, help="DV polygon alpha.")

    ap.add_argument(
    "--detector-boxes-include",
    default=r"^(TargetArea|UpstreamTagger|DecayGas|SBT|Tracker\d*|SpectrometerMagnet|TimingDetector|ECal|HCal|MuonDetector|MuonShield)$"
    )
    ap.add_argument("--detector-boxes-exclude", default=None)
    ap.add_argument("--detector-boxes-labels", action="store_true", default=False)

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

    ap.add_argument(
        "--inputdir",
        required=True,
        nargs="+",
        help="One or more directories containing job_* (or a single job directory). "
            "Example: --inputdir /path/prodA /path/prodB"
    )
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
    ap.add_argument("--dv-z-first", type=float, default=None)
    ap.add_argument("--dv-z-last",  type=float, default=None)

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
    ap.add_argument("--hexbin-gridsize", type=int, default=100)
    ap.add_argument("--hexbin-mincnt", type=int, default=1)
    ap.add_argument("--hexbin-bins", choices=["log","linear"], default="log")
    ap.add_argument("--label-mode", choices=["none","panel","plot"], default="panel",
                    help="Label mode: none, panel (info panel), or plot (staggered text labels).")
    ap.add_argument("--mark-dv-zlines", action="store_true", default=False,
                    help="Mark DV entrance/exit z as dashed lines.")
    ap.add_argument("--auto-limits", action=BooleanOptionalAction, default=True,
                    help="Auto-compute axis limits from detector boxes when no explicit limits are set.")
    ap.add_argument("--write-summary", action="store_true", default=True)
    ap.add_argument("--dump-z-transitions", action="store_true", default=False,
                    help="Scan along z at (beam_x,beam_y) and print volume transitions (helps locate muon-shield end / gap).")

    # axis limits (same defaults as surviving_xyzplots)
    ap.add_argument("--zlim", nargs=2, type=float, default=None)
    ap.add_argument("--xlim", nargs=2, type=float, default=None)
    ap.add_argument("--ylim", nargs=2, type=float, default=None)
    ap.add_argument("--xylim-x", nargs=2, type=float, default=None)
    ap.add_argument("--xylim-y", nargs=2, type=float, default=None)

    args = ap.parse_args()
    random.seed(args.seed)

    in_dirs = [Path(p).resolve() for p in args.inputdir]
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

    def collect_job_dirs(base: Path, jobs_glob: str):
        # "single job directory" if it directly contains *_rec.root
        if list(base.glob("*_rec.root")):
            return [base]
        return sorted([p for p in base.glob(jobs_glob) if p.is_dir()])

    job_dirs = []
    seen = set()
    for base in in_dirs:
        if not base.exists():
            print(f"! skip inputdir (does not exist): {base}")
            continue
        for jd in collect_job_dirs(base, args.jobs):
            if jd not in seen:
                job_dirs.append(jd)
                seen.add(jd)

    if args.max_jobs is not None:
        job_dirs = job_dirs[:max(0, args.max_jobs)]

    if not job_dirs:
        raise RuntimeError(f"No job dirs found in any inputdir {in_dirs} with pattern {args.jobs}")

    # global containers (merged across jobs)
    pts = {"zx":[[],[]], "zy":[[],[]], "xy":[[],[]]}
    counters = Counter()
    region_counts = Counter()
    dv_info_global = {
        "prefix": args.dv_prefix, "z_first": None, "z_last": None, "mode": args.geometry,
        "edgecolor": args.dv_edgecolor, "facecolor": args.dv_facecolor, "alpha": args.dv_alpha
    }
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
        if args.dv_z_first is not None and args.dv_z_last is not None:
            z_first, z_last = args.dv_z_first, args.dv_z_last
        elif args.geometry == "fixed":
            z_first, z_last = -2500.0, 2500.0
        else:
            z_first, z_last = find_z_range_for_prefix(args.dv_prefix, args.zscan_min, args.zscan_max,
                                                    args.zscan_step, args.beam_x, args.beam_y)

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
            if not hasattr(tree, "MCTrack") or tree.MCTrack.GetEntries() < 1:
                continue

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

    label_mode = args.label_mode
    if ("--label-mode" not in sys.argv) and (args.detector_boxes_labels or args.overlay_labels):
        label_mode = "plot"

    detector_boxes = None
    detector_boxes_mode = None
    detector_boxes_path = args.detector_boxes_csv
    if detector_boxes_path:
        if not Path(detector_boxes_path).exists():
            raise RuntimeError(f"Detector boxes CSV not found: {detector_boxes_path}")
    else:
        auto_csv = find_detector_boxes_csv([Path.cwd(), *in_dirs, out_dir])
        if auto_csv:
            detector_boxes_path = str(auto_csv)
            print(f"Auto-loaded detector boxes from: {detector_boxes_path}")

    if detector_boxes_path:
        detector_boxes, detector_boxes_mode = load_subsystem_boxes_csv(detector_boxes_path)
        print(f"Detector boxes CSV mode: {detector_boxes_mode}")
        if detector_boxes_mode == "z-only":
            print("! Warning: z-only detector boxes CSV loaded; XY boxes are not available.")
        detector_boxes = filter_detector_boxes(detector_boxes, args.detector_boxes_include, args.detector_boxes_exclude)
        if not detector_boxes:
            print("! Warning: detector boxes CSV loaded but no boxes matched include/exclude filters.")
    else:
        print("! Warning: no detector boxes CSV found. Target/MuonShield boxes will not be drawn. "
              "Pass --detector-boxes-csv to enable.")

    detector_boxes_cfg = {
        "include": None,
        "exclude": None,
        "labels": (label_mode == "plot"),
        "mode": detector_boxes_mode or "full-bbox",
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

    explicit_limits = any(
        v is not None for v in (args.zlim, args.xlim, args.ylim, args.xylim_x, args.xylim_y)
    )
    auto_limits = args.auto_limits and not explicit_limits
    limits_from_boxes_cfg = limits_from_boxes(detector_boxes) if auto_limits else None
    if auto_limits and limits_from_boxes_cfg:
        zlim = list(limits_from_boxes_cfg["zlim"])
        if "xlim" in limits_from_boxes_cfg:
            xlim = list(limits_from_boxes_cfg["xlim"])
            ylim = list(limits_from_boxes_cfg["ylim"])
            xylim_x = list(limits_from_boxes_cfg["xlim"])
            xylim_y = list(limits_from_boxes_cfg["ylim"])
        else:
            xlim = args.xlim if args.xlim is not None else DEFAULT_LIMITS["xlim"]
            ylim = args.ylim if args.ylim is not None else DEFAULT_LIMITS["ylim"]
            xylim_x = args.xylim_x if args.xylim_x is not None else DEFAULT_LIMITS["xylim_x"]
            xylim_y = args.xylim_y if args.xylim_y is not None else DEFAULT_LIMITS["xylim_y"]
    else:
        if auto_limits and not detector_boxes:
            print("! Warning: auto-limits requested but no detector boxes available. Using default limits.")
        zlim = args.zlim if args.zlim is not None else DEFAULT_LIMITS["zlim"]
        xlim = args.xlim if args.xlim is not None else DEFAULT_LIMITS["xlim"]
        ylim = args.ylim if args.ylim is not None else DEFAULT_LIMITS["ylim"]
        xylim_x = args.xylim_x if args.xylim_x is not None else DEFAULT_LIMITS["xylim_x"]
        xylim_y = args.xylim_y if args.xylim_y is not None else DEFAULT_LIMITS["xylim_y"]

    axes_cfg = {
        "zlim": zlim,
        "xlim": xlim,
        "ylim": ylim,
        "xylim_x": xylim_x,
        "xylim_y": xylim_y,
    }

    mark_z = []
    if args.mark_dv_zlines:
        if dv_info_global.get("z_first") is not None:
            mark_z.append((dv_info_global["z_first"], "DV entrance"))
        if dv_info_global.get("z_last") is not None:
            mark_z.append((dv_info_global["z_last"], "DV exit"))

    overlay_labels = (label_mode == "plot")
    hexbin_cfg = {"gridsize": args.hexbin_gridsize, "mincnt": args.hexbin_mincnt, "bins": args.hexbin_bins}

    panel_lines = []
    if label_mode == "panel":
        panel_lines.append("Subsystem boxes:")
        if detector_boxes:
            for name in sorted(detector_boxes.keys()):
                panel_lines.append(f"- {name}")
        else:
            panel_lines.append("- (none)")
        panel_lines.append("")
        panel_lines.append(f"z-range: {axes_cfg['zlim'][0]:.0f} .. {axes_cfg['zlim'][1]:.0f} cm")
        if detector_boxes_mode == "z-only":
            panel_lines.append("XY boxes not available (z-only CSV)")
        if mark_z:
            for zline, label in mark_z:
                panel_lines.append(f"{label}: z={zline:.0f} cm")
    panel_info = "\n".join(panel_lines) if panel_lines else None

    overlays_state = [
        f"DV polygon={args.geometry}",
        f"detector boxes={'yes' if detector_boxes else 'no'}",
        f"beamline overlays={'on' if beam_overlays else 'off'}",
        f"mark dv z-lines={'on' if args.mark_dv_zlines else 'off'}",
        f"label mode={label_mode}",
        f"auto limits={'on' if auto_limits else 'off'}",
    ]
    print("Overlays: " + ", ".join(overlays_state))
    print(f"Axis limits: z={axes_cfg['zlim']}, x={axes_cfg['xlim']}, y={axes_cfg['ylim']}, "
          f"x/y={axes_cfg['xylim_x']},{axes_cfg['xylim_y']}")

    plot_points(
    pts, str(outpath), title, dv_info_global, axes_cfg,
    show_detectors=args.show_detectors,
    mark_z=mark_z,
    plot_mode=args.plot_mode,
    beam_overlays=beam_overlays,
    overlay_labels=overlay_labels,
    detector_boxes=detector_boxes,
    detector_boxes_cfg=detector_boxes_cfg,
    label_mode=label_mode,
    panel_info=panel_info,
    hexbin_cfg=hexbin_cfg,
    point_color=args.point_color,
    point_alpha=args.point_alpha,
    point_size=args.point_size,
    )






    print(f"\nSaved: {outpath}")
    print(f"Stored points (after filters): {len(pts['zx'][0])}")
    print(f"Density counts (before filters): {dict(counters)}")
    print(f"Region counts  (before filters): {dict(region_counts)}")
    #print(f"DV scan z-range: {dv_info_global.get('z_first')} .. {dv_info_global.get('z_last')}  (prefix={args.dv_prefix})")

    if args.write_summary:
        summary = {
            "inputdirs": [str(d) for d in in_dirs],
            "jobs_used": [str(p) for p in job_dirs],
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
