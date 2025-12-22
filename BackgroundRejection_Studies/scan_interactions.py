#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import json
from collections import Counter, defaultdict
from typing import Optional

import ROOT
ROOT.gROOT.SetBatch(True)

# -----------------------------
# Helpers
# -----------------------------

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def pick_first(patterns, base_dir):
    """Return first file matching any of the patterns in base_dir, else None."""
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(base_dir, pat)))
        if hits:
            return hits[0]
    return None

def load_geom(geofile: str) -> ROOT.TFile:
    fg = ROOT.TFile.Open(geofile, "READ")
    if not fg or fg.IsZombie():
        raise RuntimeError(f"Cannot open geofile: {geofile}")

    # This usually initialises gGeoManager in FairShip geofiles.
    _ = fg.Get("FAIRGeom")
    if not ROOT.gGeoManager:
        raise RuntimeError("gGeoManager was not initialised after loading FAIRGeom.")
    return fg

def find_z_range_for_volume_prefix(prefix: str, zmin: float, zmax: float, zstep: float, x0: float, y0: float):
    """
    Scan along z at fixed (x0,y0) and find first/last z where volume name startswith(prefix).
    Returns (z_first, z_last) or (None, None) if not found.
    """
    z_first, z_last = None, None
    z = zmin
    while z <= zmax:
        node = ROOT.gGeoManager.FindNode(x0, y0, z)
        if node and not ROOT.gGeoManager.IsOutside():
            volname = node.GetVolume().GetName()
            if volname.startswith(prefix):
                if z_first is None:
                    z_first = z
                z_last = z
        z += zstep
    return z_first, z_last

def density_bin(rho: float, he_thr: float, air_thr: float):
    """
    Simple density binning:
      rho < he_thr            -> He-like (very low)
      he_thr <= rho < air_thr -> air-like (low)
      rho >= air_thr          -> solids/other (high)
    """
    if rho != rho:  # NaN
        return "nan"
    if rho < he_thr:
        return f"very_low(<{he_thr:g}) [He-like]"
    if rho < air_thr:
        return f"low({he_thr:g}..{air_thr:g}) [air-like]"
    return f"high(>={air_thr:g}) [solids/other]"

# -----------------------------
# Core scan
# -----------------------------

def scan_one_job(job_dir: str,
                 max_events: Optional[int],
                 dv_prefix: str,
                 zmin: float,
                 zmax: float,
                 zstep: float,
                 beam_x: float,
                 beam_y: float,
                 he_thr: float,
                 air_thr: float,
                 air_path_samples: int,
                 dump_path_region: str):

    """
    Scan truth interaction vertices (MCTrack[0] start vertex) for a given job directory.
    Writes no files; returns a dict summary + some bulky counters.
    """

    # Rec file / geofile discovery (robust, production-independent)
    rec = pick_first(
        patterns=[
            "*_rec.root",
            "ship.*_rec.root",
            "ship.conical*_rec.root",
        ],
        base_dir=job_dir,
    )
    geofile = pick_first(
        patterns=[
            "geofile*.root",
            "*geofile*.root",
            "geofile_full*.root",
        ],
        base_dir=job_dir,
    )

    if not rec:
        raise RuntimeError(f"[{job_dir}] No *_rec.root found.")
    if not geofile:
        raise RuntimeError(f"[{job_dir}] No geofile*.root found.")

    frec = ROOT.TFile.Open(rec, "READ")
    if not frec or frec.IsZombie():
        raise RuntimeError(f"[{job_dir}] Cannot open rec file: {rec}")
    tree = frec.Get("cbmsim")
    if not tree:
        raise RuntimeError(f"[{job_dir}] Tree 'cbmsim' not found in {rec}")

    fgeo = load_geom(geofile)

    # Determine DV entrance/exit by scanning along the beam axis (x=beam_x, y=beam_y)
    dv_z_first, dv_z_last = find_z_range_for_volume_prefix(
        prefix=dv_prefix, zmin=zmin, zmax=zmax, zstep=zstep, x0=beam_x, y0=beam_y
    )

    n_entries = int(tree.GetEntries())
    n_scan = n_entries if max_events is None else min(n_entries, int(max_events))

    vol_counts = Counter()
    mat_counts = Counter()
    dens_counts = Counter()

    # Region bookkeeping (to implement your checklist logic)
    # region by z relative to dv_z_first + whether volume is DV
    region_counts = Counter()
    airlike_region_counts = Counter()
    dv_density_counts = Counter()  # densities only for events inside DV volumes
    dv_vol_counts = Counter()      # which DV volumes are actually hit

    # Save a few representative "air-like" events with paths
    air_paths = []
    air_paths_collected = 0

    vtx = ROOT.TVector3()

    for i in range(n_scan):
        tree.GetEntry(i)

        # Truth vertex of the neutrino interaction:
        # In your analysis code this is exactly how heliumCase/vesselCase gets assigned.
        tree.MCTrack[0].GetStartVertex(vtx)
        x, y, z = float(vtx.X()), float(vtx.Y()), float(vtx.Z())

        node = ROOT.gGeoManager.FindNode(x, y, z)
        if not node or ROOT.gGeoManager.IsOutside():
            vol = "OUTSIDE"
            mat = "OUTSIDE"
            rho = float("nan")
        else:
            vol = node.GetVolume().GetName()
            material = node.GetVolume().GetMaterial()
            mat = material.GetName() if material else "NO_MATERIAL"
            rho = float(material.GetDensity()) if material else float("nan")

        vol_counts[vol] += 1
        mat_counts[mat] += 1

        dbin = density_bin(rho, he_thr=he_thr, air_thr=air_thr)
        dens_counts[dbin] += 1

        # Determine region
        in_dv_by_name = vol.startswith(dv_prefix)
        if dv_z_first is None:
            # DV not found by scan -> region only by name / fallback
            region = "DV" if in_dv_by_name else "NonDV"
        else:
            if in_dv_by_name:
                region = "DV"
            elif z < dv_z_first:
                region = "UpstreamOfDV"
            else:
                region = "DownstreamOrSide"

        region_counts[region] += 1

        # Track air-like specifically, split by region
        if "[air-like]" in dbin:
            airlike_region_counts[region] += 1

        # DV-specific checks (your “He is really He?” sanity)
        if in_dv_by_name:
            dv_vol_counts[vol] += 1
            dv_density_counts[dbin] += 1

        # Collect a few air-like paths to identify if it is the SND-hole volume path
        if air_paths_collected < air_path_samples and "[air-like]" in dbin:
            # optionally restrict to a region to make it more diagnostic
            if dump_path_region == "any" or dump_path_region == region:
                path = ROOT.gGeoManager.GetPath()
                air_paths.append({
                    "event": i,
                    "x": x, "y": y, "z": z,
                    "vol": vol,
                    "mat": mat,
                    "rho": rho,
                    "region": region,
                    "path": path
                })
                air_paths_collected += 1

    frec.Close()
    fgeo.Close()

    # Build summary numbers (your checklist as numbers)
    summary = {
        "job_dir": job_dir,
        "rec_file": rec,
        "geofile": geofile,
        "entries_total": n_entries,
        "entries_scanned": n_scan,
        "dv_prefix": dv_prefix,
        "dv_z_first": dv_z_first,
        "dv_z_last": dv_z_last,
        "top_volumes": vol_counts.most_common(20),
        "density_bins": dens_counts,
        "region_counts": region_counts,
        "airlike_region_counts": airlike_region_counts,
        "dv_hit_fraction": (sum(dv_vol_counts.values()) / n_scan) if n_scan > 0 else 0.0,
        "dv_density_bins": dv_density_counts,
        "dv_top_volumes": dv_vol_counts.most_common(10),
        "airlike_path_samples": air_paths,
    }

    return summary, vol_counts, dens_counts

# -----------------------------
# Main
# -----------------------------

def main():
    import argparse

    ap = argparse.ArgumentParser(
        description="Scan truth interaction vertices (MCTrack[0]) and classify volumes/material densities. "
                    "Designed to diagnose air vs helium effects and upstream 'holes'."
    )
    ap.add_argument("--inputdir", required=True, help="Production directory containing job_* folders (or a single job folder).")
    ap.add_argument("--outputdir", required=True, help="Where to write JSON summaries and text outputs.")

    ap.add_argument("--jobs", default="job_*", help="Glob for job folders inside inputdir (default: job_*)")
    ap.add_argument("--max-jobs", type=int, default=None, help="Limit number of jobs scanned (default: no limit).")
    ap.add_argument("--max-events", type=int, default=None, help="Max events per job (default: all).")

    ap.add_argument("--dv-prefix", default="DecayVacuum", help="Volume name prefix for the Helium decay volume (default: DecayVacuum)")
    ap.add_argument("--zmin", type=float, default=-50000.0, help="z scan start (cm) for DV entrance detection.")
    ap.add_argument("--zmax", type=float, default=150000.0, help="z scan end (cm) for DV entrance detection.")
    ap.add_argument("--zstep", type=float, default=50.0, help="z scan step (cm) for DV entrance detection.")

    ap.add_argument("--beam-x", type=float, default=0.0, help="x (cm) for beamline DV scan (default: 0).")
    ap.add_argument("--beam-y", type=float, default=0.0, help="y (cm) for beamline DV scan (default: 0).")

    ap.add_argument("--he-thr", type=float, default=3e-4, help="Density threshold (g/cm^3) below which we call it He-like.")
    ap.add_argument("--air-thr", type=float, default=3e-3, help="Density threshold (g/cm^3) above which we call it solids/other.")

    ap.add_argument("--air-path-samples", type=int, default=25, help="How many air-like event paths to dump per job.")
    ap.add_argument("--dump-path-region", default="UpstreamOfDV",
                    choices=["any", "DV", "UpstreamOfDV", "DownstreamOrSide", "NonDV"],
                    help="Restrict which region to sample for gGeoManager.GetPath() dumps (default: UpstreamOfDV).")

    args = ap.parse_args()

    in_dir = os.path.abspath(args.inputdir)
    out_dir = os.path.abspath(args.outputdir)
    ensure_dir(out_dir)

    # Determine job dirs
    # If inputdir itself looks like a job folder (contains *_rec.root), scan it directly.
    if glob.glob(os.path.join(in_dir, "*_rec.root")):
        job_dirs = [in_dir]
    else:
        job_dirs = sorted(glob.glob(os.path.join(in_dir, args.jobs)))
        job_dirs = [d for d in job_dirs if os.path.isdir(d)]

    if args.max_jobs is not None:
        job_dirs = job_dirs[:max(0, int(args.max_jobs))]

    if not job_dirs:
        raise RuntimeError(f"No job directories found in {in_dir} with pattern {args.jobs}")

    combined_vol = Counter()
    combined_dens = Counter()

    all_summaries = []

    for jd in job_dirs:
        print(f"\n=== Scanning job: {jd}")
        summary, vol_counts, dens_counts = scan_one_job(
            job_dir=jd,
            max_events=args.max_events,
            dv_prefix=args.dv_prefix,
            zmin=args.zmin,
            zmax=args.zmax,
            zstep=args.zstep,
            beam_x=args.beam_x,
            beam_y=args.beam_y,
            he_thr=args.he_thr,
            air_thr=args.air_thr,
            air_path_samples=args.air_path_samples,
            dump_path_region=args.dump_path_region,
        )

        all_summaries.append(summary)
        combined_vol.update(vol_counts)
        combined_dens.update(dens_counts)

        # Write per-job JSON
        job_name = os.path.basename(os.path.normpath(jd))
        job_out = os.path.join(out_dir, f"scan_{job_name}.json")
        with open(job_out, "w") as f:
            json.dump(summary, f, indent=2, sort_keys=True)
        print(f"  -> wrote {job_out}")

        # Human-readable quick print (your checklist)
        print(f"  entries_scanned: {summary['entries_scanned']}")
        print(f"  DV z-range (by scan @ beamline): {summary['dv_z_first']} .. {summary['dv_z_last']}")
        print(f"  DV-hit fraction (by name startswith '{args.dv_prefix}'): {summary['dv_hit_fraction']:.4f}")

        print("  Density bins:")
        for k, v in summary["density_bins"].items():
            print(f"    {v:8d}  {k}")

        print("  Regions:")
        for k, v in summary["region_counts"].items():
            print(f"    {v:8d}  {k}")

        print("  Air-like split by region:")
        for k, v in summary["airlike_region_counts"].items():
            print(f"    {v:8d}  {k}")

        # Write sample paths as a text file too (easy grep)
        path_out = os.path.join(out_dir, f"air_paths_{job_name}.txt")
        with open(path_out, "w") as f:
            for e in summary["airlike_path_samples"]:
                f.write(
                    f"event={e['event']} region={e['region']} "
                    f"x={e['x']:.2f} y={e['y']:.2f} z={e['z']:.2f} "
                    f"rho={e['rho']:.4g} vol={e['vol']} mat={e['mat']}\n"
                    f"path={e['path']}\n\n"
                )
        print(f"  -> wrote {path_out}")

        # DV density sanity: is DV really He-like?
        print("  DV density bins (events whose volume name startswith DV prefix):")
        for k, v in summary["dv_density_bins"].items():
            print(f"    {v:8d}  {k}")

        print("  DV top volumes actually hit:")
        for k, v in summary["dv_top_volumes"]:
            print(f"    {v:8d}  {k}")

    # Combined summary
    combined = {
        "n_jobs": len(job_dirs),
        "combined_top_volumes": combined_vol.most_common(30),
        "combined_density_bins": combined_dens,
        "jobs": [os.path.basename(os.path.normpath(j)) for j in job_dirs],
    }

    combined_out = os.path.join(out_dir, "scan_combined.json")
    with open(combined_out, "w") as f:
        json.dump(combined, f, indent=2, sort_keys=True)
    print(f"\n=== Combined summary written to {combined_out}")

if __name__ == "__main__":
    main()
