#!/usr/bin/env python3
"""
Combine selection-parameter ROOT files from SHiP job_* directories into one
file with up to three sub-directories: /all, /heliumCase, /vesselCase.

Works with uproot 4.1.0   (no call to identify.to_TH1).

Options
-------
  --test   merge only the first 10 ROOT files for each flavour
  -j N     N parallel workers (default 3, set 1 for serial)
  -o FILE  custom output file name
"""

# ------------------------------------------------------------------ imports
import sys, time, multiprocessing as mp
from pathlib import Path
from argparse import ArgumentParser

import numpy as np
import uproot
import tqdm

# ------------------------------------------------------------------ CLI
cli = ArgumentParser(description=__doc__)
bg = cli.add_mutually_exclusive_group(required=True)
bg.add_argument("--muonDIS", action="store_true")
bg.add_argument("--neuDIS",  action="store_true")
bg.add_argument("--mupi",    action="store_true")
bg.add_argument("--mumuv",   action="store_true")

cli.add_argument("-j", "--jobs", type=int, default=3,
                 help="parallel workers (default 3)")
cli.add_argument("-o", "--outfile", default="combined_selectionparameters.root",
                 help="output ROOT file")
cli.add_argument("--test", action="store_true",
                 help="process only first 10 files per filetype")

ARGS = cli.parse_args()

# ------------------------------------------------------------------ background paths
if   ARGS.muonDIS:
    TAG, BASES = "muonDIS", [
        "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/muonDIS/SBT",
        "/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/muonDIS/Tr",
    ]
elif ARGS.neuDIS:
    TAG, BASES = "neuDIS", ["/eos/experiment/ship/user/anupamar/BackgroundStudies/corrected/neuDIS"]
elif ARGS.mupi:
    TAG, BASES = "mupi",   ["/eos/experiment/ship/user/anupamar/BackgroundStudies/mupi_EventCalc"]
elif ARGS.mumuv:
    TAG, BASES = "mumuv",  ["/eos/experiment/ship/user/anupamar/BackgroundStudies/2muv_EventCalc"]
else:
    sys.exit("No background type selected.")

if ARGS.outfile == "combined_selectionparameters.root":
    ARGS.outfile = f"combined_selectionparameters_{TAG}.root"

# three flavours are always attempted
FLAVOURS = [
    ("all",        "selectionparameters_all.root"),
    ("heliumCase", "selectionparameters_heliumCase.root"),
    ("vesselCase", "selectionparameters_vesselCase.root"),
]
FILE_LIMIT = 10 if ARGS.test else None

# ------------------------------------------------------------------ worker
def merge_one(pair):
    label, fname = pair
    tic = time.perf_counter()

    # gather job_* files
    files = [
        p / fname
        for base in BASES
        for p in Path(base).glob("job_*")
        if (p / fname).exists()
    ]
    if not files:
        return None
    if FILE_LIMIT and len(files) > FILE_LIMIT:
        files = files[:FILE_LIMIT]
        print(f"[{label}] TEST MODE – using {len(files)} files")

    merged = {}                       # hname -> [vals, vars, edges]

    for pf in tqdm.tqdm(files, desc=f"{label:10}", unit="file", leave=False):
        with uproot.open(pf) as f:
            for hname, h in f.items():
                if not hasattr(h, "to_numpy"):
                    continue
                vals, edges = h.to_numpy()      # 4.1 returns (vals, edges)
                vars_       = h.variances()     # Sumw2 (zeros if none)

                if hname not in merged:
                    merged[hname] = [vals.copy(), vars_.copy(), edges]
                else:
                    np.add(merged[hname][0], vals,  out=merged[hname][0])
                    np.add(merged[hname][1], vars_, out=merged[hname][1])

    return label, merged, len(files), time.perf_counter() - tic

# ------------------------------------------------------------------ run workers
t0 = time.perf_counter()
with mp.Pool(processes=min(ARGS.jobs, len(FLAVOURS))) as pool:
    results = [r for r in pool.map(merge_one, FLAVOURS) if r]

if not results:
    sys.exit("❌  No selection-parameter files found.")

# ------------------------------------------------------------------ write output
out = Path(ARGS.outfile).expanduser().resolve()
with uproot.recreate(out) as fout:
    for label, merged, n_files, dt in results:
        d = fout.mkdir(label)
        for hname, (vals, vars_, edges) in merged.items():
            # uproot 4.1: build TH1 directly from numpy (values, edges, variances)
            d[hname] = uproot.behaviors.TH1.from_numpy((vals, edges, vars_))
        print(f"[{label:10}] {n_files:5} files  → {len(merged):4} histos  (Δt {dt:.1f}s)")

print(f"\n✅  Created {out}")
print(f"Total wall-clock time: {time.perf_counter() - t0:.1f} s")
