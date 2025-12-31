#!/usr/bin/env python3
import argparse
import glob
import math
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)

# --- helper: build Lorentz vector from fitted track state and mass hypothesis
def lv_from_fittrack(fittrack, mass):
    st = fittrack.getFittedState()
    mom = st.getMom()  # TVector3-like
    px, py, pz = mom.X(), mom.Y(), mom.Z()
    p2 = px*px + py*py + pz*pz
    E = math.sqrt(p2 + mass*mass)
    v = ROOT.TLorentzVector(px, py, pz, E)
    return v

def invmass_from_two_fittracks(ft1, ft2, m1, m2):
    v = lv_from_fittrack(ft1, m1) + lv_from_fittrack(ft2, m2)
    return v.M()

def get_tree(f):
    # Try common tree names
    for name in ("cbmsim", "tree", "events"):
        t = f.Get(name)
        if t and isinstance(t, ROOT.TTree):
            return t
    # fallback: first TTree in file
    for key in f.GetListOfKeys():
        obj = f.Get(key.GetName())
        if isinstance(obj, ROOT.TTree):
            return obj
    return None

def analyze_files(files, max_events, max_candidates_per_event, out_prefix):
    # masses in GeV
    m_e  = 0.00051099895
    m_mu = 0.1056583755
    m_pi = 0.13957039

    # Histograms
    h_cand = ROOT.TH1F("h_cand",  "Candidate GetMass(); m [GeV]; Entries", 200, 0.0, 5.0)
    h_pipi = ROOT.TH1F("h_pipi",  "m( #pi#pi ) from FitTracks; m [GeV]; Entries", 200, 0.0, 5.0)
    h_ee   = ROOT.TH1F("h_ee",    "m( ee ) from FitTracks; m [GeV]; Entries", 200, 0.0, 5.0)
    h_mumu = ROOT.TH1F("h_mumu",  "m( #mu#mu ) from FitTracks; m [GeV]; Entries", 200, 0.0, 5.0)

    n_events_total = 0
    n_cand_total   = 0

    min_cand = float("inf")
    min_pipi = float("inf")
    min_ee   = float("inf")
    min_mumu = float("inf")

    for path in files:
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            print(f"[WARN] Cannot open: {path}")
            continue
        t = get_tree(f)
        if not t:
            print(f"[WARN] No TTree found in: {path}")
            f.Close()
            continue

        n_entries = t.GetEntries()
        n_to_run = n_entries if max_events < 0 else min(n_entries, max_events)

        for i in range(n_to_run):
            t.GetEntry(i)
            n_events_total += 1

            # Expect: event.Particles contains ShipParticle candidates
            # Guard in case branch name differs
            if not hasattr(t, "Particles"):
                continue

            nP = len(t.Particles)
            if nP < 1:
                continue

            # analyze up to N candidates per event (usually 1)
            n_take = min(nP, max_candidates_per_event)
            for ic in range(n_take):
                cand = t.Particles[ic]
                n_cand_total += 1

                # candidate mass as stored in ShipParticle
                m_c = float(cand.GetMass())
                h_cand.Fill(m_c)
                min_cand = min(min_cand, m_c)

                # manual masses using fitted tracks + different mass hypotheses
                if not hasattr(t, "FitTracks"):
                    continue
                try:
                    d1 = cand.GetDaughter(0)
                    d2 = cand.GetDaughter(1)
                except Exception:
                    continue
                if d1 < 0 or d2 < 0:
                    continue
                if d1 >= len(t.FitTracks) or d2 >= len(t.FitTracks):
                    continue

                ft1 = t.FitTracks[d1]
                ft2 = t.FitTracks[d2]

                m_pipi = invmass_from_two_fittracks(ft1, ft2, m_pi, m_pi)
                m_ee   = invmass_from_two_fittracks(ft1, ft2, m_e,  m_e)
                m_mumu = invmass_from_two_fittracks(ft1, ft2, m_mu, m_mu)

                h_pipi.Fill(m_pipi); min_pipi = min(min_pipi, m_pipi)
                h_ee.Fill(m_ee);     min_ee   = min(min_ee, m_ee)
                h_mumu.Fill(m_mumu); min_mumu = min(min_mumu, m_mumu)

        f.Close()

    print("\n=== Summary ===")
    print(f"Files: {len(files)}")
    print(f"Events processed: {n_events_total}")
    print(f"Candidates processed: {n_cand_total}")
    print(f"Min candidate GetMass(): {min_cand:.6f} GeV")
    print(f"Min m(pipi) from tracks: {min_pipi:.6f} GeV  (2*m_pi ~ {2*0.13957039:.6f} GeV)")
    print(f"Min m(ee)   from tracks: {min_ee:.6f} GeV")
    print(f"Min m(mumu) from tracks: {min_mumu:.6f} GeV")

    # Save histograms
    out_root = out_prefix + ".root"
    fout = ROOT.TFile(out_root, "RECREATE")
    for h in (h_cand, h_pipi, h_ee, h_mumu):
        h.Write()
    fout.Close()
    print(f"[OK] Wrote histograms to: {out_root}")

    # Save quick PDF
    c = ROOT.TCanvas("c", "c", 900, 700)
    c.SetLogy(True)
    h_cand.Draw("HIST")
    h_pipi.SetLineColor(ROOT.kRed);  h_pipi.Draw("HIST SAME")
    h_ee.SetLineColor(ROOT.kBlue);   h_ee.Draw("HIST SAME")
    h_mumu.SetLineColor(ROOT.kGreen+2); h_mumu.Draw("HIST SAME")
    leg = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
    leg.AddEntry(h_cand, "candidate.GetMass()", "l")
    leg.AddEntry(h_pipi, "track-based #pi#pi", "l")
    leg.AddEntry(h_ee,   "track-based ee", "l")
    leg.AddEntry(h_mumu, "track-based #mu#mu", "l")
    leg.Draw()
    out_pdf = out_prefix + ".pdf"
    c.SaveAs(out_pdf)
    print(f"[OK] Wrote plot to: {out_pdf}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, nargs="+",
                    help="Input ROOT file(s) or glob pattern(s)")
    ap.add_argument("--max-events", type=int, default=200000,
                    help="Max events per file (-1 = all). Default: 200k")
    ap.add_argument("--max-cands", type=int, default=1,
                    help="Max candidates per event to inspect. Default: 1")
    ap.add_argument("-o", "--out", default="invmass_check",
                    help="Output prefix (without extension). Default: invmass_check")
    args = ap.parse_args()

    files = []
    for inp in args.input:
        matched = sorted(glob.glob(inp))
        if matched:
            files.extend(matched)
        elif os.path.isfile(inp):
            files.append(inp)

    files = sorted(set(files))

    if not files:
        print(f"[ERROR] No files matched: {args.input}")
        sys.exit(1)

    analyze_files(files, args.max_events, args.max_cands, args.out)

if __name__ == "__main__":
    main()
