#!/usr/bin/env python3
import ROOT, argparse, csv, re
from collections import Counter, defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--rec", required=True, help="Eingabe-ROOT (…_rec.root)")
parser.add_argument("--geo", required=True, help="Geometrie-ROOT (geofile_full.root)")
parser.add_argument("--max", type=int, default=100000, help="max Events (default: alle)")
parser.add_argument("--out", default="volumes_seen.csv", help="CSV-Output (rohe Liste)")
parser.add_argument("--explain-out", default="volumes_explained.csv",
                    help="CSV mit Material/Dichte/Medium/Pfad")
parser.add_argument("--groups-out", default="volumes_groups.csv",
                    help="CSV mit gruppierter Zusammenfassung")
parser.add_argument("--top", type=int, default=80,
                    help="Wieviele häufigste Volumes in der Konsole detailliert erklären")
parser.add_argument("--apply-filter", action="store_true",
                    help="Prefix-Filter wie in ip_category() anwenden (default: AUS)")
args = parser.parse_args()

# ----------------------------
# Geometrie laden (wichtig für gGeoManager/FindNode)
# ----------------------------
fgeo = ROOT.TFile.Open(args.geo)
sGeo = fgeo.Get("FAIRGeom")  # sorgt dafür, dass gGeoManager gefüllt ist
if not sGeo:
    raise RuntimeError("FAIRGeom nicht gefunden – Geometrie-Datei korrekt?")

# ----------------------------
# Eingabedatei
# ----------------------------
f = ROOT.TFile.Open(args.rec)
sTree = f.Get("cbmsim")
if not sTree:
    raise RuntimeError("TTree 'cbmsim' nicht gefunden")

# ----------------------------
# Einstellungen / Datenstrukturen
# ----------------------------
exclude_prefixes = ("LiSc","VetoInnerWall","VetoOuterWall",
                    "VetoVerticalRib","VetoLongitRib","DecayVacuum")

seen = Counter()                     # Volume -> Count
path_map = {}                        # Volume -> repräsentativer TGeo-Pfad
miss = 0

# ----------------------------
# Event-Loop
# ----------------------------
ip = ROOT.TVector3()
nTot = sTree.GetEntries() if args.max < 0 else min(args.max, sTree.GetEntries())
for i in range(nTot):
    sTree.GetEntry(i)
    if not hasattr(sTree,"MCTrack") or sTree.MCTrack.GetEntries() < 1:
        continue
    sTree.MCTrack[0].GetStartVertex(ip)
    try:
        node = ROOT.gGeoManager.FindNode(ip.X(), ip.Y(), ip.Z())
        vol = node.GetVolume().GetName() if node else ""
    except Exception:
        vol = ""
    if not vol:
        miss += 1
        continue

    # Optional: Filter aktivieren wie in ip_category()
    if args.apply_filter:
        if vol.startswith(exclude_prefixes):
            continue

    seen[vol] += 1

    # repräsentativen Pfad für dieses Volume einmalig sichern
    if vol not in path_map:
        nav = ROOT.gGeoManager.GetCurrentNavigator()
        # Fallback: kein Navigator → setze nur Volumen-Namen
        path_map[vol] = nav.GetPath() if nav else f"/?/{vol}"

# ----------------------------
# Ausgabe (Konsole)
# ----------------------------
print(f"Unterschiedliche Volumennamen (nach{' ' if args.apply_filter else ' ohne '}Filter): {len(seen)}")
for name, cnt in seen.most_common():
    print(f"{name:40s}  {cnt}")
print(f"Fälle ohne gültiges Volume (FindNode-Fehler/außerhalb): {miss}")

# ----------------------------
# CSV 1: rohe Liste
# ----------------------------
with open(args.out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["volume_name","count"])
    for name, cnt in seen.most_common():
        w.writerow([name, cnt])
print(f"→ CSV geschrieben: {args.out}")

# ----------------------------
# Zusatzinfos pro Volume: Material / Dichte / Medium / Pfad
# ----------------------------
def vol_material_tuple(vname: str):
    v = ROOT.gGeoManager.GetVolume(vname)
    mat_name, dens, medium = "", "", ""
    if v and v.GetMedium() and v.GetMedium().GetMaterial():
        mat = v.GetMedium().GetMaterial()
        mat_name = mat.GetName()
        dens = f"{mat.GetDensity():.6f}"  # g/cm3
        medium = v.GetMedium().GetName()
    return mat_name, dens, medium

with open(args.explain_out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["volume_name","count","material","density_g_cm3","medium","example_tgeo_path"])
    for name, cnt in seen.most_common():
        mat_name, dens, medium = vol_material_tuple(name)
        path = path_map.get(name, "")
        w.writerow([name, cnt, mat_name, dens, medium, path])
print(f"→ CSV geschrieben: {args.explain_out}")

# ----------------------------
# Gruppierte Übersicht per Präfix
# ----------------------------
groups = {
  "SBT_LiSc": re.compile(r"^LiSc"),
  "Vessel_walls": re.compile(r"^Veto(Inner|Outer)Wall"),
  "Vessel_ribs": re.compile(r"^(VetoVerticalRib|vLongitRib[XY])"),
  "UBT_glass": re.compile(r"^glass"),
  "UBT_optics": re.compile(r"^(pmma|FR4|Al)"),
  "Tracker_straw": re.compile(r"^(straw_|gas_)"),
  "Tracker_trplanes": re.compile(r"^Tr[12]"),
  "Decay_gas": re.compile(r"^DecayVacuum"),
  "Other": re.compile(r".*"),
}

agg = defaultdict(int)
for name, cnt in seen.items():
    for label, rx in groups.items():
        if rx.match(name):
            agg[label] += cnt
            break

with open(args.groups_out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["group","count"])
    for label in groups.keys():
        w.writerow([label, agg.get(label, 0)])
print(f"→ CSV geschrieben: {args.groups_out}")

# ----------------------------
# „Top N“ Detailausgabe in die Konsole
# ----------------------------
print(f"\nTop {args.top} Volumes – Material/Dichte/Medium + Beispielpfad:\n")
printed = 0
for name, cnt in seen.most_common():
    if printed >= args.top: break
    mat_name, dens, medium = vol_material_tuple(name)
    path = path_map.get(name, "")
    print(f"{name:40s}  hits={cnt:6d}  material={mat_name:20s}  density={dens:>10s}  medium={medium}")
    if path:
        print(f"  path: {path}")
    printed += 1

print("""
Tipp:
• Mit --apply-filter aktivierst du denselben Prefix-Filter wie ip_category().
• 'example_tgeo_path' zeigt dir die Eltern-Hierarchie (Subsystem) an der Stelle, wo das Volume erstmals gefunden wurde.
• Für tieferes Debugging kannst du einzelne Events mit FindNode(...) inspizieren und direkt nav.GetPath() ausgeben.
""")
