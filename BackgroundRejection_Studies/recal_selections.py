import ROOT
df = ROOT.RDataFrame("cands", "selection_eventwise.root")

expr = ("pass_preselection && pass_UBT && pass_AdvSBT45 "
        "&& pass_PID && pass_inv_mass")

# sum of expected events in 15y:
n15y = df.Filter(expr).Sum("w15y").GetValue()
print("nEvents15y =", n15y)

# by category (helium only)
n15y_he = df.Filter("(cat==2) && ("+expr+")").Sum("w15y").GetValue()
print("heliumCase 15y =", n15y_he)
