import ROOT
"""
file1=ROOT.TFile.Open("ip_plots_neuDIS_vesselCase_wdist2innerwall.root")
file2=ROOT.TFile.Open("ip_plots_muDIS_vesselCase_wdist2innerwall.root")
file3=ROOT.TFile.Open("compare_impactpar_mumuvu.root")

hist1 = file1.Get("vesselCase_wsel_dist_to_innerwall")
hist2 = file2.Get("vesselCase_wsel_dist_to_innerwall")
hist3 = file3.Get("dist_to_innerwall")

canvas = ROOT.TCanvas("c1", "", 1000, 700)
canvas.SetLogz()
canvas.SetMargin(0.13, 0.15, 0.12, 0.05)
ROOT.gStyle.SetLegendBorderSize(0)
s_leg  = ROOT.TLegend(0.18, 0.75, 0.55, 0.85)

s_leg.SetTextSize(0.035) 
s_leg.SetFillStyle(0)       # Fully transparent (no background box)

for i,hist in enumerate([hist1,hist2,hist3]):
	print(i)

	hist.SetStats(0)
	hist.SetTitle("")

	hist.GetXaxis().SetTitle("Distance to inner SBT wall [cm]")
	hist.GetYaxis().SetTitle("a.u")
	
	hist.Rebin(5)

	max_val = hist.GetMaximum()
	
	if max_val > 0:
	    hist.Scale(1.0/max_val)
	
	hist.GetXaxis().SetTitleSize(0.045)
	hist.GetYaxis().SetTitleSize(0.045)
	hist.GetXaxis().SetTitleOffset(1.2)
	hist.GetYaxis().SetTitleOffset(1.4)

	hist.SetLineWidth(3)
	
	if i==0:


		label='#nu-DIS in SBT'
		hist.SetLineColor(ROOT.kBlue+1)    # Blue
		hist.SetFillColor(ROOT.kBlue+1)    # Blue
		hist.SetFillStyle(3003) 
		hist.Draw("HIST")
	if i==1:
		label='#mu-DIS in SBT'
		hist.SetLineColor(ROOT.kGreen)    # Blue
		hist.SetFillColor(ROOT.kGreen)  
		#hist.SetFillStyle(3010) 
		hist.SetFillStyle(3003) 
		hist.Draw("HISTSAME")
	if i==2:
		label='#mu#mu#nu'

		hist.SetLineColor(ROOT.kRed+1)    # Blue
		hist.SetFillColor(ROOT.kRed+1)  
		#hist.SetFillStyle(3002) 
		hist.SetFillStyle(3003) 
		hist.Draw("HISTSAME")


	canvas.Update()
	s_leg.AddEntry(hist, label, "f")

s_leg.Draw()
canvas.Update()

canvas.SaveAs(f"compare_dist2innerwall.pdf")
canvas.SaveAs(f"compare_dist2innerwall.png")
canvas.Close()

exit()

file=ROOT.TFile.Open("compare_impactpar_mumuvu.root")
canvas = ROOT.TCanvas("c1", "", 1000, 700)
canvas.SetLogz()
canvas.SetMargin(0.13, 0.15, 0.12, 0.05)
hist = file.Get(f"impact_parameter_vs_zvtx")

hist.GetXaxis().SetTitle("Candidate vertex z [cm]")
hist.GetYaxis().SetTitle("Impact parameter [cm]")
hist.GetZaxis().SetTitle("a.u")


#hist.GetYaxis().SetRangeUser(0, 0.5)
#hist.GetYaxis().SetTitle("Energy Deposition in the SBT cell (GeV)")

hist.SetStats(0)
hist.SetTitle("")
hist.Rebin2D(400,8); 
hist.Draw("COLZ")
#hist.ProfileX().Draw("SAME"); 


hist.GetXaxis().SetTitleSize(0.045)
hist.GetYaxis().SetTitleSize(0.045)
hist.GetXaxis().SetTitleOffset(1.2)
hist.GetYaxis().SetTitleOffset(1.4)

label = ROOT.TLatex()
label.SetNDC()
label.SetTextSize(0.038)

label.DrawLatex(0.15, 0.88, f"#mu#mu#nu")

canvas.SaveAs(f"zvtx_mumunu.pdf")
canvas.SaveAs(f"zvtx_mumunu.png")
canvas.Close()


hist = file.Get(f"wsel_impact_parameter_vs_zvtx")

canvas = ROOT.TCanvas("c1", "", 1000, 700)
canvas.SetLogz()
canvas.SetMargin(0.13, 0.15, 0.12, 0.05)
hist.GetXaxis().SetTitle("Candidate vertex z [cm]")
hist.GetYaxis().SetTitle("Impact parameter [cm]")
hist.GetZaxis().SetTitle("a.u")


#hist.GetYaxis().SetRangeUser(0, 0.5)
#hist.GetYaxis().SetTitle("Energy Deposition in the SBT cell (GeV)")

hist.SetStats(0)
hist.SetTitle("")
hist.Rebin2D(400,8); 
hist.Draw("COLZ")
#hist.ProfileX().Draw("SAME"); 


hist.GetXaxis().SetTitleSize(0.045)
hist.GetYaxis().SetTitleSize(0.045)
hist.GetXaxis().SetTitleOffset(1.2)
hist.GetYaxis().SetTitleOffset(1.4)

label = ROOT.TLatex()
label.SetNDC()
label.SetTextSize(0.038)

label.DrawLatex(0.15, 0.88, f"#mu#mu#nu: wSel")

canvas.SaveAs(f"zvtx_wSel_mumunu.pdf")
canvas.SaveAs(f"zvtx_wSel_mumunu.png")
canvas.Close()


exit()
"""
filenames=["ip_plots_muDIS_heliumCase","ip_plots_muDIS_vesselCase"]#"ip_plots_neuDIS_heliumCase","ip_plots_neuDIS_vesselCase",

for filename in filenames:

	file = ROOT.TFile.Open(f"{filename}.root")

	canvas = ROOT.TCanvas("c1", "", 1000, 700)
	canvas.SetLogz()
	canvas.SetMargin(0.13, 0.15, 0.12, 0.05)

	hist = file.Get(f"{filename.split('_')[-1]}_impact_parameter_vs_z_interac_point")

	hist.GetXaxis().SetTitle("Interaction point z [cm]")
	hist.GetYaxis().SetTitle("Impact parameter [cm]")
	hist.GetZaxis().SetTitle("a.u")


	#hist.GetYaxis().SetRangeUser(0, 0.5)
	#hist.GetYaxis().SetTitle("Energy Deposition in the SBT cell (GeV)")

	hist.SetStats(0)
	hist.SetTitle("")
	hist.Rebin2D(400,8); 
	hist.Draw("COLZ")
	#hist.ProfileX().Draw("SAME"); 


	hist.GetXaxis().SetTitleSize(0.045)
	hist.GetYaxis().SetTitleSize(0.045)
	hist.GetXaxis().SetTitleOffset(1.2)
	hist.GetYaxis().SetTitleOffset(1.4)

	label = ROOT.TLatex()
	label.SetNDC()
	label.SetTextSize(0.038)
	if "neuDIS" in filename:
		tag='#nu-DIS'
	else:
		tag='#mu-DIS'
	label.DrawLatex(0.15, 0.88, f"{tag} : {filename.split('_')[-1]}")

	canvas.SaveAs(f"zIP_{filename}.pdf")
	canvas.SaveAs(f"zIP_{filename}.png")
	canvas.Close()

	canvas = ROOT.TCanvas("c2", "", 1000, 700)
	canvas.SetLogz()
	canvas.SetMargin(0.13, 0.15, 0.12, 0.05)

	hist = file.Get(f"{filename.split('_')[-1]}_impact_parameter_vs_zvtx")

	hist.GetXaxis().SetTitle("Candidate vertex z [cm]")
	hist.GetYaxis().SetTitle("Impact parameter [cm]")
	hist.GetZaxis().SetTitle("a.u")


	#hist.GetYaxis().SetRangeUser(0, 0.5)
	#hist.GetYaxis().SetTitle("Energy Deposition in the SBT cell (GeV)")

	hist.SetStats(0)
	hist.SetTitle("")
	hist.Rebin2D(400,8); 
	hist.Draw("COLZ")
	#hist.ProfileX().Draw("SAME"); 


	hist.GetXaxis().SetTitleSize(0.045)
	hist.GetYaxis().SetTitleSize(0.045)
	hist.GetXaxis().SetTitleOffset(1.2)
	hist.GetYaxis().SetTitleOffset(1.4)

	label = ROOT.TLatex()
	label.SetNDC()
	label.SetTextSize(0.038)
	if "neuDIS" in filename:
		tag='#nu-DIS'
	else:
		tag='#mu-DIS'
	label.DrawLatex(0.15, 0.88, f"{tag} : {filename.split('_')[-1]}")

	canvas.SaveAs(f"zvtx_{filename}.pdf")
	canvas.SaveAs(f"zvtx_{filename}.png")
	canvas.Close()

	#--------------

	canvas = ROOT.TCanvas("c1", "", 1000, 700)
	canvas.SetLogz()
	canvas.SetMargin(0.13, 0.15, 0.12, 0.05)

	hist = file.Get(f"{filename.split('_')[-1]}_wsel_impact_parameter_vs_z_interac_point")

	hist.GetXaxis().SetTitle("Interaction point z [cm]")
	hist.GetYaxis().SetTitle("Impact parameter [cm]")
	hist.GetZaxis().SetTitle("a.u")


	#hist.GetYaxis().SetRangeUser(0, 0.5)
	#hist.GetYaxis().SetTitle("Energy Deposition in the SBT cell (GeV)")

	hist.SetStats(0)
	hist.SetTitle("")
	hist.Rebin2D(400,8); 
	hist.Draw("COLZ")
	#hist.ProfileX().Draw("SAME"); 


	hist.GetXaxis().SetTitleSize(0.045)
	hist.GetYaxis().SetTitleSize(0.045)
	hist.GetXaxis().SetTitleOffset(1.2)
	hist.GetYaxis().SetTitleOffset(1.4)

	label = ROOT.TLatex()
	label.SetNDC()
	label.SetTextSize(0.038)
	if "neuDIS" in filename:
		tag='#nu-DIS'
	else:
		tag='#mu-DIS'
	label.DrawLatex(0.15, 0.88, f"{tag} : {filename.split('_')[-1]} wSel")

	canvas.SaveAs(f"zIP_wsel_{filename}.pdf")
	canvas.SaveAs(f"zIP_wsel_{filename}.png")
	canvas.Close()

	canvas = ROOT.TCanvas("c2", "", 1000, 700)
	canvas.SetLogz()
	canvas.SetMargin(0.13, 0.15, 0.12, 0.05)

	
	hist = file.Get(f"{filename.split('_')[-1]}_wsel_impact_parameter_vs_zvtx")

	hist.GetXaxis().SetTitle("Candidate vertex z [cm]")
	hist.GetYaxis().SetTitle("Impact parameter [cm]")
	hist.GetZaxis().SetTitle("a.u")


	#hist.GetYaxis().SetRangeUser(0, 0.5)
	#hist.GetYaxis().SetTitle("Energy Deposition in the SBT cell (GeV)")

	hist.SetStats(0)
	hist.SetTitle("")
	hist.Rebin2D(400,8); 
	hist.Draw("COLZ")
	#hist.ProfileX().Draw("SAME"); 


	hist.GetXaxis().SetTitleSize(0.045)
	hist.GetYaxis().SetTitleSize(0.045)
	hist.GetXaxis().SetTitleOffset(1.2)
	hist.GetYaxis().SetTitleOffset(1.4)

	label = ROOT.TLatex()
	label.SetNDC()
	label.SetTextSize(0.038)
	if "neuDIS" in filename:
		tag='#nu-DIS'
	else:
		tag='#mu-DIS'
	label.DrawLatex(0.15, 0.88, f"{tag} : {filename.split('_')[-1]} wSel")

	canvas.SaveAs(f"zvtx_wsel_{filename}.pdf")
	canvas.SaveAs(f"zvtx_wsel_{filename}.png")
	canvas.Close()