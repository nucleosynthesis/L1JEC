import ROOT
import sys

print "qcd from ", sys.argv[1]
print "ttbar from ", sys.argv[2]

ROOT.gROOT.SetBatch(1)
treen = "Uncalib/valid"

fout = ROOT.TFile("comparevars.root","RECREATE")

fqcd = ROOT.TFile.Open(sys.argv[1])
fttb = ROOT.TFile.Open(sys.argv[2])

tqcd = fqcd.Get(treen)
tttb = fttb.Get(treen)

vars = ["sopt","(er4+hr4)/pt","(er1+hr1)/pt","(er2+hr2)/pt","(er3+hr3)/pt","(er1+hr2+er2+hr2)/(er3+hr3+er4+hr4)"]

tttb.Draw("pt>>httb")
tpt = ROOT.gROOT.FindObject("httb")

tqcd.Draw("pt>>hqcd")
qpt = ROOT.gROOT.FindObject("hqcd")

class weight_func:
  def __init__(self,hist):
    self.hist = hist.Clone()
  def __call__(self,x):
    numerate = float(self.hist.Integral()/self.hist.GetNbinsX())
    denom    = self.hist.GetBinContent(self.hist.FindBin(x))
    if denom>0 :return numerate/denom
    else : return 1;


qcd_weight_f = weight_func( qpt )
ttb_weight_f = weight_func( tpt )

def emptyHist(h):
 for b in range(h.GetNbinsX()): 
 	h.SetBinContent(b+1,0)
 	h.SetBinError(b+1,0)

for v in vars: 

  c = ROOT.TCanvas(v,v,600,600)

  tqcd.Draw("%s>>hqcd"%v,"")
  hqcd = ROOT.gROOT.FindObject("hqcd")
  emptyHist(hqcd)
  for ev in range(tqcd.GetEntries()):
    tqcd.GetEntry(ev)
    hqcd.Fill(getattr(tqcd,v),qcd_weight_f(getattr(tqcd,"pt")))
  hqcd.SetLineColor(1)
  hqcd.SetMarkerColor(1)

  tttb.Draw("%s>>httb"%v,"")
  httb = ROOT.gROOT.FindObject("httb")
  emptyHist(httb)
  for ev in range(tttb.GetEntries()):
    tttb.GetEntry(ev)
    httb.Fill(getattr(tttb,v),ttb_weight_f(getattr(tttb,"pt")))
  httb.SetLineColor(2)
  httb.SetMarkerColor(2)

  hqcd.Scale(1./hqcd.Integral())
  httb.Scale(1./httb.Integral())

  minX = min([hqcd.GetXaxis().GetXmin(),httb.GetXaxis().GetXmin()])
  maxX = max([hqcd.GetXaxis().GetXmax(),httb.GetXaxis().GetXmax()])

  hqcd.GetXaxis().SetRangeUser(minX,maxX)

  hqcd.SetTitle("")
  hqcd.GetXaxis().SetTitle(v)

  hqcd.Draw("hist")
  httb.Draw("histsame")
  hqcd.Draw("sameP")
  httb.Draw("sameP")

  fout.WriteTObject(c)

print ".. Histos in  ", fout.GetName()
fout.Close()
