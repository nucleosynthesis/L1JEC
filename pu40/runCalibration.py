import ROOT
import sys
import array 
import numpy 

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(1111)

#what = "AllJets"

# definition of the response function to fit 
fitfcn = ROOT.TF1("fitfcn","[0] + [1]/(TMath::Log(x)*TMath::Log(x) + [2] ) + [3]*TMath::Exp(-1*[4]*TMath::Log(x-[5])*TMath::Log(x-[5]))",20,250) 
fitfcn.SetParameter(0,0.5)
fitfcn.SetParameter(1,27)
fitfcn.SetParameter(2,8.78)
fitfcn.SetParameter(3,0.05)
fitfcn.SetParameter(4,0.0)
fitfcn.SetParameter(5,-11.)

def makeResponseCurves(inputfile,outputfile,ptBins_in, absetamin,absetamax):
 
 name = "ResponseProj"
 nptbins = len(ptBins_in)-1
 coll = "Off"
 ext = 'PT'
 obj = 'P_{T}^{Gen}'
 unit = '(GeV)'
 
 nb      = 300
 min,max = 0,300

 output_f =  outputfile.mkdir('eta_%g_%g'%(absetamin,absetamax)) 
 output_f_hists =  output_f.mkdir("Histograms") 
 
 lat = ROOT.TLatex()
 lat.SetNDC()
 lat.SetTextFont(42); lat.SetTextSize(0.04)

 leg = ROOT.TLegend(0.55,0.7,0.94,0.94)
 leg.SetFillColor(0); leg.SetTextFont(42); leg.SetTextSize(0.04)

 leg2 = ROOT.TLegend(0.6,0.75,0.94,0.94)
 leg2.SetFillColor(0); leg.SetTextFont(42); leg.SetTextSize(0.04)

 tree_raw    = inputfile.Get("Uncalib/valid")

 # From the input file, we make the histograms?
 
 cstr = " TMath::Abs(eta)<%g && TMath::Abs(eta) > %g "%(absetamax,absetamin)

 tree_raw.Draw("1./rsp:rsp*pt>>h2d_raw(%d,%g,%g,200,0,6)"%(nb,min,max),cstr)
 #tree_raw.Draw("pt:rsp*pt>>h2d_raw_l1(%d,%g,%g,200,0,200)"%(nb,min,max),cstr)
 h2d_raw = ROOT.gROOT.FindObject("h2d_raw")

 
 ptBins = []
 binindeces = []
 # first run through the bins and put the actual values of the bin edges there
 for i,ptR in enumerate(ptBins_in[0:-1]):
    bin1 = h2d_raw.GetXaxis().FindBin(ptR)
    bin2 = h2d_raw.GetXaxis().FindBin(ptBins_in[i+1])-1
    xlow = h2d_raw.GetXaxis().GetBinLowEdge(bin1)
    xup  = h2d_raw.GetXaxis().GetBinLowEdge(bin2+1)
    binindeces.append([bin1,bin2])
    ptBins.append(xlow)
 ptBins.append(xup)  # only need this last one

 gr = ROOT.TGraphErrors();
 grc=0
 for i,ptR in enumerate(ptBins[0:-1]):
 
    bin1 = binindeces[i][0]#h2d_calib.GetXaxis().FindBin(ptR)
    bin2 = binindeces[i][1]#h2d_calib.GetXaxis().FindBin(ptBins[i+1])-1
    #print "Binning mis-matches", ptR, ptBins[i+1], h2d_calib.GetXaxis().GetBinLowEdge(bin1),h2d_calib.GetXaxis().GetBinLowEdge(bin2+1)
   
   ########################### CALIBRATED #############################
    hc = h2d_raw.ProjectionY("prj_%s_%sBin%d"%(name,ext,i),bin1,bin2)
    xlow  = ptR
    xhigh = ptBins[i+1]
    tree_raw.Draw("pt>>hpt",cstr+" && pt*rsp < %g && pt*rsp > %g "%(xhigh,xlow))
    hpt = ROOT.gROOT.FindObject("hpt")

    hpt.SetName("L1_pt_genpt_%g_%g"%(xlow,xhigh))
    hc.SetName("Rsp_genpt_%g_%g"%(xlow,xhigh))
    output_f_hists.WriteTObject(hpt)
    if hc.GetEntries()>0: hc.Fit("gaus","Q","R",hc.GetMean()-1.*hc.GetRMS(),hc.GetMean()+1.*hc.GetRMS())
    output_f_hists.WriteTObject(hc)

    if not hpt.GetEntries()>0: continue
    if not hc.GetEntries()>0: continue

    mean = hc.GetFunction("gaus").GetParameter(1)
    err  = hc.GetFunction("gaus").GetParError(1)
    if not err > 0 : continue
    gr.SetPoint(grc,hpt.GetMean(),1./mean)
    gr.SetPointError(grc,hpt.GetMeanError(),err)
    grc+=1

 thisfit = fitfcn.Clone()
 thisfit.SetName(fitfcn.GetName()+'eta_%g_%g'%(absetamin,absetamax))
 gr.Fit(thisfit.GetName(),"","R+",20,250)
 gr.SetName('l1corr_eta_%g_%g'%(absetamin,absetamax))
 gr.GetXaxis().SetTitle("<p_{T}^{L1}>")
 gr.GetYaxis().SetTitle("1/<p_{T}^{Gen}/p_{T}^{L1}>")

 outputfile.WriteTObject(gr)
 outputfile.WriteTObject(thisfit)

########### MAIN ########################

#inputf = ROOT.TFile('~/store/l1jec_upgrade/phase2calib/OutputJetsQcdHighPtE.root')
inputf = ROOT.TFile(sys.argv[1])
output_f = ROOT.TFile(sys.argv[2],"RECREATE")
#input_ttbar = ROOT.TFile('~/store/l1jec_upgrade/phase2calib/OutputJetsTtbarHighPtEta2p5.root')

#ptBins  = [10,18,24,32,40,48,56,70,80,90,100,120,140,160,180,200,220,240,260,280,300,325,350]
ptBins_1 = numpy.arange(5,105,5)
ptBins_2 = numpy.arange(100,300,10)
ptBins = list(ptBins_1)[:]
ptBins+=list(ptBins_2)
etaBins = [ 0.0,0.348, 0.695, 1.044, 1.392, 1.74, 2.172, 3.0]
#etaBins = [ 1.392,1.74]
# Vs PT
for i,eta in enumerate(etaBins[0:-1]):
  emin = eta
  emax = etaBins[i+1]
  makeResponseCurves(inputf,output_f,ptBins,emin,emax) # 0 = pt, 1 = eta

