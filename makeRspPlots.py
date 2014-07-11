import ROOT
import sys
import array 

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(1111)

what = "3Jets"
cstr = "pt*rsp>10"
#what = "AllJets"
def makeResponsePlots(inputfile,outputfile,ptBins_in,iseta):
 
 nptbins = len(ptBins_in)-1
 coll = "Off"
 ext = 'PT'
 obj = 'P_{T}^{Gen}'
 unit = '(GeV)'
 
 nb      = 300
 min,max = 0,300

 if iseta: 
  ext = 'Eta'
  obj = '|#eta|^{L1}'
  unit = ''
  coll = "L1"
  nb = 100
  min,max =-2.5,2.5

 output_f =  outputfile.mkdir(ext,ext) 
 output_f_hists =  output_f.mkdir("Histograms") 
 
 lat = ROOT.TLatex()
 lat.SetNDC()
 lat.SetTextFont(42); lat.SetTextSize(0.04)

 leg = ROOT.TLegend(0.55,0.7,0.94,0.94)
 leg.SetFillColor(0); leg.SetTextFont(42); leg.SetTextSize(0.04)

 leg2 = ROOT.TLegend(0.6,0.75,0.94,0.94)
 leg2.SetFillColor(0); leg.SetTextFont(42); leg.SetTextSize(0.04)

 tree_raw    = inputfile.Get("Uncalib/valid")
# tree_lpu    = inputfile.Get("LPUS/valid")
 tree_cal    = inputfile.Get("Calib/valid")

 # From the input file, we make the histograms?
 
 if iseta:
   tree_cal.Draw("1./rsp:eta>>h2d_calib(%d,%g,%g,80,0,4)"%(nb,min,max),cstr)
   h2d_calib = ROOT.gROOT.FindObject("h2d_calib")
   #tree_lpu.Draw("1./rsp:eta>>h2d_lpu(%d,%g,%g,80,0,4)"%(nb,min,max),cstr)
   #h2d_lpu = ROOT.gROOT.FindObject("h2d_lpu")
   tree_raw.Draw("1./rsp:eta>>h2d_raw(%d,%g,%g,80,0,4)"%(nb,min,max),cstr)
   h2d_raw = ROOT.gROOT.FindObject("h2d_raw")
 else: 
   tree_cal.Draw("1./rsp:rsp*pt>>h2d_calib(%d,%g,%g,80,0,4)"%(nb,min,max),cstr)
   h2d_calib = ROOT.gROOT.FindObject("h2d_calib")
   #tree_lpu.Draw("1./rsp:rsp*pt>>h2d_lpu(%d,%g,%g,80,0,4)"%(nb,min,max),cstr)
   #h2d_lpu = ROOT.gROOT.FindObject("h2d_lpu")
   tree_raw.Draw("1./rsp:rsp*pt>>h2d_raw(%d,%g,%g,80,0,4)"%(nb,min,max),cstr)
   h2d_raw = ROOT.gROOT.FindObject("h2d_raw")

 # put them in a dict  index: label,fillcol,linecol,hist2D 
 schemas = {
	    #0: ["No PU + Calib", ROOT.kBlue-8,ROOT.kBlue+3, h2d_calib,0,0] 
	    #0: ["PU + calib", ROOT.kBlue-8,ROOT.kBlue+3, h2d_calib,0] 
	    0: ["Raw"	    	     , ROOT.kRed-8,ROOT.kRed+2, h2d_raw,0,0]
	  #  ,1: ["Local PU-sub"      , ROOT.kBlue-8,ROOT.kBlue+3, h2d_lpu,0,0]    
	    ,1: ["Local PU + calib"  , ROOT.kGreen-2,ROOT.kGreen+4, h2d_calib,0,0] 
	   }                  
 
 ptBins = []
 binindeces = []
 # first run through the bins and put the actual values of the bin edges there
 for i,ptR in enumerate(ptBins_in[0:-1]):
    bin1 = h2d_calib.GetXaxis().FindBin(ptR)
    bin2 = h2d_calib.GetXaxis().FindBin(ptBins_in[i+1])-1
    xlow = h2d_calib.GetXaxis().GetBinLowEdge(bin1)
    xup  = h2d_calib.GetXaxis().GetBinLowEdge(bin2+1)
    binindeces.append([bin1,bin2])
    ptBins.append(xlow)
 ptBins.append(xup)  # only need this last one

 for index in schemas.keys():

   name      = schemas[index][0]
   fc        = schemas[index][1]
   lc        = schemas[index][2]
   h2d_calib = schemas[index][3].Clone()
   
   hcalib   = ROOT.TH1F("h_%s"%name,";%s;<P_{T}^{L1}/P_{T}^{Gen}> (GeV)"%obj,nptbins,array.array('d',ptBins))
   
   for i,ptR in enumerate(ptBins[0:-1]):
 
    bin1 = binindeces[i][0]#h2d_calib.GetXaxis().FindBin(ptR)
    bin2 = binindeces[i][1]#h2d_calib.GetXaxis().FindBin(ptBins[i+1])-1
    #print "Binning mis-matches", ptR, ptBins[i+1], h2d_calib.GetXaxis().GetBinLowEdge(bin1),h2d_calib.GetXaxis().GetBinLowEdge(bin2+1)
   
   ########################### CALIBRATED #############################
    hc = h2d_calib.ProjectionY("prj_%s_%sBin%d"%(name,ext,i),bin1,bin2)
    c = ROOT.TCanvas('c_%s'%hc.GetName(),'c_%s'%hc.GetName(),600,600)
    hc.GetYaxis().SetTitleOffset(1.2)
    hc.GetXaxis().SetTitle('<P_{T}^{L1}/P_{T}^{Gen}>')
    hc.GetYaxis().SetTitle('Arbitrary Units')
    # fit for the mean 
    #hc.Fit('gaus','','',hc.GetMean()-1.2*hc.GetRMS(),hc.GetMean()+1.2*hc.GetRMS())
    hc.SetTitle('')
    hc.SetName('%s_%sbin_%d'%(name,ext,i))
    c.SetName('rsp_%s_%sbin_%d'%(name,ext,i))
    #mean = hc.GetFunction('gaus').GetParameter(1)
    #err  = hc.GetFunction('gaus').GetParError(1)
    #print mean,err
    mean = hc.GetMean()
    err  = hc.GetMeanError()
    hcalib.SetBinContent(i+1,mean)
    hcalib.SetBinError(i+1,err)

    hc.Draw('')
    lat.DrawLatex(0.1,0.92," %.2f < %s < %.2f %s "%(ptR,obj,ptBins[i+1],unit))
    output_f_hists.WriteTObject(c)

   hcalib.SetLineColor(lc)
   hcalib.SetLineWidth(2)
   hcalib.SetFillColor(fc)
   hcalib.SetMinimum(0.35)
   if iseta : hcalib.SetMaximum(1.8)
   else : hcalib.SetMaximum(2.0)
   hcalib.GetXaxis().SetTitleOffset(1.2)
   schemas[index][4]=hcalib.Clone()

   hallproj = h2d_calib.ProjectionX("prj_%s_%s_all"%(name,ext))
   hallproj.SetLineColor(fc)
   hallproj.SetLineWidth(3)
   hallproj.SetFillColor(0)
   schemas[index][5]=hallproj.Clone()

 c_vs_pt = ROOT.TCanvas("rsp_vs_%s"%ext,"rsp_vs_%s"%ext,800,600)
 for index in schemas.keys(): 
   if index==0: schemas[index][4].Draw("E2") 
   else: schemas[index][4].Draw("E2same")
   schemas[index][4].Draw('Lsame')
   leg.AddEntry(schemas[index][4],schemas[index][0],"LEF")

 leg.Draw()
  
 lone  = ROOT.TLine(hcalib.GetXaxis().GetXmin(),1,hcalib.GetXaxis().GetXmax(),1)
 lup   = ROOT.TLine(hcalib.GetXaxis().GetXmin(),1.1,hcalib.GetXaxis().GetXmax(),1.1)
 ldown = ROOT.TLine(hcalib.GetXaxis().GetXmin(),0.9,hcalib.GetXaxis().GetXmax(),0.9)
 lone.SetLineWidth(2)

 lone.SetLineColor(1)
 lup.SetLineStyle(2)
 ldown.SetLineStyle(2)
 lone.Draw()
 lup.Draw()
 ldown.Draw()
 if not iseta:
   #c_vs_pt.SetLogx()
   lCut=ROOT.TLine(30,0.35,30,2.0)
   aCut=ROOT.TArrow(30,1.5,70,1.5)
   lCut.SetLineColor(1)
   aCut.SetLineColor(1)
   lCut.SetLineWidth(2)
   aCut.SetLineWidth(2)
   lCut.Draw()
   aCut.Draw()

 # Now just plot the projection along X (-> Matched distribution)
 c_proj = ROOT.TCanvas("dist_%s"%ext,"dist_%s"%ext,800,600)
 for index in schemas.keys(): 
   if index==0: schemas[index][5].Draw("") 
   else: schemas[index][5].Draw("same")
   leg2.AddEntry(schemas[index][5],schemas[index][0],"L")
 leg2.Draw()
  

 output_f.WriteTObject(c_vs_pt)
 output_f.WriteTObject(c_proj)
 outname = (outputfile.GetName().split('.'))[0]
      
 #c_vs_pt.SaveAs("%s_%s.pdf"%(outname,c_vs_pt.GetName()))
 #c_vs_pt.SaveAs("%s_%s.C"%(outname,c_vs_pt.GetName()))
 #c_vs_pt.SaveAs("%s_%s.png"%(outname,c_vs_pt.GetName()))

 #c_proj.SaveAs("%s_%s.pdf"%(outname,c_proj.GetName()))
 #c_proj.SaveAs("%s_%s.C"%(outname,c_proj.GetName()))
 #c_proj.SaveAs("%s_%s.png"%(outname,c_proj.GetName()))

########### MAIN ########################

#inputf = ROOT.TFile('~/store/l1jec_upgrade/phase2calib/OutputJetsQcdHighPtE.root')
inputf = ROOT.TFile(sys.argv[1])
output_f = ROOT.TFile(sys.argv[2],"RECREATE")
#input_ttbar = ROOT.TFile('~/store/l1jec_upgrade/phase2calib/OutputJetsTtbarHighPtEta2p5.root')

ptBins  = [10,12,14,16,18,20,22.5,25,30,35,40,45,50,55,60,65,70,75,80,90,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,325,350]
etaBins = [-2.5,-2.2,-2.0,-1.75,-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.2,2.5]
# Vs PT

makeResponsePlots(inputf,output_f,ptBins,0) # 0 = pt, 1 = eta
makeResponsePlots(inputf,output_f,etaBins,1) # 0 = pt, 1 = eta

