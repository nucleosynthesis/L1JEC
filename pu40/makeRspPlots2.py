import ROOT
import sys
import array 
import numpy

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(1111)

what = "3Jets"
cstr = " jrank<=1 && l1dr > 1 "
genPtCut = 20  # cut applied to other plots, NOT rsp vs pt of course
#what = "AllJets"

def makeResponsePlots(inputfile,outputfile,ptBins_in,what,ext,obj,unit,nb,min,max,cutGenPt):

 schemas = {
	    0:  ["Raw"	    	     	     , ROOT.kRed-8,ROOT.kRed+2, 0,0,0,0,"Uncalib/valid"]
	    ,1: ["Calib old style"  	     , ROOT.kGreen-2,ROOT.kGreen+4, 0,0,0,0,"Calib_param/valid"] 
	    ,2: ["Calib Regression 2 var"    , ROOT.kBlue-8,ROOT.kBlue+3, 0,0,0,0,"Calib_2var/valid"] 
	    ,3: ["Calib Regression 3 var H/E"    , ROOT.kPink-3,ROOT.kPink+1, 0,0,0,0,"Calib_3var_hoe/valid"] 
	    ,4: ["Calib Regression 3 var Seed/pT"    , ROOT.kOrange-3,ROOT.kRed, 0,0,0,0,"Calib_3var_sopt/valid"] 
	    ,5: ["Calib Regression 4 var"       , ROOT.kMagenta-3,ROOT.kMagenta, 0,0,0,0,"Calib_4var/valid"] 
 }                 
 
 nptbins = len(ptBins_in)-1

 output_f =  outputfile.mkdir(ext,ext) 
 output_f_hists =  output_f.mkdir("Histograms") 
 
 lat = ROOT.TLatex()
 lat.SetNDC()
 lat.SetTextFont(42); lat.SetTextSize(0.04)

 leg = ROOT.TLegend(0.55,0.7,0.94,0.94)
 leg.SetFillColor(0); leg.SetTextFont(42); leg.SetTextSize(0.04)

 leg2 = ROOT.TLegend(0.6,0.75,0.94,0.94)
 leg2.SetFillColor(0); leg.SetTextFont(42); leg.SetTextSize(0.04)

 leg3 = ROOT.TLegend(0.6,0.75,0.94,0.94)
 leg3.SetFillColor(0); leg.SetTextFont(42); leg.SetTextSize(0.04)


 # From the input file, we make the histograms?
 

 # put them in a dict  index: label,fillcol,linecol,hist2D 

 tdirkeys = [t.GetName() for t in list(inputfile.GetListOfKeys())]
 for index in schemas.keys():
   
  if not schemas[index][-1].split("/")[0] in tdirkeys : 
  	print "Directory -- %s not found, skipping"%schemas[index][-1]
  	schemas.pop(index)

 for index in schemas.keys():
  tree_raw    = inputfile.Get(schemas[index][-1])
  
  if cutGenPt:
   tree_raw.Draw("1./rsp:%s>>h2d_raw_%s(%d,%g,%g,120,0,6)"%(what,schemas[index][0]+what,nb,min,max),cstr+ " && matched_gpt > %g "%(genPtCut))
  else: 
   tree_raw.Draw("1./rsp:%s>>h2d_raw_%s(%d,%g,%g,120,0,6)"%(what,schemas[index][0]+what,nb,min,max),cstr)
  h2d_raw = ROOT.gROOT.FindObject("h2d_raw_%s"%(schemas[index][0]+what))
  schemas[index][3]=h2d_raw.Clone()
 
 ptBins = []
 binindeces = []

 # first run through the bins and put the actual values of the bin edges there
 h2d_calib = schemas[schemas.keys()[0]][3].Clone()
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
    if hc.GetEntries() > 0: 
    #  hc.Fit("gaus","Q","R",hc.GetMean()-hc.GetRMS(),hc.GetMean()+1*hc.GetRMS())
    #  mean = hc.GetFunction('gaus').GetParameter(1)
    #  err  = hc.GetFunction('gaus').GetParError(1)i

    #print mean,err
     mean = hc.GetMean()
     err  = hc.GetMeanError()
    else: continue

    hcalib.SetBinContent(i+1,mean)
    hcalib.SetBinError(i+1,err)

    hc.Draw('')
    lat.DrawLatex(0.1,0.92," %.2f < %s < %.2f %s "%(ptR,obj,ptBins[i+1],unit))
    output_f_hists.WriteTObject(c)

   hcalib.SetLineColor(lc)
   hcalib.SetLineWidth(2)
   hcalib.SetFillColor(fc)
   hcalib.SetMarkerStyle(24)
   hcalib.SetMarkerSize(0.6)
   hcalib.SetMarkerColor(lc)
   hcalib.SetMinimum(0.35)
   hcalib.SetMaximum(1.8)
   hcalib.GetXaxis().SetTitleOffset(1.2)
   schemas[index][4]=hcalib.Clone()

   hallproj = h2d_calib.ProjectionX("prj_%s_%s_all"%(name,ext))
   hallproj.SetLineColor(fc)
   hallproj.SetLineWidth(3)
   hallproj.SetFillColor(0)
   schemas[index][5]=hallproj.Clone()

 c_vs_pt = ROOT.TCanvas("rsp_vs_%s"%ext,"rsp_vs_%s"%ext,800,600)
 for index in schemas.keys(): 
   if index==0: schemas[index][4].Draw("P") 
   else: schemas[index][4].Draw("Psame")
   #schemas[index][4].Draw('Lsame')
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
 if not cutGenPt:
   #c_vs_pt.SetLogx()
   lCut=ROOT.TLine(genPtCut,0.35,genPtCut,2.0)
   aCut=ROOT.TArrow(genPtCut,1.5,genPtCut+30,1.5)
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
  
 c_proj_rsp = ROOT.TCanvas("response_%s"%ext,"response_%s"%ext,800,600)
 c_proj_rsp.SetTitle(c_proj_rsp.GetTitle()+'_'+cstr )
 if cutGenPt: c_proj_rsp.SetName(c_proj_rsp.GetTitle()+"_ptGen_gt_%g"%genPtCut)
 
 for index in schemas.keys(): 
   name      = schemas[index][0]
   fc        = schemas[index][1]
   lc        = schemas[index][2]
   h2d_calib = schemas[index][3].Clone()

   hallprojY = h2d_calib.ProjectionY("prjY_%s_%s_all"%(name,ext))
   hallprojY.SetLineColor(fc)
   hallprojY.SetLineWidth(3)
   hallprojY.SetFillColor(0)
   hallprojY.GetYaxis().SetTitle("p_{T}^{L1}/p_{T}^{Gen}")
   schemas[index][6]=hallprojY.Clone()

   #schemas[index][6].Fit("gaus")
   #schemas[index][6].GetFunction("gaus").SetLineColor(fc)

   #mu = schemas[index][6].GetFunction("gaus").GetParameter(1)
   #sigma = schemas[index][6].GetFunction("gaus").GetParameter(2)
   mu     = schemas[index][6].GetMean();
   sigma  = schemas[index][6].GetRMS();
   #schemas[index][6].GetFunction("gaus").SetName("gaus_%s"%hallprojY.GetName())

   leg3.AddEntry(schemas[index][6],schemas[index][0]+" <rsp>=%.2f, RMS=%.2f"%(mu,sigma),"L")

 for index in schemas.keys(): 
   if index==0: schemas[index][6].Draw("hist") 
   else: schemas[index][6].Draw("samehist")

 leg3.Draw()

 output_f.WriteTObject(c_vs_pt)
 output_f.WriteTObject(c_proj)
 output_f.WriteTObject(c_proj_rsp)

 for index in schemas.keys():
   c2d = ROOT.TCanvas("response_2d_vs_%s_%s"%(what,schemas[index][0])) 
   schemas[index][3].SetTitle(schemas[index][0])
   schemas[index][3].GetYaxis().SetTitle("p_{T}^{L1}/p_{T}^{Gen}")
   schemas[index][3].GetXaxis().SetTitle(obj+' '+unit)
   schemas[index][3].Draw("colz")
   output_f.WriteTObject(c2d)
 outname = (outputfile.GetName().split('.'))[0]
 
 # write the profile histos too 
 for index in schemas.keys():
  output_f_hists.WriteTObject(schemas[index][4])
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

#ptBins  = [4,6,8,10,12,14,16,18,20,22.5,25,30,35,40,45,50,55,60,65,70,75,80,90,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,325,350]
#etaBins = [-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0,1.2,1.4,1.6,1.8,2.0,2.1,2.2,2.3,2.4,2.5]
ptBins  = range(5,300,4)
etaBins = list(numpy.arange(-3,3,0.1))
vtxBins  = range(20,65,1)
soptBins = list(numpy.arange(0,1,0.025))
hoeBins =  list(numpy.arange(0,1,0.025))
# Vs PT

makeResponsePlots(inputf,output_f,ptBins,"matched_gpt","PT",'P_{T}^{Gen}','(GeV)',300,0,300,False) # 0 = pt, 1 = eta
makeResponsePlots(inputf,output_f,etaBins,"matched_geta","Eta",'|#eta|^{Gen}','',600,-3,3,True) # 0 = pt, 1 = eta
makeResponsePlots(inputf,output_f,vtxBins,"npv","N_pv",'N PV','',80,0,80,True) # 0 = pt, 1 = eta
makeResponsePlots(inputf,output_f,soptBins,"sopt","S_over_Pt",'Seed/pT','',400,0,1,True) # 0 = pt, 1 = eta
makeResponsePlots(inputf,output_f,hoeBins,"hoe","H_over_E",'H/(E+H)','',400,0,1,True) # 0 = pt, 1 = eta

