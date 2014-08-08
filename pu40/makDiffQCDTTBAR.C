TH1F *g_ngu_pt_hist;
TH1F *g_ttb_pt_hist;

float ngu_weight_f(float x){
  float numerate = g_ngu_pt_hist->Integral()/g_ngu_pt_hist->GetNbinsX();
  float denom	 = g_ngu_pt_hist->GetBinContent(g_ngu_pt_hist->FindBin(x));
  if (denom > 0) return numerate/denom;
  else return 1.;
}

float ttb_weight_f(float x){
  //std::cout << "weight on " << x << std::endl;
  float numerate = g_ttb_pt_hist->Integral()/g_ttb_pt_hist->GetNbinsX();
  float denom	 = g_ttb_pt_hist->GetBinContent(g_ttb_pt_hist->FindBin(x));
 // std::cout << x << " " << g_ttb_pt_hist->FindBin(x) << " " << denom << " " << numerate/denom << std::endl; 
  if (denom > 0) return numerate/denom;
  else return 1.;
}

TGraph *roc(TH1F *hist_s, TH1F *hist_b){

  TGraph *gr = new TGraph();
  gr->SetMarkerSize(0.8);
  gr->SetMarkerStyle(20);
 
  int grp = 0;
  int nb = hist_s->GetNbinsX();
  for (int i=0; i<nb ; i++){
  	gr->SetPoint(grp,hist_s->Integral(i,nb),1-hist_b->Integral(i,nb));
	grp++;
  }

  gr->SetTitle("");
  gr->GetXaxis()->SetTitle(Form("eff %s",hist_s->GetName()));
  gr->GetYaxis()->SetTitle(Form("1-eff %s",hist_b->GetName()));
  return gr;
}

void makDiffQCDTTBAR(std::string sngu, std::string sttb){
 
 gROOT->SetBatch(true);

 const int nvars = 12;
 std::string vars[nvars] = {"eta","pt","sopt","(er4+hr4)/pt","(er1+hr1)/pt","(er2+hr2)/pt","(er3+hr3)/pt","(er1+hr1+er2+hr2)/(er3+hr3+er4+hr4)","(er4+hr4)/(eseed+hseed)","rms_phi","rms_eta","hoe"};
 double xlows[nvars]     = {-3,10,0,0,0,0,0,0,0,0,0,0};
 double xhigh[nvars]     = {3,250,1,2,2,2,2,15,1,5,5,1};
 //const int nvars = 1;
 //std::string vars[nvars] = {"pt"};

 std::string treen = "Calib_2var/valid";

 TFile *fout = new TFile("comparevars.root","RECREATE");

 TFile *fngu = TFile::Open(sngu.c_str());
 TFile *fttb = TFile::Open(sttb.c_str());
 
 TTree *tngu = (TTree*)fngu->Get(treen.c_str());
 TTree *tttb = (TTree*)fttb->Get(treen.c_str());

 tttb->Draw("pt>>httb(50,10,250)","jrank<2 && pt > 10 && pt < 250 ");
 TH1F *tpt = (TH1F*)gROOT->FindObject("httb");
 tngu->Draw("pt>>hngu(50,10,250)","jrank<2 && pt > 10 && pt < 250");
 TH1F *qpt = (TH1F*)gROOT->FindObject("hngu");
 

 g_ngu_pt_hist = (TH1F*)qpt->Clone();
 g_ttb_pt_hist = (TH1F*)tpt->Clone();

 TLegend *lega = new TLegend(0.6,0.72,0.89,0.89);
 qpt->SetLineColor(1);
 tpt->SetLineColor(2);
 qpt->SetLineWidth(2);
 tpt->SetLineWidth(2);

 lega->SetFillColor(0);
 lega->SetTextFont(42);
 lega->AddEntry(qpt,"Neutrino PU-40", "L");
 lega->AddEntry(tpt,"t#bar{t} PU-40", "L");
 qpt->Scale(1./qpt->Integral());
 tpt->Scale(1./tpt->Integral());

 TCanvas * ca = new TCanvas("pt_unweighted","pt",600,600);
 qpt->SetTitle("");
 qpt->GetXaxis()->SetTitle("p_{T}");
 qpt->Draw("hist");
 tpt->Draw("histsame");
 lega->Draw();
 
 fout->WriteTObject(ca);

 // Now draw them
 for (int vi = 0 ; vi < nvars;  vi++){
 
    std::string v = vars[vi];
    std::cout << "Plotting " << v.c_str() << std::endl;
    TCanvas * c = new TCanvas(v.c_str(),v.c_str(),600,600);

    TH1F* hngu = new TH1F("hngu","hngu",50,xlows[vi],xhigh[vi]) ;
    hngu->Sumw2();
    tngu->Draw(Form("%s>>hngu",v.c_str()),"ngu_weight_f(pt)*(jrank<2 && pt > 10 && pt < 250)");
    //tngu->Draw(Form("%s>>hngu(100,%g,%g)",v.c_str(),xlows[vi],xhigh[vi]),"");
    //(TH1F*) gROOT->FindObject("hngu");
    hngu->SetLineColor(1);
    hngu->SetMarkerColor(1);

    TH1F* httb = new TH1F("httb","httb",50,xlows[vi],xhigh[vi]) ;
    httb->Sumw2();
    tttb->Draw(Form("%s>>httb",v.c_str()),"ttb_weight_f(pt)*(jrank<2 && pt>10 && pt < 250)");
    //tttb->Draw(Form("%s>>httb(100,%g,%g)",v.c_str(),xlows[vi],xhigh[vi]),"");
    //TH1F* httb = (TH1F*) gROOT->FindObject("httb");
    httb->SetLineColor(2);
    httb->SetMarkerColor(2);

    hngu->Scale(1./hngu->Integral());
    httb->Scale(1./httb->Integral());

    httb->SetLineWidth(2);
    hngu->SetLineWidth(2);
    //hngu->GetXaxis()->SetRangeUser(minX,maxX);

    hngu->SetTitle("");
    hngu->GetXaxis()->SetTitle(v.c_str());

     TLegend *leg = new TLegend(0.6,0.72,0.89,0.89);
     leg->SetFillColor(0);
     leg->SetTextFont(42);
     leg->AddEntry(hngu,"Neutrino PU-40", "L");
     leg->AddEntry(httb,"t#bar{t} PU-40", "L");
    hngu->Draw("hist");
    httb->Draw("histsame");
    hngu->Draw("sameP");
    httb->Draw("sameP");
     leg->Draw();

    fout->WriteTObject(c);

    TGraph *gr_v = (TGraph*) roc(httb,hngu);
    gr_v->SetName(Form("Graph_%s",v.c_str()));
    fout->WriteTObject(gr_v);

 }
 fout->Close();
}
