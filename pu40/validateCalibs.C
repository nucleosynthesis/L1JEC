#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TMVA/Reader.h"
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

const int nEtaBins=7;
bool matchpt = false;

const int MAXEVENTS = 200000;
const int MAXJETS = 99;

double minseed = 2.5;  // do we really want this guy?
double minpt     = 1.0;
double mingenpt  = 5.0;
double maxeta = 3.0;
double maxdr  = 0.5; 
int maxn = 2;  // 0,1,2

bool SUBDONUT=true;

int CALIBMODE=0;

/* Calib modes ...
	0 = no calibration
	1 = Old-style calibs
	2 = 2-parameter Regression
	3 = 3-parameter Regression (H/E)
	4 = 3-parameter Regression (S/PT)
	5 = 4-parameter Regression (combine above)
*/

std::string ptreweightfile = "qcdPtReweight.root";
std::string pthist	   = "qcdpt";

struct jet {
   double pt;
   double phi;
   double eta;

   double dnut;
   double es;
   double hs;
   double e1;
   double e2;
   double e3;
   double e4;
   double h1;
   double h2;
   double h3;
   double h4;

   double hoe;
   double sopt;

   double rms_phi;
   double rms_eta;
   double corr;
   int jid;
}; 

struct jpair {

   jet l1_jet;
   jet rf_jet;

   double dr;
   double deta;
   double dphi;

   int nnear;
};

double deltaphi(double phi1, double phi2){
	// do somthing to sort out delta R??
	double mpi = TMath::Pi();
	double res = phi1 - phi2;
	while (res > mpi) res-= 2*mpi;
	while (res <= -mpi) res+= 2*mpi;
	return res;

}

double deltaR(double phi1, double phi2, double eta1, double eta2){
	
	double dele = fabs(eta1-eta2);
	double delp = fabs(deltaphi(phi1,phi2));

	return sqrt(dele*dele+delp*delp);
}

bool sortl1pt(jpair i, jpair j){
	return i.l1_jet.pt > j.l1_jet.pt;
}

bool sortpt(jet i, jet j){
	return i.pt > j.pt;
}

bool ismatched(int i, std::vector<int> vint){
	for (std::vector<int>::iterator vint_it = vint.begin(); vint_it!=vint.end();vint_it++){
	   if (i==*vint_it) return true;
	}
	return false;
}
void makeTrees(std::string l1_collection ,std::string outDir, TFile *fout, TFile *fin){


	TRandom3 *globalRand = new TRandom3();

	TFile *fptQCD = TFile::Open(ptreweightfile.c_str());

	TH1F *hptQCD;

	if (fptQCD){
	   hptQCD  = (TH1F*)fptQCD->Get(pthist.c_str());
 	   hptQCD->Scale(1./hptQCD->Integral());
	   std::cout << " Found QCD PT weights " << hptQCD->GetName() << std::endl;
	}

	TDirectory *tdir = fout->mkdir(outDir.c_str());
	tdir->cd();
	TTree *evtree  = new TTree("events","events");

	int evt_n_rfjet,evt_n_l1jet,evt_ev_number;

	evtree->Branch("n_rfjet"   ,&evt_n_rfjet ,"n_rfjet/Int_t");
	evtree->Branch("n_l1jet"   ,&evt_n_l1jet ,"n_l1jet/Int_t");
	evtree->Branch("ev_number"   ,&evt_ev_number ,"ev_number/Int_t");

	TTree *outtree = new TTree("valid","valid");
	
	std::string calibtype="BDTG";
	
	double etaBins[nEtaBins+1] = { 0.0, 0.348, 0.695, 1.044, 1.392, 1.74, 2.172, 3.0};
	std::vector<TF1*> tf1Calibs;
	if (CALIBMODE==1){
		TFile *fCalib = TFile::Open("calibration.root");
		for (int ei=0;ei<nEtaBins;ei++){
			std::cout << Form("fitfcneta_%g_%g",etaBins[ei],etaBins[ei+1]) << std::endl;
			TF1 *fc = (TF1*) fCalib->Get(Form("fitfcneta_%g_%g",etaBins[ei],etaBins[ei+1]));
			std::cout << "Adding Calibration Function" << fc->GetName() << std::endl;
			tf1Calibs.push_back(fc);
		}
	}

	float tmPT, tmETA, tmHOE, tmSOPT;  // tmvaVariables to hold on to!
	std::string rgeweights;
	if (CALIBMODE==2) rgeweights =  Form("weights/TMVARegresion_%s.weights.xml",calibtype.c_str());
	if (CALIBMODE==3) rgeweights =  Form("weights/TMVARegresion_%s_3var_hoe.weights.xml",calibtype.c_str());
	if (CALIBMODE==4) rgeweights =  Form("weights/TMVARegresion_%s_3var_sopt.weights.xml",calibtype.c_str());
	if (CALIBMODE==5) rgeweights =  Form("weights/TMVARegresion_%s_4var.weights.xml",calibtype.c_str());
	TMVA::Reader *tmvaReader_ = new TMVA::Reader();
 	tmvaReader_->AddVariable("pt" ,&tmPT );
 	tmvaReader_->AddVariable("eta",&tmETA);
 	if (CALIBMODE==3)  tmvaReader_->AddVariable("hoe",&tmHOE);
 	if (CALIBMODE==4)  tmvaReader_->AddVariable("sopt",&tmSOPT);
 	if (CALIBMODE==5)  {
		tmvaReader_->AddVariable("hoe",&tmHOE);
		tmvaReader_->AddVariable("sopt",&tmSOPT);
	}
	if (CALIBMODE>=2  ) tmvaReader_->BookMVA(calibtype.c_str(),rgeweights.c_str());

	int n_l1,n_rf; 
	int l1_jrank;
	double rho;

	float out_ptweight;

	//float pt_l1[MAXJETS], eta_l1[MAXJETS], phi_l1[MAXJETS];
	//float pt_rf[MAXJETS], eta_rf[MAXJETS], phi_rf[MAXJETS];
	std::vector<float> *pt_rf = 0;
	std::vector<float> *eta_rf = 0;
	std::vector<float> *phi_rf = 0;

	std::vector<float> *pt_l1 = 0;
	std::vector<float> *eta_l1 = 0;
	std::vector<float> *phi_l1 = 0;
		
	std::vector<float> *es_l1 = 0;
	std::vector<float> *e1_l1 = 0;
	std::vector<float> *e2_l1 = 0;
	std::vector<float> *e3_l1 = 0;
	std::vector<float> *e4_l1 = 0;

	std::vector<float> *hs_l1 = 0;
	std::vector<float> *h1_l1 = 0;
	std::vector<float> *h2_l1 = 0;
	std::vector<float> *h3_l1 = 0;
	std::vector<float> *h4_l1 = 0;

	std::vector<float> *donut_l1 = 0;
	std::vector<float> *jet_area = 0;
	std::vector<float> *jet_area_strip1= 0;
	std::vector<float> *jet_area_strip2= 0;
	std::vector<float> *jet_energy_strip1= 0;
	std::vector<float> *jet_energy_strip2= 0;
	std::vector<float> *rms_phi= 0;
	std::vector<float> *rms_eta= 0;
	// output branches!
	int out_nnear;
	int out_npv;
	float out_rho, out_hoe, out_l1dr;
	float out_pt, out_eta, out_phi, out_rsp, out_dr, out_deta, out_dphi;
	float out_eseed, out_hseed, out_er1,out_er2,out_er3,out_er4,out_hr1,out_hr2,out_hr3,out_hr4;
	float out_donut, out_rms_phi, out_rms_eta;
	float out_matchedgpt, out_matchedgeta,out_l1corr, out_dnut_on_pt, out_seed_on_pt;

	outtree->Branch("pt"   ,&out_pt ,"pt/Float_t");
	outtree->Branch("eta"  ,&out_eta,"eta/Float_t");
	outtree->Branch("phi"  ,&out_phi,"phi/Float_t");
	outtree->Branch("rsp"  ,&out_rsp,"rsp/Float_t");
	outtree->Branch("dr"   ,&out_dr ,"dr/Float_t");
	outtree->Branch("rho"  ,&out_rho ,"rho/Float_t");
	outtree->Branch("deta" ,&out_deta ,"deta/Float_t");
	outtree->Branch("dphi" ,&out_dphi ,"dphi/Float_t");
	outtree->Branch("nj"   ,&n_l1 ,"nj/Int_t");
	outtree->Branch("ngen" ,&n_rf ,"ngen/Int_t");
	outtree->Branch("nnear",&out_nnear ,"nnear/Int_t");
	outtree->Branch("jrank",&l1_jrank ,"jrank/Int_t");

	// Some "Donut" jets (pre-pu removal)
	outtree->Branch("eseed",&out_eseed,"eseed/Float_t");
	outtree->Branch("er1"  ,&out_er1  ,"er1/Float_t");
	outtree->Branch("er2"  ,&out_er2  ,"er2/Float_t");
	outtree->Branch("er3"  ,&out_er3  ,"er3/Float_t");
	outtree->Branch("er4"  ,&out_er4  ,"er4/Float_t");

	outtree->Branch("hseed",&out_hseed,"hseed/Float_t");
	outtree->Branch("hr1"  ,&out_hr1  ,"hr1/Float_t");
	outtree->Branch("hr2"  ,&out_hr2  ,"hr2/Float_t");
	outtree->Branch("hr3"  ,&out_hr3  ,"hr3/Float_t");
	outtree->Branch("hr4"  ,&out_hr4  ,"hr4/Float_t");
	
	outtree->Branch("hoe"  ,&out_hoe  ,"hoe/Float_t");
	outtree->Branch("donut",&out_donut,"donut/Float_t");
	outtree->Branch("rms_phi",&out_rms_phi,"rms_phi/Float_t");
	outtree->Branch("rms_eta",&out_rms_eta,"rms_eta/Float_t");
	outtree->Branch("l1dr",&out_l1dr,"l1dr/Float_t");

	outtree->Branch("matched_gpt",&out_matchedgpt,"matched_gpt/Float_t");
	outtree->Branch("matched_geta",&out_matchedgeta,"matched_geta/Float_t");
	outtree->Branch("l1corr",&out_l1corr,"l1corr/Float_t");
	outtree->Branch("dopt",&out_dnut_on_pt,"dopt/Float_t");
	outtree->Branch("sopt",&out_seed_on_pt,"sopt/Float_t");

	outtree->Branch("ptweight",&out_ptweight,"ptweight/Float_t");

	outtree->Branch("npv",&out_npv,"npv/Int_t");
	
	// choose the L1/Ref collectiom to train
	//std::string l1_collection = "JetTree/PrePUS";
	//std::string rf_collection = "Ak5";

        TTree *l1_tree =( TTree*) fin->Get(Form("%s",l1_collection.c_str()));
        //TTree *rf_tree =( TTree*) fin->Get(Form("%s",l1_collection.c_str()));

	//if (rf_tree->GetEntries()!=l1_tree->GetEntries())  {
	//	std::cout << "Different number of entries in the trees???!?!?!? " << std::endl;
	//	assert(0);
	//}

	//l1_tree->SetBranchAddress("njets",&n_l1);
	l1_tree->SetBranchAddress("jetPt_L1_for_Nick",&pt_l1);
	l1_tree->SetBranchAddress("jetEta_L1_for_Nick",&eta_l1);
	l1_tree->SetBranchAddress("jetPhi_L1_for_Nick",&phi_l1);
	l1_tree->SetBranchAddress("jetDonut_L1_for_Nick",&donut_l1);

	l1_tree->SetBranchAddress("jetRingSumsECAL_0_L1_for_Nick",&es_l1);
	l1_tree->SetBranchAddress("jetRingSumsECAL_1_L1_for_Nick",&e1_l1);
	l1_tree->SetBranchAddress("jetRingSumsECAL_2_L1_for_Nick",&e2_l1);
	l1_tree->SetBranchAddress("jetRingSumsECAL_3_L1_for_Nick",&e3_l1);
	l1_tree->SetBranchAddress("jetRingSumsECAL_4_L1_for_Nick",&e4_l1);
	
	l1_tree->SetBranchAddress("jetRingSumsHCAL_0_L1_for_Nick",&hs_l1);
	l1_tree->SetBranchAddress("jetRingSumsHCAL_1_L1_for_Nick",&h1_l1);
	l1_tree->SetBranchAddress("jetRingSumsHCAL_2_L1_for_Nick",&h2_l1);
	l1_tree->SetBranchAddress("jetRingSumsHCAL_3_L1_for_Nick",&h3_l1);
	l1_tree->SetBranchAddress("jetRingSumsHCAL_4_L1_for_Nick",&h4_l1);

	//rf_tree->SetBranchAddress("njets",&n_rf);
	l1_tree->SetBranchAddress("jetPt_ak4_gen",&pt_rf);
	l1_tree->SetBranchAddress("jetEta_ak4_gen",&eta_rf);
	l1_tree->SetBranchAddress("jetPhi_ak4_gen",&phi_rf);

	l1_tree->SetBranchAddress("medianRho",&rho);

	l1_tree->SetBranchAddress("jetArea_L1_for_Nick",&jet_area);
	l1_tree->SetBranchAddress("jetOuterStripsArea_1_L1_for_Nick",&jet_area_strip1);
	l1_tree->SetBranchAddress("jetOuterStripsArea_2_L1_for_Nick",&jet_area_strip2);
	l1_tree->SetBranchAddress("jetOuterStripsEnergy_1_L1_for_Nick",&jet_energy_strip1);
	l1_tree->SetBranchAddress("jetOuterStripsEnergy_2_L1_for_Nick",&jet_energy_strip2);

	l1_tree->SetBranchAddress("jetSecEta_L1_for_Nick",&rms_eta);
	l1_tree->SetBranchAddress("jetSecPhi_L1_for_Nick",&rms_phi);

	l1_tree->SetBranchAddress("mNPV",&out_npv);

	int nTotalTreeEntries = l1_tree->GetEntries() ;

	// Can take a random percentage of the events 
	double fraction=1.;
	if ( (MAXEVENTS) > 0 ){
	   fraction = (double)MAXEVENTS/l1_tree->GetEntries() ;
	   if (fraction>1) fraction = 1.;
	}

	int n_processed_events=0;

	// Now BEGIN LOOP!
	std::cout << "fraction of events = " << fraction << std::endl;

	// fill a list of "jets" which pass some basic cuts
	for (unsigned int ev_i=0;ev_i<l1_tree->GetEntries();ev_i++){
	  if (( n_processed_events > MAXEVENTS) && MAXEVENTS > 0 ) break;

	  if ( (MAXEVENTS > 0) ) {
	  	// skip 1-fraction of events!
		if ( globalRand->Uniform(0,1) < 1.-fraction) continue;
	  }

	  l1_tree->GetEntry(ev_i);
	  //rf_tree->GetEntry(ev_i);

	  std::vector<jet> l1_jets;
	  std::vector<jet> rf_jets;

	  // L1 collection
	  //n_l1 = pt_l1->size(); 
	  for (unsigned int i=0;i<pt_l1->size();i++){
	  	if ( i > MAXJETS) break;
	  	if ( i > pt_l1->size()) break;
		if ( fabs((*eta_l1)[i]) > maxeta ) continue;

		tmETA = (*eta_l1)[i];
		tmPT  = (*pt_l1)[i];

	        // apply seed threshold 
		if (0.5*( (*es_l1)[i] + (*hs_l1)[i] ) < minseed) continue;

		if (SUBDONUT) {  // leave in unphysical scale 
			double stripE = (*jet_energy_strip1)[i] + (*jet_energy_strip2)[i] ;
			double dnutE = stripE*( (*jet_area)[i]/( (*jet_area_strip1)[i] + (*jet_area_strip2)[i] ) );
			//std::cout << "pt, eta, dnut " << tmPT <<", "<< tmETA << ", " << dnutE << std::endl;
			tmPT -= dnutE; 
		}
		tmPT = (float) 0.5*tmPT;  // not in physical scale yet"
		if ( tmPT < minpt ) continue;   // cut pt before calibration?

		tmSOPT = ((*hs_l1)[i]+(*es_l1)[i])/(*pt_l1)[i] ;  // before Donut
		tmHOE  = ( (*hs_l1)[i] 
                          + (*h1_l1)[i]
                          + (*h2_l1)[i]
                          + (*h3_l1)[i]
                          + (*h4_l1)[i])
                          / ((*e1_l1)[i]
                          + (*e2_l1)[i]
                          + (*e3_l1)[i]
                          + (*e4_l1)[i]
                          + (*es_l1)[i]
			  + (*hs_l1)[i]
			  + (*h1_l1)[i]
			  + (*h2_l1)[i]
			  + (*h3_l1)[i]
			  + (*h4_l1)[i]
			  );
			
		// Now give a chance to calibrate
		double corr =1.;
		if (CALIBMODE>=2){
			corr = tmvaReader_->EvaluateRegression(0,calibtype.c_str());
		//	std::cout << "pt, eta, correction " << tmPT <<", "<< tmETA << ", " << corr<< std::endl;
		} else if (CALIBMODE==1){
			//Find the Eta Bin and pick up the TF1 Function!
			double abseta = TMath::Abs(tmETA);
			int cIndex = 0;
			for (int ei=0;ei<nEtaBins;ei++){
				if (abseta > etaBins[ei] && abseta <= etaBins[ei+1]) 
				cIndex = ei;
				break;
			}
			corr = tf1Calibs[cIndex]->Eval(tmPT);
		}
		tmPT = tmPT*corr;
		
		jet l1;
		l1.pt  = tmPT;
	        l1.sopt = tmSOPT;	
		l1.corr = corr;
		l1.eta = (*eta_l1)[i];  
		l1.phi = (*phi_l1)[i];  
		l1.dnut  = 0.5*((*jet_energy_strip1)[i] + (*jet_energy_strip2)[i]);//(*donut_l1)[i];

		l1.es = (*es_l1)[i];
		l1.hs = (*hs_l1)[i];

		l1.e1    = (*e1_l1)[i] ;
		l1.e2    = (*e2_l1)[i] ;
		l1.e3    = (*e3_l1)[i] ;
		l1.e4    = (*e4_l1)[i] ;

		l1.h1    = (*h1_l1)[i] ;
		l1.h2    = (*h2_l1)[i] ;
		l1.h3    = (*h3_l1)[i] ;
		l1.h4    = (*h4_l1)[i] ;
		l1.hoe   =( (*hs_l1)[i] 
                          + (*h1_l1)[i]
                          + (*h2_l1)[i]
                          + (*h3_l1)[i]
                          + (*h4_l1)[i])
                          / ( (*e1_l1)[i]
                          + (*e2_l1)[i]
                          + (*e3_l1)[i]
                          + (*e4_l1)[i]
                          + (*es_l1)[i]
			  + (*hs_l1)[i]
			  + (*h1_l1)[i]
			  + (*h2_l1)[i]
			  + (*h3_l1)[i]
			  + (*h4_l1)[i]
			  );

		l1.rms_phi = (*rms_phi)[i];
		l1.rms_eta = (*rms_eta)[i];
		l1.jid = i;
		l1_jets.push_back(l1);
	  }
	  // Ref collection 
	  n_rf = pt_rf->size(); 
	  for (unsigned int i=0;i<n_rf;i++){
	  	if ( i > MAXJETS) break;
	  	if ( i > pt_rf->size()) break;
		if ( fabs((*eta_rf)[i]) > maxeta ) continue;
		if ( ((*pt_rf)[i]) < mingenpt ) continue;
		jet rf;
		rf.pt  = (*pt_rf)[i];  
		rf.eta = (*eta_rf)[i];  
		rf.phi = (*phi_rf)[i];  
		rf.jid=i;
		rf_jets.push_back(rf);
	  }

	  std::vector<jpair> jet_pairs;
	  std::sort(l1_jets.begin(),l1_jets.end(),sortpt);
	  std::sort(rf_jets.begin(),rf_jets.end(),sortpt);

	  // loop through each l1 jet and find closest ref jet, then make a pair 
	  std::vector<int> matched_refs;
	  for (std::vector<jet>::iterator l1_it = l1_jets.begin(); l1_it!=l1_jets.end(); l1_it++){
	    double dr = 999., ptdiff = 999.;
	    double deta = 0., dphi = 0.;
	    int rfid  = 0, i_rf = 0, nnear = 0;

	    if (matchpt){
	      for (std::vector<jet>::iterator rf_it = rf_jets.begin(); rf_it!=rf_jets.end(); rf_it++){
	    	  double ndr = deltaR(l1_it->phi,rf_it->phi,l1_it->eta,rf_it->eta);
		  if (ndr > maxdr) continue;
		  double ptd = fabs(l1_it->pt - rf_it->pt);
		  if (ptd <= ptdiff){
		     if (ptd == ptdiff) {
			if (ndr > dr) continue;
		     }
			dr = ndr; 
			ptdiff = ptd;
			deta = l1_it->eta - rf_it->eta ; 
			dphi = deltaphi(l1_it->phi, rf_it->phi) ; 
			rfid = i_rf;	
		  }
		  i_rf++;
	      }
	    } else {

	      for (std::vector<jet>::iterator rf_it = rf_jets.begin(); rf_it!=rf_jets.end(); rf_it++){
	    	  double ndr = deltaR(l1_it->phi,rf_it->phi,l1_it->eta,rf_it->eta);
	    	  if (ndr < maxdr) nnear+=1;
	    	  if (ismatched(i_rf,matched_refs)) continue;
		  if (ndr < dr){
			dr = ndr; 
			deta = l1_it->eta - rf_it->eta ; 
			dphi = deltaphi(l1_it->phi, rf_it->phi) ; 
			rfid = i_rf;
		  }
		  i_rf++;
	      }
	     }
	     if ( dr > maxdr ) continue;
	     jpair cpair;
	     cpair.l1_jet = *l1_it;
	     cpair.rf_jet = rf_jets[rfid];
	     cpair.dr     = dr;
	     cpair.deta     = deta;
	     cpair.dphi     = dphi;
	     cpair.nnear    = nnear;

	     jet_pairs.push_back(cpair);
	     matched_refs.push_back(rfid);
	  }

	  // print jets in the event (up to 10 say?)
	  /*
	  std::cout << " Gen Jets " << std::endl;
	  for (std::vector<jet>::iterator rf_it = rf_jets.begin(); rf_it!=rf_jets.end(); rf_it++){
	  	
	  	std::cout << "pt " << rf_it->pt << ", eta " << rf_it->eta << ", phi " << rf_it->phi <<std::endl;
		
	  }
	  std::cout << " L1 Jets " << std::endl;
	  for (std::vector<jet>::iterator rf_it = l1_jets.begin(); rf_it!=l1_jets.end(); rf_it++){
	  	
	  	std::cout << "pt " << rf_it->pt << ", eta " << rf_it->eta << ", phi " << rf_it->phi <<std::endl;
		
	  }
	  */
	  // Sort pairs according to L1 jet pt
	  /*
	  std::cout << "Unsorted ..... " << std::endl;
	  for (std::vector<jpair>::iterator pit = jet_pairs.begin(); pit!=jet_pairs.end(); pit++){
	  	std::cout << "pt " << pit->l1_jet.pt << ", " << pit->rf_jet.pt << std::endl;
	  	std::cout << "eta " << pit->l1_jet.eta << ", " << pit->rf_jet.eta << std::endl;
	  	std::cout << "phi " << pit->l1_jet.phi << ", " << pit->rf_jet.phi << std::endl;
	  	std::cout << "dr " << pit->dr << std::endl;
	  }
	  */
	  std::sort(jet_pairs.begin(),jet_pairs.end(),sortl1pt);
	  /*std::cout << "Sorted ..... " << std::endl;
	  for (std::vector<jpair>::iterator pit = jet_pairs.begin(); pit!=jet_pairs.end(); pit++){
	  	std::cout << "pt " << pit->l1_jet.pt << ", " << pit->rf_jet.pt << std::endl;
	  	std::cout << "eta " << pit->l1_jet.eta << ", " << pit->rf_jet.eta << std::endl;
	  	std::cout << "phi " << pit->l1_jet.phi << ", " << pit->rf_jet.phi << std::endl;
	  	std::cout << "dr " << pit->dr << std::endl;
	  }*/

	  // loop through pairs, then fill the output tree, target is ref_pt/l1_pt
	  int jpair_i=0;
	  for (std::vector<jpair>::iterator pit = jet_pairs.begin(); pit!=jet_pairs.end(); pit++){
	    if (jpair_i > maxn) break;

	    int noverlap = 0;
	    double closel1dr = 999.0;
	    // for each L1 jet we also count how many Overlapping jets there are .... / deltaR < 0.4 == overlap!
	    for (std::vector<jet>::iterator l1_it = l1_jets.begin(); l1_it!=l1_jets.end(); l1_it++){
	    	  if (pit->l1_jet.jid == l1_it->jid) continue;
	    	  double ndr = deltaR(l1_it->phi,pit->l1_jet.phi,l1_it->eta,pit->l1_jet.eta);
		  if (ndr < 0.4) noverlap++;	
		  if (ndr < closel1dr) closel1dr = ndr;

	    } 
	    noverlap--; // always count this jet so "uncount" it

	    out_pt  = pit->l1_jet.pt;
	    out_eta = pit->l1_jet.eta;
	    out_phi = pit->l1_jet.phi;
	    out_rsp = pit->rf_jet.pt/pit->l1_jet.pt;
	    out_dr  = pit->dr;
	    out_l1dr  = closel1dr;
	    out_deta = pit->deta;
	    out_dphi  = pit->dphi;
	    out_nnear = noverlap;//pit->nnear;
	    out_rho   = (float)rho;

	    out_donut =  pit->l1_jet.dnut;
	    out_eseed =  pit->l1_jet.es;
	    out_er1   =  pit->l1_jet.e1;
	    out_er2   =  pit->l1_jet.e2;
	    out_er3   =  pit->l1_jet.e3;
	    out_er4   =  pit->l1_jet.e4;
		 
	    out_hseed =  pit->l1_jet.hs;
	    out_hr1   =  pit->l1_jet.h1;
	    out_hr2   =  pit->l1_jet.h2;
	    out_hr3   =  pit->l1_jet.h3;
	    out_hr4   =  pit->l1_jet.h4;
	    out_hoe = pit->l1_jet.hoe;
	    out_rms_phi = pit->l1_jet.rms_phi;
	    out_rms_eta = pit->l1_jet.rms_eta;
	    out_matchedgpt = pit->rf_jet.pt;
	    out_matchedgeta = pit->rf_jet.eta;
	    out_l1corr     = pit->l1_jet.corr;
	    out_dnut_on_pt = (pit->l1_jet.dnut)/pit->l1_jet.pt;
	    out_seed_on_pt = pit->l1_jet.sopt;

	    n_l1 = l1_jets.size();
	    n_rf = rf_jets.size();

	    l1_jrank = jpair_i;

	    out_ptweight = 1.;
	    if (hptQCD) {
	      if (pit->rf_jet.pt < hptQCD->GetXaxis()->GetXmax()){
	         out_ptweight = 1./(hptQCD->GetBinContent(hptQCD->FindBin(pit->rf_jet.pt)));
	      }
	    }

	    outtree->Fill();
	    jpair_i++;
	  }
	 n_processed_events++;

	 //std::cout << "Processed " << ev_i << std::endl;
	 evt_n_rfjet = rf_jets.size(); 
	 evt_n_l1jet = l1_jets.size(); 
	 evt_ev_number = ev_i;

	 evtree->Fill();
	}
     tdir->cd();
     outtree->Write();
     evtree->Write();
     fptQCD->Close();
}

void validateCalibs(std::string filenam, std::string outname){
	TFile *fin  = TFile::Open(filenam.c_str());
	TFile *fout = new TFile(outname.c_str(),"RECREATE");
	CALIBMODE=0;
	makeTrees("demo/L1Tree","Uncalib",fout,fin);
	CALIBMODE=1;
	makeTrees("demo/L1Tree","Calib_param",fout,fin);	
	CALIBMODE=2;
	makeTrees("demo/L1Tree","Calib_2var",fout,fin);
	CALIBMODE=3;
	makeTrees("demo/L1Tree","Calib_3var_hoe",fout,fin);
	CALIBMODE=4;
	makeTrees("demo/L1Tree","Calib_3var_sopt",fout,fin);
	CALIBMODE=5;
	makeTrees("demo/L1Tree","Calib_4var",fout,fin);
	fout->Close();
	fin->Close();	
//	makeTrees("LPUS",fout,fin);
//	makeTrees("CalibLPUS",fout,fin);
}

