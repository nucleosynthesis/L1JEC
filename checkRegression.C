{

	float opt, oeta;

	TMVA::Reader *reader2 = new TMVA::Reader();
	reader2->AddVariable("pt",&opt);
	reader2->AddVariable("eta",&oeta);
	reader2->BookMVA("KNN","weights/TMVARegresion_KNN.weights.xml");

	TH2F *h22 = new TH2F("h22","h22",100,5,200,50,0,2.5);

	for (int i=0;i<100;i++){
	 for (int j=0;j<50;j++){
		opt = h22->GetXaxis()->GetBinCenter(i+1);
		oeta = h22->GetYaxis()->GetBinCenter(j+1);
		
		h22->SetBinContent(i+1,j+1,reader2->EvaluateRegression(0,"KNN"));
	 }
	}	
	TCanvas *can2 = new TCanvas();
	h22->SetTitle("BDTG");
	h22->Draw("colz");
	
}
