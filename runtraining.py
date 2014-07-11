# TMVA regression training for the response

import ROOT as r 
fin     = r.TFile.Open("validation.root")
fout     = r.TFile("jetresponsetrained.root","RECREATE")
fout2    = r.TFile("jetresponsetrained.root","RECREATE")
factory = r.TMVA.Factory ("TMVARegresion",fout,"!V:!Silent:Color:DrawProgressBar")
factory.AddVariable("pt","L1 Jet PT","GeV", 'F')
factory.AddVariable("eta","L1 Jet Eta","", 'F')
#factory.AddVariable("rho","Jet Rho (Global)","", 'F')
factory.AddTarget("rsp")
tree = fin.Get("Uncalib/valid")
factory.AddRegressionTree(tree,1.)
factory.PrepareTrainingAndTestTree(r.TCut("jrank<3"),'nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V')
factory.BookMethod(r.TMVA.Types.kBDT,"BDTG","!H:!V:NTrees=800:BoostType=Grad:MaxDepth=5");
#factory.BookMethod(r.TMVA.Types.kBDT,"BDTA","!H:!V:BoostType=AdaBoost");
factory.BookMethod(r.TMVA.Types.kKNN, "KNN", "nkNN=100");
#factory.BookMethod(r.TMVA.Types.kMLP, "MLP", "");
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()


factory2 = r.TMVA.Factory ("TMVARegresion",fout,"!V:!Silent:Color:DrawProgressBar")
factory2.AddVariable("pt","L1 Jet PT","GeV", 'F')
factory2.AddVariable("eta","L1 Jet Eta","", 'F')
factory2.AddVariable("dopt","L1 Donut / L1 PT","", 'F')
factory2.AddTarget("rsp")
factory2.AddRegressionTree(tree,1.)
factory2.PrepareTrainingAndTestTree(r.TCut("jrank<3"),'nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V')
factory2.BookMethod(r.TMVA.Types.kBDT,"BDTG_3var","!H:!V:NTrees=500:BoostType=Grad:MaxDepth=8");
#factory2.BookMethod(r.TMVA.Types.kBDT,"BDTA_3var","!H:!V:BoostType=AdaBoost");
factory2.BookMethod(r.TMVA.Types.kKNN, "KNN_3var", "nkNN=100");
#factory.BookMethod(r.TMVA.Types.kMLP, "MLP", "");
factory2.TrainAllMethods()
factory2.TestAllMethods()
factory2.EvaluateAllMethods()
