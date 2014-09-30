# TMVA regression training for the response
import sys
import ROOT as r 
fin     = r.TFile.Open("qcd_training_globalrho.root")
tree = fin.Get("Uncalib/valid")

REWEIGHT=1

which = int(sys.argv[1])

if which==0:
  fout     = r.TFile("jetresponsetrained.root","RECREATE")
  factory = r.TMVA.Factory ("TMVARegresion",fout,"!V:!Silent:Color:DrawProgressBar")
  factory.AddVariable("pt","L1 Jet PT","GeV", 'F')
  factory.AddVariable("eta","L1 Jet Eta","", 'F')
  factory.AddTarget("rsp")
  factory.SetWeightExpression("ptweight","Regression")
  factory.AddRegressionTree(tree,1.,r.TMVA.Types.kTraining)
  factory.PrepareTrainingAndTestTree(r.TCut("jrank<=1"),'!V:nTrain_regression=500000:nTest_regression=500000:SplitMode=Random')
  factory.BookMethod(r.TMVA.Types.kBDT,"BDTG","!H:!V:NTrees=1000:BoostType=Grad:MaxDepth=5:Shrinkage=0.1:UseBaggedGrad:GradBaggingFraction=0.5");
  factory.TrainAllMethods()

if which==1:
  fout     = r.TFile("jetresponsetrained_3var_sopt.root","RECREATE")
  factory2 = r.TMVA.Factory ("TMVARegresion",fout,"!V:!Silent:Color:DrawProgressBar")
  factory2.AddVariable("pt","L1 Jet PT","GeV", 'F')
  factory2.AddVariable("eta","L1 Jet Eta","", 'F')
  factory2.AddVariable("sopt","L1 Seed on Pt","", 'F')
  factory2.AddTarget("rsp")
  factory2.SetWeightExpression("ptweight","Regression")
  factory2.AddRegressionTree(tree,1.,r.TMVA.Types.kTraining)
  factory2.PrepareTrainingAndTestTree(r.TCut("jrank<=1"),'!V:nTrain_regression=500000:nTest_regression=500000:SplitMode=Random')
  factory2.BookMethod(r.TMVA.Types.kBDT,"BDTG_3var_sopt","!H:!V:NTrees=1000:BoostType=Grad:MaxDepth=5:Shrinkage=0.1:UseBaggedGrad:GradBaggingFraction=0.5");
  factory2.TrainAllMethods()

if which==2:
  fout     = r.TFile("jetresponsetrained_3var_hoe.root","RECREATE")
  factory3 = r.TMVA.Factory ("TMVARegresion",fout,"!V:!Silent:Color:DrawProgressBar")
  factory3.AddVariable("pt","L1 Jet PT","GeV", 'F')
  factory3.AddVariable("eta","L1 Jet Eta","", 'F')
  factory3.AddVariable("hoe","L1 H/E","", 'F')
  factory3.AddTarget("rsp")
  factory3.SetWeightExpression("ptweight","Regression")
  factory3.AddRegressionTree(tree,1.,r.TMVA.Types.kTraining)
  factory3.PrepareTrainingAndTestTree(r.TCut("jrank<=1"),'!V:nTrain_regression=500000:nTest_regression=500000:SplitMode=Random')
  factory3.BookMethod(r.TMVA.Types.kBDT,"BDTG_3var_hoe","!H:!V:NTrees=1000:BoostType=Grad:MaxDepth=5:Shrinkage=0.1:UseBaggedGrad:GradBaggingFraction=0.5");
  factory3.TrainAllMethods()

if which==3:
  fout     = r.TFile("jetresponsetrained_4var.root","RECREATE")
  factory4 = r.TMVA.Factory ("TMVARegresion",fout,"!V:!Silent:Color:DrawProgressBar")
  factory4.AddVariable("pt","L1 Jet PT","GeV", 'F')
  factory4.AddVariable("eta","L1 Jet Eta","", 'F')
  factory4.AddVariable("hoe","L1 H/E","", 'F')
  factory4.AddVariable("sopt","L1 Seed on Pt","", 'F')
  factory4.AddTarget("rsp")
  factory4.SetWeightExpression("ptweight","Regression")
  factory4.AddRegressionTree(tree,1.,r.TMVA.Types.kTraining)
  factory4.PrepareTrainingAndTestTree(r.TCut("jrank<=1"),'!V:nTrain_regression=500000:nTest_regression=500000:SplitMode=Random')
  factory4.BookMethod(r.TMVA.Types.kBDT,"BDTG_4var","!H:!V:NTrees=1000:BoostType=Grad:MaxDepth=5:Shrinkage=0.1:UseBaggedGrad:GradBaggingFraction=0.5");
  factory4.TrainAllMethods()
