void projectPPMB()
{
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* infData = new TFile("../../DNtuple/skim_ncand_D0Dntuple_crab_pp_ALLMinimumBias_AOD_D0_tkpt0p5_Ds_01212016.root");
  
  TTree* tData = (TTree*)infData->Get("ntDkpi");
  
  TFile* outf = new TFile("bFeedDownPPMB.hist.root","recreate");
  
  TCut cutpp ="Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>1.0&&Dtrk2Pt>1.0&&Dtrk1PtErr/Dtrk1Pt<0.3&&Dtrk2PtErr/Dtrk2Pt<0.3&&abs(Dtrk1Eta)<2.0&&abs(Dtrk2Eta)<2.0&&(DsvpvDistance/DsvpvDisErr)>3.5&&(DlxyBS/DlxyBSErr)>2.5&&Dchi2cl>0.05&&Dalpha<0.12";

  TCut cutpt = "Dpt<15";
  TCut cutSignal = "abs(Dmass-1.8649)<0.025";
  TCut cutSideband = "abs(Dmass-1.8649)>0.075&&abs(Dmass-1.8649)<0.1";

  const int nBinX = 9;
  Float_t binsX[nBinX+1] = {2,6,8,10,15,20,25,35,50,100};

  const int nBinY = 20;
  Float_t binsY[nBinY+1];
  float firstBinYWidth = 0.001;
  float binYWidthRatio = 1.27;
  binsY[0]=0;
  for(int i=1; i<=nBinY; i++)
    binsY[i] = binsY[i-1]+firstBinYWidth*pow(binYWidthRatio,i-1);
  cout<<"last y bin: "<<binsY[nBinY]<<endl;

  const int nBinZ = 20;
  Float_t binsZ[nBinZ+1];
  float firstBinZ = 3.5;
  float binWidthZ = 5;
  for(int i=0; i<=nBinZ; i++)
    binsZ[i] = firstBinZ+binWidthZ*i;
  
  TH3D* hSideband = new TH3D("hSideband",";p_{T} (GeV/c};D0 DCA (cm);flight distance significance",nBinX,binsX,nBinY,binsY,nBinZ,binsZ);
  TH3D* hData = new TH3D("hData",";p_{T} (GeV/c};D0 DCA (cm);flight distance significance",nBinX,binsX,nBinY,binsY,nBinZ,binsZ);
  
  hSideband->Sumw2();
  hData->Sumw2();
  
  tData->Draw("DsvpvDistance/DsvpvDisErr:DsvpvDistance*sin(Dalpha):Dpt>>hSideband",cutpt&&cutpp&&cutSideband);
  tData->Draw("DsvpvDistance/DsvpvDisErr:DsvpvDistance*sin(Dalpha):Dpt>>hData",cutpt&&cutpp&&cutSignal);
  
  outf->Write();
}


