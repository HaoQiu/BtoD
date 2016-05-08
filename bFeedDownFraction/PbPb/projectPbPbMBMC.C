void projectPbPbMBMC()
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

  TFile* infMCP = new TFile("../../DNtuple/ntD_EvtBase_20160303_Dfinder_20160302_PbPb_Pythia8_prompt_D0_dPt1tkPt1_pthatweight.root");
  TFile* infMCNP = new TFile("../../DNtuple/ntD_EvtBase_20160303_Dfinder_20160302_PbPb_Pythia8_nonprompt_D0_dPt1tkPt1_pthatweight.root");
  
  TTree* tMCP = (TTree*)infMCP->Get("ntDkpi");
  TTree* tMCNP = (TTree*)infMCNP->Get("ntDkpi");
  TTree* tMCPHI = (TTree*)infMCP->Get("ntHi");
  TTree* tMCNPHI = (TTree*)infMCNP->Get("ntHi");
  tMCP->AddFriend(tMCPHI);
  tMCNP->AddFriend(tMCNPHI);
  
  TFile* outf = new TFile("bFeedDownPbPbMBMC.hist.root","recreate");
  
  TCut cutPbPbMB = "Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>1.0&&Dtrk2Pt>1.0&&Dtrk1PtErr/Dtrk1Pt<0.3&&Dtrk2PtErr/Dtrk2Pt<0.3&&abs(Dtrk1Eta)<1.5&&abs(Dtrk2Eta)<1.5&&((DlxyBS/DlxyBSErr)>1.5&&Dalpha<0.12&&((Dpt>1&&Dpt<2&&(DsvpvDistance/DsvpvDisErr)>6.0&&Dchi2cl>0.25)||(Dpt>2&&Dpt<4&&(DsvpvDistance/DsvpvDisErr)>5.86&&Dchi2cl>0.224)||(Dpt>4&&Dpt<5&&(DsvpvDistance/DsvpvDisErr)>5.46&&Dchi2cl>0.196)||(Dpt>5&&Dpt<6&&(DsvpvDistance/DsvpvDisErr)>4.86&&Dchi2cl>0.170)||(Dpt>6&&Dpt<8&&(DsvpvDistance/DsvpvDisErr)>4.54&&Dchi2cl>0.125)||(Dpt>8&&Dpt<10&&(DsvpvDistance/DsvpvDisErr)>4.42&&Dchi2cl>0.091)||(Dpt>10&&Dpt<15&&(DsvpvDistance/DsvpvDisErr)>4.06&&Dchi2cl>0.069)||(Dpt>15&&Dpt<20&&(DsvpvDistance/DsvpvDisErr)>3.71&&Dchi2cl>0.056)||(Dpt>20&&Dpt<25&&(DsvpvDistance/DsvpvDisErr)>3.25&&Dchi2cl>0.054)||(Dpt>25&&(DsvpvDistance/DsvpvDisErr)>2.97&&Dchi2cl>0.050)))";

  TCut cutmc = "(Dgen==23333||Dgen==23344)";
  TCut cutmcSignal = "Dgen==23333";
  TCut cutmcSwapped = "Dgen==23344";
  TCut cutpt = "Dpt<20";
  TCut cutSignal = "abs(Dmass-1.8649)<0.025";
  TCut cutSideband = "abs(Dmass-1.8649)>0.075&&abs(Dmass-1.8649)<0.1";
  TCut weightfunctionreco = "pthatweight";

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
  
  const int nBinM = 60;
  Float_t binsM[nBinM];
  float minMassBin = 1.7;
  float massBinWidth = 0.005;
  for(int i=0; i<=nBinM; i++)
    binsM[i] = minMassBin + massBinWidth*i;

  TH3D* hMCPSignal = new TH3D("hMCPSignal",";p_{T} (GeV/c);D0 DCA (cm);flight distance significance",nBinX,binsX,nBinY,binsY,nBinZ,binsZ);
  TH3D* hMCNPSignal = new TH3D("hMCNPSignal",";p_{T} (GeV/c);D0 DCA (cm);flight distance significance",nBinX,binsX,nBinY,binsY,nBinZ,binsZ);
  TH2D* hPtMMCPSignal = new TH2D("hPtMMCPSignal",";p_{T} (GeV/c);m (GeV/c^{2})",nBinX,binsX,nBinM,binsM);
  TH2D* hPtMMCPSwapped = new TH2D("hPtMMCPSwapped",";p_{T} (GeV/c);m (GeV/c^{2})",nBinX,binsX,nBinM,binsM);
  
  hMCPSignal->Sumw2();
  hMCNPSignal->Sumw2();
  hPtMMCPSignal->Sumw2();
  hPtMMCPSwapped->Sumw2();
  
  tMCP->Draw("DsvpvDistance/DsvpvDisErr:DsvpvDistance*sin(Dalpha):Dpt>>hMCPSignal",weightfunctionreco*(cutpt&&cutPbPbMB&&cutSignal&&cutmc));
  tMCNP->Draw("DsvpvDistance/DsvpvDisErr:DsvpvDistance*sin(Dalpha):Dpt>>hMCNPSignal",weightfunctionreco*(cutpt&&cutPbPbMB&&cutSignal&&cutmc));
  tMCP->Draw("Dmass:Dpt>>hPtMMCPSignal",weightfunctionreco*(cutpt&&cutPbPbMB&&cutmcSignal));
  tMCP->Draw("Dmass:Dpt>>hPtMMCPSwapped",weightfunctionreco*(cutpt&&cutPbPbMB&&cutmcSwapped));
  
  outf->Write();
}

