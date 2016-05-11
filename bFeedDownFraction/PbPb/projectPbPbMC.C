void projectPbPbMC()
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
  
  TFile* outf = new TFile("bFeedDownPbPbMC.hist.root","recreate");
  
  TCut cutPbPb = "Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>8.5&&Dtrk2Pt>8.5&&Dtrk1PtErr/Dtrk1Pt<0.1&&Dtrk2PtErr/Dtrk2Pt<0.1&&abs(Dtrk1Eta)<2.0&&abs(Dtrk2Eta)<2.0&&(Dtrk1PixelHit+Dtrk1StripHit)>=11&&(Dtrk2PixelHit+Dtrk2StripHit)>=11&&(Dtrk1Chi2ndf/(Dtrk1nStripLayer+Dtrk1nPixelLayer)<0.15)&&(Dtrk2Chi2ndf/(Dtrk2nStripLayer+Dtrk2nPixelLayer)<0.15)&&((DlxyBS/DlxyBSErr)>1.5&&Dalpha<0.12&&((Dpt>2&&Dpt<4&&(DsvpvDistance/DsvpvDisErr)>5.86&&Dchi2cl>0.224)||(Dpt>4&&Dpt<5&&(DsvpvDistance/DsvpvDisErr)>5.46&&Dchi2cl>0.196)||(Dpt>5&&Dpt<6&&(DsvpvDistance/DsvpvDisErr)>4.86&&Dchi2cl>0.170)||(Dpt>6&&Dpt<8&&(DsvpvDistance/DsvpvDisErr)>4.54&&Dchi2cl>0.125)||(Dpt>8&&Dpt<10&&(DsvpvDistance/DsvpvDisErr)>4.42&&Dchi2cl>0.091)||(Dpt>10&&Dpt<15&&(DsvpvDistance/DsvpvDisErr)>4.06&&Dchi2cl>0.069)||(Dpt>15&&Dpt<20&&(DsvpvDistance/DsvpvDisErr)>3.71&&Dchi2cl>0.056)||(Dpt>20&&Dpt<25&&(DsvpvDistance/DsvpvDisErr)>3.25&&Dchi2cl>0.054)||(Dpt>25&&(DsvpvDistance/DsvpvDisErr)>2.97&&Dchi2cl>0.050)))";

  TCut cutmc = "(Dgen==23333||Dgen==23344)&&pthat>(1.5*Dpt)";
  TCut cutmcSignal = "Dgen==23333&&pthat>(1.5*Dpt)";
  TCut cutmcSwapped = "Dgen==23344&&pthat>(1.5*Dpt)";
  TCut cutpt = "Dpt>20";
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
  TH3D* hPtMD0DcaMCPSignal = new TH3D("hPtMD0DcaMCPSignal",";p_{T} (GeV/c);m (GeV/c^{2});D^{0} DCA (cm)",nBinX,binsX,nBinM,binsM,nBinY,binsY);
  TH3D* hPtMD0DcaMCPSwapped = new TH3D("hPtMD0DcaMCPSwapped",";p_{T} (GeV/c);m (GeV/c^{2});D^{0} DCA (cm)",nBinX,binsX,nBinM,binsM,nBinY,binsY);
    
  hMCPSignal->Sumw2();
  hMCNPSignal->Sumw2();
  hPtMD0DcaMCPSignal->Sumw2();
  hPtMD0DcaMCPSwapped->Sumw2();
  
  tMCP->Draw("DsvpvDistance/DsvpvDisErr:DsvpvDistance*sin(Dalpha):Dpt>>hMCPSignal",weightfunctionreco*(cutpt&&cutPbPb&&cutSignal&&cutmc));
  tMCNP->Draw("DsvpvDistance/DsvpvDisErr:DsvpvDistance*sin(Dalpha):Dpt>>hMCNPSignal",weightfunctionreco*(cutpt&&cutPbPb&&cutSignal&&cutmc));
  tMCP->Draw("DsvpvDistance*sin(Dalpha):Dmass:Dpt>>hPtMD0DcaMCPSignal",weightfunctionreco*(cutpt&&cutPbPb&&cutmcSignal));
  tMCP->Draw("DsvpvDistance*sin(Dalpha):Dmass:Dpt>>hPtMD0DcaMCPSwapped",weightfunctionreco*(cutpt&&cutPbPb&&cutmcSwapped));

  outf->Write();
}

