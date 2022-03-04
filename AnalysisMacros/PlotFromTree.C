#include <TFile.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TF1.h>
#include <TProfile.h>
#include <TTree.h>
#include <TCanvas.h>

const int maxFiles=5;
TH1F* hFallTime[maxFiles];
TH1F* hAmpl[maxFiles];
TH1F* hRecoTime[maxFiles];
TH2F* hFallTimeTrigChanVsAmplTrigChan[maxFiles];
TH2F* hFallTimeTrigChanVsAmplTrigChanCluSiz1[maxFiles];
TH2F* hFallTimeTrigChanVsAmplTrigChanCluSizGt1[maxFiles];
TProfile* profFallTime[maxFiles];
TProfile* profFallTimeCluSiz1[maxFiles];
TProfile* profFallTimeCluSizGt1[maxFiles];
TH1F* hTotAmpl[maxFiles];
TH1F* hCluSiz[maxFiles];
TH1F* hCluTyp[maxFiles];
TH2F* hAmplTrigChanVsCluSiz[maxFiles];
TH1F* hFallTimeTrigChanCluSiz1[maxFiles];
TProfile* profAmplTrigChan[maxFiles];
TH1F* hAmplCh1FastResp[maxFiles];
TH1F* hAmplCh2FastResp[maxFiles];
TH1F* hAmplCh3FastResp[maxFiles];
TH1F* hAmplCh4FastResp[maxFiles];
TH1F* hCluSizFastResp[maxFiles];
TH1F* hDeltaTime12[maxFiles];
TH1F* hDeltaTime13[maxFiles];
TH1F* hDeltaTime14[maxFiles];
TH1F* hDeltaTime21[maxFiles];
TH1F* hDeltaTime31[maxFiles];
TH1F* hDeltaTime41[maxFiles];


void FormatHistos(TObjArray* arrHisto, int theColor){
  int nHist=arrHisto->GetEntries();
  for(int jh=0; jh<nHist; jh++){
    TString cname=((TObject*)arrHisto->At(jh))->ClassName();
    if(cname.Contains("TH1")){
      TH1* h=(TH1*)arrHisto->At(jh);
      h->SetLineWidth(2);
      h->SetLineColor(theColor);
      h->GetYaxis()->SetTitleOffset(1.25);
    }else if(cname.Contains("TH2")){
      TH2* h=(TH2*)arrHisto->At(jh);
      h->SetLineColor(theColor);
    }    
  }
}

void FillHistosFromTree(TFile* f, int jfil, int iTrigChan){
  printf("Fill histos wit iTrigChan = %d\n",iTrigChan);
  TTree* tree=(TTree*)f->Get("treeParams");
  int ev;
  double t10Vec[4];
  double t50Vec[4];
  double t90Vec[4];
  double amplVec[4];
  double baselVec[4];
  double recoTimeVec[4];
  double fallTimeVec[4];
  tree->SetBranchAddress("Event",&ev);
  for(int k=0; k<4; k++){
    tree->SetBranchAddress(Form("t10Ch%d",k+1),&t10Vec[k]);
    tree->SetBranchAddress(Form("t50Ch%d",k+1),&t50Vec[k]);
    tree->SetBranchAddress(Form("t90Ch%d",k+1),&t90Vec[k]);
    tree->SetBranchAddress(Form("FallTimeCh%d",k+1),&fallTimeVec[k]);
    tree->SetBranchAddress(Form("RecovTimeExpoCh%d",k+1),&recoTimeVec[k]);
    tree->SetBranchAddress(Form("SignalAmplCh%d",k+1),&amplVec[k]);
  }

  for(int ient=0; ient<tree->GetEntriesFast(); ient++){
    tree->GetEvent(ient);
    if(amplVec[iTrigChan]>=0){
      hFallTime[jfil]->Fill(fallTimeVec[iTrigChan]/1000.);
      hAmpl[jfil]->Fill(amplVec[iTrigChan]*1000.);
      hRecoTime[jfil]->Fill(recoTimeVec[iTrigChan]/1e6);
      hFallTimeTrigChanVsAmplTrigChan[jfil]->Fill(amplVec[iTrigChan]*1000.,fallTimeVec[iTrigChan]/1000.);
    }
    int clusiz=0;
    int clutyp=0;
    double totampl=0;
    for(int k=0; k<4; k++){
      if(amplVec[k]>0){
  	clutyp+=1<<k;
  	clusiz+=1;
  	totampl+=amplVec[k]*1000;
      }
    }
    
    if(amplVec[iTrigChan]>0){
      hCluSiz[jfil]->Fill(clusiz);
      hCluTyp[jfil]->Fill(clutyp);
      hTotAmpl[jfil]->Fill(totampl);
      hAmplTrigChanVsCluSiz[jfil]->Fill(clusiz,amplVec[iTrigChan]*1000);
    }
    if(amplVec[iTrigChan]>0 && clusiz==1){
      hFallTimeTrigChanCluSiz1[jfil]->Fill(fallTimeVec[iTrigChan]/1000.);
      hFallTimeTrigChanVsAmplTrigChanCluSiz1[jfil]->Fill(amplVec[iTrigChan]*1000.,fallTimeVec[iTrigChan]/1000.);
    }
    if(amplVec[iTrigChan]>0 && clusiz>1) hFallTimeTrigChanVsAmplTrigChanCluSizGt1[jfil]->Fill(amplVec[iTrigChan]*1000.,fallTimeVec[iTrigChan]/1000.);
    if(amplVec[iTrigChan]>0 && fallTimeVec[iTrigChan]<400.){
      hCluSizFastResp[jfil]->Fill(clusiz);
      hAmplCh1FastResp[jfil]->Fill(TMath::Max(0.,amplVec[0]*1000.));
      hAmplCh2FastResp[jfil]->Fill(TMath::Max(0.,amplVec[1]*1000.));
      hAmplCh3FastResp[jfil]->Fill(TMath::Max(0.,amplVec[2]*1000.));
      hAmplCh4FastResp[jfil]->Fill(TMath::Max(0.,amplVec[3]*1000.));
    }
    if(amplVec[iTrigChan]>0){
      if(amplVec[1]>0 && amplVec[1]<amplVec[0]) hDeltaTime12[jfil]->Fill((t50Vec[1]-t50Vec[0])/1000.);
      else if(amplVec[1]>0 && amplVec[1]>=amplVec[0]) hDeltaTime21[jfil]->Fill((t50Vec[1]-t50Vec[0])/1000.);
      if(amplVec[2]>0 && amplVec[2]<amplVec[0]) hDeltaTime13[jfil]->Fill((t50Vec[2]-t50Vec[0])/1000.);
      else if(amplVec[2]>0 && amplVec[2]>=amplVec[0]) hDeltaTime31[jfil]->Fill((t50Vec[2]-t50Vec[0])/1000.);
      if(amplVec[3]>0 && amplVec[3]<amplVec[0]) hDeltaTime14[jfil]->Fill((t50Vec[3]-t50Vec[0])/1000.);
      else if(amplVec[3]>0 && amplVec[3]>=amplVec[0]) hDeltaTime41[jfil]->Fill((t50Vec[3]-t50Vec[0])/1000.);
    }
  }
  for(int k=1; k<=hCluSiz[jfil]->GetNbinsX(); k++){
    hCluSiz[jfil]->SetBinContent(k, hCluSiz[jfil]->GetBinContent(k)/hCluSiz[jfil]->GetEntries());
  }
  hCluSizFastResp[jfil]->Scale(1/hAmplCh1FastResp[jfil]->GetEntries());
  for(int k=1; k<=hCluTyp[jfil]->GetNbinsX(); k++){
    hCluTyp[jfil]->SetBinContent(k, hCluTyp[jfil]->GetBinContent(k)/hCluTyp[jfil]->GetEntries());
  }
  delete tree;
  return;
}


void PlotFromTree(TString configFile="configuration.txt", bool normalizeToArea=kTRUE){

  int nFiles=0;
  TString fileNames[maxFiles];
  TString trigChan[maxFiles];
  int jTrigChan[maxFiles];
  int cols[maxFiles]={kMagenta+1,1};
  TString legEntry[maxFiles];
  FILE* cFile=fopen(configFile.Data(),"r");
  char suff[100];
  char line[200];
  fgets(line,200,cFile);
  sscanf(line,"%d %s",&nFiles,suff);
  if(nFiles>maxFiles){
    printf("ERROR: maximum number of files is %d\n",maxFiles);
    return;
  }
  TString suffix=suff;
  int readFiles=0;
  for(int jf=0; jf<nFiles; jf++){
    fgets(line,200,cFile);
    TString theLine(line);
    TObjArray* arrEnt=theLine.Tokenize(";");
    int nEnt=arrEnt->GetEntries();
    if(nEnt!=4){
      printf("ERROR: expect filename ; trigchan ; color ; legendtext\n");
      break;
    }
    for(int k=0; k<nEnt; k++){
      TObjString* str=(TObjString*)arrEnt->At(k);
      TString theStr=str->GetString();
      theStr.ReplaceAll("\n","");
      if(k==0) fileNames[jf]=theStr.Data();
      if(k==1) trigChan[jf]=theStr.Data();
      else if(k==2) cols[jf]=theStr.Atoi();
      else if(k==3) legEntry[jf]=theStr.Data();
      jTrigChan[jf]=1;
      if(trigChan[jf].Contains("J5")) jTrigChan[jf]=1;
      else if(trigChan[jf].Contains("J10")) jTrigChan[jf]=3;
    }
    //    arrEnt->Delete();
    delete arrEnt;
    readFiles++;
    if(feof(cFile)) break;
  }
  fclose(cFile);
  if(readFiles!=nFiles){
    printf("ERROR: mismatch between number of expected files (%d) and number of read lines (%d)\n",nFiles,readFiles);
    return;
  }
  printf("Number of files to be analyzed = %d suffix for plots = %s\n",nFiles,suffix.Data());
  for(int jf=0; jf<nFiles; jf++){
    printf("File %d = %s  trigger channel = %s(%d)  Color = %d  legend Entry = %s\n",jf, fileNames[jf].Data(),trigChan[jf].Data(),jTrigChan[jf],cols[jf],legEntry[jf].Data());
  }
  if(nFiles==0) return;
 
  TObjArray* arrHisto = new TObjArray();
  double cnt04[4];
  double cntall[4];
  double cntFT04[4];
  double cntFTall[4];
  for(int j=0; j<nFiles; j++){
    hFallTime[j]=new TH1F(Form("hFallTime%d",j)," All clusters ; Fall Time TrigChan (ns) ; Entries",100,0.,10.);
    hAmpl[j]=new TH1F(Form("hAmpl%d",j)," ; Signal Amplitude TrigChan (mV) ; Entries",100,0.,100.);
    hRecoTime[j]=new TH1F(Form("hRecoTime%d",j)," ; Recovery Time (#mus) ; Entries",100,0.,3.);
    hFallTimeTrigChanVsAmplTrigChan[j]=new TH2F(Form("hFallTimeTrigChanVsAmplTrigChan%d",j)," ; Signal Amplitude TrigChan (mV) ; Fall Time TrigChan (ns) ; Entries",50,0.,100.,100,0.,10.);
    hFallTimeTrigChanVsAmplTrigChanCluSiz1[j]=new TH2F(Form("hFallTimeTrigChanVsAmplTrigChanCluSiz1%d",j)," Cluster size = 1  ; Signal Amplitude TrigChan (mV) ; Fall Time TrigChan (ns) ; Entries",50,0.,100.,100,0.,10.);
    hFallTimeTrigChanVsAmplTrigChanCluSiz1[j]->SetTitle(Form("%s Cluster size = 1",legEntry[j].Data()));
    hFallTimeTrigChanVsAmplTrigChanCluSizGt1[j]=new TH2F(Form("hFallTimeTrigChanVsAmplTrigChanCluSizGt1%d",j)," Cluster size > 1  ; Signal Amplitude TrigChan (mV) ; Fall Time TrigChan (ns) ; Entries",50,0.,100.,100,0.,10.);
    hFallTimeTrigChanVsAmplTrigChanCluSizGt1[j]->SetTitle(Form("%s Cluster size > 1",legEntry[j].Data()));
    hTotAmpl[j]=new TH1F(Form("hTotAmpl%d",j)," ; Signal Amplitude 4-pixels (mV) ; Entries",100,0.,100.);
    hCluSiz[j]=new TH1F(Form("hCluSiz%d",j)," ; Cluster Size ; Fraction of events",5,-0.5,4.5);
    hCluTyp[j]=new TH1F(Form("hCluTyp%d",j)," ; Cluster Shape ; Fraction of events",16,-0.5,15.5);
    hAmplTrigChanVsCluSiz[j]=new TH2F(Form("hAmplTrigChanVsCluSiz%d",j)," ; Cluster Size ; Signal Amplitude TrigChan (mV) ; Entries",5,-0.5,4.5,100,0.,100.);    
    hFallTimeTrigChanCluSiz1[j]=new TH1F(Form("hFallTimeTrigChanCluSiz1%d",j)," Cluster size = 1 ; Fall Time TrigChan (ns) ; Entries",100,0.,10.);    
    hAmplCh1FastResp[j]=new TH1F(Form("hAmplCh1FastResp%d",j)," Signals with fall time < 0.4 ns ; Signal Amplitude Ch1 (mV) ; Entries",100,0.,100.);
    hAmplCh2FastResp[j]=new TH1F(Form("hAmplCh2FastResp%d",j)," Signals with fall time < 0.4 ns ; Signal Amplitude Ch2 (mV) ; Entries",100,0.,100.);
    hAmplCh3FastResp[j]=new TH1F(Form("hAmplCh3FastResp%d",j)," Signals with fall time < 0.4 ns ; Signal Amplitude Ch3 (mV) ; Entries",100,0.,100.);
    hAmplCh4FastResp[j]=new TH1F(Form("hAmplCh4FastResp%d",j)," Signals with fall time < 0.4 ns ; Signal Amplitude Ch4 (mV) ; Entries",100,0.,100.);
    hCluSizFastResp[j]=new TH1F(Form("hCluSizFastResp%d",j)," Signals with fall time < 0.4 ns ; Cluster Size ; Fraction of events",5,-0.5,4.5);
    hDeltaTime12[j]=new TH1F(Form("hDeltaTime12%d",j)," Events with Ampl1 > Ampl2 ; #Deltat_{50} (Ch2 - Ch1) (ns) ; Entries",100,-6,6.);
    hDeltaTime13[j]=new TH1F(Form("hDeltaTime13%d",j)," Events with Ampl1 > Ampl3 ; #Deltat_{50} (Ch3 - Ch1) (ns) ; Entries",100,-6,6.);
    hDeltaTime14[j]=new TH1F(Form("hDeltaTime14%d",j)," Events with Ampl1 > Ampl4 ; #Deltat_{50} (Ch4 - Ch1) (ns) ; Entries",100,-6,6.);
    hDeltaTime21[j]=new TH1F(Form("hDeltaTime21%d",j)," Events with Ampl1 < Ampl2 ; #Deltat_{50} (Ch2 - Ch1) (ns) ; Entries",100,-6,6.);
    hDeltaTime31[j]=new TH1F(Form("hDeltaTime31%d",j)," Events with Ampl1 < Ampl3 ; #Deltat_{50} (Ch3 - Ch1) (ns) ; Entries",100,-6,6.);
    hDeltaTime41[j]=new TH1F(Form("hDeltaTime41%d",j)," Events with Ampl1 < Ampl4 ; #Deltat_{50} (Ch4 - Ch1) (ns) ; Entries",100,-6,6.);
    for(int k=0; k<16; k++){
      TString lab="";
      for(int jbit=0; jbit<4;jbit++){	
	if(k&(1<<jbit)){
	  if(lab!="") lab+="+";
	  lab+=jbit+1;
	}
      }
      hCluTyp[j]->GetXaxis()->SetBinLabel(k+1,lab.Data());
    }

    int indexh=0;
    arrHisto->AddAtAndExpand(hFallTime[j],indexh++);
    arrHisto->AddAtAndExpand(hAmpl[j],indexh++);
    arrHisto->AddAtAndExpand(hRecoTime[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeTrigChanVsAmplTrigChan[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeTrigChanVsAmplTrigChanCluSiz1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeTrigChanVsAmplTrigChanCluSizGt1[j],indexh++);
    arrHisto->AddAtAndExpand(hTotAmpl[j],indexh++);
    arrHisto->AddAtAndExpand(hCluSiz[j],indexh++);
    arrHisto->AddAtAndExpand(hCluTyp[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplTrigChanVsCluSiz[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeTrigChanCluSiz1[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh1FastResp[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh2FastResp[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh3FastResp[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh4FastResp[j],indexh++);
    arrHisto->AddAtAndExpand(hCluSizFastResp[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime12[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime13[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime14[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime21[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime31[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime41[j],indexh++);

    TFile* f=new TFile(fileNames[j].Data());
    fileNames[j].ReplaceAll("_TTree.root","");
    FillHistosFromTree(f,j,jTrigChan[j]-1);
    FormatHistos(arrHisto,cols[j]);
    arrHisto->Clear();
    profFallTime[j]=hFallTimeTrigChanVsAmplTrigChan[j]->ProfileX(Form("profFallTimeVsAmpl%d",j));
    profFallTime[j]->GetYaxis()->SetTitle("<Fall Time TrigChan> (ns)");
    hFallTimeTrigChanVsAmplTrigChan[j]->SetStats(0);
    profFallTime[j]->SetLineColor(cols[j]);
    profFallTime[j]->SetLineWidth(2);
    profFallTime[j]->SetStats(0);
    profFallTimeCluSiz1[j]=hFallTimeTrigChanVsAmplTrigChanCluSiz1[j]->ProfileX(Form("profFallTimeVsAmplCluSiz1%d",j));
    profFallTimeCluSiz1[j]->GetYaxis()->SetTitle("<Fall Time TrigChan> (ns)");
    profFallTimeCluSiz1[j]->SetLineColor(cols[j]);
    profFallTimeCluSiz1[j]->SetMarkerColor(cols[j]);
    profFallTimeCluSiz1[j]->SetMarkerStyle(20);
    profFallTimeCluSiz1[j]->SetLineWidth(2);
    profFallTimeCluSiz1[j]->SetStats(0);
    profFallTimeCluSizGt1[j]=hFallTimeTrigChanVsAmplTrigChanCluSizGt1[j]->ProfileX(Form("profFallTimeVsAmplCluSizGt1%d",j));
    profFallTimeCluSizGt1[j]->GetYaxis()->SetTitle("<Fall Time TrigChan> (ns)");
    profFallTimeCluSizGt1[j]->SetLineColor(cols[j]);
    profFallTimeCluSizGt1[j]->SetMarkerColor(cols[j]);
    profFallTimeCluSizGt1[j]->SetMarkerStyle(25);
    profFallTimeCluSizGt1[j]->SetLineWidth(2);
    profFallTimeCluSizGt1[j]->SetStats(0);
    
    cnt04[j]=hFallTimeTrigChanCluSiz1[j]->Integral(1,hFallTimeTrigChanCluSiz1[j]->GetXaxis()->FindBin(0.3999));
    cntall[j]=hFallTimeTrigChanCluSiz1[j]->GetEntries();
    cntFT04[j]=hFallTime[j]->Integral(1,hFallTime[j]->GetXaxis()->FindBin(0.3999));
    cntFTall[j]=hFallTime[j]->GetEntries();
    if(normalizeToArea){
      hTotAmpl[j]->Scale(1./hTotAmpl[j]->Integral());
      hTotAmpl[j]->GetYaxis()->SetTitle("Enrites (a.u.)");
      hAmpl[j]->Scale(1./hAmpl[j]->Integral());
      hAmpl[j]->GetYaxis()->SetTitle("Enrites (a.u.)");
      hFallTime[j]->Scale(1./hFallTime[j]->Integral());
      hFallTime[j]->GetYaxis()->SetTitle("Enrites (a.u.)");
      hFallTimeTrigChanCluSiz1[j]->Scale(1./hFallTimeTrigChanCluSiz1[j]->Integral());
      hFallTimeTrigChanCluSiz1[j]->GetYaxis()->SetTitle("Enrites (a.u.)");
      hRecoTime[j]->Scale(1./hRecoTime[j]->Integral());
      hRecoTime[j]->GetYaxis()->SetTitle("Enrites (a.u.)");
      hDeltaTime12[j]->Scale(1./hDeltaTime12[j]->Integral());
      hDeltaTime12[j]->GetYaxis()->SetTitle("Enrites (a.u.)");
      hDeltaTime13[j]->Scale(1./hDeltaTime13[j]->Integral());
      hDeltaTime13[j]->GetYaxis()->SetTitle("Enrites (a.u.)");
      hDeltaTime14[j]->Scale(1./hDeltaTime14[j]->Integral());
      hDeltaTime14[j]->GetYaxis()->SetTitle("Enrites (a.u.)");
      hDeltaTime21[j]->Scale(1./hDeltaTime21[j]->Integral());
      hDeltaTime21[j]->GetYaxis()->SetTitle("Enrites (a.u.)");
      hDeltaTime31[j]->Scale(1./hDeltaTime31[j]->Integral());
      hDeltaTime31[j]->GetYaxis()->SetTitle("Enrites (a.u.)");
      hDeltaTime41[j]->Scale(1./hDeltaTime41[j]->Integral());
      hDeltaTime41[j]->GetYaxis()->SetTitle("Enrites (a.u.)");
    }

    if(hFallTime[j]->GetMaximum()>hFallTime[0]->GetMaximum()) hFallTime[0]->SetMaximum(hFallTime[j]->GetMaximum()*1.05);
    if(hFallTimeTrigChanCluSiz1[j]->GetMaximum()>hFallTimeTrigChanCluSiz1[0]->GetMaximum()) hFallTimeTrigChanCluSiz1[0]->SetMaximum(hFallTimeTrigChanCluSiz1[j]->GetMaximum()*1.05);
    if(hRecoTime[j]->GetMaximum()>hRecoTime[0]->GetMaximum()) hRecoTime[0]->SetMaximum(hRecoTime[j]->GetMaximum()*1.05);
    if(hAmpl[j]->GetMaximum()>hAmpl[0]->GetMaximum()) hAmpl[0]->SetMaximum(hAmpl[j]->GetMaximum()*1.05);
    if(hAmplCh1FastResp[j]->GetMaximum()>hAmplCh1FastResp[0]->GetMaximum()) hAmplCh1FastResp[0]->SetMaximum(hAmplCh1FastResp[j]->GetMaximum()*1.05);
    if(hAmplCh2FastResp[j]->GetMaximum()>hAmplCh2FastResp[0]->GetMaximum()) hAmplCh2FastResp[0]->SetMaximum(hAmplCh2FastResp[j]->GetMaximum()*1.05);
    if(hAmplCh3FastResp[j]->GetMaximum()>hAmplCh3FastResp[0]->GetMaximum()) hAmplCh3FastResp[0]->SetMaximum(hAmplCh3FastResp[j]->GetMaximum()*1.05);
    if(hAmplCh4FastResp[j]->GetMaximum()>hAmplCh4FastResp[0]->GetMaximum()) hAmplCh4FastResp[0]->SetMaximum(hAmplCh4FastResp[j]->GetMaximum()*1.05);
    if(hDeltaTime12[j]->GetMaximum()>hDeltaTime12[0]->GetMaximum()) hDeltaTime12[0]->SetMaximum(hDeltaTime12[j]->GetMaximum()*1.05);
    if(hDeltaTime13[j]->GetMaximum()>hDeltaTime13[0]->GetMaximum()) hDeltaTime13[0]->SetMaximum(hDeltaTime13[j]->GetMaximum()*1.05);
    if(hDeltaTime14[j]->GetMaximum()>hDeltaTime14[0]->GetMaximum()) hDeltaTime14[0]->SetMaximum(hDeltaTime14[j]->GetMaximum()*1.05);
    if(hDeltaTime21[j]->GetMaximum()>hDeltaTime21[0]->GetMaximum()) hDeltaTime21[0]->SetMaximum(hDeltaTime21[j]->GetMaximum()*1.05);
    if(hDeltaTime31[j]->GetMaximum()>hDeltaTime31[0]->GetMaximum()) hDeltaTime31[0]->SetMaximum(hDeltaTime31[j]->GetMaximum()*1.05);
    if(hDeltaTime41[j]->GetMaximum()>hDeltaTime41[0]->GetMaximum()) hDeltaTime41[0]->SetMaximum(hDeltaTime41[j]->GetMaximum()*1.05);
    if(profFallTime[j]->GetMaximum()>profFallTime[0]->GetMaximum()*0.9) profFallTime[0]->SetMaximum(profFallTime[j]->GetMaximum()*1.2);
    if(profFallTimeCluSiz1[j]->GetMaximum()>profFallTimeCluSiz1[0]->GetMaximum()*0.9) profFallTimeCluSiz1[0]->SetMaximum(profFallTimeCluSiz1[j]->GetMaximum()*1.2);
    if(profFallTimeCluSizGt1[j]->GetMaximum()>profFallTimeCluSizGt1[0]->GetMaximum()*0.9) profFallTimeCluSizGt1[0]->SetMaximum(profFallTimeCluSizGt1[j]->GetMaximum()*1.2);
  
    profAmplTrigChan[j]=hAmplTrigChanVsCluSiz[j]->ProfileX(Form("profAmplTrigChan%d",j));
    profAmplTrigChan[j]->GetYaxis()->SetTitle("<Amplitude TrigChan> (mV)");
    hCluTyp[j]->SetStats(0);
    hAmplTrigChanVsCluSiz[j]->SetStats(0);
    profAmplTrigChan[j]->SetLineColor(cols[j]);
    profAmplTrigChan[j]->SetLineWidth(2);
    profAmplTrigChan[j]->SetStats(0);
    if(hTotAmpl[j]->GetMaximum()>hTotAmpl[0]->GetMaximum()) hTotAmpl[0]->SetMaximum(hTotAmpl[j]->GetMaximum()*1.05);
    if(hCluTyp[j]->GetMaximum()>hCluTyp[0]->GetMaximum()) hCluTyp[0]->SetMaximum(hCluTyp[j]->GetMaximum()*1.05);
    if(hCluSiz[j]->GetMaximum()>hCluSiz[0]->GetMaximum()) hCluSiz[0]->SetMaximum(hCluSiz[j]->GetMaximum()*1.05);
    if(hCluSizFastResp[j]->GetMaximum()>hCluSizFastResp[0]->GetMaximum()) hCluSizFastResp[0]->SetMaximum(hCluSizFastResp[j]->GetMaximum()*1.05);
    if(profAmplTrigChan[j]->GetMaximum()>profAmplTrigChan[0]->GetMaximum()*0.9) profAmplTrigChan[0]->SetMaximum(profAmplTrigChan[j]->GetMaximum()*1.1);

  }

  TCanvas* c1 = new TCanvas("c1","",1650,900);
  c1->Divide(3,2);
  c1->cd(1);
  TLegend* leg= new TLegend(0.12,0.6,0.5,0.89);
  leg->SetMargin(0.1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime[j]->Draw("histo");
    else hFallTime[j]->Draw("histo,sames");
    leg->AddEntry(hFallTime[j],legEntry[j].Data(),"L")->SetTextColor(cols[j]);
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime[j]->GetLineColor());
    gPad->Modified();
  } 
  c1->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hRecoTime[j]->Draw("histo");
    else hRecoTime[j]->Draw("histo,sames");
    gPad->Update();
    TPaveStats *st=(TPaveStats*)hRecoTime[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hRecoTime[j]->GetLineColor());
    gPad->Modified();
  }
  c1->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmpl[j]->Draw("histo");
    else hAmpl[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmpl[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmpl[j]->GetLineColor());
    gPad->Modified();
  }
  c1->cd(4);
  gPad->SetLogz();
  hFallTimeTrigChanVsAmplTrigChan[0]->Draw("box");
  for(int j=1; j<nFiles; j++) hFallTimeTrigChanVsAmplTrigChan[j]->Draw("same,box");
  c1->cd(5);
  profFallTime[0]->Draw("");
  for(int j=1; j<nFiles; j++) profFallTime[j]->Draw("same");
  c1->cd(6);
  for(int j=0; j<nFiles; j++){
    if(j==0) hTotAmpl[j]->Draw("histo");
    else hTotAmpl[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hTotAmpl[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hTotAmpl[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
  c1->SaveAs(Form("Histos_%s.png",suffix.Data()));
  
  TCanvas* c2 = new TCanvas("c2","",1650,900);
  c2->Divide(3,2);
  c2->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hCluSiz[j]->Draw();
    else hCluSiz[j]->Draw("sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hCluSiz[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hCluSiz[j]->GetLineColor());
    gPad->Modified();
  }
  c2->cd(2);
  hCluTyp[0]->Draw();
  for(int j=1; j<nFiles; j++) hCluTyp[j]->Draw("same");
  c2->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTime[j]->Draw("histo");
    else hFallTime[j]->Draw("histo,sames");
    TLatex* tfrach=new TLatex(0.5,0.72-0.08*j,Form("%.0f/%.0f=%.3f",cntFT04[j],cntFTall[j],cntFT04[j]/cntFTall[j]));
    tfrach->SetNDC();
    tfrach->SetTextColor(cols[j]);
    tfrach->SetTextFont(43);
    tfrach->SetTextSize(20);
    tfrach->Draw();
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTime[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTime[j]->GetLineColor());
    gPad->Modified();
  } 
  TLatex* tfrac=new TLatex(0.48,0.8,"Fall time < 0.4 ns");
  tfrac->SetNDC();
  tfrac->SetTextFont(43);
  tfrac->SetTextSize(20);
  tfrac->Draw();
  c2->cd(4);
  hAmplTrigChanVsCluSiz[0]->Draw("box");
  for(int j=1; j<nFiles; j++) hAmplTrigChanVsCluSiz[j]->Draw("same,box");
  c2->cd(5);
  profAmplTrigChan[0]->Draw("");
  for(int j=1; j<nFiles; j++) profAmplTrigChan[j]->Draw("same");
  c2->cd(6);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeTrigChanCluSiz1[j]->Draw("histo");
    else hFallTimeTrigChanCluSiz1[j]->Draw("histo,sames");
    TLatex* tfrach=new TLatex(0.5,0.72-0.08*j,Form("%.0f/%.0f=%.3f",cnt04[j],cntall[j],cnt04[j]/cntall[j]));
    tfrach->SetNDC();
    tfrach->SetTextColor(cols[j]);
    tfrach->SetTextFont(43);
    tfrach->SetTextSize(20);
    tfrach->Draw();
    gPad->Update();
    TPaveStats *st=(TPaveStats*)hFallTimeTrigChanCluSiz1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hTotAmpl[j]->GetLineColor());
    gPad->Modified();
  }
  tfrac->Draw();
  leg->Draw();
  c2->SaveAs(Form("CluSiz_%s.png",suffix.Data()));
  
  TCanvas* c3 = new TCanvas("c3","",1650,900);
  c3->Divide(3,2);
  c3->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hCluSizFastResp[j]->Draw("histo");
    else hCluSizFastResp[j]->Draw("histosames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hCluSizFastResp[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hCluSizFastResp[j]->GetLineColor());
    gPad->Modified();
  }
  c3->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplCh1FastResp[j]->Draw();
    else hAmplCh1FastResp[j]->Draw("sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplCh1FastResp[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplCh1FastResp[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
  c3->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplCh2FastResp[j]->Draw();
    else hAmplCh2FastResp[j]->Draw("sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplCh2FastResp[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplCh2FastResp[j]->GetLineColor());
    gPad->Modified();
  }
  c3->cd(5);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplCh3FastResp[j]->Draw();
    else hAmplCh3FastResp[j]->Draw("sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplCh3FastResp[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplCh3FastResp[j]->GetLineColor());
    gPad->Modified();
  }
  c3->cd(6);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplCh4FastResp[j]->Draw();
    else hAmplCh4FastResp[j]->Draw("sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplCh4FastResp[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplCh4FastResp[j]->GetLineColor());
    gPad->Modified();
  }
  c3->cd(4);
  for(int j=0; j<nFiles; j++){
    TLatex* t1=new TLatex(0.05,0.8-0.25*j,fileNames[j].Data());
    t1->SetTextFont(43);
    t1->SetTextSize(18);
    t1->SetNDC();
    t1->SetTextColor(cols[j]);
    t1->Draw();
    TLatex* t2=new TLatex(0.05,0.72-0.25*j,Form("Fraction with fast response = %.0f / %.0f = %.3f",hAmplCh1FastResp[j]->GetEntries(),hAmpl[j]->GetEntries(),hAmplCh1FastResp[j]->GetEntries()/hAmpl[j]->GetEntries()));
    t2->SetTextFont(43);
    t2->SetTextSize(18);
    t2->SetNDC();
    t2->SetTextColor(cols[j]);
    t2->Draw();
  }
  c3->SaveAs(Form("FastResp_%s.png",suffix.Data()));
  
  TCanvas* c4 = new TCanvas("c4","",1650,900);
  c4->Divide(3,2);
  c4->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime12[j]->Draw("histo");
    else hDeltaTime12[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hDeltaTime12[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime12[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
  c4->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime13[j]->Draw("histo");
    else hDeltaTime13[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hDeltaTime13[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime13[j]->GetLineColor());
    gPad->Modified();
  }
  c4->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime14[j]->Draw("histo");
    else hDeltaTime14[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hDeltaTime14[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime14[j]->GetLineColor());
    gPad->Modified();
  }
  c4->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime21[j]->Draw("histo");
    else hDeltaTime21[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hDeltaTime21[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime21[j]->GetLineColor());
    gPad->Modified();
  }
  c4->cd(5);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime31[j]->Draw("histo");
    else hDeltaTime31[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hDeltaTime31[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime31[j]->GetLineColor());
    gPad->Modified();
  }
  c4->cd(6);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime41[j]->Draw("histo");
    else hDeltaTime41[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hDeltaTime41[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime41[j]->GetLineColor());
    gPad->Modified();
  }
  // c4->cd(8);
  // for(int j=0; j<nFiles; j++){
  //   TLatex* t1=new TLatex(0.05,0.8-0.1*j,fileNames[j].Data());
  //   t1->SetTextFont(43);
  //   t1->SetTextSize(18);
  //   t1->SetNDC();
  //   t1->SetTextColor(cols[j]);
  //   t1->Draw();
  // }


  c4->SaveAs(Form("DeltaTim_%s.png",suffix.Data()));

  TCanvas* c5=new TCanvas("c5","",1650,900);
  c5->Divide(nFiles+1,2);
  for(int j=0; j<nFiles; j++){
    c5->cd(j+1);
    hFallTimeTrigChanVsAmplTrigChanCluSiz1[j]->Draw("colz");
    c5->cd(j+2+nFiles);
    hFallTimeTrigChanVsAmplTrigChanCluSizGt1[j]->Draw("colz");
  }
  c5->cd(nFiles+1);
  profFallTimeCluSiz1[0]->SetTitle("");
  profFallTimeCluSiz1[0]->Draw("");
  for(int j=1; j<nFiles; j++) profFallTimeCluSiz1[j]->Draw("same");
  for(int j=0; j<nFiles; j++) profFallTimeCluSizGt1[j]->Draw("same");
  TLegend* lpr=new TLegend(0.5,0.69,0.89,0.89);
  for(int j=0; j<nFiles; j++){
    lpr->AddEntry(profFallTimeCluSiz1[j],Form("%s cluster size = 1 ",legEntry[j].Data()),"P");
    lpr->AddEntry(profFallTimeCluSizGt1[j],Form("%s cluster size > 1 ",legEntry[j].Data()),"P");
  }
  lpr->Draw();
  c5->SaveAs(Form("FallTimeVsAmpl_%s.png",suffix.Data()));
}
