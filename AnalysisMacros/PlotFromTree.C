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
TH1F* hFallTimeTrigChan[maxFiles];
TH1F* hFallTimeCh1[maxFiles];
TH1F* hFallTimeCh2[maxFiles];
TH1F* hFallTimeCh3[maxFiles];
TH1F* hFallTimeCh4[maxFiles];
TH1F* hFallTimeMaxSigPix[maxFiles];
TH1F* hFallTimeMaxSigPixCluSiz1[maxFiles];
TH1F* hAmplTrigChan[maxFiles];
TH1F* hAmplCh1[maxFiles];
TH1F* hAmplCh2[maxFiles];
TH1F* hAmplCh3[maxFiles];
TH1F* hAmplCh4[maxFiles];
TH1F* hRecoTimeTrigChan[maxFiles];
TH1F* hRecoTimeCh1[maxFiles];
TH1F* hRecoTimeCh2[maxFiles];
TH1F* hRecoTimeCh3[maxFiles];
TH1F* hRecoTimeCh4[maxFiles];
TH1F* hRecoTimeMaxSigPix[maxFiles];
TH1F* hRecoTimeMaxSigPixCluSiz1[maxFiles];
TH2F* hFallTimeTrigChanVsAmplTrigChan[maxFiles];
TH2F* hFallTimeTrigChanVsAmplTrigChanCluSiz1[maxFiles];
TH2F* hFallTimeTrigChanVsAmplTrigChanCluSizGt1[maxFiles];
TH2F* hFallTimeMaxSigPixVsAmplMaxSigPix[maxFiles];
TH2F* hFallTimeMaxSigPixVsAmplMaxSigPixCluSiz1[maxFiles];
TH2F* hFallTimeMaxSigPixVsAmplMaxSigPixCluSizGt1[maxFiles];
TProfile* profFallTime[maxFiles];
TProfile* profFallTimeCluSiz1[maxFiles];
TProfile* profFallTimeCluSizGt1[maxFiles];
TH1F* hTotAmpl[maxFiles];
TH1F* hCluSiz[maxFiles];
TH1F* hCluTyp[maxFiles];
TH2F* hAmplTrigChanVsCluSiz[maxFiles];
TH2F* hAmplMaxSigPixVsCluSiz[maxFiles];
TH1F* hFallTimeTrigChanCluSiz1[maxFiles];
TProfile* profAmplTrigChan[maxFiles];
TProfile* profAmplMaxSigPix[maxFiles];
TH1F* hAmplCh1FastResp[maxFiles];
TH1F* hAmplCh2FastResp[maxFiles];
TH1F* hAmplCh3FastResp[maxFiles];
TH1F* hAmplCh4FastResp[maxFiles];
TH1F* hCluSizFastResp[maxFiles];
TH1F* hDeltaTime12[maxFiles];
TH1F* hDeltaTime13[maxFiles];
TH1F* hDeltaTime14[maxFiles];
TH1F* hDeltaTime23[maxFiles];
TH1F* hDeltaTime24[maxFiles];
TH1F* hDeltaTime34[maxFiles];
TH1F* hDeltaTime21[maxFiles];
TH1F* hDeltaTime31[maxFiles];
TH1F* hDeltaTime41[maxFiles];
TH1F* hDeltaTime32[maxFiles];
TH1F* hDeltaTime42[maxFiles];
TH1F* hDeltaTime43[maxFiles];


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

void NormalizeHistos(TObjArray* arrHisto){
  int nHist=arrHisto->GetEntries();
  for(int jh=0; jh<nHist; jh++){
    TString cname=((TObject*)arrHisto->At(jh))->ClassName();
    if(cname.Contains("TH1")){
      TH1* h=(TH1*)arrHisto->At(jh);
      TString hname=h->GetName();
      if(hname.Contains("CluSiz") || hname.Contains("CluTyp") ) continue;
      double tot=h->Integral();
      if(tot>0){
	h->Scale(1./tot);
	h->GetYaxis()->SetTitle("Entries (a.u.)");
      }
    }
  }
}


void SetHistoMaximum(TH1* hcur, TH1* href, double scal=1.05){
  if(hcur->GetMaximum() > href->GetMaximum()) href->SetMaximum(hcur->GetMaximum()*scal);
}

void SetProfMaximum(TProfile* pcur, TProfile* pref, double scal1, double scal2){
  if(pcur->GetMaximum() > pref->GetMaximum()*scal1) pref->SetMaximum(pcur->GetMaximum()*scal2);
}

void SetMaxima(int nfils){

  for(int j=0; j<nfils; j++){
    SetHistoMaximum(hFallTimeTrigChan[j],hFallTimeTrigChan[0]);
    SetHistoMaximum(hFallTimeCh1[j],hFallTimeCh1[0]);
    SetHistoMaximum(hFallTimeCh2[j],hFallTimeCh2[0]);
    SetHistoMaximum(hFallTimeCh3[j],hFallTimeCh3[0]);
    SetHistoMaximum(hFallTimeCh4[j],hFallTimeCh4[0]);
    SetHistoMaximum(hFallTimeMaxSigPix[j],hFallTimeMaxSigPix[0]);
    SetHistoMaximum(hFallTimeMaxSigPixCluSiz1[j],hFallTimeMaxSigPixCluSiz1[0]);    
    SetHistoMaximum(hFallTimeTrigChanCluSiz1[j],hFallTimeTrigChanCluSiz1[0]);
    SetHistoMaximum(hRecoTimeTrigChan[j],hRecoTimeTrigChan[0]);
    SetHistoMaximum(hRecoTimeCh1[j],hRecoTimeCh1[0]);
    SetHistoMaximum(hRecoTimeCh2[j],hRecoTimeCh2[0]);
    SetHistoMaximum(hRecoTimeCh3[j],hRecoTimeCh3[0]);
    SetHistoMaximum(hRecoTimeCh4[j],hRecoTimeCh4[0]);
    SetHistoMaximum(hRecoTimeMaxSigPix[j],hRecoTimeMaxSigPix[0]);
    SetHistoMaximum(hRecoTimeMaxSigPixCluSiz1[j],hRecoTimeMaxSigPixCluSiz1[0]);    
    SetHistoMaximum(hAmplTrigChan[j],hAmplTrigChan[0]);
    SetHistoMaximum(hAmplCh1[j],hAmplCh1[0]);
    SetHistoMaximum(hAmplCh2[j],hAmplCh2[0]);
    SetHistoMaximum(hAmplCh3[j],hAmplCh3[0]);
    SetHistoMaximum(hAmplCh4[j],hAmplCh4[0]);
    SetHistoMaximum(hAmplCh1FastResp[j],hAmplCh1FastResp[0]);
    SetHistoMaximum(hAmplCh2FastResp[j],hAmplCh2FastResp[0]);
    SetHistoMaximum(hAmplCh3FastResp[j],hAmplCh3FastResp[0]);
    SetHistoMaximum(hAmplCh4FastResp[j],hAmplCh4FastResp[0]);
    SetHistoMaximum(hDeltaTime12[j],hDeltaTime12[0]);
    SetHistoMaximum(hDeltaTime13[j],hDeltaTime13[0]);
    SetHistoMaximum(hDeltaTime14[j],hDeltaTime14[0]);
    SetHistoMaximum(hDeltaTime23[j],hDeltaTime23[0]);
    SetHistoMaximum(hDeltaTime24[j],hDeltaTime24[0]);
    SetHistoMaximum(hDeltaTime34[j],hDeltaTime34[0]);
    SetHistoMaximum(hDeltaTime21[j],hDeltaTime21[0]);
    SetHistoMaximum(hDeltaTime31[j],hDeltaTime31[0]);
    SetHistoMaximum(hDeltaTime41[j],hDeltaTime41[0]);
    SetHistoMaximum(hDeltaTime32[j],hDeltaTime32[0]);
    SetHistoMaximum(hDeltaTime42[j],hDeltaTime42[0]);
    SetHistoMaximum(hDeltaTime43[j],hDeltaTime43[0]);
    SetHistoMaximum(hTotAmpl[j],hTotAmpl[0]);
    SetHistoMaximum(hCluTyp[j],hCluTyp[0]);
    SetHistoMaximum(hCluSiz[j],hCluSiz[0]);
    SetHistoMaximum(hCluSizFastResp[j],hCluSizFastResp[0]);
    SetProfMaximum(profFallTime[j],profFallTime[0],0.9,1.2);
    SetProfMaximum(profFallTimeCluSiz1[j],profFallTimeCluSiz1[0],0.9,1.2);
    SetProfMaximum(profFallTimeCluSizGt1[j],profFallTimeCluSizGt1[0],0.9,1.2);
    SetProfMaximum(profAmplTrigChan[j],profAmplTrigChan[0],0.9,1.1);
    SetProfMaximum(profAmplMaxSigPix[j],profAmplMaxSigPix[0],0.9,1.1);
  }

}

void FillHistosFromTree(TFile* f, int jfil, int iTrigChan){
  printf("Fill histos with iTrigChan = %d\n",iTrigChan);
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
    if(iTrigChan>=0 && amplVec[iTrigChan]>=0){
      hFallTimeTrigChan[jfil]->Fill(fallTimeVec[iTrigChan]/1000.);
      hFallTimeCh1[jfil]->Fill(fallTimeVec[0]/1000.);
      hFallTimeCh2[jfil]->Fill(fallTimeVec[1]/1000.);
      hFallTimeCh3[jfil]->Fill(fallTimeVec[2]/1000.);
      hFallTimeCh4[jfil]->Fill(fallTimeVec[3]/1000.);
      hAmplTrigChan[jfil]->Fill(amplVec[iTrigChan]*1000.);
      hAmplCh1[jfil]->Fill(amplVec[0]*1000.);
      hAmplCh2[jfil]->Fill(amplVec[1]*1000.);
      hAmplCh3[jfil]->Fill(amplVec[2]*1000.);
      hAmplCh4[jfil]->Fill(amplVec[3]*1000.);
      hRecoTimeTrigChan[jfil]->Fill(recoTimeVec[iTrigChan]/1e6);
      hRecoTimeCh1[jfil]->Fill(recoTimeVec[0]/1e6);
      hRecoTimeCh2[jfil]->Fill(recoTimeVec[1]/1e6);
      hRecoTimeCh3[jfil]->Fill(recoTimeVec[2]/1e6);
      hRecoTimeCh4[jfil]->Fill(recoTimeVec[3]/1e6);
      hFallTimeTrigChanVsAmplTrigChan[jfil]->Fill(amplVec[iTrigChan]*1000.,fallTimeVec[iTrigChan]/1000.);
    }else if(iTrigChan<0){
      if(amplVec[0]>=0 && fallTimeVec[0]>=0){
	hFallTimeCh1[jfil]->Fill(fallTimeVec[0]/1000.);
	hAmplCh1[jfil]->Fill(amplVec[0]*1000.);
	hRecoTimeCh1[jfil]->Fill(recoTimeVec[0]/1e6);
      }
      if(amplVec[1]>=0 && fallTimeVec[1]>=0){
	hFallTimeCh2[jfil]->Fill(fallTimeVec[1]/1000.);
	hAmplCh2[jfil]->Fill(amplVec[1]*1000.);
	hRecoTimeCh2[jfil]->Fill(recoTimeVec[1]/1e6);
      }     
      if(amplVec[2]>=0 && fallTimeVec[2]>=0){
	hFallTimeCh3[jfil]->Fill(fallTimeVec[2]/1000.);
	hAmplCh3[jfil]->Fill(amplVec[2]*1000.);
	hRecoTimeCh3[jfil]->Fill(recoTimeVec[2]/1e6);
      }     
      if(amplVec[3]>=0 && fallTimeVec[3]>=0){
	hFallTimeCh4[jfil]->Fill(fallTimeVec[3]/1000.);
	hAmplCh4[jfil]->Fill(amplVec[3]*1000.);
	hRecoTimeCh4[jfil]->Fill(recoTimeVec[3]/1e6);
      }
    }
    int clusiz=0;
    int clutyp=0;
    double totampl=0;
    int maxsigpix=-1;
    double maxsig=-999.;
    for(int k=0; k<4; k++){
      if(amplVec[k]>0){
  	clutyp+=1<<k;
  	clusiz+=1;
  	totampl+=amplVec[k]*1000;
	if(amplVec[k]>maxsig){
	  maxsig=amplVec[k];
	  maxsigpix=k;
	}
      }
    }    
    if((iTrigChan>=0 && amplVec[iTrigChan]>0) || (iTrigChan<0 && clusiz>0)){
      hCluSiz[jfil]->Fill(clusiz);
      hCluTyp[jfil]->Fill(clutyp);
      hTotAmpl[jfil]->Fill(totampl);
      hAmplMaxSigPixVsCluSiz[jfil]->Fill(clusiz,amplVec[maxsigpix]*1000);
      if(iTrigChan>=0) hAmplTrigChanVsCluSiz[jfil]->Fill(clusiz,amplVec[iTrigChan]*1000);
      hFallTimeMaxSigPix[jfil]->Fill(fallTimeVec[maxsigpix]/1000.);
      hRecoTimeMaxSigPix[jfil]->Fill(recoTimeVec[maxsigpix]/1e6);
      hFallTimeMaxSigPixVsAmplMaxSigPix[jfil]->Fill(amplVec[maxsigpix]*1000.,fallTimeVec[maxsigpix]/1000.);
      if(clusiz==1){
	hFallTimeMaxSigPixCluSiz1[jfil]->Fill(fallTimeVec[maxsigpix]/1000.);
	hRecoTimeMaxSigPixCluSiz1[jfil]->Fill(recoTimeVec[maxsigpix]/1e6);
 	hFallTimeMaxSigPixVsAmplMaxSigPixCluSiz1[jfil]->Fill(amplVec[maxsigpix]*1000.,fallTimeVec[maxsigpix]/1000.);
      }
      if(clusiz>1) hFallTimeMaxSigPixVsAmplMaxSigPixCluSizGt1[jfil]->Fill(amplVec[maxsigpix]*1000.,fallTimeVec[maxsigpix]/1000.);
    }
    if(iTrigChan>=0){
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
    }else{
      if(fallTimeVec[maxsigpix]<400.){
	hCluSizFastResp[jfil]->Fill(clusiz);
	hAmplCh1FastResp[jfil]->Fill(TMath::Max(0.,amplVec[0]*1000.));
	hAmplCh2FastResp[jfil]->Fill(TMath::Max(0.,amplVec[1]*1000.));
	hAmplCh3FastResp[jfil]->Fill(TMath::Max(0.,amplVec[2]*1000.));
	hAmplCh4FastResp[jfil]->Fill(TMath::Max(0.,amplVec[3]*1000.));	
      }
    }
    if(clusiz>1){
      if(amplVec[1]>0 && amplVec[0]>0 && amplVec[1]<amplVec[0]) hDeltaTime12[jfil]->Fill((t50Vec[1]-t50Vec[0])/1000.);
      else if(amplVec[1]>0  && amplVec[0]>0 && amplVec[1]>=amplVec[0]) hDeltaTime21[jfil]->Fill((t50Vec[1]-t50Vec[0])/1000.);
      if(amplVec[2]>0  && amplVec[0]>0 && amplVec[2]<amplVec[0]) hDeltaTime13[jfil]->Fill((t50Vec[2]-t50Vec[0])/1000.);
      else if(amplVec[2]>0  && amplVec[0]>0 && amplVec[2]>=amplVec[0]) hDeltaTime31[jfil]->Fill((t50Vec[2]-t50Vec[0])/1000.);
      if(amplVec[3]>0 && amplVec[0]>0  && amplVec[3]<amplVec[0]) hDeltaTime14[jfil]->Fill((t50Vec[3]-t50Vec[0])/1000.);
      else if(amplVec[3]>0  && amplVec[0]>0 && amplVec[3]>=amplVec[0]) hDeltaTime41[jfil]->Fill((t50Vec[3]-t50Vec[0])/1000.);
      if(amplVec[2]>0 && amplVec[1]>0  && amplVec[2]<amplVec[1]) hDeltaTime23[jfil]->Fill((t50Vec[2]-t50Vec[1])/1000.);
      else if(amplVec[2]>0  && amplVec[1]>0 && amplVec[2]>=amplVec[1]) hDeltaTime32[jfil]->Fill((t50Vec[2]-t50Vec[1])/1000.);
      if(amplVec[3]>0 && amplVec[1]>0  && amplVec[3]<amplVec[1]) hDeltaTime24[jfil]->Fill((t50Vec[3]-t50Vec[1])/1000.);
      else if(amplVec[3]>0  && amplVec[1]>0 && amplVec[3]>=amplVec[1]) hDeltaTime42[jfil]->Fill((t50Vec[3]-t50Vec[1])/1000.);
      if(amplVec[3]>0 && amplVec[2]>0  && amplVec[3]<amplVec[2]) hDeltaTime34[jfil]->Fill((t50Vec[3]-t50Vec[2])/1000.);
      else if(amplVec[3]>0  && amplVec[2]>0 && amplVec[3]>=amplVec[2]) hDeltaTime43[jfil]->Fill((t50Vec[3]-t50Vec[2])/1000.);
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
      jTrigChan[jf]=0;
      if(trigChan[jf].Contains("J5")) jTrigChan[jf]=1;
      else if(trigChan[jf].Contains("J10")) jTrigChan[jf]=3;
      else if(trigChan[jf].Contains("OR")) jTrigChan[jf]=0;
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
  double cntMaxFT04[4];
  double cntMax[4];
  for(int j=0; j<nFiles; j++){
    hFallTimeTrigChan[j]=new TH1F(Form("hFallTimeTrigChan_%d",j)," All clusters ; Fall Time TrigChan (ns) ; Entries",100,0.,10.);
    hFallTimeCh1[j]=new TH1F(Form("hFallTimeCh1_%d",j)," All clusters ; Fall Time Ch1 (ns) ; Entries",100,0.,10.);
    hFallTimeCh2[j]=new TH1F(Form("hFallTimeCh2_%d",j)," All clusters ; Fall Time Ch2 (ns) ; Entries",100,0.,10.);
    hFallTimeCh3[j]=new TH1F(Form("hFallTimeCh3_%d",j)," All clusters ; Fall Time Ch3 (ns) ; Entries",100,0.,10.);
    hFallTimeCh4[j]=new TH1F(Form("hFallTimeCh4_%d",j)," All clusters ; Fall Time Ch4 (ns) ; Entries",100,0.,10.);
    hFallTimeMaxSigPix[j]=new TH1F(Form("hFallTimeMaxSigPix_%d",j)," All clusters ; Fall Time MaxSigPix (ns) ; Entries",100,0.,10.);
    hFallTimeMaxSigPixCluSiz1[j]=new TH1F(Form("hFallTimeMaxSigPixCluSiz1_%d",j)," Cluster size = 1 ; Fall Time MaxSigPix (ns) ; Entries",100,0.,10.);
    hAmplTrigChan[j]=new TH1F(Form("hAmplTrigChan_%d",j)," ; Signal Amplitude TrigChan (mV) ; Entries",100,0.,100.);
    hAmplCh1[j]=new TH1F(Form("hAmplCh1_%d",j)," ; Signal Amplitude Ch1 (mV) ; Entries",100,0.,100.);
    hAmplCh2[j]=new TH1F(Form("hAmplCh2_%d",j)," ; Signal Amplitude Ch2 (mV) ; Entries",100,0.,100.);
    hAmplCh3[j]=new TH1F(Form("hAmplCh3_%d",j)," ; Signal Amplitude Ch3 (mV) ; Entries",100,0.,100.);
    hAmplCh4[j]=new TH1F(Form("hAmplCh4_%d",j)," ; Signal Amplitude Ch4 (mV) ; Entries",100,0.,100.);
    hRecoTimeTrigChan[j]=new TH1F(Form("hRecoTimeTrigChan_%d",j)," ; Recovery Time TrigChan (#mus) ; Entries",100,0.,3.);
    hRecoTimeCh1[j]=new TH1F(Form("hRecoTimeCh1_%d",j)," All clusters ; Recovery Time Ch1 (#mus) ; Entries",100,0.,3.);
    hRecoTimeCh2[j]=new TH1F(Form("hRecoTimeCh2_%d",j)," All clusters ; Recovery Time Ch2 (#mus) ; Entries",100,0.,3.);
    hRecoTimeCh3[j]=new TH1F(Form("hRecoTimeCh3_%d",j)," All clusters ; Recovery Time Ch3 (#mus) ; Entries",100,0.,3.);
    hRecoTimeCh4[j]=new TH1F(Form("hRecoTimeCh4_%d",j)," All clusters ; Recovery Time Ch4 (#mus) ; Entries",100,0.,3.);
    hRecoTimeMaxSigPix[j]=new TH1F(Form("hRecoTimeMaxSigPix_%d",j)," All clusters ; Recovery Time MaxSigPix (#mus) ; Entries",100,0.,3.);
    hRecoTimeMaxSigPixCluSiz1[j]=new TH1F(Form("hRecoTimeMaxSigPixCluSiz1_%d",j)," Cluster size = 1 ; Recovery Time MaxSigPix (#mus) ; Entries",100,0.,3.);
    hFallTimeTrigChanVsAmplTrigChan[j]=new TH2F(Form("hFallTimeTrigChanVsAmplTrigChan_%d",j)," ; Signal Amplitude TrigChan (mV) ; Fall Time TrigChan (ns) ; Entries",50,0.,100.,100,0.,10.);
    hFallTimeTrigChanVsAmplTrigChanCluSiz1[j]=new TH2F(Form("hFallTimeTrigChanVsAmplTrigChanCluSiz1_%d",j)," Cluster size = 1  ; Signal Amplitude TrigChan (mV) ; Fall Time TrigChan (ns) ; Entries",50,0.,100.,100,0.,10.);
    hFallTimeTrigChanVsAmplTrigChanCluSiz1[j]->SetTitle(Form("%s Cluster size = 1",legEntry[j].Data()));
    hFallTimeTrigChanVsAmplTrigChanCluSizGt1[j]=new TH2F(Form("hFallTimeTrigChanVsAmplTrigChanCluSizGt1_%d",j)," Cluster size > 1  ; Signal Amplitude TrigChan (mV) ; Fall Time TrigChan (ns) ; Entries",50,0.,100.,100,0.,10.);
    hFallTimeTrigChanVsAmplTrigChanCluSizGt1[j]->SetTitle(Form("%s Cluster size > 1",legEntry[j].Data()));
    hFallTimeMaxSigPixVsAmplMaxSigPix[j]=new TH2F(Form("hFallTimeMaxSigPixVsAmplMaxSigPix_%d",j)," ; Signal Amplitude MaxSigPix (mV) ; Fall Time MaxSigPix (ns) ; Entries",50,0.,100.,100,0.,10.);
    hFallTimeMaxSigPixVsAmplMaxSigPixCluSiz1[j]=new TH2F(Form("hFallTimeMaxSigPixVsAmplMaxSigPixCluSiz1_%d",j)," Cluster size = 1  ; Signal Amplitude MaxSigPix (mV) ; Fall Time MaxSigPix (ns) ; Entries",50,0.,100.,100,0.,10.);
    hFallTimeMaxSigPixVsAmplMaxSigPixCluSiz1[j]->SetTitle(Form("%s Cluster size = 1",legEntry[j].Data()));
    hFallTimeMaxSigPixVsAmplMaxSigPixCluSizGt1[j]=new TH2F(Form("hFallTimeMaxSigPixVsAmplMaxSigPixCluSizGt1_%d",j)," Cluster size > 1  ; Signal Amplitude MaxSigPix (mV) ; Fall Time MaxSigPix (ns) ; Entries",50,0.,100.,100,0.,10.);
    hFallTimeMaxSigPixVsAmplMaxSigPixCluSizGt1[j]->SetTitle(Form("%s Cluster size > 1",legEntry[j].Data()));
    hTotAmpl[j]=new TH1F(Form("hTotAmpl_%d",j)," ; Signal Amplitude 4-pixels (mV) ; Entries",100,0.,100.);
    hCluSiz[j]=new TH1F(Form("hCluSiz_%d",j)," ; Cluster Size ; Fraction of events",5,-0.5,4.5);
    hCluTyp[j]=new TH1F(Form("hCluTyp_%d",j)," ; Cluster Shape ; Fraction of events",16,-0.5,15.5);
    hAmplTrigChanVsCluSiz[j]=new TH2F(Form("hAmplTrigChanVsCluSiz_%d",j)," ; Cluster Size ; Signal Amplitude TrigChan (mV) ; Entries",5,-0.5,4.5,100,0.,100.);    
    hAmplMaxSigPixVsCluSiz[j]=new TH2F(Form("hAmplMaxSigPixVsCluSiz_%d",j)," ; Cluster Size ; Signal Amplitude MaxSigPix (mV) ; Entries",5,-0.5,4.5,100,0.,100.);    
    hFallTimeTrigChanCluSiz1[j]=new TH1F(Form("hFallTimeTrigChanCluSiz1_%d",j)," Cluster size = 1 ; Fall Time TrigChan (ns) ; Entries",100,0.,10.);    
    hAmplCh1FastResp[j]=new TH1F(Form("hAmplCh1FastResp_%d",j)," Signals with fall time < 0.4 ns ; Signal Amplitude Ch1 (mV) ; Entries",100,0.,100.);
    hAmplCh2FastResp[j]=new TH1F(Form("hAmplCh2FastResp_%d",j)," Signals with fall time < 0.4 ns ; Signal Amplitude Ch2 (mV) ; Entries",100,0.,100.);
    hAmplCh3FastResp[j]=new TH1F(Form("hAmplCh3FastResp_%d",j)," Signals with fall time < 0.4 ns ; Signal Amplitude Ch3 (mV) ; Entries",100,0.,100.);
    hAmplCh4FastResp[j]=new TH1F(Form("hAmplCh4FastResp_%d",j)," Signals with fall time < 0.4 ns ; Signal Amplitude Ch4 (mV) ; Entries",100,0.,100.);
    hCluSizFastResp[j]=new TH1F(Form("hCluSizFastResp_%d",j)," Signals with fall time < 0.4 ns ; Cluster Size ; Fraction of events",5,-0.5,4.5);
    hDeltaTime12[j]=new TH1F(Form("hDeltaTime12_%d",j)," Events with Ampl1 > Ampl2 ; #Deltat_{50} (Ch2 - Ch1) (ns) ; Entries",100,-6,6.);
    hDeltaTime13[j]=new TH1F(Form("hDeltaTime13_%d",j)," Events with Ampl1 > Ampl3 ; #Deltat_{50} (Ch3 - Ch1) (ns) ; Entries",100,-6,6.);
    hDeltaTime14[j]=new TH1F(Form("hDeltaTime14_%d",j)," Events with Ampl1 > Ampl4 ; #Deltat_{50} (Ch4 - Ch1) (ns) ; Entries",100,-6,6.);
    hDeltaTime23[j]=new TH1F(Form("hDeltaTime23_%d",j)," Events with Ampl2 > Ampl3 ; #Deltat_{50} (Ch3 - Ch2) (ns) ; Entries",100,-6,6.);
    hDeltaTime24[j]=new TH1F(Form("hDeltaTime24_%d",j)," Events with Ampl2 > Ampl4 ; #Deltat_{50} (Ch4 - Ch2) (ns) ; Entries",100,-6,6.);
    hDeltaTime34[j]=new TH1F(Form("hDeltaTime34_%d",j)," Events with Ampl3 > Ampl4 ; #Deltat_{50} (Ch4 - Ch3) (ns) ; Entries",100,-6,6.);
    hDeltaTime21[j]=new TH1F(Form("hDeltaTime21_%d",j)," Events with Ampl1 < Ampl2 ; #Deltat_{50} (Ch2 - Ch1) (ns) ; Entries",100,-6,6.);
    hDeltaTime31[j]=new TH1F(Form("hDeltaTime31_%d",j)," Events with Ampl1 < Ampl3 ; #Deltat_{50} (Ch3 - Ch1) (ns) ; Entries",100,-6,6.);
    hDeltaTime41[j]=new TH1F(Form("hDeltaTime41_%d",j)," Events with Ampl1 < Ampl4 ; #Deltat_{50} (Ch4 - Ch1) (ns) ; Entries",100,-6,6.);
    hDeltaTime32[j]=new TH1F(Form("hDeltaTime32_%d",j)," Events with Ampl2 < Ampl3 ; #Deltat_{50} (Ch3 - Ch2) (ns) ; Entries",100,-6,6.);
    hDeltaTime42[j]=new TH1F(Form("hDeltaTime42_%d",j)," Events with Ampl2 < Ampl4 ; #Deltat_{50} (Ch4 - Ch2) (ns) ; Entries",100,-6,6.);
    hDeltaTime43[j]=new TH1F(Form("hDeltaTime43_%d",j)," Events with Ampl3 < Ampl4 ; #Deltat_{50} (Ch4 - Ch3) (ns) ; Entries",100,-6,6.);
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
    arrHisto->AddAtAndExpand(hFallTimeTrigChan[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeCh1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeCh2[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeCh3[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeCh4[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeMaxSigPix[j],indexh++);    
    arrHisto->AddAtAndExpand(hFallTimeMaxSigPixCluSiz1[j],indexh++);    
    arrHisto->AddAtAndExpand(hAmplTrigChan[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh1[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh2[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh3[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh4[j],indexh++);
    arrHisto->AddAtAndExpand(hRecoTimeTrigChan[j],indexh++);
    arrHisto->AddAtAndExpand(hRecoTimeCh1[j],indexh++);
    arrHisto->AddAtAndExpand(hRecoTimeCh2[j],indexh++);
    arrHisto->AddAtAndExpand(hRecoTimeCh3[j],indexh++);
    arrHisto->AddAtAndExpand(hRecoTimeCh4[j],indexh++);
    arrHisto->AddAtAndExpand(hRecoTimeMaxSigPix[j],indexh++);    
    arrHisto->AddAtAndExpand(hRecoTimeMaxSigPixCluSiz1[j],indexh++);    
    arrHisto->AddAtAndExpand(hFallTimeTrigChanVsAmplTrigChan[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeTrigChanVsAmplTrigChanCluSiz1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeTrigChanVsAmplTrigChanCluSizGt1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeMaxSigPixVsAmplMaxSigPix[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeMaxSigPixVsAmplMaxSigPixCluSiz1[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeMaxSigPixVsAmplMaxSigPixCluSizGt1[j],indexh++);
    arrHisto->AddAtAndExpand(hTotAmpl[j],indexh++);
    arrHisto->AddAtAndExpand(hCluSiz[j],indexh++);
    arrHisto->AddAtAndExpand(hCluTyp[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplTrigChanVsCluSiz[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplMaxSigPixVsCluSiz[j],indexh++);
    arrHisto->AddAtAndExpand(hFallTimeTrigChanCluSiz1[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh1FastResp[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh2FastResp[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh3FastResp[j],indexh++);
    arrHisto->AddAtAndExpand(hAmplCh4FastResp[j],indexh++);
    arrHisto->AddAtAndExpand(hCluSizFastResp[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime12[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime13[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime14[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime23[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime24[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime34[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime21[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime31[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime41[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime32[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime42[j],indexh++);
    arrHisto->AddAtAndExpand(hDeltaTime43[j],indexh++);
 
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
    profFallTimeCluSiz1[j]=hFallTimeMaxSigPixVsAmplMaxSigPixCluSiz1[j]->ProfileX(Form("profFallTimeVsAmplCluSiz1%d",j));
    profFallTimeCluSiz1[j]->GetYaxis()->SetTitle("<Fall Time MaxSigPix> (ns)");
    profFallTimeCluSiz1[j]->SetLineColor(cols[j]);
    profFallTimeCluSiz1[j]->SetMarkerColor(cols[j]);
    profFallTimeCluSiz1[j]->SetMarkerStyle(20);
    profFallTimeCluSiz1[j]->SetLineWidth(2);
    profFallTimeCluSiz1[j]->SetStats(0);
    profFallTimeCluSizGt1[j]=hFallTimeMaxSigPixVsAmplMaxSigPixCluSizGt1[j]->ProfileX(Form("profFallTimeVsAmplCluSizGt1%d",j));
    profFallTimeCluSizGt1[j]->GetYaxis()->SetTitle("<Fall Time MaxSigPix> (ns)");
    profFallTimeCluSizGt1[j]->SetLineColor(cols[j]);
    profFallTimeCluSizGt1[j]->SetMarkerColor(cols[j]);
    profFallTimeCluSizGt1[j]->SetMarkerStyle(25);
    profFallTimeCluSizGt1[j]->SetLineWidth(2);
    profFallTimeCluSizGt1[j]->SetStats(0);
    
    cnt04[j]=hFallTimeMaxSigPixCluSiz1[j]->Integral(1,hFallTimeMaxSigPixCluSiz1[j]->GetXaxis()->FindBin(0.3999));
    cntall[j]=hFallTimeMaxSigPixCluSiz1[j]->GetEntries();
    cntFT04[j]=hFallTimeMaxSigPix[j]->Integral(1,hFallTimeMaxSigPix[j]->GetXaxis()->FindBin(0.3999));
    cntFTall[j]=hFallTimeMaxSigPix[j]->GetEntries();
    cntMaxFT04[j]=hCluSizFastResp[j]->GetEntries();
    cntMax[j]=hCluSiz[j]->GetEntries();


    if(normalizeToArea) NormalizeHistos(arrHisto);

    profAmplTrigChan[j]=hAmplTrigChanVsCluSiz[j]->ProfileX(Form("profAmplTrigChan%d",j));
    profAmplTrigChan[j]->GetYaxis()->SetTitle("<Amplitude TrigChan> (mV)");
    profAmplMaxSigPix[j]=hAmplMaxSigPixVsCluSiz[j]->ProfileX(Form("profAmplMaxSigPix%d",j));
    profAmplMaxSigPix[j]->GetYaxis()->SetTitle("<Amplitude MaxSigPix> (mV)");
    hCluTyp[j]->SetStats(0);
    hAmplTrigChanVsCluSiz[j]->SetStats(0);
    hAmplMaxSigPixVsCluSiz[j]->SetStats(0);
    profAmplTrigChan[j]->SetLineColor(cols[j]);
    profAmplTrigChan[j]->SetLineWidth(2);
    profAmplTrigChan[j]->SetStats(0);
    profAmplMaxSigPix[j]->SetLineColor(cols[j]);
    profAmplMaxSigPix[j]->SetLineWidth(2);
    profAmplMaxSigPix[j]->SetStats(0);
  }
  SetMaxima(nFiles);

  TCanvas* cA = new TCanvas("cA","Amplitudes",1650,900);
  cA->Divide(3,2);
  cA->cd(1);
  TLegend* leg= new TLegend(0.12,0.6,0.5,0.89);
  leg->SetMargin(0.1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplCh1[j]->Draw("histo");
    else hAmplCh1[j]->Draw("histo,sames");
    leg->AddEntry(hAmplCh1[j],legEntry[j].Data(),"L")->SetTextColor(cols[j]);
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplCh1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplCh1[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
  cA->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplCh2[j]->Draw("histo");
    else hAmplCh2[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplCh2[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplCh2[j]->GetLineColor());
    gPad->Modified();
  }
  cA->cd(5);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplCh3[j]->Draw("histo");
    else hAmplCh3[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplCh3[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplCh3[j]->GetLineColor());
    gPad->Modified();
  } 
  cA->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplCh4[j]->Draw("histo");
    else hAmplCh4[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplCh4[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplCh4[j]->GetLineColor());
    gPad->Modified();
  }
  cA->cd(6);
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
  cA->SaveAs(Form("Amplitudes_%s.png",suffix.Data()));

  TCanvas* cF = new TCanvas("cF","FallTimes",1650,900);
  cF->Divide(3,2);
  cF->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeCh1[j]->Draw("histo");
    else hFallTimeCh1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeCh1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeCh1[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
  cF->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeCh2[j]->Draw("histo");
    else hFallTimeCh2[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeCh2[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeCh2[j]->GetLineColor());
    gPad->Modified();
  }
  cF->cd(5);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeCh3[j]->Draw("histo");
    else hFallTimeCh3[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeCh3[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeCh3[j]->GetLineColor());
    gPad->Modified();
  } 
  cF->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeCh4[j]->Draw("histo");
    else hFallTimeCh4[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeCh4[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeCh4[j]->GetLineColor());
    gPad->Modified();
  }
  cF->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeMaxSigPix[j]->Draw("histo");
    else hFallTimeMaxSigPix[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeMaxSigPix[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeMaxSigPix[j]->GetLineColor());
    gPad->Modified();
  }  
  cF->cd(6);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeMaxSigPixCluSiz1[j]->Draw("histo");
    else hFallTimeMaxSigPixCluSiz1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeMaxSigPixCluSiz1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeMaxSigPixCluSiz1[j]->GetLineColor());
    gPad->Modified();
  }
  cF->SaveAs(Form("FallTime_%s.png",suffix.Data()));

  TCanvas* cR = new TCanvas("cR","RecoTimes",1650,900);
  cR->Divide(3,2);
  cR->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hRecoTimeCh1[j]->Draw("histo");
    else hRecoTimeCh1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hRecoTimeCh1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hRecoTimeCh1[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
  cR->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hRecoTimeCh2[j]->Draw("histo");
    else hRecoTimeCh2[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hRecoTimeCh2[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hRecoTimeCh2[j]->GetLineColor());
    gPad->Modified();
  }
  cR->cd(5);
  for(int j=0; j<nFiles; j++){
    if(j==0) hRecoTimeCh3[j]->Draw("histo");
    else hRecoTimeCh3[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hRecoTimeCh3[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hRecoTimeCh3[j]->GetLineColor());
    gPad->Modified();
  } 
  cR->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hRecoTimeCh4[j]->Draw("histo");
    else hRecoTimeCh4[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hRecoTimeCh4[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hRecoTimeCh4[j]->GetLineColor());
    gPad->Modified();
  }
  cR->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hRecoTimeMaxSigPix[j]->Draw("histo");
    else hRecoTimeMaxSigPix[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hRecoTimeMaxSigPix[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hRecoTimeMaxSigPix[j]->GetLineColor());
    gPad->Modified();
  }  
  cR->cd(6);
  for(int j=0; j<nFiles; j++){
    if(j==0) hRecoTimeMaxSigPixCluSiz1[j]->Draw("histo");
    else hRecoTimeMaxSigPixCluSiz1[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hRecoTimeMaxSigPixCluSiz1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hRecoTimeMaxSigPixCluSiz1[j]->GetLineColor());
    gPad->Modified();
  }
  cR->SaveAs(Form("RecoTime_%s.png",suffix.Data()));
  
  TCanvas* c1 = new TCanvas("c1","Histos TrigChan",1650,900);
  c1->Divide(3,2);
  c1->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeTrigChan[j]->Draw("histo");
    else hFallTimeTrigChan[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeTrigChan[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeTrigChan[j]->GetLineColor());
    gPad->Modified();
  } 
  c1->cd(2);
  for(int j=0; j<nFiles; j++){
    if(j==0) hRecoTimeTrigChan[j]->Draw("histo");
    else hRecoTimeTrigChan[j]->Draw("histo,sames");
    gPad->Update();
    TPaveStats *st=(TPaveStats*)hRecoTimeTrigChan[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hRecoTimeTrigChan[j]->GetLineColor());
    gPad->Modified();
  }
  c1->cd(3);
  for(int j=0; j<nFiles; j++){
    if(j==0) hAmplTrigChan[j]->Draw("histo");
    else hAmplTrigChan[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hAmplTrigChan[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hAmplTrigChan[j]->GetLineColor());
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
  c1->SaveAs(Form("TrigChan_%s.png",suffix.Data()));
  
  TCanvas* c2 = new TCanvas("c2","Histos CluSiz",1650,900);
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
    if(j==0) hFallTimeMaxSigPix[j]->Draw("histo");
    else hFallTimeMaxSigPix[j]->Draw("histo,sames");
    TLatex* tfrach=new TLatex(0.5,0.72-0.08*j,Form("%.0f/%.0f=%.3f",cntFT04[j],cntFTall[j],cntFT04[j]/cntFTall[j]));
    tfrach->SetNDC();
    tfrach->SetTextColor(cols[j]);
    tfrach->SetTextFont(43);
    tfrach->SetTextSize(20);
    tfrach->Draw();
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hFallTimeMaxSigPix[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hFallTimeMaxSigPix[j]->GetLineColor());
    gPad->Modified();
  } 
  TLatex* tfrac=new TLatex(0.48,0.8,"Fall time < 0.4 ns");
  tfrac->SetNDC();
  tfrac->SetTextFont(43);
  tfrac->SetTextSize(20);
  tfrac->Draw();
  c2->cd(4);
  hAmplMaxSigPixVsCluSiz[0]->Draw("box");
  for(int j=1; j<nFiles; j++) hAmplMaxSigPixVsCluSiz[j]->Draw("same,box");
  c2->cd(5);
  profAmplMaxSigPix[0]->Draw("");
  for(int j=1; j<nFiles; j++) profAmplMaxSigPix[j]->Draw("same");
  c2->cd(6);
  for(int j=0; j<nFiles; j++){
    if(j==0) hFallTimeMaxSigPixCluSiz1[j]->Draw("histo");
    else hFallTimeMaxSigPixCluSiz1[j]->Draw("histo,sames");
    TLatex* tfrach=new TLatex(0.5,0.72-0.08*j,Form("%.0f/%.0f=%.3f",cnt04[j],cntall[j],cnt04[j]/cntall[j]));
    tfrach->SetNDC();
    tfrach->SetTextColor(cols[j]);
    tfrach->SetTextFont(43);
    tfrach->SetTextSize(20);
    tfrach->Draw();
    gPad->Update();
    TPaveStats *st=(TPaveStats*)hFallTimeMaxSigPixCluSiz1[j]->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hTotAmpl[j]->GetLineColor());
    gPad->Modified();
  }
  tfrac->Draw();
  leg->Draw();
  c2->SaveAs(Form("CluSiz_%s.png",suffix.Data()));
  
  TCanvas* c3 = new TCanvas("c3","FastResponse",1650,900);
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
    TLatex* t2=new TLatex(0.05,0.72-0.25*j,Form("Fraction with FT<0.4 ns on peak pixel = %.0f / %.0f = %.3f",cntMaxFT04[j],cntMax[j],cntMaxFT04[j]/cntMax[j]));
    t2->SetTextFont(43);
    t2->SetTextSize(18);
    t2->SetNDC();
    t2->SetTextColor(cols[j]);
    t2->Draw();
  }
  c3->SaveAs(Form("FastResp_%s.png",suffix.Data()));
  
  TCanvas* c4 = new TCanvas("c4","Delta t50",1650,900);
  c4->Divide(6,2);
  c4->cd(1);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime12[j]->Draw("histo");
    else hDeltaTime12[j]->Draw("histo,sames");
    gPad->Update();    
    TPaveStats *st=(TPaveStats*)hDeltaTime12[j]->GetListOfFunctions()->FindObject("stats");
    st->SetX1NDC(0.6);
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
    st->SetX1NDC(0.6);
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
    st->SetX1NDC(0.6);
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime14[j]->GetLineColor());
    gPad->Modified();
  }
  c4->cd(4);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime23[j]->Draw("histo");
    else hDeltaTime23[j]->Draw("histo,sames");
    gPad->Update();
    TPaveStats *st=(TPaveStats*)hDeltaTime23[j]->GetListOfFunctions()->FindObject("stats");
    st->SetX1NDC(0.6);
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime23[j]->GetLineColor());
    gPad->Modified();
  }
  leg->Draw();
  c4->cd(5);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime24[j]->Draw("histo");
    else hDeltaTime24[j]->Draw("histo,sames");
    gPad->Update();
    TPaveStats *st=(TPaveStats*)hDeltaTime24[j]->GetListOfFunctions()->FindObject("stats");
    st->SetX1NDC(0.6);
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime24[j]->GetLineColor());
    gPad->Modified();
  }
  c4->cd(6);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime34[j]->Draw("histo");
    else hDeltaTime34[j]->Draw("histo,sames");
    gPad->Update();
    TPaveStats *st=(TPaveStats*)hDeltaTime34[j]->GetListOfFunctions()->FindObject("stats");
    st->SetX1NDC(0.6);
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime34[j]->GetLineColor());
    gPad->Modified();
  }
  c4->cd(7);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime21[j]->Draw("histo");
    else hDeltaTime21[j]->Draw("histo,sames");
    gPad->Update();
    TPaveStats *st=(TPaveStats*)hDeltaTime21[j]->GetListOfFunctions()->FindObject("stats");
    st->SetX1NDC(0.6);
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime21[j]->GetLineColor());
    gPad->Modified();
  }
  c4->cd(8);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime31[j]->Draw("histo");
    else hDeltaTime31[j]->Draw("histo,sames");
    gPad->Update();
    TPaveStats *st=(TPaveStats*)hDeltaTime31[j]->GetListOfFunctions()->FindObject("stats");
    st->SetX1NDC(0.6);
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime31[j]->GetLineColor());
    gPad->Modified();
  }
  c4->cd(9);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime41[j]->Draw("histo");
    else hDeltaTime41[j]->Draw("histo,sames");
    gPad->Update();
    TPaveStats *st=(TPaveStats*)hDeltaTime41[j]->GetListOfFunctions()->FindObject("stats");
    st->SetX1NDC(0.6);
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime41[j]->GetLineColor());
    gPad->Modified();
  }
  c4->cd(10);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime32[j]->Draw("histo");
    else hDeltaTime32[j]->Draw("histo,sames");
    gPad->Update();
    TPaveStats *st=(TPaveStats*)hDeltaTime32[j]->GetListOfFunctions()->FindObject("stats");
    st->SetX1NDC(0.6);
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime32[j]->GetLineColor());
    gPad->Modified();
  }
  c4->cd(11);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime42[j]->Draw("histo");
    else hDeltaTime42[j]->Draw("histo,sames");
    gPad->Update();
    TPaveStats *st=(TPaveStats*)hDeltaTime42[j]->GetListOfFunctions()->FindObject("stats");
    st->SetX1NDC(0.6);
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime42[j]->GetLineColor());
    gPad->Modified();
  }
  c4->cd(12);
  for(int j=0; j<nFiles; j++){
    if(j==0) hDeltaTime43[j]->Draw("histo");
    else hDeltaTime43[j]->Draw("histo,sames");
    gPad->Update();
    TPaveStats *st=(TPaveStats*)hDeltaTime43[j]->GetListOfFunctions()->FindObject("stats");
    st->SetX1NDC(0.6);
    st->SetY1NDC(0.77-0.2*j);
    st->SetY2NDC(0.96-0.2*j);
    st->SetTextColor(hDeltaTime43[j]->GetLineColor());
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

  TCanvas* c5=new TCanvas("c5","FallTime vs Ampl",1650,900);
  c5->Divide(nFiles+1,2);
  for(int j=0; j<nFiles; j++){
    c5->cd(j+1);
    hFallTimeMaxSigPixVsAmplMaxSigPixCluSiz1[j]->Draw("colz");
    c5->cd(j+2+nFiles);
    hFallTimeMaxSigPixVsAmplMaxSigPixCluSizGt1[j]->Draw("colz");
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
