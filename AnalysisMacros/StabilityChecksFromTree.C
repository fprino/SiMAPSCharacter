#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TTimeStamp.h>
#include <TLegendEntry.h>
#include <TCanvas.h>


void StabilityChecksFromTree(TString treeFileName="APTS10_Vbb-4V_WalerConfig_2000trg_TTree.root", int trigchan=3){
  int ev;
  ULong64_t timest;
  double basel[4],ampl[4],falltime[4];
  TFile* f=new TFile(treeFileName.Data());
  TTree* tree=(TTree*)f->Get("treeParams");
  tree->SetBranchAddress("Event",&ev);
  tree->SetBranchAddress("Timestamp",&timest);
  for(int jpix=0; jpix<4; jpix++){
    tree->SetBranchAddress(Form("BaselineCh%d",jpix+1),&basel[jpix]);
    tree->SetBranchAddress(Form("SignalAmplCh%d",jpix+1),&ampl[jpix]);
    tree->SetBranchAddress(Form("FallTimeCh%d",jpix+1),&falltime[jpix]);
  }
  TGraph* gTimeFromStart=new TGraph(0);
  TGraph* gDeltaTim=new TGraph(0);
  gDeltaTim->SetTitle("");
  gDeltaTim->GetXaxis()->SetTitle("Event");
  gDeltaTim->GetYaxis()->SetTitle("Time to next event (s)");
  gTimeFromStart->SetTitle("");
  gTimeFromStart->GetXaxis()->SetTitle("Event");
  gTimeFromStart->GetYaxis()->SetTitle("Time from start of run (s)");
  TH1F* hDeltaTim=new TH1F("hDeltaTim"," ; Time to next event (s) ; Entries",100,0.,400.);
  TGraph** gBasel=new TGraph*[4];
  TGraph** gBaselSig=new TGraph*[4];
  TGraph** gFallTime=new TGraph*[4];
  TGraph** gAmpl=new TGraph*[4];
  for(int jpix=0; jpix<4; jpix++){
    gBasel[jpix]=new TGraph(0);
    gBaselSig[jpix]=new TGraph(0);
    gFallTime[jpix]=new TGraph(0);
    gAmpl[jpix]=new TGraph(0);
    gBasel[jpix]->SetTitle(Form("Channel %d",jpix+1));
    gBasel[jpix]->GetXaxis()->SetTitle("Event");
    gBasel[jpix]->GetYaxis()->SetTitle("Baseline");
    gBaselSig[jpix]->SetTitle(Form("Channel %d, events with signal",jpix+1));
    gBaselSig[jpix]->GetXaxis()->SetTitle("Event");
    gBaselSig[jpix]->GetYaxis()->SetTitle("Baseline");
    gFallTime[jpix]->SetTitle(Form("Channel %d",jpix+1));
    gFallTime[jpix]->GetXaxis()->SetTitle("Event");
    gFallTime[jpix]->GetYaxis()->SetTitle("Fall Time (ns)");
    gAmpl[jpix]->SetTitle(Form("Channel %d",jpix+1));
    gAmpl[jpix]->GetXaxis()->SetTitle("Event");
    gAmpl[jpix]->GetYaxis()->SetTitle("Signal Amplitude (mV)");
  }
  ULong64_t timestold=0;
  ULong64_t tstart=0;
  double totEv=tree->GetEntriesFast();
  double evWithSignal[4]={0.,0.,0.,0.};
  double averampl[4]={0.,0.,0.,0.};
  for(int ient=0; ient<tree->GetEntriesFast(); ient++){
    tree->GetEvent(ient);
    TTimeStamp tcur(timest);
    if(ient==0) tstart=timest;
    ULong64_t timFromStart=timest-tstart;
    gTimeFromStart->SetPoint(ient,ev,(double)timFromStart);
    if(timestold>0){
      ULong64_t dtim=timest-timestold;
      if(timestold<=timest){
	gDeltaTim->SetPoint(gDeltaTim->GetN(),ev,(double)dtim);
	hDeltaTim->Fill((double)dtim);
      }else{
	if(timestold>0){
	  TTimeStamp tprev(timestold);
	  printf("Reverted time at event %d? Current %d %d   Previous %d %d\n",ev,tcur.GetDate(),tcur.GetTime(),tprev.GetDate(),tprev.GetTime());
	}
      }
    }
    timestold=timest;
    if(ampl[trigchan-1]<0) printf("Event %d: no signal in ch %d\n",ient,trigchan);
    for(int jpix=0; jpix<4; jpix++){
      if(ampl[jpix]>0){
	evWithSignal[jpix]+=1;
	averampl[jpix]+=ampl[jpix]*1000.;
      }
      gBasel[jpix]->SetPoint(ient,ev,basel[jpix]);
      if(ampl[jpix]>0){
	gBaselSig[jpix]->SetPoint(ient,ev,basel[jpix]);
	gAmpl[jpix]->SetPoint(ient,ev,ampl[jpix]*1000.);
      }
      if(ampl[jpix]>0 && falltime[jpix]>=0) gFallTime[jpix]->SetPoint(ient,ev,falltime[jpix]/1000.);
    }
  }
  for(int jpix=0; jpix<4; jpix++){
    if(evWithSignal[jpix]>0) averampl[jpix]/=evWithSignal[jpix];
    gBasel[jpix]->SetMarkerStyle(7);
    gBaselSig[jpix]->SetMarkerStyle(7);
    gBaselSig[jpix]->SetMarkerColor(kGreen+1);
    gFallTime[jpix]->SetMarkerStyle(7);
    gAmpl[jpix]->SetMarkerStyle(7);
  }
  
  TCanvas* cT=new TCanvas("cT","",1600,800);
  cT->Divide(2,2);
  cT->cd(1);
  gTimeFromStart->SetMarkerStyle(7);
  gTimeFromStart->Draw("AP");
  cT->cd(2);
  gDeltaTim->SetMarkerStyle(7);
  gDeltaTim->Draw("AP");
  cT->cd(3);
  gPad->SetLogy();
  hDeltaTim->Draw();

  TCanvas* cB=new TCanvas("cB","",1600,800);
  cB->Divide(2,2);
  for(int jpix=0; jpix<4; jpix++){
    cB->cd(jpix+1);
    gBasel[jpix]->Draw("AP");
    gBaselSig[jpix]->Draw("PSAME");
    TLatex* t1=new TLatex(0.15,0.85,Form("Fraction of events with signal in ch%d= %.0f / %.0f = %.3f",jpix+1,evWithSignal[jpix],totEv,evWithSignal[jpix]/totEv));
    t1->SetTextColor(4);
    t1->SetNDC();
    t1->Draw();
    if(jpix==1){
      TLegend* leg=new TLegend(0.2,0.2,0.6,0.4);
      leg->AddEntry(gBasel[jpix],"All Events","P");
      leg->AddEntry(gBaselSig[jpix],"Events with signal in ch1","P")->SetTextColor(gBaselSig[jpix]->GetMarkerColor());
      leg->Draw();
    }
  }

  TCanvas* cA=new TCanvas("cA","",1600,800);
  cA->Divide(2,2);
  for(int jpix=0; jpix<4; jpix++){
    cA->cd(jpix+1);
    gAmpl[jpix]->Draw("AP");
    TLatex* tF=new TLatex(0.15,0.85,Form("Fraction of events with signal in ch %d=%.0f / %.0f = %.3f",jpix+1,evWithSignal[jpix],totEv,evWithSignal[jpix]/totEv));
    tF->SetTextColor(4);
    tF->SetNDC();
    tF->Draw();
    TLatex* tA=new TLatex(0.15,0.78,Form("<Amplitude> Ch %d = %.2f mV",jpix+1,averampl[jpix]));
    tA->SetTextColor(4);
    tA->SetNDC();
    tA->Draw();
  }

  TCanvas* cFT=new TCanvas("cFT","",1600,800);
  cFT->Divide(2,2);
  for(int jpix=0; jpix<4; jpix++){
    cFT->cd(jpix+1);
    gFallTime[jpix]->Draw("AP");
    TLatex* tF=new TLatex(0.15,0.85,Form("Fraction of events with signal in ch%d= %.0f / %.0f = %.3f",jpix+1,evWithSignal[jpix],totEv,evWithSignal[jpix]/totEv));
    tF->SetTextColor(4);
    tF->SetNDC();
    tF->Draw();
  }
}
