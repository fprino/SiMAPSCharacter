#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TH1.h>
#include <TLegend.h>
#include <TTimeStamp.h>
#include <TLegendEntry.h>
#include <TCanvas.h>


void StabilityChecksFromTree(TString treeFileName="APTS10_Vbb-4V_WalerConfig_2000trg_TTree.root"){
  int ev;
  ULong64_t timest;
  double basel1,basel2,basel3,basel4,ampl1;
  TFile* f=new TFile(treeFileName.Data());
  TTree* tree=(TTree*)f->Get("treeParams");
  tree->SetBranchAddress("Event",&ev);
  tree->SetBranchAddress("Timestamp",&timest);
  tree->SetBranchAddress("BaselineCh1",&basel1);
  tree->SetBranchAddress("BaselineCh2",&basel2);
  tree->SetBranchAddress("BaselineCh3",&basel3);
  tree->SetBranchAddress("BaselineCh4",&basel4);
  tree->SetBranchAddress("SignalAmplCh1",&ampl1);
  TGraph* gTimeFromStart=new TGraph(0);
  TGraph* gDeltaTim=new TGraph(0);
  TGraph* gBasel1=new TGraph(0);
  TGraph* gBasel2=new TGraph(0);
  TGraph* gBasel3=new TGraph(0);
  TGraph* gBasel4=new TGraph(0);
  TGraph* gBasel1Sig=new TGraph(0);
  TGraph* gBasel2Sig=new TGraph(0);
  TGraph* gBasel3Sig=new TGraph(0);
  TGraph* gBasel4Sig=new TGraph(0);
  gBasel1->SetTitle("Channel 1");
  gBasel2->SetTitle("Channel 2");
  gBasel3->SetTitle("Channel 3");
  gBasel4->SetTitle("Channel 4");
  gBasel1->GetXaxis()->SetTitle("Event");
  gBasel2->GetXaxis()->SetTitle("Event");
  gBasel3->GetXaxis()->SetTitle("Event");
  gBasel4->GetXaxis()->SetTitle("Event");
  gBasel1->GetYaxis()->SetTitle("Baseline");
  gBasel2->GetYaxis()->SetTitle("Baseline");
  gBasel3->GetYaxis()->SetTitle("Baseline");
  gBasel4->GetYaxis()->SetTitle("Baseline");
  gBasel1Sig->SetTitle("Channel 1, events with signal");
  gBasel2Sig->SetTitle("Channel 2, events with signal");
  gBasel3Sig->SetTitle("Channel 3, events with signal");
  gBasel4Sig->SetTitle("Channel 4, events with signal");
  gBasel1Sig->GetXaxis()->SetTitle("Event");
  gBasel2Sig->GetXaxis()->SetTitle("Event");
  gBasel3Sig->GetXaxis()->SetTitle("Event");
  gBasel4Sig->GetXaxis()->SetTitle("Event");
  gBasel1Sig->GetYaxis()->SetTitle("Baseline");
  gBasel2Sig->GetYaxis()->SetTitle("Baseline");
  gBasel3Sig->GetYaxis()->SetTitle("Baseline");
  gBasel4Sig->GetYaxis()->SetTitle("Baseline");
  gDeltaTim->SetTitle("");
  gDeltaTim->GetXaxis()->SetTitle("Event");
  gDeltaTim->GetYaxis()->SetTitle("Time to next event (s)");
  gTimeFromStart->SetTitle("");
  gTimeFromStart->GetXaxis()->SetTitle("Event");
  gTimeFromStart->GetYaxis()->SetTitle("Time from start of run (s)");
  TH1F* hDeltaTim=new TH1F("hDeltaTim"," ; Time to next event (s) ; Entries",100,0.,200.);
  ULong64_t timestold=-9;
  ULong64_t tstart=-9;
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
	TTimeStamp tprev(timestold);
	printf("Reverted time at event %d? Current %d %d   Previous %d %d\n",ev,tcur.GetDate(),tcur.GetTime(),tprev.GetDate(),tprev.GetTime());
      }
    }
    timestold=timest;
    gBasel1->SetPoint(ient,ev,basel1);
    gBasel2->SetPoint(ient,ev,basel2);
    gBasel3->SetPoint(ient,ev,basel3);
    gBasel4->SetPoint(ient,ev,basel4);
    if(ampl1>0){
      gBasel1Sig->SetPoint(ient,ev,basel1);
      gBasel2Sig->SetPoint(ient,ev,basel2);
      gBasel3Sig->SetPoint(ient,ev,basel3);
      gBasel4Sig->SetPoint(ient,ev,basel4);
    }
  }
  gBasel1->SetMarkerStyle(7);
  gBasel2->SetMarkerStyle(7);
  gBasel3->SetMarkerStyle(7);
  gBasel4->SetMarkerStyle(7);
  gBasel1Sig->SetMarkerStyle(7);
  gBasel2Sig->SetMarkerStyle(7);
  gBasel3Sig->SetMarkerStyle(7);
  gBasel4Sig->SetMarkerStyle(7);
  gBasel1Sig->SetMarkerColor(kGreen+1);
  gBasel2Sig->SetMarkerColor(kGreen+1);
  gBasel3Sig->SetMarkerColor(kGreen+1);
  gBasel4Sig->SetMarkerColor(kGreen+1);
  TCanvas* cB=new TCanvas("cB","",1600,800);
  cB->Divide(2,2);
  cB->cd(1);
  gBasel1->Draw("AP");
  gBasel1Sig->Draw("PSAME");
  cB->cd(2);
  gBasel2->Draw("AP");
  gBasel2Sig->Draw("PSAME");
  TLegend* leg=new TLegend(0.2,0.2,0.5,0.4);
  leg->AddEntry(gBasel2,"All Events","P");
  leg->AddEntry(gBasel2Sig,"Events with signal","P")->SetTextColor(kGreen+1);
  leg->Draw();
  cB->cd(3);
  gBasel3->Draw("AP");
  gBasel3Sig->Draw("PSAME");
  cB->cd(4);
  gBasel4->Draw("AP");
  gBasel4Sig->Draw("PSAME");

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
}
