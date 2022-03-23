#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TCanvas.h>


void DoFit(TH1F* h, TF1* funcKa, TF1* funcKb){
  double roughmean=h->GetBinCenter(h->GetMaximumBin());
  double roughpeak=h->GetBinContent(h->GetMaximumBin());
  double roughbase=h->GetBinContent(h->GetXaxis()->FindBin(65.));
  funcKa->SetParameters(roughbase,roughpeak,roughmean,2);
  h->Fit(funcKa);
  int peakb=-1;
  double cntpeakb=0;
  for(int jb=1; jb<h->GetNbinsX(); jb++){
    double x=h->GetBinCenter(jb);
    double y=h->GetBinContent(jb);
    if(x>funcKa->GetParameter(2)+3*funcKa->GetParameter(3) && y>cntpeakb){
      peakb=jb;
      cntpeakb=y;
    }
  }
  funcKb->SetParameters(h->GetBinContent(peakb),h->GetBinCenter(peakb),funcKa->GetParameter(3));
  h->Fit(funcKb,"","",funcKa->GetParameter(2)+3*funcKa->GetParameter(3),100);
  return;
}

void ComputeLinearParams(TGraphErrors* g, double& slop, double& cnst){
  double x1,y1,x2,y2;
  g->GetPoint(0,x1,y1);
  g->GetPoint(1,x2,y2);
  slop=(y2-y1)/(x2-x1);
  cnst=y1-slop*x1;
  return;
}

void FitKab(TString filename="APTS10_Vbb-4V_WaltConf_trgOR_20220311_TTree.root", double maxFallTime=400){

  double eA=5.89875;
  double eB=6.49045;
  
  TH1F* hAmplChan[4];
  TH1F* hAmplChanCluSiz1[4];
  TH1F* hAmplChanFastCluSiz1[4];

  TFile* f=new TFile(filename.Data());
  TTree* tree=(TTree*)f->Get("treeParams");
  double amplVec[4];
  double fallTimeVec[4];
  for(int k=0; k<4; k++){
    tree->SetBranchAddress(Form("FallTimeCh%d",k+1),&fallTimeVec[k]);
    tree->SetBranchAddress(Form("SignalAmplCh%d",k+1),&amplVec[k]);
    hAmplChan[k]=new TH1F(Form("hAmplChan%d",k+1),Form(" All clusters ; Signal Amplitude Ch%d (mV) ; Entries",k+1),80,60.,100.);
    hAmplChanCluSiz1[k]=new TH1F(Form("hAmplChanCluSiz1%d",k+1),Form(" Cluster Size = 1 ; Signal Amplitude Ch%d (mV) ; Entries",k+1),80,60.,100.);
    hAmplChanFastCluSiz1[k]=new TH1F(Form("hAmplChanFastCluSiz1%d",k+1),Form("Cluster Size = 1, FallTime < %.0f; Signal Amplitude Ch%d (mV) ; Entries",maxFallTime,k+1),80,60.,100.);
  }
  for(int ient=0; ient<tree->GetEntriesFast(); ient++){
    tree->GetEvent(ient);
    int clusiz=0;
    double maxsig=-999.;
    int maxsigpix=-1;
    for(int k=0; k<4; k++){
      if(amplVec[k]>0){
  	clusiz+=1;
	hAmplChan[k]->Fill(amplVec[k]*1000.);
	if(amplVec[k]>maxsig){
	  maxsig=amplVec[k];
	  maxsigpix=k;
	}
      }
    }
    if(clusiz==1){
      hAmplChanCluSiz1[maxsigpix]->Fill(amplVec[maxsigpix]*1000.);
      if(fallTimeVec[maxsigpix]<maxFallTime) hAmplChanFastCluSiz1[maxsigpix]->Fill(amplVec[maxsigpix]*1000.);
    }
  }

  TF1* funcKa=new TF1("funcKa","[0]*(1-TMath::Erf((x-[2])/[3]))+[1]*TMath::Exp(-0.5*(x-[2])*(x-[2])/[3]/[3])",60,90.);
  TF1* funcKb=new TF1("funcKb","[0]*TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])",80,100.);
  funcKa->SetLineWidth(2);
  funcKa->SetLineColor(kRed+1);
  funcKb->SetLineWidth(2);
  funcKb->SetLineColor(4);

  TGraphErrors* gAall[4];
  TGraphErrors* gAclu1[4];
  TGraphErrors* gAfast[4];
  int cols[4]={1,kRed+1,kGreen+2,kBlue};
  for(int k=0; k<4; k++){
    gAall[k]=new TGraphErrors(0);
    gAclu1[k]=new TGraphErrors(0);
    gAfast[k]=new TGraphErrors(0);
    gAall[k]->SetMarkerColor(cols[k]);
    gAall[k]->SetLineColor(cols[k]);
    gAall[k]->SetMarkerStyle(25);
    gAall[k]->SetMarkerSize(1.2);
    gAclu1[k]->SetMarkerColor(cols[k]);
    gAclu1[k]->SetLineColor(cols[k]);
    gAclu1[k]->SetMarkerStyle(20);
    gAfast[k]->SetMarkerColor(cols[k]);
    gAfast[k]->SetLineColor(cols[k]);
    gAfast[k]->SetMarkerStyle(47);
  }
  
  double slopeAll[4];  
  double constAll[4];
  double slopeClu1[4];  
  double constClu1[4];
  double slopeFast[4];  
  double constFast[4];
  
  TCanvas* cA = new TCanvas("cA","Amplitudes",1650,900);
  cA->Divide(2,2);
  for(int k=0; k<4; k++){
    cA->cd(k+1);
    DoFit(hAmplChan[k],funcKa,funcKb);
    hAmplChan[k]->Draw();
    funcKa->DrawClone("same");
    funcKb->DrawClone("same");
    gAall[k]->SetPoint(0,eA,funcKa->GetParameter(2));
    gAall[k]->SetPointError(0,0.,funcKa->GetParError(2));
    gAall[k]->SetPoint(1,eB,funcKb->GetParameter(1));
    gAall[k]->SetPointError(1,0.,funcKb->GetParError(1));
    ComputeLinearParams(gAall[k],slopeAll[k],constAll[k]);
  }
  cA->SaveAs("Fit-Kab-AllClu.png");
  
  TCanvas* c1 = new TCanvas("c1","Amplitudes CluSiz1",1650,900);
  c1->Divide(2,2);
  for(int k=0; k<4; k++){
    c1->cd(k+1);
    DoFit(hAmplChanCluSiz1[k],funcKa,funcKb);
    hAmplChanCluSiz1[k]->Draw();
    funcKa->DrawClone("same");
    funcKb->DrawClone("same");
    gAclu1[k]->SetPoint(0,eA,funcKa->GetParameter(2));
    gAclu1[k]->SetPointError(0,0.,funcKa->GetParError(2));
    gAclu1[k]->SetPoint(1,eB,funcKb->GetParameter(1));
    gAclu1[k]->SetPointError(1,0.,funcKb->GetParError(1));
    ComputeLinearParams(gAclu1[k],slopeClu1[k],constClu1[k]);
  }
  c1->SaveAs("Fit-Kab-CluSiz1.png");

  TCanvas* cF = new TCanvas("cF","Amplitudes Fast CluSiz1",1650,900);
  cF->Divide(2,2);
  for(int k=0; k<4; k++){
    cF->cd(k+1);
    DoFit(hAmplChanFastCluSiz1[k],funcKa,funcKb);
    hAmplChanFastCluSiz1[k]->Draw();
    funcKa->DrawClone("same");
    funcKb->DrawClone("same");
    gAfast[k]->SetPoint(0,eA,funcKa->GetParameter(2));
    gAfast[k]->SetPointError(0,0.,funcKa->GetParError(2));
    gAfast[k]->SetPoint(1,eB,funcKb->GetParameter(1));
    gAfast[k]->SetPointError(1,0.,funcKb->GetParError(1));
    ComputeLinearParams(gAfast[k],slopeFast[k],constFast[k]);
  }
  cF->SaveAs("Fit-Kab-Fast.png");

  TH2F* hFrame=new TH2F("hFrame"," ; Energy (keV) ; Signal (mV)",100,5.,7.5,100,70.,95.);
  hFrame->SetStats(0);

  TCanvas* cen=new TCanvas("cen","",800,600);
  hFrame->Draw();
  TLegend* leg=new TLegend(0.65,0.2,0.89,0.8);
  for(int k=0; k<4; k++){
    gAall[k]->Draw("PSAME");
    gAclu1[k]->Draw("PSAME");
    gAfast[k]->Draw("PSAME");
    leg->AddEntry(gAall[k],Form("Pixel %d all clu",k+1),"P")->SetTextColor(gAall[k]->GetMarkerColor());
    leg->AddEntry(gAclu1[k],Form("Pixel %d clusiz=1",k+1),"P")->SetTextColor(gAclu1[k]->GetMarkerColor());
    leg->AddEntry(gAfast[k],Form("Pixel %d fast",k+1),"P")->SetTextColor(gAfast[k]->GetMarkerColor());
  }
  leg->Draw();
  cen->SaveAs("PeakCalib.png");

  printf("Slopes (mV/keV)\n");
  for(int k=0; k<4; k++){
    printf("Pixel %d:  all clu=%f    clusiz1=%f    fast=%f\n",k+1,slopeAll[k],slopeClu1[k],slopeFast[k]);
  }
}
