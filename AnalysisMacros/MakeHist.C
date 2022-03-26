#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMath.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <Riostream.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TMultiGraph.h>
using namespace std;
#endif

void MakeHist(TString filein="APTS10_Vbb-4V_WaltConf_trgOR_20220311_TTree.root", TString fileout="APTS10_Fast_noGC.root")
{
// Project (from TTree) Amplitudes for the "fast" cluster size 1 events
// Histograms are saved in a Canvas onto a Root file
int nbins = 80; // numer of bins
float low = 0.060; // lower limit (V)
float high = 0.100; // higher limit (V)
TH1F* h1 = (TH1F*)(new TH1F("h1","h1",nbins,low,high));
TH1F* h2 = (TH1F*)(new TH1F("h2","h2",nbins,low,high));
TH1F* h3 = (TH1F*)(new TH1F("h3","h3",nbins,low,high));
TH1F* h4 = (TH1F*)(new TH1F("h4","h4",nbins,low,high));
TFile* f = new TFile(filein);
gDirectory->cd("Rint://");
TTree* myTree=(TTree*)f->Get("treeParams");
TCanvas* c1 = new TCanvas("c1","APTS10 Vbb -4V",0,0,1200,800);
 
myTree->Draw("SignalAmplCh1>>h1","SignalAmplCh1>0 && FallTimeCh1>0 && FallTimeCh1<400 && SignalAmplCh2<0 && SignalAmplCh3<0 && SignalAmplCh4<0","box");
myTree->Draw("SignalAmplCh2>>h2","SignalAmplCh2>0 && FallTimeCh2>0 && FallTimeCh2<400 && SignalAmplCh1<0 && SignalAmplCh3<0 && SignalAmplCh4<0","box");
myTree->Draw("SignalAmplCh3>>h3","SignalAmplCh3>0 && FallTimeCh3>0 && FallTimeCh3<400 && SignalAmplCh1<0 && SignalAmplCh2<0 && SignalAmplCh4<0","box");
myTree->Draw("SignalAmplCh4>>h4","SignalAmplCh4>0 && FallTimeCh4>0 && FallTimeCh4<400 && SignalAmplCh1<0 && SignalAmplCh2<0 && SignalAmplCh3<0","box");

h1->SetLineColor(kBlack);
h2->SetLineColor(kRed);
h3->SetLineColor(kBlue);
h4->SetLineColor(kGreen);

c1->cd();
c1->Draw();
h1->DrawClone("");
h2->DrawClone("SAME");
h3->DrawClone("SAME");
h4->DrawClone("SAME");
c1->Modified();
c1->Update();
 
c1->SaveAs(fileout);
}
