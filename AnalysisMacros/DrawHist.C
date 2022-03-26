#include <TROOT.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TCanvas.h>

void DrawHist(TString filename="APTS10_Fast.root", double myhismax=180.){

  gStyle->SetLegendTextSize(0.03);
  
  const char *namhist[4] = {"h1", "h2", "h3", "h4"}; // histogram names:
  //                         h1-h4 (no GC), h11-h14 (GC1)
  // cout << namhist[0] <<endl;
  
  // Open file with histograms stored in a canvas:
  
  TFile *file = new TFile(filename);
  file->ls();  
  file->GetListOfKeys()->Print();
  int nobj = file->GetNkeys();
  TList* objects = file->GetListOfKeys();
  
  for(int i=0; i<nobj; ++i) {
    TKey* obj = (TKey*) objects->At(i);
    TCanvas* c = (TCanvas*) file->Get(obj->GetName());
    c->GetListOfPrimitives()->Print();

    TH1F* gHfast[4];
    
    // Get individual histograms from the canvas:
    
    TH1F *h1;
    h1 = (TH1F*) c->GetPrimitive(namhist[0]);
    gHfast[0] = h1;
    h1->SetMaximum(myhismax);
    h1->SetLineWidth(2);
    h1->GetXaxis()->SetTitle("Amplitude (V)");
    h1->GetYaxis()->SetTitle("Counts");
    h1->Draw();

    TH1F *h2;
    h2 = (TH1F*) c->GetPrimitive(namhist[1]);
    gHfast[1] = h2;
    h2->SetMaximum(myhismax);
    h2->SetLineWidth(2);
    h2->SetFillColor(42);
    h2->Draw("SAME");

    TH1F *h3;
    h3 = (TH1F*) c->GetPrimitive(namhist[2]);
    gHfast[2] = h3;
    h3->SetMaximum(myhismax);
    h3->SetLineWidth(2);
    h3->Draw("SAME");

    TH1F *h4;
    h4 = (TH1F*) c->GetPrimitive(namhist[3]);
    gHfast[3] = h4;
    h4->SetMaximum(myhismax);
    h4->SetLineWidth(2);
    h4->Draw("SAME");

    // Add a legend:
    TLegend* leg = new TLegend(0.1,0.55,0.3,0.85);
    
    for(int k=0; k<4; k++) {
       leg->AddEntry(gHfast[k],Form("Pixel %d",k+1),"L")->SetTextColor(gHfast[k]->GetLineColor());
    }
    
    c->Draw();
    leg->Draw();
  }

}
