#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLatex.h>
#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>

TGraph* Smooth(TGraph* g, int nsm=2){
  TGraph* gsm=new TGraph(0);
  int npts=0;
  for(int j=nsm; j<g->GetN()-nsm; j++){
    double x,y,x1,y1;
    g->GetPoint(j,x,y);
    double sum=y;
    for(int k=1; k<=nsm;k++){
      g->GetPoint(j-k,x1,y1);
      sum+=y1;
      g->GetPoint(j+k,x1,y1);
      sum+=y1;
    }
    sum/=(2*nsm+1);
    gsm->SetPoint(npts++,x,sum);
  }
  return gsm;
}

double GetMaxX(TGraph* g){
  double xmax=0.;
  for(int j=0; j<g->GetN(); j++){
    double x,y;
    g->GetPoint(j,x,y);
    if(x>xmax) xmax=x;
  }
  return xmax;
}

double ComputeDerivative(TGraph* g, int j, int npts=5){
  if(npts==3){
    double xm1,ym1,xp1,yp1;
    g->GetPoint(j-1,xm1,ym1);
    g->GetPoint(j+1,xp1,yp1);
    double der=(yp1-ym1)/(xp1-xm1);
    return der;
  }else{
    double xm2,ym2,xm1,ym1,xp1,yp1,xp2,yp2;
    g->GetPoint(j-2,xm2,ym2);
    g->GetPoint(j-1,xm1,ym1);
    g->GetPoint(j+1,xp1,yp1);
    g->GetPoint(j+2,xp2,yp2);
    double der=(ym2-8*ym1+8*yp1-yp2)/(xp2-xp1)/12.;
    return der;
  }
}

TGraph* GetDerivative(TGraph* g, int nsm=4){
  TGraph* gder=new TGraph(0);
  int npts=0;
  for(int j=2; j<g->GetN()-nsm-2; j++){
    double x,y;
    g->GetPoint(j,x,y);
    double der=0;
    double nnn=0;
    for(int k=0; k<=nsm;k++){
      der+=ComputeDerivative(g,j+k);
      nnn+=1.;
    }
    if(nnn>0){
      der/=nnn;
      gder->SetPoint(npts++,x,der);
    }
  }
  return gder;
}

TGraph* CountNextNegativeDer(TGraph* g){
  TGraph* gn=new TGraph(0);
  for(int j=0; j<g->GetN()-1; j++){
    double x,y;
    g->GetPoint(j,x,y);
    int cntneg=0;
    for(Int_t k=j; k<g->GetN()-1; k++){
      double der=ComputeDerivative(g,k);
      if(der>=0) break;
      else cntneg++;
    }
    gn->SetPoint(j,x,cntneg);
  }
  return gn;
}

void GetMeanAndRMSCounts(TGraph* g, double xmin, double xmax, double& mean, double& rms){
  double sum=0,sum2=0,cnts=0;
  for(int j=0; j<g->GetN(); j++){
    double x,c;
    g->GetPoint(j,x,c);
    if(x>xmin && x<xmax){
      cnts+=1.;
      sum+=c;
      sum2+=(c*c);
    }
  }
  if(cnts>0){
    mean=sum/cnts;
    rms=TMath::Sqrt(sum2/cnts-mean*mean);
  }else{
    mean=0;
    rms=0;
  }
  return;
}
double FindOnGraph(TGraph* gcount, double y, double xmin, double xmax, int interpolate, bool backw=kFALSE){
  int jfirst=0;
  int dstep=1;
  if(backw){
    jfirst=gcount->GetN();
    dstep=-1;
  }
  for(int jstep=0; jstep<gcount->GetN(); jstep++){
    int j=jfirst+dstep*jstep;
    double x,c,xbef,cbef,xaft,caft,xaft2,caft2;
    gcount->GetPoint(j,x,c);
    gcount->GetPoint(j-dstep,xbef,cbef);
    gcount->GetPoint(j+dstep,xaft,caft);
    gcount->GetPoint(j+2*dstep,xaft2,caft2);
    if((dstep==1 && c<y && cbef>y && caft<y) || (dstep==-1 && c>y && cbef<y && caft>y) ){
      if(interpolate==0) return x;
      else{
	double sumx=0,sumx2=0,sumy=0,sumxy=0,npts=0;
	for(int k=j-interpolate; k<=j+interpolate; k++){
	  double xP,yP;
	  gcount->GetPoint(k,xP,yP);
	  sumx+=xP;
	  sumy+=yP;
	  sumxy+=(xP*yP);
	  sumx2+=(xP*xP);
	  npts+=1;
	}
	double m=(npts*sumxy-sumx*sumy)/(npts*sumx2-sumx*sumx);
	double q=(sumy*sumx2-sumx*sumxy)/(npts*sumx2-sumx*sumx);
	double xinterp=(y-q)/m;
	if(xinterp<xmin || xinterp>xmax || TMath::Abs(xinterp-x)>1000.) continue;
	return xinterp;
      }
    }
  }
  return -999.;
}



void FindEdge(TGraph* gcount,TGraph* gnegder, TGraph* gder, double& endplateau, double& edgeleft, double& edgeright){
  
  // first very rough: compute flat levels on the left and on the right and check their difference
  double maxTime=GetMaxX(gcount);
  double levleft,rmsleft,levright,rmsright;
  GetMeanAndRMSCounts(gcount,0.,2000.,levleft,rmsleft);
  GetMeanAndRMSCounts(gcount,maxTime-2000,maxTime,levright,rmsright);
  double y50=0.5*(levleft+levright);
  double t50fromleft=FindOnGraph(gcount,y50,0.,maxTime,4);
  double t50fromright=FindOnGraph(gcount,y50,0.,maxTime,4,kTRUE);
  double roughsig=levleft-levright;
  printf("Rough signal = %f Rough edge position = %f %f\n",roughsig,t50fromleft,t50fromright);
  double minSearchWindow=0;
  double maxSearchWindow=maxTime;
  if(roughsig>0.0005){
    minSearchWindow=TMath::Min(t50fromleft,t50fromright)-6000.;
    if(minSearchWindow<0) minSearchWindow=0;
    maxSearchWindow=TMath::Max(t50fromleft,t50fromright)+6000.;
    if(maxSearchWindow>maxTime) maxSearchWindow=maxTime;
  }
  printf("Search window = %f %f\n",minSearchWindow,maxSearchWindow);
  
  // second step: search for accumulation of adjacent points with negative derivative
  double xmaxn=-1;
  double cmaxn=-1;
  int jmaxn=-1;
  if(gnegder){
    for(int j=0; j<gnegder->GetN(); j++){
      double x,c;
      gnegder->GetPoint(j,x,c);
      if(x<minSearchWindow || x>maxSearchWindow) continue;
      if(c>cmaxn){
	cmaxn=c;
	xmaxn=x;
	jmaxn=j;
      }
      if(c==cmaxn){
	int sum0=0;
	int sum1=0;
	for(int k=1; k<20; k++){
	  double xk,ck;
	  gnegder->GetPoint(jmaxn+k,xk,ck);
	  sum0+=ck;
	  gnegder->GetPoint(j+k,xk,ck);
	  sum1+=ck;
	}
	if(sum1>sum0){
	  cmaxn=c;
	  xmaxn=x;
	  jmaxn=j;
	}
      }
    }
    printf("Maximum adjacent points with negative derivative: t_maxn=%f   n_neg=%f\n",xmaxn,cmaxn);
  }
  
  // third step: search for minimum of derivative and range where derivative differs from 0
  double xminder=-1;
  double dermin=99999.;
  int jminder=-1;
  for(int j=0; j<gder->GetN(); j++){
    double x,d;
    gder->GetPoint(j,x,d);
    if(x<minSearchWindow || x>maxSearchWindow) continue;
    if(d<dermin){
      dermin=d;
      xminder=x;
      jminder=j;
    }
  }
  if(jminder<0){
    endplateau=0;
    edgeleft=0;
    edgeright=0;
    return;
  }
  printf("Minimum of derivative: xminder=%f   dermin=%f\n",xminder,dermin);
  int jleft=-1;
  double dthresh=-1e-7;
  for(int j=jminder; j>0; j--){
    double x,d;
    gder->GetPoint(j,x,d);
    if(d>dthresh){
      jleft=j;
      break;
    }
  }
  int jright=-1;
  for(int j=jminder; j<gder->GetN(); j++){
    double x,d;
    gder->GetPoint(j,x,d);
    if(d>dthresh){
      jright=j;
      break;
    }
  }
  double xleft,xright,dum;
  gder->GetPoint(jleft,xleft,dum);
  gder->GetPoint(jright,xright,dum);
  printf("Region of negative derivative: xleft=%f   xright=%f\n",xleft,xright);
  if(xmaxn>0 && TMath::Abs(xmaxn-xleft)<5000 && TMath::Abs(xmaxn-xright)<5000){
    if(xleft>xmaxn) xleft=xmaxn;
    if(xright<xmaxn) xright=xmaxn;
  }
  printf("Edge range after analysis of derivative: xleft=%f   xright=%f\n",xleft,xright);

  // Fourth step: start from left and seach for N points with couns < baseline-3sigma
  double cmean,crms;
  GetMeanAndRMSCounts(gcount,0.,xleft,cmean,crms);
  printf("Mean before edge = %f rms = %f\n",cmean,crms);
  double thresh=cmean-3*crms;
  int threshbin=TMath::Nint(TMath::Max(10.,cmaxn/3.));
  if(cmaxn<0) threshbin=3;
  double xleft2=-1.;
  for(int j=0; j<gcount->GetN(); j++){
    double x,c;
    gcount->GetPoint(j,x,c);
    int nbelow=0;
    for(Int_t k=j; k<gcount->GetN()-1; k++){
      double x2,c2;
      gcount->GetPoint(k,x2,c2);
      if(c2<thresh) nbelow++;
      else break;
    }
    if(nbelow>threshbin){
      xleft2=x;
      break;
    }
  }
  printf("Left Edge from baseline-N*rms = %f\n",xleft2);
  if(xleft2>0){
    endplateau=TMath::Min(xleft,xleft2);
    edgeleft=TMath::Max(xleft,xleft2);
    edgeright=xright;
    printf("Edge range after all steps: endplateau=%f   edgeleft=%f   edgeright=%f\n",endplateau,edgeleft,edgeright);
    return;
  }
  endplateau=0;
  edgeleft=0;
  edgeright=0;
  return;
}


void ProcessEvent(TGraph* g, TGraph* glong, double params[20], bool plot){

  TGraph* gs=Smooth(g,10);
  gs->GetXaxis()->SetTitle(g->GetXaxis()->GetTitle());
  gs->SetTitle(g->GetTitle());
  gs->GetYaxis()->SetTitle("Amplitude (smoothened)");
  TGraph* gnegd=CountNextNegativeDer(g);
  gnegd->GetXaxis()->SetTitle(g->GetXaxis()->GetTitle());
  gnegd->SetTitle(g->GetTitle());
  gnegd->GetYaxis()->SetTitle("N. adjacent samplings with negative derivative");
  TGraph* gsnegd=CountNextNegativeDer(gs);
  TGraph* gd=GetDerivative(gs,40);
  gd->GetXaxis()->SetTitle(g->GetXaxis()->GetTitle());
  gd->SetTitle(g->GetTitle());
  gd->GetYaxis()->SetTitle("Amplitude derivative (smoothened)");
  double endpl,edgeleft,edgeright;
  printf("--- Edge finding on fine graph ---\n");
  FindEdge(g,gnegd,gd,endpl,edgeleft,edgeright);
  TF1* fbas=new TF1("fbas","[0]");
  TF1* flow=new TF1("flow","[0]");
  fbas->SetLineColor(2);
  flow->SetLineColor(4);
  double maxTime=GetMaxX(g);
  if(endpl<0.0001) endpl=maxTime;
  
  printf("--- Edge parameters ---\n");
  g->Fit(flow,"Q+","",edgeright,maxTime);
  g->Fit(fbas,"Q+","",0.,endpl);
  double basel=fbas->GetParameter(0);
  double baselsig=fbas->GetParError(0);
  double baselc,baselsigc;
  GetMeanAndRMSCounts(g,0.,endpl,baselc,baselsigc);
  printf("Baseline from fit = %f+-%f    from counts=%f rms=%f\n",basel,baselsig,baselc,baselsigc);
  double sign=flow->GetParameter(0);
  double signsig=flow->GetParError(0);
  double signc,signsigc;
  GetMeanAndRMSCounts(g,edgeright,25000,signc,signsigc);
  printf("Signal level from fit = %f+-%f    from counts=%f rms %f\n",sign,signsig,signc,signsigc);

  bool isSignal=kFALSE;
  if(edgeright>edgeleft && TMath::Abs(sign-basel)>0.001) isSignal=kTRUE;
  printf("IsSignal = %d  edgeleft=%f  edgeright=%f\n",isSignal,edgeleft,edgeright);
  double amplitude=-999.;
  double t10=-999.;
  double t50=-999.;
  double t90=-999.;
  if(isSignal){
    amplitude=TMath::Abs(sign-basel);
    t10=FindOnGraph(g,basel-0.1*amplitude,0.,maxTime,4,kTRUE);
    t50=FindOnGraph(g,basel-0.5*amplitude,0.,maxTime,4);
    t90=FindOnGraph(g,basel-0.9*amplitude,0.,maxTime,4);
  }
  printf("t10 = %f  t50=%f   t90 = %f\n",t10,t50,t90);
  
  params[0]=basel;
  params[1]=sign;
  params[2]=amplitude;
  params[3]=t10;
  params[4]=t50;
  params[5]=t90;
  params[6]=t90-t10;

  TGraph* gdlong=GetDerivative(glong,2);
  double MaxTimeLong=GetMaxX(glong);
  gdlong->GetXaxis()->SetTitle(g->GetXaxis()->GetTitle());
  gdlong->SetTitle(g->GetTitle());
  gdlong->GetYaxis()->SetTitle("Amplitude derivative");
  double endpl2,edgeleft2,edgeright2;
  printf("--- Edge finding on coarse graph ---\n");
  FindEdge(glong,0x0,gdlong,endpl2,edgeleft2,edgeright2);
  if(!isSignal){
    endpl2=MaxTimeLong;
    edgeleft2=0;
    edgeright2=0;
  }
  TF1* fbas2=new TF1("fbas2","[0]");
  fbas2->SetLineColor(2);
  glong->Fit(fbas2,"Q+","",0.,endpl2);
  double basel2=fbas2->GetParameter(0);
  TF1* fexp2=new TF1("fexp2","[0]+[1]*TMath::Exp(-x/[2])");
  fexp2->SetParameters(-basel2,-basel2,1000000.);
  fexp2->SetLineColor(kGreen+2);
  double expConst=-999;
  if(isSignal){
    glong->Fit(fexp2,"Q+","",edgeright2,MaxTimeLong);
    expConst=fexp2->GetParameter(2);
  }
  TF1* flin2=new TF1("flin2","[0]+[1]*x",edgeright2,MaxTimeLong);
  flin2->SetLineColor(kMagenta);
  flin2->SetLineWidth(1);
  double xcross=-999;
  if(isSignal){  
    glong->Fit(flin2,"Q+","",edgeright2,edgeright2+expConst);
    double slope=flin2->GetParameter(1);
    double q=flin2->GetParameter(0);
    xcross=(basel2-q)/slope;
  }

  params[7]=basel2;
  params[8]=expConst;
  params[9]=xcross;
  
  if(!plot){
    delete gs;
    delete gnegd;
    delete gsnegd;
    delete gd;
    delete gdlong;
    delete fbas;
    delete flow;
    delete fexp2;
    delete flin2;
    return;
  }
  
  TCanvas* c1= new TCanvas("c1","",1600,800);
  c1->Divide(3,2);
  c1->cd(1);
  gPad->SetRightMargin(0.04);
  g->Draw();
  flow->Draw("same");
  fbas->Draw("same");
  TLine* lp=new TLine(endpl,g->GetYaxis()->GetXmin(),endpl,g->GetYaxis()->GetXmax());
  lp->SetLineColor(2);
  lp->SetLineStyle(2);
  lp->Draw();
  TLine* ll=new TLine(edgeleft,g->GetYaxis()->GetXmin(),edgeleft,g->GetYaxis()->GetXmax());
  ll->SetLineColor(kGreen+1);
  ll->SetLineStyle(9);
  ll->Draw();
  TLine* lr=new TLine(edgeright,g->GetYaxis()->GetXmin(),edgeright,g->GetYaxis()->GetXmax());
  lr->SetLineColor(kGreen+1);
  lr->SetLineStyle(9);
  lr->Draw();
  TLine* l10=new TLine(t10,g->GetYaxis()->GetXmin(),t10,g->GetYaxis()->GetXmax());
  l10->SetLineColor(kMagenta+1);
  l10->SetLineStyle(7);
  l10->Draw();
  TLine* l90=new TLine(t90,g->GetYaxis()->GetXmin(),t90,g->GetYaxis()->GetXmax());
  l90->SetLineColor(kMagenta+1);
  l90->SetLineStyle(7);
  l90->Draw();
  TLatex* t1=new TLatex(edgeright+1400,basel,Form("Baseline = %.2f mV\n",basel*1000.));
  t1->SetTextColor(2);
  t1->Draw();
  TLatex* t2=new TLatex(edgeright+1400,basel-0.01*basel,Form("Signal ampl = %.2f mV\n",amplitude*1000.));
  t2->SetTextColor(4);
  t2->Draw();
  TLatex* t3=new TLatex(edgeright+1400,basel-0.02*basel,Form("Fall time = %.1f ps\n",t90-t10));
  t3->SetTextColor(kMagenta+1);
  t3->Draw();
  c1->cd(2);
  gPad->SetRightMargin(0.04);
  gs->Draw();
  c1->cd(4);
  gPad->SetRightMargin(0.04);
  gnegd->Draw();
  c1->cd(5);
  gPad->SetRightMargin(0.04);
  gd->Draw();
  c1->cd(3);
  gPad->SetRightMargin(0.04);
  glong->Draw();
  TLine* lp2=new TLine(endpl2,glong->GetYaxis()->GetXmin(),endpl2,glong->GetYaxis()->GetXmax());
  lp2->SetLineColor(2);
  lp2->SetLineStyle(2);
  lp2->Draw();
  TLine* lr2=new TLine(edgeright2,glong->GetYaxis()->GetXmin(),edgeright2,glong->GetYaxis()->GetXmax());
  lr2->SetLineColor(kGreen+1);
  lr2->SetLineStyle(9);
  lr2->Draw();
  TLine* b3=new TLine(endpl2,basel2,xcross,basel2);
  b3->SetLineColor(kRed-7);
  b3->SetLineStyle(7);
  b3->Draw();
  TLine* lcr=new TLine(xcross,glong->GetYaxis()->GetXmin(),xcross,basel2);
  lcr->SetLineColor(kMagenta);
  lcr->SetLineStyle(7);
  lcr->Draw();
  fbas2->Draw("same");
  if(isSignal){
    flin2->Draw("same");
    fexp2->Draw("same");
    TLatex* tlin=new TLatex(edgeright2+1e6,flin2->Eval(edgeright2),Form("Recovery time (line) = %.2f ns",(xcross-edgeright2)/1000.));
    tlin->SetTextColor(kMagenta);
    tlin->Draw();
    TLatex* texp=new TLatex(edgeright2+1e6,flin2->Eval(edgeright2)-amplitude*0.1,Form("Recovery time (expo) = %.2f ns",fexp2->GetParameter(2)/1000.));
    texp->SetTextColor(kGreen+2);
    texp->Draw();
  }
  c1->cd(6);
  gPad->SetRightMargin(0.04);
  gdlong->Draw();

}

void ProcessEvent(TString filnam, int ev, int chan){
  if(chan<1 || chan >4){
    printf("ERROR: channel number should be 1, 2, 3 or 4\n");
    return;
  }
  TFile* f=new TFile(filnam.Data());
  double params[10];
  TGraph* g=(TGraph*)f->Get(Form("grEv%dChanC%dsamp25",ev,chan));
  TGraph* glong=(TGraph*)f->Get(Form("grEv%dChanC%dsamp10000",ev,chan));
  if(g==0x0){
    printf("TGraph %s not found in root file\n",Form("grEv%dChanC%dsamp25",ev,chan));
    return;
  }
  if(glong==0x0){
    printf("TGraph %s not found in root file\n",Form("grEv%dChanC%dsamp10000",ev,chan));
    return;
  }
  ProcessEvent(g,glong,params,kTRUE);
}

void ProcessFile(TString filnam){

  const int nParsPerChan=10;
  TString varNames[nParsPerChan]={"Baseline","MinLevel","SignalAmpl","t10","t50","t90",
				  "FallTime","BaselCoarse","RecovTimeExpo","RecovTimeLin"};
  double paramschan[20];
  double params[4*nParsPerChan];
  int ev;
  ulong timest;

  TString outfilnam=filnam.Data();
  outfilnam.ReplaceAll(".root","_TTree.root");
  TFile* outFile=new TFile(outfilnam.Data(),"recreate");
  TTree* outTree=new TTree("treeParams","tree of parameters");
  outTree->Branch("Event",&ev,"Event/I");
  outTree->Branch("Timestamp",&timest,"Timestamp/l");
  for(int ichan=0; ichan<4; ichan++){
    for(int ivar=0; ivar<nParsPerChan; ivar++){
      outTree->Branch(Form("%sCh%d",varNames[ivar].Data(),ichan+1),&params[ivar+ichan*nParsPerChan],Form("%sCh%d/D",varNames[ivar].Data(),ichan+1));
    }
  }

  TFile* f=new TFile(filnam.Data());
  int totEv=0;
  int nkeys=f->GetNkeys();
  TList* lkeys=f->GetListOfKeys();
  int period;
  int lastev=0;
  int chan;
  for(Int_t j=0; j<nkeys; j++){
    TKey* k=(TKey*)lkeys->At(j);
    TString cname=k->GetClassName();
    TString oname=k->GetName();
    if(cname=="TGraph"){
      sscanf(oname.Data(),"grEv%dChanC%dsamp%d",&ev,&chan,&period);
      if(ev>lastev) lastev=ev;
    }
  }

  int dum;
  char ddum[2];
  for(int iev=0; iev<=lastev; iev++){
    printf("----- Event %d -----\n",iev);
    for(int ichan=1; ichan<=4; ichan++){
      printf("--- Channel %d ---\n",ichan);
      TGraph* g=(TGraph*)f->Get(Form("grEv%dChanC%dsamp25",iev,ichan));
      TGraph* glong=(TGraph*)f->Get(Form("grEv%dChanC%dsamp10000",iev,ichan));
      if(!g ||!glong) continue;
      const char* grtit=g->GetTitle();
      sscanf(grtit,"Event %d Channel %2s Time %ld",&dum,ddum,&timest);
      ev=iev;
      ProcessEvent(g,glong,paramschan,kFALSE);
      for(int ivar=0; ivar<nParsPerChan; ivar++) params[ivar+(ichan-1)*nParsPerChan]=paramschan[ivar];
    }
    outTree->Fill();
  }
  //  outTree->Scan();
  
  outFile->cd();
  outTree->Write();
  outFile->Close();
  
}
