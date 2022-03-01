#include <TGraph.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TTimeStamp.h>
#include <TMath.h>
#include <TAxis.h>

void ConvertTxtToRoot(TString filname="APTS03_Vbb_0V_92pt.txt"){
  FILE* infil=fopen(filname,"r");
  char line[500];
  char channel[5];
  char strmonth[5];
  int iEv=-1;
  int iEvOld=-1;
  int jSampling;
  float timeBin;
  int day,month,year,hour,mins,sec,nanosec;
  float secfloat;
  TDirectory* curDir=gDirectory;
  TString outRootFilName=filname.ReplaceAll(".txt",".root");
  TFile* outRoot=new TFile(outRootFilName.Data(),"recreate");
  TTimeStamp tst;
  
  while(!feof(infil)){
    char* rc=fgets(line,500,infil);
    if(rc && strstr(line,"Trigger n.")) break;
  }
  while(!feof(infil)){
    sscanf(line,"Trigger n.   %d",&iEv);
    if(iEv!=iEvOld){
      printf("New Event: %d\n",iEv);
      iEvOld=iEv;
      char* rc=fgets(line,500,infil);
      if(!rc || !strstr(line,"TRIGGER_TIME")){
	printf("ERROR! Line with TRIGGER_TIME not found .... Exiting\n");
	return;
      }else{
	// add lines of code in case we want to extract info on date and time
	printf("Read time for event %d\n",iEv);
	sscanf(line,"\"TRIGGER_TIME       : Date = %3s %d, %d, Time = %d:%d:%f\"",strmonth,&day,&year,&hour,&mins,&secfloat);
	sec=(int)secfloat;
	nanosec=0;//(int)(secfloat*1e9);
	if(strcmp(strmonth,"FEB")==0) month=2;
	else if(strcmp(strmonth,"MAR")==0) month=3;
	tst.Set(year,month,day,hour,mins,sec,nanosec,kFALSE,0);
	printf("Event %d Time: %d%02d%02d %02d:%02d:%02d --> %ld \n",iEv,year,month,day,hour,mins,sec,tst.GetSec());	
      }
      int nWaveforms=0;
      while(1){
	char* rc=fgets(line,500,infil);
	if(feof(infil)) break;
	if(rc && strstr(line,"Waveform")){
	  sscanf(line,"Waveform with SP of %d, Horiz_interval of %f, Channel %2s:",&jSampling,&timeBin,channel);
	  timeBin=TMath::Nint(timeBin*1e12); // convert to an integer number of ps
	  printf(" Read Waveform for Event %d channel %s  time sampling = %f ps\n",iEv,channel,timeBin*jSampling);
	  ++nWaveforms;
	  char* rc=fgets(line,500,infil);
	  if(!rc || !strstr(line,"\"")){
	    printf("ERROR: wrong format, expecting \"");
	    return;
	  }
	  //  printf("start reading\n");
	  int nRead=0;
	  TGraph* g=new TGraph(0);
	  g->SetName(Form("grEv%dChan%ssamp%d",iEv,channel,TMath::Nint(timeBin*jSampling)));
	  g->SetTitle(Form("Event %d Channel %s Time %ld",iEv,channel,tst.GetSec()));
	  g->GetYaxis()->SetTitle("Amplitude");
	  g->GetXaxis()->SetTitle("Time (ps)");
	  while(1){
	    char* rc=fgets(line,500,infil);
	    if(rc && strstr(line,"\"")){
	      //	      printf("stop reading\n");
	      break;
	    }else{
	      while(1){
		char* remaining;
		double sig=strtod(line,&remaining);
		if(remaining == line) break;
		g->SetPoint(nRead,nRead*timeBin*jSampling,sig);
		nRead++;
		strcpy(line,remaining);
	      }
	    }
	  }
	  printf("   Read %d values\n",nRead);
	  outRoot->cd();
	  g->Write();
	  curDir->cd();
	}
	if(rc && strstr(line,"Trigger n.")) break;
      }
      printf(" --> Number of read waveforms = %d\n",nWaveforms);
      // possibly add a check that all waverforms were found?
    }
  }
  return;
}
  
