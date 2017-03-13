/* Program: example.cxx
 *       Created  by Ken Teh,            Aug. 2005
 */
 
// Header Files
using namespace std; //used to eliminate deprecated header file error message

//  NOTE - the advantage of having useful histogram names is the ability to find them in the code!!!


#include <sys/stat.h>
#include <cstdio>
#include <iostream>
#include "daphuserfunc.h"
#include "ScarletEvnt.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TRandom.h"
#include "TMath.h"
#include "TDirectory.h"
#include <fstream>

#define NSCALERS 12

//TFile *f; //used to create ROOT file

Float_t totals[NSCALERS];
Int_t stopped;

// Declaration of Histograms

/* 1-D histograms: */
//TH1 *hE[24];
TH1 *hELUM[6];

/* 2-D histograms */
TH2F *hADC[6];
TH2F *hTDC;
TH2F *hXFXN[24];
TH2F *hRDT[4];
TH2F *hEDE0;
TH2F *hDE0_RF;
TH2F *hELUM_RF[6];

/* function to count the number of set bits in a 16 bit word */
Int_t cntbit(Int_t word)
{
  Int_t nbits=0;
  for (Int_t ibit=0; ibit<16; ibit++) {
    if (word & (Int_t) TMath::Power(2,ibit)) {nbits++;}
  }
  return nbits;
}

/* The userentry() function:  
 */
int userentry()
{
/* stopped flag for Elliot's scaler program */
  for (Int_t i=0; i<NSCALERS; i++) totals[i]=0;
  stopped = 1; 

 //Open ROOT file
  //f = new TFile("H007_online.root", "recreate");
 // 1d - 2d histograms
 hTDC=new TH2F("hTDC","hTDC",512,0,4096,17,0,17);
 //E0 DE0
 hEDE0=new TH2F("hEDE0","hEDE0",512,0,4096,512,0,4096);
 hDE0_RF=new TH2F("hDE0_RF","hDE0_RF",512,0,4096,512,0,4096);

 for(int a=0;a<7;++a){
   TString name="hADC";
   TString title="raw ADC";
   name+=(a+1);
   title+=(a+1);
   hADC[a]=new TH2F(name,title,1024,0,4096,17,0,17);
 }
 
 for(int a=0;a<24;++a){
   TString name="hXFXN";
   TString title="XF vs. XN detector ";
   name+=(a+1);
   title+=(a+1);
   hXFXN[a]=new TH2F(name,title,512,0,4096,512,0,4096);
 

   //   TString nameE="hE";
   //  TString titleE="hE";
   //   nameE+=(a+1);
   //   titleE+=(a+1);
   //   hE[a]=new TH1F(nameE,titleE,1024,0,4096);
 }

for(int a=0;a<4;++a){
   TString name="hRDT";
   TString title="raw RDT";
   name+=(a+1);
   title+=(a+1);
   hRDT[a]=new TH2F(name,title,512,0,4096,512,0,4096);
 }
 
 for(int a=0;a<6;++a){
   TString name="hELUM";
   TString title="raw ELUM";
   name+=(a+1);
   title+=(a+1);
   hELUM[a]=new TH1F(name,title,1024,0,4096);

   TString name1="hELUM_RF";
   TString title1="raw ELUM_RF";
   name1+=(a+1);
   title1+=(a+1);
   hELUM_RF[a]=new TH2F(name1,title1,512,0,4096,512,0,4096);
 }
 
 return 0;
}

/* function to deal with scalers, adapted from Elliot's program */
void scalers(ScarletEvnt &e)
{
    /* Adapted from Kanter's scaler program
     */
    unsigned int *p = reinterpret_cast<unsigned int*>(e.body());
    unsigned int ttotal, tdiff, ithscaler, ithrate;
    FILE *sf;

    /* On the first sync event after a stop, the following tests for the
     * existence of a file called scalers.zap.  If it exists, the scaler
     * totals are reset.
     */
    if (stopped) {
        struct stat st;
        stopped = 0;
        if (stat("scalers.zap", &st) == 0) {
            for (int i = 0; i < NSCALERS; ++i) totals[i] = 0.0;
        }
    }


    if ((sf = fopen("scalers.dat", "w")) == 0) return;
    ttotal = *p++;
    tdiff = *p++;
    fprintf(sf, "%u %u\n", ttotal, tdiff);
    for (int i = 0; i < NSCALERS; ++i) {
        ithscaler = *p++ & 0x00ffffff;
        totals[i] += ithscaler;

        ithrate = tdiff != 0 ? ithscaler/tdiff : 0;
        fprintf(sf, "%.0f %u\n", totals[i], ithrate);
    }
    fclose(sf);
}


int userdecode(ScarletEvnt &event) {
  
  ScarletEvnt subevent1;
  Int_t dataword;

  subevent1=event[1];
  int *p1 = reinterpret_cast<int*>(subevent1.body());

  /*  The data are parsed into Channel ID and Channel Data and then output 
    into a raw array scheme (5x16) and then mapped to an array scheme (24x3).  
  */

  Int_t nhits[5]={0,0,0,0,0};
  Int_t Chan[5][16];
  Int_t RawData[5][16];
  Int_t Data[24][3];
  Int_t RawAux[16];
  Int_t RawTDC[16];

  Int_t MapDet[5][16]={{ 1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2},
  		       {11,10, 9, 8, 7, 6, 5, 4, 7, 6,11,10, 9, 8, 7, 6},
 		       {15,14,13,12,11,10, 9, 8,17,16,15,14,13,12,17,16},
                       {20,21,22,23,18,19,20,21,12,13,14,15,16,17,18,19},
                       {-1,-1,-1,-1,-1,-1,-1,-1,22,23,18,19,20,21,22,23}}; //for ADC5 only channels 0-7 
                                                                           //are used for Si signals

  Int_t MapSig[5][16]={{1,1,0,0,0,0,0,0,2,2,2,2,1,1,1,1}, //0->E, 1->XF, 2->XN
                       {0,0,0,0,0,0,2,2,2,2,1,1,1,1,1,1},
                       {0,0,0,0,2,2,2,2,1,1,1,1,1,1,0,0},
                       {0,0,0,0,1,1,1,1,2,2,2,2,2,2,0,0},
		       {-1,-1,-1,-1,-1,-1,-1,-1,1,1,2,2,2,2,2,2}};

  for(Int_t ii=0; ii<24; ii++)
    {
      Data[ii][0]=0;
      Data[ii][1]=0;
      Data[ii][2]=0;
    }
    
  for(int i=0;i<16;++i) {
    RawAux[i]=0;
    RawTDC[i]=0;
  }
  //Beginning Event Unpacking...   

  //Read In Aux Detectors
  // ch - 1 de0
  // ch - 2 e0
  // ch - 3 elum1
  // ch - 4 elum2
  // ch - 5 elum3
  // ch - 6 elum4
  // ch - 7 elum5
  // ch - 8 elum6  
  // ch - 9 recoil de1
  // ch - 10 recoil e1
  // ch - 11 recoil de2
  // ch - 12 recoil e2
  // ch - 13 recoil de3
  // ch - 14 recoil e3
  // ch - 15 recoil de4
  // ch - 16 recoil e4
  for(int i=0;i<16;++i) {
    dataword=*p1++;
    RawAux[i]=(dataword & 0x00000fff);
    hADC[5]->Fill((dataword & 0x00000fff),i); 
  }

  //Read in TDC
  // DE0 - RF
  // ELUM - RF
  // RDT - RF
  // ARRAY - RF

  for(int z=0;z<6;++z) { //read in TDC data
    dataword=*p1++;
    RawTDC[z]=(dataword & 0x00000fff);
    hTDC->Fill(RawTDC[z],z);
  }

 // Fill Aux histograms
  for(int i=0;i<6;i++) {
    hELUM[i]->Fill(RawAux[i+2]); 
    if(RawAux[i+3]>0){
      hELUM_RF[i]->Fill(RawAux[i+3],RawTDC[1]);
    }
  }
  
  //Recoils
 for(int i=0;i<4;++i) {
    hRDT[i]->Fill(RawAux[i*2+8],RawAux[i*2+9]);
  }

 // Fill raw e0 de0  
 hEDE0->Fill(RawAux[0],RawAux[1]);
 hDE0_RF->Fill(RawAux[0],RawTDC[0]);
 
  // Read in ARRAY
  for(Int_t nadc=0;nadc<5;++nadc){ //loop over ADCs 1-5    
    dataword=*p1++;
    nhits[nadc]=cntbit(dataword); //determines number of hits in ADC by counting set bits (ones) in hit register
    for(Int_t i=0;i<nhits[nadc];++i){ //loop over ADC hits, if any
      dataword=*p1++;
      Chan[nadc][i]=((dataword & 0x0000f000)>>12);
      RawData[nadc][i]=(dataword & 0x00000fff);
      hADC[nadc]->Fill(RawData[nadc][i],Chan[nadc][i]);
      if(nadc>-1&&nadc<5&&i>-1&&i<16){ // Tests for sensible ADC and Channel numbers
	if((MapDet[nadc][Chan[nadc][i]])>-1&&(MapDet[nadc][Chan[nadc][i]])<24&&
	   (MapSig[nadc][Chan[nadc][i]])>-1&&(MapSig[nadc][Chan[nadc][i]])<3){ //Tests for sensible values in 
	                                                      
                 //remap matrices.  Note -1 is excluded.
	  Data[MapDet[nadc][Chan[nadc][i]]][MapSig[nadc][Chan[nadc][i]]]=RawData[nadc][i]; //Remaps RawData array
	  
	}
      }
    }
  }


  /****Done unpacking event, filling raw histograms, and remapping data.*****/      
     
  //Filling histograms with (24x3) detector mapping  
  Float_t e=0,xf=0,xn=0;
                       
  for(Int_t i=0;i<24;++i){
    e=Data[i][0];
    xf=Data[i][1];
    xn=Data[i][2];
    if (xn>0 && xf>0) hXFXN[i]->Fill(xn,xf);
    //  if (xn>0 && xf>0 && e>0) hE[i]->Fill(e);

  }       
  return 0;
}
/* The userfunc() function
:  This function is called per event.  The event
 * is supplied by daphne.  Unpack the event and fill your histograms here.
 */
int userfunc(const struct ScarletEvntHdr* h)
{
    ScarletEvnt event, subevent;
    event = h;
 
  switch (event.eventtype()) {

  case SE_TYPE_TRIGGERED:
    userdecode(event);
    break;

  case SE_TYPE_SYNC:
    subevent=event[1];
    scalers(subevent);
    break;

  case SE_TYPE_STOP: 
    stopped=1;
    printf("received stop signal.\n");
    break;
  }
    return 0;
}

/* The userexit() function:  This function is called when the sort thread is
 * stopped.  The sort thread is stopped either by explicitly stopping it or
 * when a new sort is started.  It is not stopped if the sorting completes
 * or terminates.  By not stopping, the user retains access to the ROOT
 * objects even after a sort finishes.  Typically, the userexit() function
 * is used to close whatever root files were opened in userentry().
 */
int userexit()
{
  //    f->Write();
  //  f->Close();
  //    delete f;
    return 0;
}
