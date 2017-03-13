/* Program: example.cxx
 *       Created  by Ken Teh,            Aug. 2005
 *       Modified by Xiaodong Tang,      Jul. 2006
 *       Modified by Masahiro Notani,    Aug. 2006 (146Sm M.Paul's exp)
 *       Modified by Hyeyoung Lee,       Jan. 2008 
 *       Modified by Batman & The Joker  Jul. 2008
 *       Modified by Dr. Oesterman       Feb. 2009 (Helios 12B)
 *       Modified by Scott Marley        Aug. 2009 (Helios 14C runs)
 * Purpose: 
 *       SCARLET Data Acquisition & Histograming for experiments
 *       with Helical Orbit Spectrometer (HELIOS) at ATLAS/ANL
 *
 * ROOT-daphne example sort program.  This program to be used with the
 * fakebldr program in this directory.
 *
 * The user must define userfunc() which is called for each event.  In
 * addition, the user may define userentry() and userexit() which are called
 * once at the start of a sort and when it terminates if defined.  The
 * functions return an int and should return 0.  A non-zero return value
 * terminates the sort.
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

TFile *f,*cutfile; //used to create ROOT file
Float_t totals[NSCALERS];
Int_t stopped;

Bool_t bOldCal=1;
Bool_t bPrintCal=0;
Int_t DoCal[4]={0,0,0,0};
/*Set Calibration level: 
[calibrate E]                    [calibrate X]                    
0 - No calibration.              0 - No calibration.              
1 - Correct position dependence  1 - "gain match" XF & XN         
2 - Energy calibration (MeV)     
*/
Bool_t bOffline=1; //Turn on/off "offline" histograms//
Bool_t bDiag=1; //Turn on/off "diagnostic" histograms


//Constants
Float_t pi=4.0*atan(1.0); 
Float_t MeV=1.602176487E-13;    // J/Mev
Float_t amu=1.660538782E-27;    // kg/amu
Float_t light=299792485;        // m/s
Float_t MeV_amu=amu/MeV*light*light;// (MeV/c^2)/amu

//Array Specifications
Int_t offset=-350; //Position (in mm) of the leading edge of the active area of the detector array relative to the target
Float_t active=50.5; //Length of active area in mm
Float_t positions[7]={offset-active/2+positions[1], //Position of Ta slits
		      66.76,124.12,182.48,241.11,299.87,358.68};//Detector-Center Positions (from schematic)

//Reaction Properties
Float_t mass=(2*MeV_amu+13.1357)/MeV_amu*amu; //Mass of detected particle in kg
Float_t Vcm=2.380E7;      //Center-of-mass velocity in m/s
Float_t Tcyc=45.898;    //cyclotron period in ns
Float_t slopeEcm=((mass*Vcm)/(Tcyc*1E-9))/MeV/1000; //Slope of kinematic curves in 
                                                    //hEZ plot in Mev/mm
Float_t intercept=11.672;  // ground state  b=(1/2.0)*mass*(V0^2-Vcm^2)

//Histograms Set-up
Int_t maxX=3000; //Default value of histogram maximum for uncalibrated Xf,XN plots
Int_t maxE=maxX; //Default value of histogram maximum for uncalibrated energy plots 

Float_t scaleX=0.1; //+/- expansion factor for X-position plots

Float_t minT=0;    //Default value of histogram minimum for uncalibrated time plots
Float_t maxT=1500; //Default value of histogram maximum for time plots 

Float_t minZ=(-positions[6]-active/2+positions[0]-10);
Float_t maxZ=(-positions[1]-active/2+positions[0]+active+10);

Float_t minEc=floor(0-maxZ*slopeEcm);
Float_t maxEc=ceil(maxE-minZ*slopeEcm);

enum {re1,re2,re3,re4,zdg,tgt,wbd};

Float_t XCal[24][10];
Float_t ECal[24][10];
Float_t TCal[24][10];
Float_t _Cal[24][10];


enum {CSI1,CSI2,CSI3,CSI4};
// Declaration of Histograms

/* 1-D histograms: */
TH1F *hRF_AR,*hRF_REC,*hRF_WBD,*hAR_REC;
TH1F *hTAC;
TH1 *hELUM[6];
//TH2F *hECSI;


/* 2-D histograms */
TH2F *hADC[6];
TH2F *hECSIall;
TH2F *hXFXN[25];
TH2F *hEDiff[24];
TH2F *hEX[24];
TH2F *hEcX[24];
TH2F *hESum[24];
TH2F *hDiffX[24];
TH2F *hETAC[4];
TH2F *hETACg[4];
TH2F *hRDT[4];
TH2F *hEDE0;
TH2F *hDE0_RF;
TH2F *hELUM_RF[6];

//TH2F *hEDE[8];
TH2F *hEDE1g, *hEDE2g, *hEDE3g, *hEDE4g;
TH2F *hEDEg2[4];
TH2F *hE, *hXN, *hXF,*hEZ,*hEZg,*hEZgg, *hEZg1,*hEZg2,*hEZg3, *hEZg4;
TH2F *hRF_AR_REC,*hRF_AR_RECg,*hARREC_RFREC,*hARREC_RFAR, *hARREC_RFARg;
TH2F *hTARREC1_RF,*hTARREC2_RF,*hTARREC3_RF,*hTARREC4_RF,*hTARRF_REC1,*hTARRF_REC2,*hTARRF_REC3,*hTARRF_REC4,*hTRECRF_AR;
TH2F *hTARREC1_RFg,*hTARREC2_RFg,*hTARREC3_RFg,*hTARREC4_RFg;
TH2F *hTARC,*hTDC;
TH2F *hETAC_ALL,*hECSISI,*hETCSI;
TH2F *hETACg_ALL,*hECSISIg,*hETCSIg;
TH2F *hEarrESi;

//Graphical Cuts//
TCutG *cEDE2_b12,*cEDE3_b12;
TCutG *cB_EDE1,*cB_EDE2,*cB_EDE3,*cB_EDE4;
TCutG *cEZ, *cEZ_rough,*TAC_p, *TAC_p2, *TAC_ran, *TAC_ran2;
TCutG *cEDE1,*cEDE2,*cEDE3,*cEDE4;

Int_t readcuts(Char_t *cfn="3alpha_cuts.root"){
  cutfile=new TFile(cfn);
  cutfile->ls();
  /*  cB_EDE1=(TCutG *) gROOT->FindObject("cB_EDE1");
  cB_EDE2=(TCutG *) gROOT->FindObject("cB_EDE2");
  cB_EDE3=(TCutG *) gROOT->FindObject("cB_EDE3");
  cB_EDE4=(TCutG *) gROOT->FindObject("cB_EDE4");*/
  cEZ_rough=(TCutG *) gROOT->FindObject("cEZ_rough");
  cEDE1=(TCutG *) gROOT->FindObject("cEDE1");
  cEDE2=(TCutG *) gROOT->FindObject("cEDE2");
  cEDE3=(TCutG *) gROOT->FindObject("cEDE3");
  cEDE4=(TCutG *) gROOT->FindObject("cEDE4");
  /*
  TAC_p=(TCutG *) gROOT->FindObject("TAC_p");
  TAC_p2=(TCutG *) gROOT->FindObject("TAC_p2");
  TAC_ran=(TCutG *) gROOT->FindObject("TAC_ran");
  TAC_ran2=(TCutG *) gROOT->FindObject("TAC_ran2");*/

  //  cEDE2_b12=(TCutG*)gROOT->FindObject("cEDE2_b12");
  //  cEDE3_b12=(TCutG*)gROOT->FindObject("cEDE3_b12");
  cutfile->Close();
  return 0;
}

/* function to count the number of set bits in a 16 bit word */
Int_t cntbit(Int_t word)
{
  Int_t nbits=0;
  for (Int_t ibit=0; ibit<16; ibit++) {
    if (word & (Int_t) TMath::Power(2,ibit)) {nbits++;}
  }
  return nbits;
}

Bool_t cexists(Char_t *cutname)
{

   Bool_t returnvalue=kFALSE;
   if (gROOT->FindObjectClassName(cutname)) returnvalue=kTRUE;
   return returnvalue;
}


/* function to see if x and y are within the boundries of a TCugG */
Bool_t checkcutg(Char_t *cutname,Float_t x, Float_t y)
{
   Bool_t returnvalue=kFALSE;
   if (cexists(cutname)) {
      TCutG *gcut=(TCutG *) gROOT->GetListOfSpecials()->FindObject(cutname);
      if (gcut->IsInside(x,y)==1) returnvalue=kTRUE;
      return returnvalue;

   } else {
     return returnvalue;
   }
}
int readcal(Char_t *calfile1="position.cal",Char_t *calfile2="energy.cal",Char_t *calfile3="ecal.cal")
{
  Bool_t showtest=0;
  Bool_t showcontent=bPrintCal;
  
  Float_t param[24][50]; 
  Int_t errorline=-1;
  cout<<"Reading in calibration files..."<<endl;

  Int_t size=sizeof(param[0])/sizeof(param[0][0]);
  if(showtest)printf("param array size is: [24][%d].\n",size); 
  Bool_t fit=kFALSE;
  for(Int_t i=0;i<24;i++)
    for(Int_t j=0;j<size;j++)
      param[i][j]=0;//initializes all array elements to zero

  FILE * infile;
  Int_t k=1;
   
  while(!fit&&k<=(size)){
    infile = fopen (calfile1,"r");
    if(calfile1!=NULL){
    for(Int_t i=0;i<24;i++){
      for(Int_t j=0;j<k;j++){
	fscanf(infile,"%f",&param[i][j]);
	}
      }
    fclose(infile);
    }
    else
      printf("File %s is NULL!\n",calfile1);
      if(showtest) printf("Testing array length %d:\n",k);

      if (param[0][0]==1)fit=kTRUE;
      else fit=kFALSE;
      for(Int_t i=0;i<24;i++){
	fit=(fit&&(param[i][0]==(i+1)));	
	if(fit){
	  if(showtest) printf("%2.0f ",param[i][0]);
	  for(Int_t j=1;j<k;j++){
	    if(showtest) printf("%7.2g ",param[i][j]);
	    }
	 if(showtest) printf("\n");
	  if((i+1)>errorline)errorline=i+1;
	}
      }
           
     if(fit)printf("File \"%s\" has %d elements per line.\n",calfile1,k);
     
    k++;
}
  
  if(fit){   
    // FILE * outfile;
    //outfile=fopen("output.txt","w");
    if(showcontent) printf("The contents of \"%s\" are:\n",calfile1);
  for(Int_t i=0;i<24;i++){
    //fprintf(outfile,"%2.0f ",param[i][0]);
    if(showcontent) printf("%2.0f ",param[i][0]);
    for(Int_t j=1;j<k-1;j++){
      //fprintf(outfile,"%g ",param[i][j]);
      if(showcontent)      printf("%7.2f ",param[i][j]);
      XCal[i][j-1]=param[i][j];
    }
    //fprintf(outfile,"\n");
    if(showcontent)    printf("\n");
  }
  //fclose(outfile);
  }
  else
    printf("File \"%s\" has more than %d elements per line, or there is an error on line %d.\n",calfile1,size,errorline);


  fit=kFALSE;
  errorline=-1;
  for(Int_t i=0;i<24;i++)
    for(Int_t j=0;j<size;j++)
      param[i][j]=0;//initializes all array elements to zero

  //FILE * infile;
  k=1;
   
  while(!fit&&k<=(size)){
    infile = fopen (calfile2,"r");
    
    if(calfile2!=NULL){
    for(Int_t i=0;i<24;i++){
      for(Int_t j=0;j<k;j++){
	fscanf(infile,"%f",&param[i][j]);
	}
      }
    fclose(infile);
}
    else
      printf("File %s is NULL!\n",calfile2);

      if(showtest) printf("Testing array length %d:\n",k);

      if (param[0][0]==1)fit=kTRUE;
      else fit=kFALSE;
      for(Int_t i=0;i<24;i++){
	fit=(fit&&(param[i][0]==(i+1)));	
	if(fit){
	  if(showtest) printf("%2.0f ",param[i][0]);
	  for(Int_t j=1;j<k;j++){
	    if(showtest) printf("%7.2f ",param[i][j]);
	    }
	 if(showtest) printf("\n");
	  if((i+1)>errorline)errorline=i+1;
	}
      }
           
     if(fit)printf("File \"%s\" has %d elements per line.\n",calfile2,k);
     
    k++;
}
  
  if(fit){   
    if(showcontent)  printf("The contents of \"%s\" are:\n",calfile2);
    for(Int_t i=0;i<24;i++){
      if(showcontent)    printf("%2.0f ",param[i][0]);
      for(Int_t j=1;j<k-1;j++){
	if(showcontent)      printf("%7.2f ",param[i][j]);
	ECal[i][j-1]=param[i][j];
      }
      if(showcontent)    printf("\n");
    }
  }
  else
    printf("File \"%s\" has more than %d elements per line, or there is an error on line %d.\n",calfile2,size,errorline);
  
  fit=kFALSE;
  errorline=-1;
  for(Int_t i=0;i<24;i++)
    for(Int_t j=0;j<size;j++)
      param[i][j]=0;//initializes all array elements to zero

  k=1;
   
  while(!fit&&k<=(size)){
    infile = fopen (calfile3,"r");
    
    if(calfile3!=NULL){
    for(Int_t i=0;i<24;i++){
      for(Int_t j=0;j<k;j++){
	fscanf(infile,"%f",&param[i][j]);
	}
      }
    fclose(infile);
}
    else
      printf("File %s is NULL!\n",calfile3);

      if(showtest) printf("Testing array length %d:\n",k);

      if (param[0][0]==1)fit=kTRUE;
      else fit=kFALSE;
      for(Int_t i=0;i<24;i++){
	fit=(fit&&(param[i][0]==(i+1)));	
	if(fit){
	  if(showtest) printf("%2.0f ",param[i][0]);
	  for(Int_t j=1;j<k;j++){
	    if(showtest) printf("%7.2f ",param[i][j]);
	    }
	 if(showtest) printf("\n");
	  if((i+1)>errorline)errorline=i+1;
	}
      }
           
     if(fit)printf("File \"%s\" has %d elements per line.\n",calfile3,k);
     
    k++;
  }
  
  if(fit){   
    if(showcontent)  printf("The contents of \"%s\" are:\n",calfile3);
    for(Int_t i=0;i<24;i++){
      if(showcontent)    printf("%2.0f ",param[i][0]);
      for(Int_t j=1;j<k-1;j++){
	if(showcontent)      printf("%7.2f ",param[i][j]);
	_Cal[i][j-1]=param[i][j];
      }
      if(showcontent)    printf("\n");
    }
  }
  else
    printf("File \"%s\" has more than %d elements per line, or there is an error on line %d.\n",calfile3,size,errorline);

  Int_t i=0;
  for(int i=0;i<24;i++){
    if(i==0)printf("Energy:   ");
    switch(DoCal[0]) //Set energy calibration level
      {
      case 0://No energy calibration 
	if(i==0) printf("[%1d] No Energy calibration\n",DoCal[0]);
	_Cal[i][0]=1;
	_Cal[i][1]=0;
	break;
      case 1: //Calibrate energy.
	if(i==0)  printf("[%1d] Correct position dependence\n",DoCal[0]);
	if(bOldCal){
	  
	}
	if(_Cal[i][0]==0)_Cal[i][0]= 1; //slope of peakfit
	break;
      case 2:
	if(i==0)  printf("[%1d] Calibrate energy in MeV\n",DoCal[0]);
	break;
	
      default:break;  
      }
    
    if(i==0)printf("Position: ");
    switch(DoCal[1]) //Set position calibration constants
      {
      case 0: //No position calibration.
	if(i==0) printf("[%1d] No position calibration\n",DoCal[1]);
	_Cal[i][2]=-1; //slope of hXFXN for a fixed energy
	_Cal[i][3]= 1; //slope of hESum
	_Cal[i][4]= 0; //intercept of hESum
	break;
      case 1: //Match XF and XN.
	if(i==0)  printf("[%1d] ""Gain match"" XF and XN\n",DoCal[1]);
	if(bOldCal){
	  _Cal[i][2]=XCal[i][1];  
	}

	if(_Cal[i][2]==2)_Cal[i][2]=-1;   //slope of hXFXN for a fixed energy
	_Cal[i][3]=1; //slope of hESum
	_Cal[i][4]=0; //intercept of hESum
	break;
      case 2: //Match (XF+XN) to E.
	if(i==0)  printf("[%1d] ""Gain match"" (XF+XN) to E\n",DoCal[1]);
	if(_Cal[i][2]==2)_Cal[i][2]=-1;
	if(_Cal[i][3]==3)_Cal[i][3]= 1; //slope of hESum
	if(_Cal[i][4]==4)_Cal[i][4]= 0; //intercept of hESum
	break;
	
      default:break;
      }
  }

  if(DoCal[0]||DoCal[1]) 
    cout<<"Applying calibration constants..."<<endl;
  if(showcontent){
    if(DoCal[0]>1){
      printf("Energy Constants:\n");
      printf("             E slope | E offset\n");
      for(Int_t i=0;i<24;i++){ //print out calibration constants
	printf("Detector %2d: %7.3f | %8.3f \n",i+1,_Cal[i][0],_Cal[i][1]);
      }
    }
    if(DoCal[0]&&bOldCal){
      printf("Position Constants:\n");
      printf("            hXFXN slope | hESum Slope | hESum offset\n");
      for(Int_t i=0;i<24;i++){ //print out calibration constants
	printf("Detector %2d:     %6.3f |      %6.3f |      %7.3f\n",
	       i+1,_Cal[i][2],_Cal[i][3],_Cal[i][4]);
      }  
    }
  }
  return 0;
}

 /* The userentry() function:  Create your ROOT objects here.  ROOT objects
  * should always be created on the heap.  That is, always allocate the
  * objects via the new operator.  If you intend to save your histograms to a
  * root file, create the file in userentry().  You can also use this
  * function to load gates and conditions from other root files.
  */
int userentry()
{
  /* stopped flag for Elliot's scaler program */
  for (Int_t i=0; i<NSCALERS; i++) totals[i]=0;
  stopped = 1; 

  //File commands
  // readcuts("3alpha_cuts.root"); 
  f = new TFile("offline.root", "recreate");
  
  if(DoCal[0]||DoCal[1]||DoCal[2]||DoCal[3]){
    if(bOldCal)
      //readcal("oldposition.cal","oldenergy.cal","oldecal.cal");
      readcal("new_position.cal","flat_cal.cal","flat_energy.cal");
    else
      readcal();
  }
  
  //Set histogram ranges, depending on calibration
  TString titleX="1/2[1+(X1-X2)/(X1+X2)]";
  
  if(DoCal[0]>1){
    if(maxZ>0)
      //Set histogram maximum for calibrated energy plots (MeV)
      maxE=(Int_t)ceil(slopeEcm*maxZ+intercept); 
    else 
      maxE=(Int_t)ceil(intercept);
    maxE=(Int_t)6;
    //    maxEc=ceil(maxE-minZ*slopeEcm);
    // minQ=floor(intercept-maxEc);      
  } 
  
  if(DoCal[2] >3){
    minT=-Tcyc;
    maxT=Tcyc+(Tcyc-minT); //Time calibration centers time plots on Tcyc.
  }

  if(((maxZ-minZ)<382)&&(DoCal[1]>2)){//Shifts Z-plots by offset if range defined by array size
    //    maxZ-=_Cal[0][14];
    //minZ-=_Cal[0][14];
  }

  Int_t bin1=256;//Sets number of bins on most histograms to conveniently reduce memory load
  
  // 1d histograms
  
  // 2d histograms  
  hTDC=new TH2F("hTDC","hTDC",512,0,4096,17,0,17);
  //E0 DE0
  hEDE0=new TH2F("hEDE0","hEDE0",512,0,4096,512,0,4096);
  hDE0_RF=new TH2F("hDE0_RF","hDE0_RF",512,0,4096,512,0,4096);
  
  hECSIall=new TH2F("hECSIall","CsI Det. vs CsI energy",bin1,0,4096,4,0,4);
  hEarrESi=new TH2F("hEarrESi","E(silicon) vs E(Array)",bin1,0,maxE,bin1,0,4096);
  
  hETAC_ALL=new TH2F("hETAC_ALL","TAC vs E(Si)",bin1,0,maxE,bin1,0,4096);
  hECSISI=new TH2F("hECSISI","Esum(CSI) vs E(Si)",bin1,0,maxE,bin1,0,16384);
  hETCSI=new TH2F("hETCSI","TAC vs Esum(CsI)",bin1,0,4096,bin1,0,16384);

  hETACg_ALL=new TH2F("hETACg_ALL","TAC vs E(Si) (goodEZ gated)",bin1,0,maxE,bin1,0,4096);
  hETCSIg=new TH2F("hETCSIg","TAC vs Esum(CsI) (goodEZ gated) ",bin1,0,4096,bin1,0,16384);
  hECSISIg=new TH2F("hECSISIg","Esum(CSI) vs E(Si) (goodEZ gated)",bin1,0,maxE,bin1,0,16384);

  hETAC[0]=new TH2F("hETAC1","TAC CsI1 vs E(Si)",bin1,0,maxE,bin1,0,4096);
  hETAC[1]=new TH2F("hETAC2","TAC CsI2 vs E(Si)",bin1,0,maxE,bin1,0,4096);
  hETAC[2]=new TH2F("hETAC3","TAC CsI3 vs E(Si)",bin1,0,maxE,bin1,0,4096);
  hETAC[3]=new TH2F("hETAC4","TAC CsI4 vs E(Si)",bin1,0,maxE,bin1,0,4096);

  hETACg[0]=new TH2F("hETAC1g","TAC CsI1 vs E(Si) (goodEZ)",bin1,0,maxE,bin1,0,4096);
  hETACg[1]=new TH2F("hETAC2g","TAC CsI2 vs E(Si) (goodEZ)",bin1,0,maxE,bin1,0,4096);
  hETACg[2]=new TH2F("hETAC3g","TAC CsI3 vs E(Si) (goodEZ)",bin1,0,maxE,bin1,0,4096);
  hETACg[3]=new TH2F("hETAC4g","TAC CsI4 vs E(Si) (goodEZ)",bin1,0,maxE,bin1,0,4096);
  
  //  hECSI=new TH2F("hECSI1","CsI1 vs E(Si)",bin1,0,maxE,bin1,0,16384);

  
  //hE=new TH2F("hE","Detector Energies (1-24)",1024,0,maxE,25,0,25);
  //hXF=new TH2F("hXF","Detector Position (far)",1024,0,maxX,25,0,25);
  //hXN=new TH2F("hXN","Detector Position (near)",1024,0,maxX,25,0,25);

  
  // hTDC=new TH2F("hTDC","hTDC",512,0,4096,16,0,16);

  hEZ=new TH2F("hEZ","Energy vs. Position",512,minZ,maxZ,512,0,maxE);
  hEZg=new TH2F("hEZg","Energy vs. Position (gated:TAC)",512,minZ,maxZ,512,0,maxE);
  hEZgg=new TH2F("hEZgg","Energy vs. Position (gated:TAC & CSI-OR)",512,minZ,maxZ,512,0,maxE);
  hEZg1=new TH2F("hEZg1","Energy vs. Position (gated: CSI1)",512,minZ,maxZ,512,0,maxE);
  hEZg2=new TH2F("hEZg2","Energy vs. Position (gated: CSI2)",512,minZ,maxZ,512,0,maxE);
  hEZg3=new TH2F("hEZg3","Energy vs. Position (gated: CSI3)",512,minZ,maxZ,512,0,maxE);
  hEZg4=new TH2F("hEZg4","Energy vs. Position (gated: CSI4)",512,minZ,maxZ,512,0,maxE);
  
  for(int a=0;a<7;++a){
    TString name="hADC";
    TString title="raw ADC";
    name+=(a+1);
    title+=(a+1);
    hADC[a]=new TH2F(name,title,1024,0,4095,17,0,17);
  }

  for(int a=0;a<25;++a){
    TString name="hXFXN";
    TString title="XF vs. XN detector ";
    if(a<24){
      name+=(a+1);
      title+=(a+1);
    }
    else
      title+="all";
    hXFXN[a]=new TH2F(name,title,512,0,maxX,512,0,maxX);
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
 

    if(bOffline){//build "offline" histograms
    for(int a=0;a<24;++a){
      TString name="hEX";
      TString title="E vs. 1/2*{1+[(XF-XN)/(XF+XN)]} det. ";
      name+=(a+1);
      title+=(a+1);
      hEX[a]=new TH2F(name,title,512,0-scaleX,1+scaleX,512,0,maxE);
    }
    for(int a=0;a<24;++a){
      TString name="hEcX";
      TString title="E-(XF+XN) vs. X det. ";
      name+=(a+1);
      title+=(a+1);
      hEcX[a]=new TH2F(name,title,bin1,0-scaleX,1+scaleX,bin1,minEc,maxEc);
    }

  }//end offline

  if(bDiag){//build "diagnostic" histograms
    hE=new TH2F("hE","Detector Energies (1-24)",1024,0,maxE,25,0,25);
    hXF=new TH2F("hXF","Detector Position (far)",1024,0,maxX,25,0,25);
    hXN=new TH2F("hXN","Detector Position (near)",1024,0,maxX,25,0,25);

 for(int a=0;a<24;++a){
      TString name="hEDiff";
      TString title="E vs.(XF-XN) detector ";
      name+=(a+1);
      title+=(a+1);
      hEDiff[a]=new TH2F(name,title,bin1,-maxX,maxX,bin1,0,maxE);
    }
  
    for(int a=0;a<24;++a){
      TString name="hESum";
      TString title="E vs. (XN+XF) det. ";
      name+=(a+1);
      title+=(a+1);
      hESum[a]=new TH2F(name,title,bin1,0,maxX,bin1,0,maxE);
  }
  
    for(int a=0;a<24;++a){
      TString name="hDiffX";
      TString title="E-(XF+XN) vs. X det. ";
      name+=(a+1);
      title+=(a+1);
      hDiffX[a]=new TH2F(name,title,bin1,0-scaleX,1+scaleX,bin1,-500,500);
    }
    
    
  }//end diag


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

  //unpack event
  /*-------
   ADC6:(see below)
   
   hitpattern ADC1
   data
   
   hitpattern ADC2
   data
   ...
   hitpattern ADC3
   data
   ...
   hitpattern ADC4
   data
   ...
   hitpattern ADC5
   data
   ...
   0x0000dead
   ---------*/
  // enum {arrf,arrec,recrf,tdc3,tdc4,tdc5,tdc6,tdc7,tdc8,e0rf,etrf,edrf,e1rf,e2rf,e3rf,e4rf};
  
  /*  The data are parsed into Channel ID and Channel Data and then output 
    into a raw array scheme (5x16) and then mapped to an array scheme (24x3).  
  */
 
  Int_t EDE[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  //array for CsI energies;  (0-3) {CsI1, CsI2,...}
  Int_t TAC=0;
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

 Int_t Toffset[24]={1959,2033,2003,1924,2267,1912,2327,2223,
		     2168,2096,2162,   0,2029,1910,2007,1918,
		     1841,1796,1569,1822,2034,1675,1732,1602};

 Float_t ecsisum=0.;
 Float_t eSi=0;


  for(Int_t ii=0; ii<24; ii++)
    {
      Data[ii][0]=0;
      Data[ii][1]=0;
      Data[ii][2]=0;
    }
  
  for(int a=0;a<16;++a)EDE[a]=0;
  for(int i=0;i<16;++i) {
    RawAux[i]=0;
    RawTDC[i]=0;
  }
 
  //Beginning Event Unpacking...   

  //read in ADC 6
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
  
    //read in ADC 1-5 (Array)
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
  

  // Define Flags

  Bool_t goodEZ=kFALSE;
  Bool_t goodT=kFALSE;
  Bool_t goodECSI=kFALSE;
  Bool_t goodECSI1=kFALSE;
  Bool_t goodECSI2=kFALSE;
  Bool_t goodECSI3=kFALSE;
  Bool_t goodECSI4=kFALSE;

  //Filling histograms with (24x3) detector mapping  
  Float_t e=0,xf=0,xn=0,x=0,z=0;

  //  Float_t sum=0;
  Int_t lowthr=0;
  // Int_t tarc=0;
  for(Int_t i=0;i<24;++i){
      e=Data[i][0];
      xf=Data[i][1];
      xn=Data[i][2];
       
      if((e>lowthr)&&(xf>lowthr)&&(xn>lowthr)) //Tests all three signals against "lowthr"
{
	  // tarc=TAC[8]+2000-Toffset[i];
	  //hTARREC_RF->Fill(tarc,TAC[recrf]);
	  // if (goodT) hTARC->Fill(tarc,i);
	  //Position Calibration Level [1] - Matches XF to XN

	  if(DoCal[1]){
	    if(bOldCal){
	      if(_Cal[i][2]<-1){ 
		xn=(-_Cal[i][2])*xn;}
	      else{
		xf=(-1/_Cal[i][2])*xf;}
	    }
	    else{//this doesn't quite work... yet
	      Float_t weight=0;
	      for(Int_t k=1;k<9;k++){
		weight+=XCal[i][k]*pow(xn,k);
	      }
	      xf+=(-1*xn-weight)/2;
	      xn+=(-1*xn-weight)/2;
	    }
	  }


	  //Position Calibration Level [2] - Matches (XF+XN) to E
	  // xf=(xf*_Cal[i][3]+_Cal[i][4]/2);
	  //xn=(xn*_Cal[i][3]+_Cal[i][4]/2); 
	  if(e>lowthr)//setting a "high" low-threshold here allows a single energy line to be fit in the hXFXN plots!  See hXFXN
	    hXFXN[i]->Fill(xn,xf);
	  hXFXN[24]->Fill(xn,xf);
	  
	  
	  x=(1/2.)*(1+((xf-xn)/(xf+xn))); //Position on detector with XN@x=0 and XF@x=1.  
	                                  //Please note at this point that the array PCBs are wired backwards, 
	                                  //so "X-Far" is closest to the target and "X-Near" is further away.

	  //Energy Calibration:	  
	  if(DoCal[0]){
	    Float_t weight=0;
	    for(Int_t k=1;k<9;k++){
	      weight+=ECal[i][k]*pow(x,k);
	    }
	    e-=weight;
	    if(DoCal[0]>1)
	    e =(( e- _Cal[i][1])   /_Cal[i][0]); //Energy in MeV
	  }
	  
	  //	  sum=e-(xf+xn);
	  // if(bDiag)
	  //  hDiffX[i]->Fill(x,sum);
	  
 	  z=-positions[(6-(i%6))]-active/2+positions[0]+(active*x); //position in magnet in mm	  
	  // printf("%f %f\n",z,e);
	  if (checkcutg("cEZ",z,e)) goodEZ=kTRUE;

	hEZ->Fill(z,e);
	if(bOffline)
	  hEX[i]->Fill(x,e);
	
	
	
	//Define conditions (gate) 
	
	if (checkcutg("cEZ_rough",z,e)) {goodEZ=kTRUE;}
	if((TAC>=140)&&(TAC<=2500)) goodT=kTRUE;
	if((EDE[0]>=100)&&(EDE[0]<=4000)) goodECSI1=kTRUE;
	if((EDE[1]>=100)&&(EDE[1]<=4000)) goodECSI2=kTRUE;
	if((EDE[2]>=100)&&(EDE[2]<=4000)) goodECSI3=kTRUE;
	if((EDE[3]>=100)&&(EDE[3]<=4000)) goodECSI4=kTRUE;
	if (goodECSI1||goodECSI2||goodECSI3||goodECSI4) goodECSI=kTRUE;
	if(TAC>350&&TAC<550) goodT=kTRUE;
	

	hEZ->Fill(z,e);
	if (TAC>50) { 
	  
	  hETAC_ALL->Fill(e,TAC);
	  
	  if(goodECSI1) hETAC[0]->Fill(e,TAC);
	  if(goodECSI2) hETAC[1]->Fill(e,TAC);
	  if(goodECSI3) hETAC[2]->Fill(e,TAC);
	  if(goodECSI4) hETAC[3]->Fill(e,TAC);

	  if(goodEZ) {
	    hETACg_ALL->Fill(e,TAC);
	    if(goodECSI1) hETACg[0]->Fill(e,TAC);
	    if(goodECSI2) hETACg[1]->Fill(e,TAC);
	    if(goodECSI3) hETACg[2]->Fill(e,TAC);
	    if(goodECSI4) hETACg[3]->Fill(e,TAC);
	  }

	}
	if (ecsisum>200) hECSISI->Fill(e,ecsisum);
	if (ecsisum>200) hETCSI->Fill(ecsisum,TAC);
	
	//Histograms with conditions (gated)
	
	if (TAC>1510 && TAC<1560) {hEarrESi->Fill(e,eSi);}
	
	
	if(ecsisum>200&&goodT) {
	  hETCSIg->Fill(ecsisum,TAC);
	  hECSISIg->Fill(e,ecsisum);
	}
	
	if(goodEZ&&goodT&&goodECSI) hEZgg->Fill(z,e);
	if(goodEZ&&goodT&&goodECSI1) hEZg1->Fill(z,e);
	if(goodEZ&&goodT&&goodECSI2) hEZg2->Fill(z,e);
	if(goodEZ&&goodT&&goodECSI3) hEZg3->Fill(z,e);
	if(goodEZ&&goodT&&goodECSI4) hEZg4->Fill(z,e);
	
	
	}
  }
  
  
  
   
   
    return 0;
}

/* The userfunc() function:  This function is called per event.  The event
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
  cout<<"Exiting sort..."<<endl;    
  f->Write();
  f->Close();
  delete f;
  printf("\a"); //"Default Beep" at sort exit.
  return 0;
}
