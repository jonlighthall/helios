/*----------------------------PostScript "pretty-print" page width-----------------------------*/
/* Program: example.cxx
 *       Created  by Ken Teh,            Aug. 2005
 *       Modified by Xiaodong Tang,      Jul. 2006
 *       Modified by Masahiro Notani,    Aug. 2006 (146Sm M.Paul's exp)
 *       Modified by Hyeyoung Lee,       Jan. 2008 
 *       Modified by Batman & The Joker, Jul. 2008 (28Si Lighthall)
 * Purpose: 
 *       SCARLET Data Acquisition & Histograming for experiments with 
 *       Helical Orbit Spectrometer (HELIOS) at ATLAS/ANL (no recoil detector)
 *
 * File Compatibility:
 *       ROC File: helios_time_roc1.c
 *       ROOT File: [separation in mm].root
 *       Calibration File: [separation in mm].cal
 *       Weight Function File: calibration.cal
 *
 * Version Information and Functionality:
 *       sort            - overall restructuring of the helios sort code
 *       sort_Si28       - histogram instanciation, etc, changed to match B12 sort files
 *        
 * ROOT-daphne example sort program.  This program to be used with the
 * fakebldr program in this directory.
 *
 * The user must define userfunc() which is called for each event.  In addition, the user may
 * define userentry() and userexit() which are called once at the start of a sort and when it 
 * terminates if defined.  The functions return an int and should return 0.  A non-zero return 
 * value terminates the sort.
 */
 
// Header Files
using namespace std; //used to eliminate deprecated header file error message
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

TFile *f; //used to create ROOT file
TFile *cutfile;
Float_t totals[NSCALERS];
Int_t stopped;

//Structures and Physical Constants
Float_t pi=4.0*atan(1.0); 
Float_t MeV=1.602E-13; //J/Mev
Float_t active=50.5; //Length of active area in mm
char buffer [50];
Int_t iter=0;     //used for debugging (set number of print-to-screen occurrences)
Int_t Counts[24]; //used to monitor the counts in each detector for each run

//Experimental Setup
TString deltaZ="500_algor"; // <--------Enter nominal target-detector separation here (in mm)
Int_t separation=atoi(deltaZ.Data());
Int_t offset=-separation; //Distance in mm between active detector area and target

Float_t mass=1*1.673E-27; //Mass of detected particle in kg
Float_t Vcm=3.174E7;      //Center-of-mass velocity in m/s
Float_t Tcyc=34.246;      //cyclotron period in ns
Float_t slopeEcm=((mass*Vcm)/(Tcyc*1E-9))/MeV/1000; //Slope of kinematic curves in 
                                                    //hEZ plot in Mev/mm
Float_t intercept=11.672;  //[0] ground state  b=(1/2.0)*mass*(V0^2-Vcm^2)
Float_t positions[7]={offset-active/2+positions[1],             //Position of Ta slits
                      66.76,124.12,182.48,241.11,299.87,358.68};//Detector-Center Positions 
                                                                //in mm (taken from schematic)
Float_t slopeT=-18.01; //time dispersion in ns/channel

//Turn channels on and off here
    Int_t include[24]={ 1, 1, 1, 1, 1, 1,
          //Detector No 1, 2, 3, 4, 5, 6, Note: 5 is not biased
	 	        1, 1, 1, 1, 1, 1,
   		     // 7, 8, 9,10,11,12, Note:  
		        1, 1, 1, 1, 1, 1,
   		     //13,14,15,16,17,18,
		        1, 1, 1, 1, 1, 1};
   		     //19,20,21,22,23,24

//Calibration
Int_t DoCal[4]={0,0,0,0};
/*Set Calibration level: 
[calibrate E]                  [calibrate X]                     [calibrate T]                  
0 - No calibration.            0 - No calibration.               0 - No calibration.            
1 - Flatten hEX plots (quad)   1 - "gain match" XF & XN          1 - Piecewise walk correction  
2 - Energy calibration (MeV)   2 - "gain match" (XF+XN) & E.     2 - Linear walk correction
                               3 - Position Calibration (offset) 3 - Flatten hTZ plots
                               4 - Position Calibration (slope)  4 - Calibrate time (ns)
[calibrate Ecm]
0 - No calibration.
1 - Q-Value calibration
*/
Float_t ECal[24][21];  //Stores calibration constants read in with readcal()

//Gating & Cuts Set-up
Float_t DoCut[5]={1, //Energy cut ON/OFF
                  1, //Position cut ON/OFF    
		  1, //Time cut ON/OFF
		  0, //TOF cut ON/OFF
		  1};//e=(xf+xn) cut ON/OFF

Float_t   cutE=1024;
Float_t sigmaE=13; //peak-fitted sigma of energy in MeV.
Float_t widthE=1; // Energy gate width in +/- sigma

Float_t   cutX=0.5; //Center of detector
Float_t sigmaX=0.5; //1/2 detector length
Float_t widthX=1.0; //Position gate in +/- sigma

Float_t   cutT=Tcyc*3.0/2.0; 
Float_t sigmaT=4.278;//3.313; //peak-fitted sigma of time in ns.
Float_t widthT=2.5; // Time gate width in +/- sigma

Float_t   cutTOF=34.02;
Float_t sigmaTOF=.548; //peak-fitted sigma of energy in MeV.; 
Float_t widthTOF=2; // Energy gate width in +/- sigma

Float_t sigmaSum=8.94;//peak-fit sigma of e-(xf+xn) [258 chan per MeV -> ~4keV per chan]
Float_t widthSum=2; // 

Float_t sigmaDiff=sigmaSum;
Float_t widthDiff=16;

//Histograms Set-up
Int_t maxX=4096; //Default value of histogram maximum for uncalibrated Xf,XN plots
Int_t maxE=maxX; //Default value of histogram maximum for uncalibrated energy plots 

Float_t scaleX=0.1; //+/- expansion factor for X-position plots

Float_t minT=0;    //Default value of histogram minimum for uncalibrated time plots
Float_t maxT=1500; //Default value of histogram maximum for time plots 

//Float_t minZ=(-positions[6]-active/2+positions[0]-10);
//Float_t maxZ=(-positions[1]-active/2+positions[0]+active+10);

Float_t minZ=-1000;
Float_t maxZ=0;

Float_t minq=0;
Float_t maxq=60;

Float_t minEc=floor(0-maxZ*slopeEcm);
Float_t maxEc=ceil(maxE-minZ*slopeEcm);

Float_t minQ=floor(intercept-maxEc);      
Float_t maxQ=ceil (intercept-minEc);

// Declaration of Histograms
TString name;
//TH2F *hist;

/* 2-D histograms */
TH2F *hADC[7];

TH2F *hE,*hXN,*hXF,*hT;

TH2F *hXFXN[24];
TH2F *hEX[24];
TH2F *hEZ;

TCutG *cTime2D; 
TCutG *cScat;

/* function to count the number of set bits inTH2F a 16 bit word */
Int_t cntbit(Int_t word)
{
  Int_t nbits=0;
  for (Int_t ibit=0; ibit<16; ibit++) {
    if (word & (Int_t) TMath::Power(2,ibit)) {nbits++;}
  }
  return nbits;
}

//Int_t readcuts(Char_t *cfn="time_cuts.root"){
Int_t readcuts(Char_t *cfn) {
  cutfile=new TFile(cfn);
  cutfile->ls();
  cTime2D=(TCutG *) gROOT->FindObject("time2d");


  cutfile->Close();
  return 0;
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


int readcal(Char_t *calfile)
{
  ifstream infile(calfile);
  Float_t detno=0;
  /*
    Calibration File structured as follows: <---out of date
    [X Detector][0 XF gain][1 XN gain][2 ESum slope][3 ESum offset][4 E slope][5 E offset]
    [6 T offset][7 XFXN slope][8 XFXN offset][9 E walk max][10 p2][11 p1]
  */
  printf("ECal array length is: %d.\n",sizeof(ECal[0])/sizeof(ECal[0][0]));
  printf("Reading in Calibration File \"%s\" with Calibration Levels:\n",calfile);
  for(Int_t i=0;i<24;i++)
    {
      infile>>detno;       //First number in each row is detector number
      if(!(detno==(i+1))){
	printf("Calibration File Corrupt on line %2d!\n",i);
	infile.close();	
	userexit();
      }
      else{ //Only read in row if first element is corrent detector number
	for(Int_t j=0;j<sizeof(ECal[0])/sizeof(ECal[0][0]);j++){
	  infile>>ECal[i][j];
	}
      }
      if(i==0)printf("Energy:   ");
      switch(DoCal[0]) //Set energy calibration level
	{
	case 0://No energy calibration 
	  if(i==0) printf("[%1d] No Energy calibration\n",DoCal[0]);
	  ECal[i][0] =1; //Energy, peakfit slope
	  ECal[i][1] =0; //Energy, peakfit offest
	  ECal[i][19]=0; //hEX X-Profile fit, p1
	  ECal[i][20]=0; //hEX X-Profile fit, p2
	  break;
	case 1: //Calibrate energy.
	  if(i==0)  printf("[%1d] Correct Energy position-dependance\n",DoCal[0]);
	  ECal[i][0]=1;                     //Energy, peakfit slope
	  ECal[i][1]=0;	                    //Energy, peakfit offest
	  if(ECal[i][19]==19)ECal[i][19]=0; //hEX X-Profile fit, p1
	  if(ECal[i][20]==20)ECal[i][20]=0; //hEX X-Profile fit, p2
	  break;
	case 2:
	  if(i==0)  printf("[%1d] Calibrate energy in MeV\n",DoCal[0]);
	  if(ECal[i][0]==0)ECal[i][0]=1;    //Energy, peakfit slope
	  if(ECal[i][1]==1)ECal[i][1]=0;    //Energy, peakfit offest
	  if(ECal[i][19]==19)ECal[i][19]=0; //hEX X-Profile fit, p1
	  if(ECal[i][20]==20)ECal[i][20]=0; //hEX X-Profile fit, p2
	  break;
	default:break;  
	}
      
      if(i==0)printf("Position: ");
      switch(DoCal[1]) //Set position calibration constants
	{
	case 0: //No position calibration
	  if(i==0) printf("[%1d] No position calibration\n",DoCal[1]);
	  ECal[i][2]=-1;                    //slope of hXFXN for fixed energy
	  ECal[i][15]=1;                    //slope of hESum
	  ECal[i][16]=0;                    //intercept of hESum
	  ECal[i][13]=(Float_t)1;                    //expansion factor (about x=0.5)
	  ECal[i][14]=0;                    //offset in mm
	  break;
	case 1: //Read in calibration for XF and XN. Overwrite others.
	  if(i==0)  printf("[%1d] ""Gain match"" XF and XN\n",DoCal[1]);
	  if(ECal[i][2]==2)ECal[i][2]=-1;   //slope of hXFXN for fix energy
	  ECal[i][15]=1;                    //slope of hESum
	  ECal[i][16]=0;                    //intercept of hESum
	  ECal[i][13]=1;                    //expansion factor (about x=0.5)
	  ECal[i][14]=0;                    //offset in mm
	  break;
	case 2: //Match (XF+XN) to E
	  if(i==0)  printf("[%1d] ""Gain match"" (XF+XN) to E\n",DoCal[1]);
	  if(ECal[i][2]  ==2)ECal[i][2]=-1; //slope of hXFXN for fix energy
	  if(ECal[i][15]==15)ECal[i][15]=1; //slope of hESum
	  if(ECal[i][16]==16)ECal[i][16]=0; //intercept of hESum
	  ECal[i][13]=1;                    //expansion factor (about x=0.5)
	  ECal[i][14]=0;                    //offset in mm
	  break;
	case 3: //calibrate array position (Z) in mm
	  if(i==0)  printf("[%1d] Adjust overall Z offset\n",DoCal[1]);
	  if(ECal[i][2]  ==2)ECal[i][2]=-1; //slope of hXFXN for fix energy
	  if(ECal[i][15]==15)ECal[i][15]=1; //slope of hESum
	  if(ECal[i][16]==16)ECal[i][16]=0; //intercept of hESum
	  ECal[i][13]=1;                    //expansion factor (about x=0.5)
	  if(ECal[i][14]==14) ECal[i][14]=0;//offset in mm
	  break;
	case 4: //calibrate X to fit simulation 
	  if(i==0)  printf("[%1d] Adjust X to correct slope\n",DoCal[1]);
	  if(ECal[i][2]  ==2)ECal[i][2]=-1; //slope of hXFXN for fix energy
	  if(ECal[i][15]==15)ECal[i][15]=1; //slope of hESum
	  if(ECal[i][16]==16)ECal[i][16]=0; //intercept of hESum
	  if(ECal[i][13]==13)ECal[i][13]=1; //expansion factor (about x=0.5)
	  if(ECal[i][14]==14)ECal[i][14]=0; //offset in mm

	  /* "Scissoring" Calibration Constants
	    if(ECal[i][7]==1&&ECal[i][8]==0){ //Enter default values on detectors 
                                              //without calibration
	    ECal[i][7]=1.051;                 //Set "Type 1" detectors
	    ECal[i][8]=-2.446;
	    if(i==(14-1)||i==(16-1)){         //Set "Type 2" detectors
	      printf("i = %2d, Type 2\n",i);
	      ECal[i][7]=1.201;
	      ECal[i][8]=-9.977;
	      }
	  }
	  */
	  break;
	default:break;
	}
     
      if(i==0)printf("Time:     ");
      
      switch(DoCal[2]) //Set time calibration level
	{
	case 0: //No time calibration
	  if(i==0) printf("[%1d] No time calibration\n",DoCal[2]);
	  ECal[i][3]=0; //hTZ ProjectionX pol4 fit, p1 
	  ECal[i][4]=0; //p2
	  ECal[i][5]=0; //p3
	  ECal[i][6]=0; //p4
	  ECal[i][7]=1; //slopeT (time dispersion)
	  ECal[i][8]=0; //proton peak
	  ECal[i][9]=0; //E max for piece-wise quadratic fit (walk correction)
	  ECal[i][10]=0; //hET ProjectionY p1
	  ECal[i][11]=0; //p2
	  ECal[i][12]=0; //hET ProjectionY slope (p1)
	  break;
	case 1: //Read in piecewise-quadratic walk E vs. T correction parameters
	  if(i==0) printf("[%1d] Walk Correction (hET#)\n",DoCal[2]);
	  ECal[i][3]=0; //hTZ pol4 fit, p1 
	  ECal[i][4]=0; //p2
	  ECal[i][5]=0; //p3
	  ECal[i][6]=0; //p4
	  ECal[i][7]=1;//slopeT (time dispersion)
	  ECal[i][8]=0;//proton peak
	  if( ECal[i][9]==9)ECal[i][9]=0;//E max for piece-wise quadratic fit (walk correction)
	  if( ECal[i][10]==10)ECal[i][10]=0;//hET ProjectionY p1
	  if( ECal[i][11]==11)ECal[i][11]=0;//p2
	  ECal[i][12]=0; //linear walk correction slope
	  break;

	case 2: //Read in linear walk E vs. T correction parameters
	  if(i==0) printf("[%1d] Walk Correction (hET#)\n",DoCal[2]);
	  ECal[i][3]=0; //hTZ pol4 fit, p1 
	  ECal[i][4]=0; //p2
	  ECal[i][5]=0; //p3
	  ECal[i][6]=0; //p4
	  ECal[i][7]=1;//slopeT (time dispersion in channels per ns)
	  ECal[i][8]=0;//proton peak
	  if( ECal[i][9]==9)ECal[i][9]=0;
	  if( ECal[i][10]==10)ECal[i][10]=0;
	  if( ECal[i][11]==11)ECal[i][11]=0;
	  if( ECal[i][12]==12)ECal[i][12]=0;//linear walk correction slope
	  break;
	case 3: //Read in T vs. Z correction parameters
	  if(i==0) printf("[%1d] Flatten hTZ# Plots\n",DoCal[2]);
	  if( ECal[i][3]==3)ECal[i][3]=0; //hTZ pol4 fit, p1 
	  if( ECal[i][4]==4)ECal[i][4]=0; //p2
	  if( ECal[i][5]==5)ECal[i][5]=0; //p3
	  if( ECal[i][6]==6)ECal[i][6]=0; //p4
	 ECal[i][7]=1;//slopeT (time dispersion in channels per ns)
	  ECal[i][8]=0;//peakfit offset, or proton peak
	  if( ECal[i][9]==9)ECal[i][9]=0;
	  if( ECal[i][10]==10)ECal[i][10]=0;
	  if( ECal[i][11]==11)ECal[i][11]=0;
	  if( ECal[i][12]==12)ECal[i][12]=0;//linear walk correction slope
	  break;
	case 4: // Read in slope and peak
	  if(i==0) printf("[%1d] Calibrate Time (ns)\n",DoCal[2]);
	  if( ECal[i][3]==3)ECal[i][3]=0;
	  if( ECal[i][4]==4)ECal[i][4]=0;
	  if( ECal[i][5]==5)ECal[i][5]=0;
	  if( ECal[i][6]==6)ECal[i][6]=0;
	  if( ECal[i][7]==7) ECal[i][7]=slopeT;	 
	  if( ECal[i][8]==8)ECal[i][8]=0; //proton peak in hET
	  if( ECal[i][9]==9)ECal[i][9]=0; //max E for piece-wise fit
	  if( ECal[i][10]==10)ECal[i][10]=0; //piece-wise p1
	  if( ECal[i][11]==11)ECal[i][11]=0;//piece-wise p1
	  if( ECal[i][12]==12)ECal[i][12]=0;//linear walk correction slope
	  break;
	default:break;  
	}

      if(i==0)printf("Q-Value:  ");
      switch(DoCal[3]) //Set Q-Value calibration level
	{
	case 0://No Q-Value calibration 
	  if(i==0) printf("[%1d] No Q-Value calibration\n",DoCal[3]);
	  ECal[i][17]=1;
	  ECal[i][18]=0;
	  break;
	case 1: //Calibrate Q-Value
	  if(i==0)  printf("[%1d] Calibrate Q-Value in MeV\n",DoCal[3]);
	  if( ECal[i][17]==17)ECal[i][17]=1;
	  if( ECal[i][18]==18)ECal[i][18]=0;
	}
    }
  printf("Calibration file successfully read.\n");
  infile.close();	
  return 0;
}


/* The userentry() function:  Create your ROOT objects here.  ROOT objects should always be 
 * created on the heap.  That is, always allocate the objects via the new operator.  If you 
 * intend to save your histograms to a root file, create the file in userentry().  You can also 
 * use this function to load gates and conditions from other root files.
 */
int userentry()
{
  sprintf(buffer,"%d_cuts.root",separation);
  // readcuts((Char_t*)(buffer));

  //readcuts("time_cuts.root"); //must be called before ROOT file is defined!(?)

  //Open ROOT file
  f = new TFile((deltaZ+".root"), "recreate");
  printf("Output ROOT file is %s\n",(deltaZ+".root").Data()); 
  
  for(Int_t i=0;i<24;i++){//
    Counts[i]=0;
  }

/* stopped flag for Elliot's scaler program */
  for (Int_t i=0; i<NSCALERS; i++) totals[i]=0;
  stopped = 1;
//Read in calibration constants
  printf("Separation = %d mm\n",separation);
  printf("Physical Constants:\n");  
  printf(" Detected Particle Mass: %g kg\n",mass);  
  printf("Center of Mass Velocity: %g m/s\n",Vcm);  
  printf("       Cyclotron Period: %g ns\n",Tcyc);  
  printf("               slopeEcm: %f MeV/mm\n",slopeEcm);  

  if(DoCal[0]==0&&DoCal[1]==0&&DoCal[2]==0){
    cout<<"Histograms being filled with RAW data"<<endl;
  }
  else{
    sprintf(buffer,"%d.cal",separation);
    readcal((Char_t*)(buffer));
    if(DoCal[0]){
    cout<<"Applying calibration constants:"<<endl;
    printf("Energy Constants:\n");
    printf("             E slope | E offset |    EX p1 |    EX p2 \n");
    for(Int_t i=0;i<24;i++){ //print out calibration constants
      printf("Detector %2d: %7.3f | %8.3f | %8.3f | %8.3f \n",i+1,ECal[i][0],ECal[i][1],ECal[i][19],ECal[i][20]);
    }}
    if(DoCal[1]){
    printf("Position Constants:\n");
    printf("Overall Offset is %5.2f mm\n",ECal[0][14]);//Previous values: 8.0@100, 22.0@350, 22.1@500 [2/4/09]
    printf("            hXFXN slope | hESum Slope | hESum offset | hEcZm slope\n");
    for(Int_t i=0;i<24;i++){ //print out calibration constants
      printf("Detector %2d:     %6.3f |      %6.3f |     %8.3f | %7.5f\n",
	     i+1,ECal[i][2],ECal[i][15],ECal[i][16], ECal[i][13]);
    }}
    
    if(DoCal[2]){
    printf("Time Constants:\n");
    printf("             p1 TZ |  p2 TZ |  p3TZ |  p4 TZ |  Tslp | peak |Emx | p1ET |p2ET | p1 ET\n");
    for(Int_t i=0;i<24;i++){ //print out calibration constants
      printf("Detector %2d: %5.0f | %6.0f |%6.0f | %6.0f | %5.1f | %4.0f | %2.0f |  %3.0f | %3.0f |%3.0f\n",
	     i+1,ECal[i][3],ECal[i][4],ECal[i][5],ECal[i][6],ECal[i][7],ECal[i][8],ECal[i][9],ECal[i][10],ECal[i][11],ECal[i][12]);
    }}
    /*
    printf("Q-Value Constants:\n");
    printf("       Q-Value slope | Q-Value offset\n");
    for(Int_t i=0;i<24;i++){ //print out calibration constants
      printf("Detector %2d: %7.3f | %8.3f \n",i+1,ECal[i][17],ECal[i][18]);
    }
    */
  }
  cout<<endl;
  printf("Gating Information:\n");
  if(DoCut[0]==1)printf("Energy gate applied at %6.1f +/- %f6.1\n",cutE,widthE*sigmaE);else printf("No energy gate applied.\n");
  if(DoCut[1]==1)printf("Position gate applied at x = %5.3f +/- %5.3f\n",cutX,widthX*sigmaX);else printf("No position gate applied.\n");
  if(DoCut[2]==1&&DoCal[2]>3)printf("Time gate applied at %f ns +/- %f ns\n",cutT,widthT*sigmaT);else printf("No time gate applied.\n");
  if(DoCut[3]==1&&DoCal[0]>1)printf("TOF gate applied at %4.1f ns +/- %4.1f ns\n",cutTOF,widthTOF*sigmaTOF);else printf("No TOF gate applied.\n");
  if(DoCut[4]==1&&DoCal[1]>1)printf("ESum gate applied on E=(XF+XN) +/- %4.1f chan\n",widthSum*sigmaSum);else printf("No ESum gate applied.\n");
  if(DoCut[4]==1&&DoCal[1]>1)printf("EDiff gate applied on E>|(XF-XN)| + %4.1f chan\n",widthDiff*sigmaDiff);else printf("No ESum gate applied.\n");
  
  cout<<endl;

//Set histogram ranges, depending on calibration
  if(DoCal[0]>1){
    if(maxZ>0)
    maxE=(Int_t)ceil(slopeEcm*maxZ+intercept); //Set histogram maximum for calibrated 
                                                   //energy plots (MeV)
    else maxE=(Int_t)ceil(intercept);
    maxEc=ceil(maxE-minZ*slopeEcm);
    minQ=floor(intercept-maxEc);      

    //    if(DoCal[1]>1)  maxX=maxE; //Set position scale to "energy" if positions are 
                               //calibrated in MeV
  } 
  
  if(DoCal[2] >3){
    minT=-Tcyc;
    maxT=Tcyc+(Tcyc-minT); //Time calibration centers time plots on Tcyc.
  }

  if(((maxZ-minZ)<382)&&(DoCal[1]>2)){//Shifts Z-plots by offset if range defined by array size
    maxZ-=ECal[0][14];
    minZ-=ECal[0][14];
  }

Int_t bin1=256;//Sets number of bins on most histograms to conveniently reduce memory load

// 2d histograms
  for(int a=0;a<5;++a){
    TString name="hADC";
    TString title="Raw ADC";
    name+=(a+1);
    title+=(a+1);
    hADC[a]=new TH2F(name,title,1024,0,4095,16,0,16);
  }
 
  hE=new  TH2F("hE","Detector Energy (1-24), ungated",1024,0,maxE,24,1,25);
  hXF=new TH2F("hXF","Detector Position (XF), ungated",1024,0,maxX,24,1,25);
  hXN=new TH2F("hXN","Detector Position (XN), ungated",1024,0,maxX,24,1,25);
  hT=new  TH2F("hT","Detector vs. Time, ungated",       1024,minT,maxT,24,1,25);

  hEZ=new TH2F("hEZ","Energy (MeV)  vs. Position (mm), ungated",2048,minZ,maxZ,1024,0,maxE);
    
  for(int a=0;a<24;++a){
    TString name="hXFXN";
    TString title="XF vs. XN detector ";
    name+=(a+1);
    title+=(a+1);
    hXFXN[a]=new TH2F(name,title,bin1,0,maxX,bin1,-maxX/8,maxX);
  }
  
  for(int a=0;a<24;++a){
    TString name="hEX";
    TString title="E vs. 1/2{1+[(XF-XN)/(XF+XN)]} det. ";
    name+=(a+1);
    title+=(a+1);
    hEX[a]=new TH2F(name,title,bin1,-scaleX,1+scaleX,3*bin1,0,maxE);
  }
  return 0;
}

/* function to deal with scalers, adapted from Elliot's program */
void scalers(ScarletEvnt &e)
{// Adapted from Kanter's scaler program
  unsigned int *p = reinterpret_cast<unsigned int*>(e.body());
  unsigned int ttotal, tdiff, ithscaler, ithrate;
    FILE *sf;

    /* On the first sync event after a stop, the following tests for the existence of a file
     * called scalers.zap.  If it exists, the scaler totals are reset.
     */
    if(stopped){
        struct stat st;
        stopped=0;
        if(stat("scalers.zap",&st)==0){
	  for(Int_t i=0;i<NSCALERS;++i)
	    totals[i]=0.0;
        }
    }

    if((sf=fopen("scalers.dat","w"))==0) return;
    ttotal=*p++;
    tdiff=*p++;
    fprintf(sf,"%u %u\n",ttotal,tdiff);
    for(int i=0;i<NSCALERS;++i){
      ithscaler=*p++ & 0x00ffffff;
      totals[i]+=ithscaler;
      ithrate=tdiff!=0 ? ithscaler/tdiff : 0;
      fprintf(sf,"%.0f %u\n",totals[i],ithrate);
    }
    fclose(sf);
}

int userdecode(ScarletEvnt &event){
  ScarletEvnt subevent1;
  Int_t dataword;
  subevent1=event[1];
  int *p1 = reinterpret_cast<int*>(subevent1.body());

  /* The online events will have the form:
   * 
   * TAC data, ADC1 hitpattern, ADC1 data,  ADC2 hitpattern, ADC2 data, ... , 
   * ADC5 data, 0x0000dead

   * The data is parsed into Channel ID and Channel Data and then output into a raw array 
   * scheme (5x16) and then mapped to an array scheme (24x3).  
   */
  Int_t time;
  Int_t nhits[5]={0,0,0,0,0};
  Int_t Chan[5][16];
  Int_t RawData[5][16];
  Int_t Data[24][3];
  Int_t MapDet[5][16]={{ 1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2},
		       {11,10, 9, 8, 7, 6, 5, 4, 7, 6,11,10, 9, 8, 7, 6},
 		       {15,14,13,12,11,10, 9, 8,17,16,15,12,14,13,16,17},
//		       {15,14,13,12,11,10, 9, 8,17,16,15,14,13,12,17,16},  //Note difference
// from straight-cable wiring.
                       {20,21,22,23,18,19,20,21,12,13,14,15,16,17,18,19},
                       {-1,-1,-1,-1,-1,-1,-1,-1,22,23,18,19,20,21,22,23}}; //on ADC5, no 0-7
  Int_t MapSig[5][16]={{ 1, 1, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1},  //0->E, 1->XF, 2->XN
                       { 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1},
                       { 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0},
                       { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0},
		       {-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 2, 2, 2, 2, 2, 2}};

  for(Int_t i=0;i<24;i++){
    Data[i][0]=0;
    Data[i][1]=0;
    Data[i][2]=0;
  }

  //Read in time 
  dataword=*p1++;
  time=(dataword & 0x00000fff);
   
  for(Int_t nadc=0;nadc<5;nadc++){ //loop over ADCs 1-5
    dataword=*p1++;
    nhits[nadc]=cntbit(dataword); //determines number of hits in ADC 
                                  //by counting set bits (ones) in hit register
    for(Int_t i=0;i<nhits[nadc];i++){ //loop over number of ADC hits, if any
      dataword=*p1++;
      
      /*Read in raw ADC data and fill ADC histograms*/
      Chan[nadc][i]=((dataword & 0x0000f000)>>12);
      RawData[nadc][i]=(dataword & 0x00000fff);
      
      hADC[nadc]->Fill(RawData[nadc][i],Chan[nadc][i]);
      
      /*Re-map data from raw (5x16) configuration to detector (24x3) configuration*/
      if(nadc>-1 && // Tests for sensible ADC and Channel numbers
	 nadc<5  &&
	    i>-1 &&
	    i<16 ){
        if((MapDet[nadc][Chan[nadc][i]])>-1&& //Tests for sensible values in re-map matrices.
	   (MapDet[nadc][Chan[nadc][i]])<24&& //Note -1 is excluded.
	   (MapSig[nadc][Chan[nadc][i]])>-1&&
	   (MapSig[nadc][Chan[nadc][i]])<3){  //Re-maps "RawData" array into "Data" array
	  Data[MapDet[nadc][Chan[nadc][i]]][MapSig[nadc][Chan[nadc][i]]]=RawData[nadc][i]; 
	}
      }
    }
 }

  //Done unpacking event, filling raw histograms, and remapping data.
  //Filling histograms with (24x3) detector mapping  
  Float_t t=time;//required to change the raw (integer) time to foating point for calibration
  Float_t e=0,xf=0,xn=0,x=0,z=0;
  Float_t E=0,XF=0,XN=0,X=0,Z=0;
  Float_t Q=0,ch=0;
  Float_t Ecm=0;
  Float_t sum=0,theta=0;
  Float_t V=0,V0=0,Z0=0;
  Float_t TOF=Tcyc;
  Float_t theta2=0;
 
  //Define software thresholds
  Int_t lowthr=75;  //Sets cut-off channel number in detector spectra
  Int_t minTime=28; //Sets cut-off channel number in time spectra 
 
  //Define tags
  Bool_t goodESum=kFALSE;
  Bool_t goodEDiff=kFALSE;
  Bool_t GoodTime=kFALSE;
  Bool_t GoodScat=kFALSE;
  
  for(Int_t i=0;i<24;i++){//Start Calibration and Histogram Fill  
    e=Data[i][0];
    xf=Data[i][1];
    xn=Data[i][2];
    
    if(i==(13-1)){
      xn=1.368*xn;//Eneter slope and intercept of the left-hand side of hEdiff 
      //plot with XF=0, i.e., slope of E=-XN line.
      xf=e-xn;
    }
    
    //Begin Histogram Fill   
    if((e>lowthr)&&(xf>lowthr)&&(xn>lowthr)&&(t>minTime)&&include[i]){ //Tests all three signals against "lowthr"
      Counts[i]=Counts[i]+1; //Stores counts per detector
      
      //Begin Calibration
      //Position Calibration	  
      //Position Calibration Level [1] - Matches XF to XN
      if(DoCal[1]){
	if(ECal[i][2]<-1){ 
	  xn=(-ECal[i][2])*xn;}
	else{
	  xf=(-1/ECal[i][2])*xf;}
        
      //Position Calibration Level [2] - Matches (XF+XN) to E
	xf=(xf*ECal[i][15]+ECal[i][16]/2);
	xn=(xn*ECal[i][15]+ECal[i][16]/2); 
      }
      
      //x=(1/2.)*(1+((2*xf-e)/e)); //position without xn
      //x=(1/2.)*(1+((e-2*xn)/e)); //position without xf
      
      x=(1/2.)*(1+((xf-xn)/(xf+xn))); //Position on detector with XN@x=0 and XF@x=1.  
      //Please note at this point that the array PCBs are 
                                          //wired backwards, so "X-Far" is closest to the 
                                          //target and "X-Near" is further away.
	  

  //Energy Calibration
      //Energy Calibration Level [1] - Correct position-dependance of energy
      if(DoCal[0]&&ECal[i][20]){
	e=e-ECal[i][20]*pow(x+(ECal[i][19]/(2*ECal[i][20])),2);
      }
	 
      
      sum=e-(xf+xn);

      
      /*Fill histograms with energy gating*/
      if((e>(cutE-widthE*sigmaE)&&e<(cutE+widthE*sigmaE))||(DoCut[0]==0)){ //Tests energy is in range OR no energy calibration applied
	//     if(e>(-(xf-xn)+(5*widthDiff*sigmaDiff))&&e>((xf-xn)+(5*widthDiff*sigmaDiff)))
	hXFXN[i]->Fill(xn,xf);
      }
  
    //Energy Calibration
      //Energy Calibration Level [2] - Calibrate Energy in MeV
      ch=e;
      if(e>0&&DoCal[0]>1){
	e =(( e- ECal[i][1])   /ECal[i][0]); //Energy in MeV
	/*
	  if(DoCal[1]>1){
	  xf=((xf-(ECal[i][1]/2))/ECal[i][0]);
	  xn=((xn-(ECal[i][1]/2))/ECal[i][0]);
	  }
	*/
      }
	  
     
      z=(active*x); //position on detector in mm
      Z=-positions[(6-(i%6))]-active/2+positions[0]+z; //position in magnet in mm
      //Position Calibration Level [4] - Offset Correction (Absolute Calibration)
      Z=Z-ECal[0][14]; //Since relative positions are fixed, only one global correction is 
      //needed.  First row value is used.
      if(iter==1){
	printf("Overall Offset is %5.2f mm\n",ECal[0][14]);
	printf("Ta Slits are located at: %7.2f\n", positions[0]);	   
	for(Int_t i=0;i<6;i++){ 
	  printf("Detector %2d Zero Position: %7.2f ",i+1,-positions[(6-(i%6))]-active/2+positions[0]+active-ECal[0][14]);
	
	}
	
      }	  
     
      E=e-slopeEcm*Z; //particle energy in MeV at 90deg in lab
      Q=(intercept-E)*(29.984/28.976); //excitation energy in MeV
      //Q-Value Calibration	  
      if(DoCal[3]){
	Q=((Q-ECal[i][18])/ECal[i][17]); //Q-Value in MeV
      }
      
      V=sqrt(2*e*MeV/mass);//Laboratory Velocity in m/s
      Z0=(e-intercept)/slopeEcm; //beam-axis intercept for given excitation energy
      TOF=Tcyc*Z/Z0; //calculated time-of-flight (TOF)
      
      V0 =sqrt((V*V)+(Vcm*Vcm)-(2*Vcm*(Z/1000)/(Tcyc*1E-9)));//Center of Mass Velocity in m/s
      
      Ecm=(1/2.0*mass*(V0*V0-Vcm*Vcm))/MeV;
      
      theta =180-(acos((V*V-V0*V0-Vcm*Vcm)/(2*V0*Vcm)) )/pi*180;//Center of mass angle in degrees, non-recursive
      theta2=180-(acos(((Z /1000)/(Tcyc*1E-9)-Vcm)/V0 ))/pi*180;//Center of mass angle in degrees, recursive
      
      /*Fill histograms without gating*/
      hE->Fill(e,i+1); 
      hXF->Fill(xf,i+1);
      hXN->Fill(xn,i+1);
      hT->Fill(t,i+1);
      
      if (checkcutg("cTime2D",t,e)) GoodTime=kTRUE;
      //      if (checkcutg("cScat",t,e)) GoodScat=kTRUE;
      
         if((goodESum&&goodEDiff)||DoCut[4]==0||DoCal[1]<2){
	   //  if((goodEDiff)||DoCut[4]==0||DoCal[1]<2){
	
	/*Fill histograms with position gating*/
	if((x>(cutX-widthX*sigmaX)&&(x<cutX+widthX*sigmaX))||(DoCut[1]==0)){
	  
	  hEX[i]->Fill(x,e);
	  hEZ->Fill(Z,e);
	  
	  
	}//end position gate
	
      }
      iter++;
    }//end fill histogram
  }//end calibration
  return 0;
}//end userdecode()

/* The userfunc() function:  This function is called per event.  The event
 * is supplied by daphne.  Unpack the event and fill your histograms here.
 */
int userfunc(const struct ScarletEvntHdr* h)
{
  Float_t CountsSum=0;
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
    printf("Received stop signal.  ");
    for(Int_t i=0;i<24;i++)CountsSum+=Counts[i];
    printf("Run sorted.  Total counts: %1.0f\n",CountsSum);
    //for(Int_t i=0;i<24;i++) printf("Detector %2d: Counts = %10d (%5.2f%%)\n",i+1,Counts[i],(Float_t)((Counts[i]/CountsSum)*100));
    break;
  }
  return 0;
}

/* The userexit() function:  This function is called when the sort thread is stopped.  The sort 
 * thread is stopped either by explicitly stopping it or when a new sort is started.  It is not
 * stopped if the sorting completes or terminates.  By not stopping, the user retains access to
 * the ROOT objects even after a sort finishes.  Typically, the userexit() function is used to 
 * close whatever root files were opened in userentry().
 */
int userexit()
{
  cout<<"Exiting sort..."<<endl;
  //  f->ls();
  f->Write();
  f->Close();
  delete f;
  printf("\a"); //"Default Beep" at sort exit.
  return 0;
}
