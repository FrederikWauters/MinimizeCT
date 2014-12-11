//Modify CT C'T to a CS CT scan (C'S(T) == CS(T))

//Frederik Wauters, August 2012
//Took an example from the cern/ROOT website, plan is to use it as an alternative for Alejandro's FORTRAN based "amoebe" method
// Addapt it to calculate limits on CT and C'T  for arbitrary input

//Nov 2014
//Float Vud and introduce quadratic terms for CS in neutron lifetime and Ft values of 0+ 0+

// Example on how to use the new Minimizer class in ROOT
//  Show usage with all the possible minimizers. 
// Minimize the Rosenbrock function (a 2D -function)
// This example is described also in 
// http://root.cern.ch/drupal/content/numerical-minimization#multidim_minim
// input : minimizer name + algorithm name
// randomSeed: = <0 : fixed value: 0 random with seed 0; >0 random with given seed 
//
//Author: L. Moneta Dec 2010

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TGraph.h"
#include "TMath.h"
#include "TGraph2D.h"
#include "TFile.h"
#include<TGraph.h>
#include "TH1.h"
#include "TH2.h"
#include "TDatime.h"
#include "TTimeStamp.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"


using namespace std;



//////////////////////////////////////////////////////////////////
//data structure to collect all info of a correlation coefficient
/////////////////////////////////////////////////////////////////

typedef struct correlationCoefficient_t{
  //const char* isotope;
  string isotope;
  //const char* parameter;
  string parameter;
  double jInit;
  double jFinal;
  int zDaughter;
  int sign; //is + for beta- and - for beta+!
  double MF;
  double MGT;
  double mE;
  double expValue;
  double error;
  bool use;
};

//////////////////////////////////
//global variables
//////////////////////////////////


//constants
const Double_t cv = 1.;// C'V = CV
const Double_t ca_fixed = -1.2763;//set C'A = CA , PDG2012 value
const Double_t alpha = 0.0072973525698;
const Double_t vUD = 0.97425; //I S Towner and J C Hardy 2010 Rep. Prog. Phys. 73 046301, Now floating!

//Includes GF, m_e, f = 1.6887 from Nilkinson NPA 377 (1982) 474 and a factor (1+RC) = 1.03886 from Abele, Prog. P&N Phys. 60 (2008) 1-81  
const Double_t nDecayStrength = 4908.5; //GF, f, and rad corrections together.
const Double_t errorNDecayStrength = 1.9; 
const Double_t superDecayStrength = 2915.64; //Towner and Hardy, ArXiv 2013, check
const Double_t errorSuperDecayStrength = 1.08;
const Int_t nDF = 2;

//variables
Double_t ctMIN, ctMAX, csMIN, csMAX, ctStep, cT,cS;
bool floatCA = 1;

Int_t nInput; //number of experimental points 

//vector of all correlation coefficients
vector<correlationCoefficient_t> vInputData;

//the input data will be saved in a Tree
TTree *outTree = new TTree("tExpInputdata","The input data used for the CT C'T scan");


TGraph2D *g = new TGraph2D();
TH2F *hChiSqr;
TH2F *hPChiSqr;
TH2F *hPNormChiSqr;
TH2F *hPDFChiSqr;
TH2F *hCA;
TH2F *hVud;

 const Int_t nPoints = 200;
 const Double_t step = 0.2;
 Double_t xValues[nPoints];
 Double_t yValues[nPoints];

 TGraph *gCL; 

//******************************************************************************************
//****************functions to calculate correlation coefficient ***************************
//******************************************************************************************

//void SetOutput(double min)

double gamma(int zDaughter)
{
   return pow(1-alpha*alpha*zDaughter*zDaughter,0.5); 
}

double lambda(double jInit, double jFinal)
{
  if( jFinal == jInit-1 )  return 1.; 
  if( jFinal == jInit)  return (1./(1.+(double)jInit));
  else return  (-1./(1.+(double)jInit));
}

double delta(double jInit, double jFinal)
{
  if( jFinal == jInit) return 1.;
  else return 0.;
}


//*******************************************************************************************
//*****************************Correlation coefficients in terms of CS and CT****************
//*******************************************************************************************

double Xi(double MF, double MGT, double ca,  double cs, double cps, double ct, double cpt, double jInit, double jFinal, int zDaughter, int sign)
{
  return pow(MF,2)*( cs*cs + cv*cv + cps*cps + cv*cv) + pow(MGT,2)*( ct*ct + ca*ca + cpt*cpt + ca*ca ); 
}


double aXi(double MF,double MGT, double ca, double cs, double cps, double ct, double cpt, double jInit, double jFinal, int zDaughter, int sign)
{
  return pow(MF,2)*( -cs*cs + cv*cv + -cps*cps + cv*cv ) + (pow(MGT,2)/3)*( ct*ct - ca*ca + cpt*cpt - ca*ca );
}


double bXi(double MF,double MGT, double ca, double cs, double cps, double ct, double cpt, double jInit, double jFinal, int zDaughter, int sign)
{
  //cout << " sign " << sign << " zDaughter " << zDaughter << "  gamma" << gamma(zDaughter) << endl;
  return sign*2*gamma(zDaughter)*(MF*MF*( cs*cv + cps*cv ) + MGT*MGT*( ct*ca + cpt*ca ));
}


double AXi(double MF,double MGT, double ca, double cs, double cps, double ct, double cpt, double jInit, double jFinal, int zDaughter, int sign)
{
  return MGT*MGT*lambda(jInit,jFinal)*( 2*sign*( ct*cpt - ca*ca) ) + delta(jInit,jFinal)*MF*MGT*sqrt((double)jInit/(1.+(double)jInit))*( 2*( cs*cpt + cps*ct - 2*ca*cv ));
}


double a(double MF,double MGT, double ca, double cs, double cps, double ct, double cpt, double jInit, double jFinal, int zDaughter, int sign)
{
  return aXi(MF,MGT,ca,cs,cps,ct,cpt,jInit,jFinal,zDaughter,sign)/Xi(MF,MGT,ca,cs,cps,ct,cpt,jInit,jFinal,zDaughter,sign);
}


double b(double MF,double MGT, double ca, double cs, double cps, double ct, double cpt, double jInit, double jFinal, int zDaughter, int sign)
{
  return bXi(MF,MGT,ca,cs,cps,ct,cpt,jInit,jFinal,zDaughter,sign)/Xi(MF,MGT,ca,cs,cps,ct,cpt,jInit,jFinal,zDaughter,sign);
}


double A(double MF,double MGT, double ca, double cs, double cps, double ct, double cpt, double jInit, double jFinal, int zDaughter, int sign)
{
  //cout << " AXi " << AXi(MF, MGT,cs,cps,ct,cpt, jInit, jFinal, zDaughter, sign) << " Xi " << Xi(MF, MGT,cs,cps,ct,cpt, jInit, jFinal, zDaughter, sign) << endl;  
  return AXi(MF, MGT,ca,cs,cps,ct,cpt, jInit, jFinal, zDaughter, sign)/Xi(MF, MGT,ca,cs,cps,ct,cpt, jInit, jFinal, zDaughter, sign);
}

vector<double> tN(double MF,double MGT, double ca, double cs, double cps, double ct, double cpt, double jInit, double jFinal, int zDaughter, int sign, double vud) //neutron lifetime value, not fierz corrected yet
{
  vector<double> vOut; //contains value and error
  double tInv = vud*vud*( MF*MF*(1+cs*cs) + ( MGT*MGT*((ct*ct+ca*ca)/(cv*cv)) ) );
  vOut.push_back(nDecayStrength/tInv);
  vOut.push_back(errorNDecayStrength/vOut.at(0));
  return vOut;
}

vector<double> tSuper(double MF,double MGT, double ca, double cs, double cps, double ct, double cpt, double jInit, double jFinal, int zDaughter, int sign, double vud)
{
  vector<double> vOut; //contains value and error
  double tInv = vud*vud*MF*MF*( 1 + (cs*cs)/(cv*cv) );
  //double tInv = vud*vud*MF*MF;//*( 1 + (cs*cs)/(cv*cv) );
  vOut.push_back(superDecayStrength/tInv);
  vOut.push_back( errorSuperDecayStrength/vOut.at(0) );
  return vOut;
}

double PDFValue(double x,int n)
{
  return ( TMath::Exp(-x/2.) * pow(x,(n/2.-1)) ) / ( TMath::Gamma(n/2.) * pow(2,n/2.) );
}

//**********************************
//fill histograms
//***********************************

void SetHisto(double min, double ca, double vud)
{
   int nPoints = g->GetN();
   //cout << "N " << nPoints << " CT :" << cT << "  CS : " << cS << "  Chi Min : " << min << endl;   
   g->SetPoint(nPoints, cT, cS, min);

  //Fill Ct vs Ctp histo
  Int_t bin = hChiSqr->FindBin(cT/ca_fixed,cS/cv);
  hChiSqr->SetBinContent(bin,min);
  hPChiSqr->SetBinContent(bin,TMath::Prob(min,nDF)); 
  hPDFChiSqr->SetBinContent(bin,PDFValue(min,nDF));
  hCA->SetBinContent(bin,ca);
  hVud->SetBinContent(bin,vud);

  

}


//**********************************
//Minimization *********************
//**********************************

//The function you want to minimize
double ChiSqrFunction(const double *xx )
{

  double  result = 0.;
  double value, error, expectation;
  //const char* a = "a";
  string littleA ("a");
  string bigA ("A");
  string bFierz ("b");
  string relPol ("R");
  string nLifeTime("tn");
  string FTSupperAllowed("Ft");


  for(Int_t i = 0; i < vInputData.size(); i++)
  {
    if(vInputData.at(i).use==1)
    {
      if(littleA.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        error = vInputData.at(i).error;
        expectation = a(vInputData.at(i).MF,vInputData.at(i).MGT, xx[0],cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
        expectation = expectation / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT, xx[0],cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
        result = result + pow((value-expectation)/error,2);
        
      }

      if(bigA.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        error = vInputData.at(i).error;
        expectation = A(vInputData.at(i).MF,vInputData.at(i).MGT, xx[0],cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
        expectation = expectation / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT, xx[0],cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
        result = result + pow((value-expectation)/error,2);
        
      }

      if(bFierz.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        error = vInputData.at(i).error;
        expectation = b(vInputData.at(i).MF,vInputData.at(i).MGT, xx[0],cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
        result = result + pow((value-expectation)/error,2);
      }

      if(relPol.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        error = vInputData.at(i).error;
        expectation = vInputData.at(i).mE*( b(1,0, xx[0],cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, -1 ) + b( 0,1, xx[0],cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, 1 )) + 1.;
        result = result + pow((value-expectation)/error,2);
      }

      if(nLifeTime.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        error = vInputData.at(i).error;
        vector<double> expectationsANDerror = tN(vInputData.at(i).MF,vInputData.at(i).MGT,xx[0],cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign,xx[1]);
        expectation =  expectationsANDerror.at(0) / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT,xx[0],cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
        double errorSq = pow(error,2) +  pow(expectationsANDerror.at(1)/expectation,2);
        result = result + pow((value-expectation),2)/errorSq;
      }

      if(FTSupperAllowed.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        error = vInputData.at(i).error;
        vector<double> expectationsANDerror = tSuper(vInputData.at(i).MF,vInputData.at(i).MGT,xx[0],cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign,xx[1]);
        expectation = expectationsANDerror.at(0) / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT,xx[0],cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
        double errorSq = pow(error,2) +  pow(expectationsANDerror.at(1)/expectation,2);
        result = result + pow((value-expectation),2)/errorSq;
      }

    }
  }

  return result;
}
 

 int NumericalMinimization(const char * minName = "Minuit2", const char *algoName = "" , int randomSeed = -1)
{ 
   // create minimizer giving a name and a name (optionally) for the specific
   // algorithm
   // possible choices are: 
   //     minName                  algoName
   // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
   //  Minuit2                     Fumili2
   //  Fumili
   //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS, 
   //                              BFGS2, SteepestDescent
   //  GSLMultiFit
   //   GSLSimAn
   //   Genetic

   ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

   // set tolerance , etc...
   min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
   min->SetMaxIterations(10000);  // for GSL 
   min->SetTolerance(0.0001);
   min->SetPrintLevel(-1);

   // create funciton wrapper for minmizer
   // a IMultiGenFunction type 
   ROOT::Math::Functor f(&ChiSqrFunction,2); 
   double step[2] = {0.0001,0.00001};


   // starting point
   double variable[2];// = { -1.,1.2};
   if (randomSeed >= 0)
   { 
      TRandom2 r(randomSeed);
      variable[0] = r.Uniform(-1.26,-1.28);
      variable[1] = r.Uniform(0.95,1.05);
   }

   min->SetFunction(f);
   // Set the free variables to be minimized!
   if(floatCA){min->SetLimitedVariable(0,"Ca",variable[0], step[0],-1.5,-1.);}
   else{min->SetFixedVariable(0,"Ca",ca_fixed);}
   min->SetLimitedVariable(1,"Vud",variable[1],step[1],0.7,1.3);

   // do the minimization
   min->Minimize(); 
 
   const double *xs = min->X();
   //std::cout << "Minimum: f(" << xs[0] << xs[1] << " ): " << min->MinValue()  << std::endl;
   
   SetHisto((double)min->MinValue(),xs[0],xs[1]);
   return 0;
}




//*************************
//******read input*********
//*************************

void ProcessInputFile()
{
  
  ifstream inputFile;
  
  inputFile.open("inputdata.dat");
  
  if (!inputFile.is_open())  cout << "error while opening file" << endl;

  cout << "Start reading input file ... " << endl;
  
  string line;
  int testint;

  //Read bool whether to float or fix CA
  while (std::getline(inputFile,line)) 
  {
    if(line[0]=='#' || line.size()<1) continue;
    else{	std::istringstream(line) >> floatCA;	break;      }
  }

  //Read in CS grid : ctMin, ctMAX, csMin, csMax, Step
  while (std::getline(inputFile,line)) 
  {
    if(line[0]=='#' || line.size()<1) continue;
    else{	std::istringstream(line) >> ctMIN  >> ctMAX  >> csMIN  >> csMAX  >> ctStep;	break;      }
  }

  cout << ctMIN << "  " << ctMAX << "  " <<csMIN << "  " <<csMAX <<"  " << ctStep << endl;
  
  //Get correlation coefficients and put them in a vector. Read data until eof
  struct correlationCoefficient_t input;
//put the input data struct in a tree, which will be written to the output file


  char iso[10];// = "Ca";
  char par[10];
  double jI = -1.0;
  double jF = -1.0;
  int zD = 0;
  //char plmn[1];
  string plmn;
  double mFermi;
  double mGamTel;
  double avME = 0.0;
  double expValue = 0.0;
  double error = 0.0;
  int useIt = 0;

  string plus("+");
  string minus("-");

  nInput = 0;

  //outTree->Branch("input",&input,"isotope/C:parameter/C:jInit/D:jFinal/D:zDaughter/I:sign/I:MF/D:MGT/D:mE/D:expValue/D:error/D:use/O");
  outTree->Branch("isotope",&iso,"iso/C");
  outTree->Branch("parameter",&par,"par/C");
  //for the first two, i write char, not strings, tree's seem not to like it
  outTree->Branch("jInit",&input.jInit,"input.jInit/D");
  outTree->Branch("jFinal",&input.jFinal,"input.jFinal/D");
  outTree->Branch("zDaughter",&input.zDaughter,"input.zDaughter/I");
  outTree->Branch("sign",&input.sign,"input.sign/I");
  outTree->Branch("MF",&input.MF,"input.MF/D");
  outTree->Branch("MGT",&input.MGT,"input.MGT/D");
  outTree->Branch("mE",&input.mE,"input.mE/D");
  outTree->Branch("expValue",&input.expValue,"input.expValue/D");
  outTree->Branch("error",&input.error,"input.error/D");
  outTree->Branch("use",&input.use,"input.use/O");



  while(std::getline(inputFile,line))
  {
    if(line[0]=='#' || line.size()<2) continue;
    
    else  
    {
      stringstream ss;
      ss << line;
      ss >> iso >> par >> jI  >> jF  >> zD >> plmn >> mFermi >> mGamTel >> avME >> expValue >>  error >> useIt;
       

      input.isotope = (char*)iso;
      input.parameter = (char*)par;
      input.jInit = jI;
      input.jFinal = jF;
      input.zDaughter = zD;
      if(plus.compare(plmn)==0) input.sign = -1; else input.sign = +1;
      input.MF = mFermi;
      input.MGT = mGamTel;
      input.mE = avME;
      input.expValue = expValue;
      input.error = error;
      input.use = useIt;

      vInputData.push_back(input);

      cout <<  vInputData.at(nInput).isotope << "  "<< vInputData.at(nInput).parameter  << "  "<< vInputData.at(nInput).jInit  << "  "<< vInputData.at(nInput).jFinal << "  " << vInputData.at(nInput).zDaughter << "  " << vInputData.at(nInput).sign << "  "<< vInputData.at(nInput).MF << " " <<vInputData.at(nInput).MGT << "  "<< vInputData.at(nInput).mE << "  "<< vInputData.at(nInput).expValue<< "  " << vInputData.at(nInput).error << " " << vInputData.at(nInput).use << endl;

      //Put experimental data in the tree which goes in the output file
      if(input.use)
      {
        outTree->Fill();
      }
      nInput++;
    }

  }

}




void testChiSqr(double ca_test)
{

  double  result = 0.;
  double value, error, expectation;
  //const char* a = "a";
  string littleA ("a");
  string bigA ("A");
  string bFierz ("b");
  string relPol ("R");
  string nLifeTime("tn");
  string FTSupperAllowed("Ft");

  cT = 0.002;
  cS = 0.0011;


  for(Int_t i = 0; i< vInputData.size(); i++)
  {
    if(vInputData.at(i).use==1)
    {
      if(littleA.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        cout << "-----value----" << value << endl;
        error = vInputData.at(i).error;
        cout << "-----error----" << error << endl;
        expectation = a(vInputData.at(i).MF,vInputData.at(i).MGT, ca_test, cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
        cout << "-----expectation----" << expectation << endl;
        expectation = expectation / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT, ca_test, cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
        cout << "-----expectation----" << expectation << endl;
        result = result + pow((value-expectation)/error,2);
        cout << "******" << vInputData.at(i).parameter << "  " <<vInputData.at(i).isotope << "  " << value << " " << expectation << endl;
      }

      if(bigA.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        cout << "-----value----  " << value << endl;
        error = vInputData.at(i).error;
        cout << "-----error----  " << error << endl;
        expectation = A(vInputData.at(i).MF,vInputData.at(i).MGT, ca_test, cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
        cout << "-----expectation----  " << expectation << endl;
        expectation = expectation / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT,ca_test, cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
        cout << "-----expectation----  " << expectation << endl;
        result = result + pow((value-expectation)/error,2);
        cout << "******" << vInputData.at(i).parameter << "  " <<vInputData.at(i).isotope << "  " << value << " " << expectation << endl;
      }

      if(bFierz.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        cout << "-----value----  " << value << endl;
        error = vInputData.at(i).error;
        cout << "-----error----  " << error << endl;
        expectation = b(vInputData.at(i).MF,vInputData.at(i).MGT, ca_test, cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
        cout << "-----expectation----  " << expectation << endl;
        result = result + pow((value-expectation)/error,2);
        cout << "******" << vInputData.at(i).parameter << "  " <<vInputData.at(i).isotope << "  " << value << " " << expectation << endl;
      }

      if(relPol.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        cout << "-----value----  " << value << endl;
        error = vInputData.at(i).error;
        cout << "-----error----  " << error << endl;
        expectation = vInputData.at(i).mE*( b(1,0, -1.27, cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, -1) + b(0,1, -1.27, cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, 1) )  +1.;
        cout << "-----expectation----  F " << vInputData.at(i).mE*( b(1,0, ca_test, cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, -1)) << endl;
         cout << "-----expectation----  GT" << vInputData.at(i).mE*( b(0,1, ca_test, cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, 1)) << endl;
        result = result + pow((value-expectation)/error,2);
        cout << "******" << vInputData.at(i).parameter << "  " <<vInputData.at(i).isotope << "  " << value << " " << expectation << endl;
      }

      if(nLifeTime.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        cout << "-----value----  " << value << endl;
        error = vInputData.at(i).error;
        cout << "-----error----  " << error << endl;
        vector<double> expectationsANDerror = tN(vInputData.at(i).MF,vInputData.at(i).MGT, ca_test, cS, cS, cT,cT,  vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign,vUD);
        expectation =  expectationsANDerror.at(0) / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT,ca_test, cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
        cout << "-----expectation----  " << expectation << endl;
        double errorSq = pow(error,2) +  pow(expectationsANDerror.at(1)/expectation,2);
         result = result + pow((value-expectation),2)/errorSq;
        cout << "******" << vInputData.at(i).parameter << "  " <<vInputData.at(i).isotope << "  " << value << " " << expectation << endl;       
      }

      if( FTSupperAllowed.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        cout << "-----value----  " << value << endl;
        error = vInputData.at(i).error;
        cout << "-----error----  " << error << endl;
        vector<double> expectationsANDerror = tSuper(vInputData.at(i).MF,vInputData.at(i).MGT, ca_test, cS, cS, cT,cT,  vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign,vUD);
        expectation = expectationsANDerror.at(0) / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT,ca_test, cS, cS, cT,cT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
        cout << "-----expectation----  " << expectation << endl;
        double errorSq = pow(error,2) +  pow(expectationsANDerror.at(1)/expectation,2);
        result = result + pow((value-expectation),2)/errorSq;
        cout << "******" << vInputData.at(i).parameter << "  " <<vInputData.at(i).isotope << "  " << value << " " << expectation << endl;
      }


    }
    cout << i << "---------  " << result << endl;
  }

}




void CalculateConfidenceLevel()
{

  //Normalize P value histogram so the max is 1
  double maxPValue = hPChiSqr->GetBinContent(hPChiSqr->GetMaximumBin());
  double binContent;
  cout << "maximum p value = " << maxPValue << endl;


  cout << "Calculate confidence levels ... " << endl;
  //make clone to mess with
  TH2F* hClone = (TH2F*)hPDFChiSqr->Clone();

  double totalInt = hClone->Integral();
  cout << "total Int : " << totalInt <<endl;

  double newInt;
  double chiSqr;

  for(Int_t i = 1; i <= nPoints ;i++)
  {
    chiSqr = i*step;
    xValues[i-1] = chiSqr;

    for(Int_t j = 1; j <= hClone->GetNbinsX(); j++)
    {
      for(Int_t k = 1; k <= hClone->GetNbinsY(); k ++)
      {
        if(hChiSqr->GetBinContent(j,k) < chiSqr) hClone->SetBinContent(j,k,0);
      }
    }
    newInt = hClone->Integral();
    yValues[i-1] = 100-newInt/totalInt*100.;
  }


  gCL = new TGraph(nPoints,xValues,yValues);
}





//****************************
//***********Set Output ******
//****************************
void SetOutput(char* h_name)
{

  //Draw stuff



  TCanvas *c2 = new TCanvas("c2","ChiSqr Limits left Vs right handed",1000,1000);
  c2->Divide(2,2);
  c2->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogz();
  hChiSqr->Draw("COLZCONT3");
  c2->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hPChiSqr->Draw("COLZCONT3");

  //TCanvas *c3 = new TCanvas("c3","Confidence limits",800,800);
  //c3->cd();
  c2->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();

  //gCL->Draw("AC");
  hPDFChiSqr->Draw("COLZCONT3");
  

 // Int_t bin = hPNormChiSqr->GetYaxis()->FindBin(0.0);
 // TH1D* hPChiSqr_Projected = hPNormChiSqr->ProjectionX("hPNormChiSqr_Projected",bin,bin);
  c2->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  //hPChiSqr_Projected->SetTitle("Limits on left-handed currents (Norm P value)");
  //hPChiSqr_Projected->SetLineWidth(2); 
  //hPChiSqr_Projected->SetFillColor(45);
  //hPChiSqrBis_Projected->Draw("");
  //hPNormChiSqrBis->Draw("COLZCONT3");
  gCL->SetTitle("Confidence levels");
  gCL->GetYaxis()->SetTitle("Confidence level (%)");
  gCL->GetXaxis()->SetTitle("Chi Sqr value");
  gCL->SetLineWidth(3);
  gCL->Draw("AC");

  //output to screen
  cout << endl << endl << "****** Input data: ********" << endl;
  for(Int_t i = 0; i< vInputData.size(); i++)
  {
    cout <<  vInputData.at(i).isotope << "  "<< vInputData.at(i).parameter  << "  "<< vInputData.at(i).jInit  << "  "<< vInputData.at(i).jFinal << "  " << vInputData.at(i).zDaughter << "  " << vInputData.at(i).sign << "  "<< vInputData.at(i).MF << " " <<vInputData.at(i).MGT << "  "<< vInputData.at(i).mE << "  "<< vInputData.at(i).expValue<< "  " << vInputData.at(i).error << " " << vInputData.at(i).use << endl;
  }

  //***************************************************
  //***********Put histograms in ROOT file*************
  //*************************************************** 

  printf("Creating File\n");

  TFile *fin = new TFile ( h_name, "ReCreate");

  if ( fin->IsZombie() ) {
    cerr << "**ERROR! Cannot open file [" << h_name << "]" << endl;
    return;
  }

  fin->cd();

  hChiSqr->Write();
  hPChiSqr->Write();
  hPDFChiSqr->Write();
  gCL->Write();
  hVud->Write();
  hCA->Write();

  //put the input data in root 
  outTree->Write();

  delete fin;
}


void GenerateTextOutput()
{
  //***********************************************
  //*************text output file *****************
  //**********************************************

  //ofstream *os = new ofstream("output.txt");
  cout << "Set textfile output " << endl;
  ofstream outfile ( "output.txt", ios::trunc );
  cout << "#It.  CT  CS   PDF " << endl;

  for(Int_t j = 1; j <= hPDFChiSqr->GetNbinsX(); j++)
  {
    for(Int_t i = 1; i <= hPDFChiSqr->GetNbinsY(); i ++)
    {
      outfile << i*j << " " << hPDFChiSqr->GetXaxis()->GetBinCenter(j) << " " << hPDFChiSqr->GetYaxis()->GetBinCenter(i) << " " <<  hPDFChiSqr->GetBinContent(j,i) << endl;
    }
    
  }

  cout << hPDFChiSqr->GetTitle() << " written to output file" << endl;

}

//*********************************************
//*******************main program**************
//*********************************************

int Minimize()
{
  
  gStyle->SetPalette(1);


  ProcessInputFile();

char h_name[1024];
  cout << "Ready to scan Cs and CT ..." << endl;
  cout << "Give output file name: " << endl;
cin >> h_name;
 

  

  
  TTimeStamp *time = new TTimeStamp;
  int seconds = time->GetSec();


  //Sca_fixednGrid
  int nSteps = (int)((ctMAX - ctMIN)/ctStep);
  int nPSteps = (int)((csMAX - csMIN)/ctStep);

  //Initialize output histogram
  hChiSqr = new TH2F("ChiSqr","Chi Sqr surface",nSteps+1,-ctMIN/ca_fixed,-ctMAX/ca_fixed,nPSteps+1,csMIN/cv,csMAX/cv);
  hChiSqr->GetXaxis()->SetTitle("C_T/C_A");
  hChiSqr->GetYaxis()->SetTitle("C_S/C_V");

  hPChiSqr = new TH2F("PChiSqr","Chi Sqr P value surface",nSteps+1,-ctMIN/ca_fixed,-ctMAX/ca_fixed,nPSteps+1,csMIN/cv,csMAX/cv);
  hPChiSqr->GetXaxis()->SetTitle("C_T/C_A");
  hPChiSqr->GetYaxis()->SetTitle("C_S/C_V");

  hPDFChiSqr = new TH2F("PDFChiSqr","Chi Sqr probability surface",nSteps+1,-ctMIN/ca_fixed,-ctMAX/ca_fixed,nPSteps+1,csMIN/cv,csMAX/cv);
  hPDFChiSqr->GetXaxis()->SetTitle("C_T/C_A");
  hPDFChiSqr->GetYaxis()->SetTitle("C_S/C_V");

  hCA = new TH2F("hCA","Minimized CA values",nSteps+1,-ctMIN/ca_fixed,-ctMAX/ca_fixed,nPSteps+1,csMIN/cv,csMAX/cv);
  hCA->GetXaxis()->SetTitle("C_T/C_A");
  hCA->GetYaxis()->SetTitle("C_S/C_V");


  hVud = new TH2F("hVud","Minimized Vud values",nSteps+1,-ctMIN/ca_fixed,-ctMAX/ca_fixed,nPSteps+1,csMIN/cv,csMAX/cv);
  hVud->GetXaxis()->SetTitle("C_T/C_A");
  hVud->GetYaxis()->SetTitle("C_S/C_V");


  testChiSqr(-1.2753);
  sleep(10);
  for(int i = 0;i<=nSteps+1;i++)
  {
    for(int j = 0; j <= nSteps+1; j++)
    {
      //set ct and spr
      cT = ctMIN+i*ctStep;
      cS = csMIN+j*ctStep;
      //minimize chisquare with cs and cps floating
      NumericalMinimization("Minuit","Combined",seconds);
    }
 }

  //cout << "seconds " << seconds << endl;

  //NumericalMinimization("Minuit","Migrad",seconds)
  
  cout << "Minimizing is finished " << endl;
  
  CalculateConfidenceLevel();

  SetOutput(h_name);

  //Minimizer options:
  //     minName                  algoName
  // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
  //  Minuit2                     Fumili2
  //  Fumili
  //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS, 
  //                              BFGS2, SteepestDescent
  //  GSLMultiFit
  //   GSLSimAn
  //   Genetic



  return 1;
}
