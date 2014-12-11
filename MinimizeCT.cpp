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
const Double_t ca_fixed = -1.2701;//set C'A = CA , PDG2012 value
const Double_t alpha = 0.0072973525698;
const Double_t vUD = 0.97425; //I S Towner and J C Hardy 2010 Rep. Prog. Phys. 73 046301

//Includes GF, m_e, f = 1.6887 from Nilkinson NPA 377 (1982) 474 and a factor (1+RC) = 1.03886 from Abele, Prog. P&N Phys. 60 (2008) 1-81
const Double_t nDecayStrength = 4908.5; //GF, f, and rad corrections together. 

const Int_t nDF = 2;

//variables
Double_t ctMIN, ctMAX, cptMIN, cptMAX, ctStep, cT,cpT;
bool floatCA = 1;

Int_t nInput; //number of experimental points 

//vector of all correlation coefficients
vector<correlationCoefficient_t> vInputData;

//the input data will be saved in a Tree
TTree *outTree = new TTree("tExpInputdata","The input data used for the CT C'T scan");


TGraph2D *g = new TGraph2D();
TH2F *hChiSqr;
TH2F *hPChiSqr;
TH2F *hChiSqrBis;
TH2F *hPChiSqrBis;
TH2F *hPNormChiSqrBis;
TH2F *hPDFChiSqrBis;

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

double tSM(double MF,double MGT, double ca, double cs, double cps, double ct, double cpt, double jInit, double jFinal, int zDaughter, int sign) //neutron SM lifetime value
{
  double tInv = vUD*vUD*(MF*MF+(MGT*MGT*ca*ca)/(cv*cv));
  return nDecayStrength/tInv;
}

double PDFValue(double x,int n)
{
  return ( TMath::Exp(-x/2.) * pow(x,(n/2.-1)) ) / ( TMath::Gamma(n/2.) * pow(2,n/2.) );
}

//**********************************
//fill histograms
//***********************************

void SetHisto(double min)
{
   int nPoints = g->GetN();
   cout << "N " << nPoints << " CT :" << cT << "  CpT : " << cpT << "  Chi Min : " << min << endl;   
   g->SetPoint(nPoints, cT, cpT, min);

  //Fill Ct vs Ctp histo
  Int_t bin = hChiSqr->FindBin(cT/ca_fixed,cpT/ca_fixed);
  hChiSqr->SetBinContent(bin,min);
  hPChiSqr->SetBinContent(bin,TMath::Prob(min,nDF)); 
  
  //fill ct+cpt vs ct-cpt histo
  bin = hChiSqrBis->FindBin((cT+cpT)/ca_fixed,(cT-cpT)/ca_fixed);
  hChiSqrBis->SetBinContent(bin,min);
  hPChiSqrBis->SetBinContent(bin,TMath::Prob(min,nDF));
  hPDFChiSqrBis->SetBinContent(bin,PDFValue(min,nDF));

  

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
  string nbFierz("bf");


  for(Int_t i = 0; i < vInputData.size(); i++)
  {
    if(vInputData.at(i).use==1)
    {
      if(littleA.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        error = vInputData.at(i).error;
        expectation = a(vInputData.at(i).MF,vInputData.at(i).MGT, xx[2],xx[0], xx[1], cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
        expectation = expectation / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT, xx[2], xx[0], xx[1], cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
        result = result + pow((value-expectation)/error,2);
        
      }

      if(bigA.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        error = vInputData.at(i).error;
        expectation = A(vInputData.at(i).MF,vInputData.at(i).MGT, xx[2], xx[0], xx[1], cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
        expectation = expectation / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT, xx[2], xx[0], xx[1], cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
        result = result + pow((value-expectation)/error,2);
        
      }

      if(bFierz.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        error = vInputData.at(i).error;
        expectation = b(vInputData.at(i).MF,vInputData.at(i).MGT, xx[2], xx[0], xx[1], cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
        result = result + pow((value-expectation)/error,2);
      }

      if(relPol.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        error = vInputData.at(i).error;
        expectation = vInputData.at(i).mE*( b(1,0, xx[2], xx[0], xx[1], cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, -1 ) + b( 0,1, xx[2], xx[0], xx[1], cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, 1 )) + 1.;
        result = result + pow((value-expectation)/error,2);
      }

      if(nLifeTime.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        error = vInputData.at(i).error;
        expectation = tSM(vInputData.at(i).MF,vInputData.at(i).MGT, xx[2], xx[0], xx[1], cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign) / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT, xx[2], xx[0], xx[1], cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
        result = result + pow((value-expectation)/error,2);
      }

      if(nbFierz.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        error = vInputData.at(i).error;
        expectation = b(vInputData.at(i).MF,vInputData.at(i).MGT, xx[2], xx[0], xx[1], cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
        result = result + pow((value-expectation)/error,2);
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
   ROOT::Math::Functor f(&ChiSqrFunction,3); 
   double step[3] = {0.0001,0.0001,0.0001};

   // starting point
   double variable[3];// = { -1.,1.2};
   if (randomSeed >= 0)
   { 
      TRandom2 r(randomSeed);
      variable[0] = r.Uniform(-0.01,0.01);
      variable[1] = r.Uniform(-0.01,0.01);
      variable[2] = r.Uniform(-1.26,-1.28);
   }
 
   min->SetFunction(f);
 
   // Set the free variables to be minimized!
   min->SetVariable(0,"CS",variable[0], step[0]);
   min->SetVariable(1,"CpS",variable[1], step[1]);
   if(floatCA){min->SetLimitedVariable(2,"Ca",variable[2], step[2],-1.5,-1.);}
   else{min->SetFixedVariable(2,"Ca",ca_fixed);}
   
   // do the minimization
   min->Minimize(); 
 
   const double *xs = min->X();
   std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] << "): " << min->MinValue()  << std::endl;
   
   SetHisto((double)min->MinValue());
   return 0;
}




//*************************
//******read input*********
//*************************

void ProcessInputFile()
{
  
  ifstream inputFile;
  
  inputFile.open("inputdata_floatCA.dat");
  
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

  //Read in CT grid : ctMin, ctMAX, ctpMin, ctpMax, ctStep
  while (std::getline(inputFile,line)) 
  {
    if(line[0]=='#' || line.size()<1) continue;
    else{	std::istringstream(line) >> ctMIN  >> ctMAX  >> cptMIN  >> cptMAX  >> ctStep;	break;      }
  }

  cout << ctMIN << "  " << ctMAX << "  " <<cptMIN << "  " <<cptMAX <<"  " << ctStep << endl;
  
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




void testChiSqr()
{

  double  result = 0.;
  double value, error, expectation;
  //const char* a = "a";
  string littleA ("a");
  string bigA ("A");
  string bFierz ("b");
  string relPol ("R");
  string nLifeTime("tn");
  string nbFierz("bf");

  cT = 0.0001;
  cpT = 0.0001;

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
        expectation = a(vInputData.at(i).MF,vInputData.at(i).MGT, -1.27, 0.000100, 0.0001, cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
        cout << "-----expectation----" << expectation << endl;
        expectation = expectation / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT, -1.27, 0.0001, 0.0001, cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
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
        expectation = A(vInputData.at(i).MF,vInputData.at(i).MGT, -1.27, 0.000100, 0.0001, cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
        cout << "-----expectation----  " << expectation << endl;
        expectation = expectation / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT, -1.27, 0.0001, 0.0001, cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
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
        expectation = b(vInputData.at(i).MF,vInputData.at(i).MGT, -1.27, 0.000100, 0.0001, cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
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
        expectation = vInputData.at(i).mE*( b(1,0, -1.27, 0.000100, 0.0001, cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, -1) + b(0,1,-1.27, 0.000100, 0.0001, cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, 1) )  +1.;
        cout << "-----expectation----  F " << vInputData.at(i).mE*( b(1,0, -1.27, 0.000100, 0.0001, cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, -1)) << endl;
         cout << "-----expectation----  GT" << vInputData.at(i).mE*( b(0,1, -1.27, 0.000100, 0.0001, cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, 1)) << endl;
        result = result + pow((value-expectation)/error,2);
        cout << "******" << vInputData.at(i).parameter << "  " <<vInputData.at(i).isotope << "  " << value << " " << expectation << endl;
      }

      if(nLifeTime.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        cout << "-----value----  " << value << endl;
        error = vInputData.at(i).error;
        cout << "-----error----  " << error << endl;
        expectation = tSM(vInputData.at(i).MF,vInputData.at(i).MGT, -1.27, 0.000100, 0.0001, cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign) / ( 1. + vInputData.at(i).mE*b(vInputData.at(i).MF,vInputData.at(i).MGT, -1.27, 0.0001, 0.0001, cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign));
        cout << "-----expectation----  " << expectation << endl;
        result = result + pow((value-expectation)/error,2);
        cout << "******" << vInputData.at(i).parameter << "  " <<vInputData.at(i).isotope << "  " << value << " " << expectation << endl;
      }

      if(nbFierz.compare(vInputData.at(i).parameter) == 0)
      {
        value = vInputData.at(i).expValue;
        cout << "-----value----  " << value << endl;
        error = vInputData.at(i).error;
        cout << "-----error----  " << error << endl;
        expectation = b(vInputData.at(i).MF,vInputData.at(i).MGT, -1.27, 0.0001, 0.0001, cT,cpT, vInputData.at(i).jInit, vInputData.at(i).jFinal, vInputData.at(i).zDaughter, vInputData.at(i).sign);
        cout << "-----expectation----  " << expectation << endl;
        result = result + pow((value-expectation)/error,2);
        cout << "******" << vInputData.at(i).parameter << "  " <<vInputData.at(i).isotope << "  " << value << " " << expectation << endl;
      }

    }
    cout << i << "---------  " << result << endl;
  }

}




void CalculateConfidenceLevel()
{

  //Normalize P value histogram so the max is 1
  double maxPValue = hPChiSqrBis->GetBinContent(hPChiSqrBis->GetMaximumBin());
  double binContent;
  cout << "maximum p value = " << maxPValue << endl;


  for(Int_t j = 1; j <= hPChiSqrBis->GetNbinsX(); j++)
  {
    for(Int_t k = 1; k <= hPChiSqrBis->GetNbinsY(); k ++)
    {
      binContent = hPChiSqrBis->GetBinContent(j,k);
      hPNormChiSqrBis->SetBinContent(j,k,binContent/maxPValue);
    }
  }

  cout << "Calculate confidence levels ... " << endl;
  //make clone to mess with
  TH2F* hClone = (TH2F*)hPDFChiSqrBis->Clone();

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
        if(hChiSqrBis->GetBinContent(j,k) < chiSqr) hClone->SetBinContent(j,k,0);
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
  TCanvas *c1 = new TCanvas("c1","ChiSqr Limits",800,1200);
  c1->Divide(1,2);
  c1->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogz();
  hChiSqr->Draw("COLZCONT3");
  c1->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hPChiSqr->Draw("COLZCONT3");


  TCanvas *c2 = new TCanvas("c2","ChiSqr Limits left Vs right handed",1000,1000);
  c2->Divide(2,2);
  c2->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogz();
  hChiSqrBis->Draw("COLZCONT3");
  c2->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hPChiSqrBis->Draw("COLZCONT3");

  //TCanvas *c3 = new TCanvas("c3","Confidence limits",800,800);
  //c3->cd();
  c2->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  gCL->SetTitle("Confidence levels");
  gCL->GetYaxis()->SetTitle("Confidence level (%)");
  gCL->GetXaxis()->SetTitle("Chi Sqr value");
  gCL->SetLineWidth(3);
  //gCL->Draw("AC");
  hPDFChiSqrBis->Draw("COLZCONT3");
  

  Int_t bin = hPNormChiSqrBis->GetYaxis()->FindBin(0.0);
  TH1D* hPChiSqrBis_Projected = hPNormChiSqrBis->ProjectionX("hPNormChiSqrBis_Projected",bin,bin);
  c2->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  hPChiSqrBis_Projected->SetTitle("Limits on left-handed currents (Norm P value)");
  hPChiSqrBis_Projected->SetLineWidth(2); 
  hPChiSqrBis_Projected->SetFillColor(45);
  //hPChiSqrBis_Projected->Draw("");
  //hPNormChiSqrBis->Draw("COLZCONT3");
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
  hChiSqrBis->Write();
  hPChiSqr->Write();
  hPChiSqrBis->Write();
  hPDFChiSqrBis->Write();
  gCL->Write();

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

  for(Int_t j = 1; j <= hPDFChiSqrBis->GetNbinsX(); j++)
  {
    for(Int_t i = 1; i <= hPDFChiSqrBis->GetNbinsY(); i ++)
    {
      outfile << i*j << " " << hPDFChiSqrBis->GetXaxis()->GetBinCenter(j) << " " << hPDFChiSqrBis->GetYaxis()->GetBinCenter(i) << " " <<  hPDFChiSqrBis->GetBinContent(j,i) << endl;
    }
    
  }

  cout << hPDFChiSqrBis->GetTitle() << " writen to output file" << endl;

}

//*********************************************
//*******************main program**************
//*********************************************

int Minimize()
{
  
  gStyle->SetPalette(1);


  ProcessInputFile();

char h_name[1024];
  cout << "Ready to make output file ..." << endl;
  cout << "Give output file name: " << endl;
cin >> h_name;
 
  double CT_PLUS_CTp;
  double CT_MINUS_CTp;
  double CT_PLUS_CTpMAX = ctMAX+cptMAX;
  double CT_MINUS_CTpMAX = ctMAX-cptMIN;
  double CT_PLUS_CTpMIN = ctMIN+cptMIN;
  double CT_MINUS_CTpMIN = ctMIN-cptMAX;
  
  int nPlusSteps = (int)((CT_PLUS_CTpMAX-CT_PLUS_CTpMIN)/ctStep);
  int nMinusSteps = (int)((CT_MINUS_CTpMAX-CT_MINUS_CTpMIN)/ctStep);

  

  
  TTimeStamp *time = new TTimeStamp;
  int seconds = time->GetSec();

  int nSteps = (int)(((ctMAX+cptMAX)-(ctMIN+cptMIN))/ctStep);
  int nPSteps = (int)(((ctMAX-cptMIN)-(ctMIN-cptMAX))/ctStep);

  hChiSqrBis = new TH2F("ChiSqrBis","Chi Sqr surface",nPlusSteps+1,-CT_PLUS_CTpMIN/ca_fixed,-CT_PLUS_CTpMAX/ca_fixed,nMinusSteps+1,-CT_MINUS_CTpMIN/ca_fixed,-CT_MINUS_CTpMAX/ca_fixed);
  hChiSqrBis->GetXaxis()->SetTitle("(C_T+C'_T)/C_A");
  hChiSqrBis->GetYaxis()->SetTitle("(C_T-C'_T)/C_A");

  hPChiSqrBis = new TH2F("PChiSqrBis","Chi Sqr P value",nPlusSteps+1,-CT_PLUS_CTpMIN/ca_fixed,-CT_PLUS_CTpMAX/ca_fixed,nMinusSteps+1,-CT_MINUS_CTpMIN/ca_fixed,-CT_MINUS_CTpMAX/ca_fixed);
  hPChiSqrBis->GetXaxis()->SetTitle("(C_T+C'_T)/C_A");
  hPChiSqrBis->GetYaxis()->SetTitle("(C_T-C'_T)/C_A");

  hPNormChiSqrBis = new TH2F("PNormChiSqrBis","Chi Sqr P value (Normalized)",nPlusSteps+1,-CT_PLUS_CTpMIN/ca_fixed,-CT_PLUS_CTpMAX/ca_fixed,nMinusSteps+1,-CT_MINUS_CTpMIN/ca_fixed,-CT_MINUS_CTpMAX/ca_fixed);
  hPNormChiSqrBis->GetXaxis()->SetTitle("(C_T+C'_T)/C_A");
  hPNormChiSqrBis->GetYaxis()->SetTitle("(C_T-C'_T)/C_A");

  hPDFChiSqrBis = new TH2F("PDFChiSqrBis","pdf(#chi^{2})",nPlusSteps+1,-CT_PLUS_CTpMIN/ca_fixed,-CT_PLUS_CTpMAX/ca_fixed,nMinusSteps+1,-CT_MINUS_CTpMIN/ca_fixed,-CT_MINUS_CTpMAX/ca_fixed);
  hPDFChiSqrBis->GetXaxis()->SetTitle("(C_T+C'_T)/C_A");
  hPDFChiSqrBis->GetYaxis()->SetTitle("(C_T-C'_T)/C_A");

  cout << "Plus MAX " << -CT_PLUS_CTpMAX/ca_fixed <<" Plus Min" << -CT_PLUS_CTpMIN/ca_fixed <<endl;
  cout << "Minus MAX " << -CT_MINUS_CTpMAX/ca_fixed <<" Minus Min" << -CT_MINUS_CTpMIN/ca_fixed <<endl;

  //Sca_fixednGrid
  nSteps = (int)((ctMAX - ctMIN)/ctStep);
  nPSteps = (int)((cptMAX - cptMIN)/ctStep);

  //Initialize output histogram
  hChiSqr = new TH2F("ChiSqr","Chi Sqr surface",nSteps+1,-ctMIN/ca_fixed,-ctMAX/ca_fixed,nPSteps+1,-cptMIN/ca_fixed,-cptMAX/ca_fixed);
  hChiSqr->GetXaxis()->SetTitle("C_T/C_A");
  hChiSqr->GetYaxis()->SetTitle("C'_T/C_A");

  hPChiSqr = new TH2F("PChiSqr","Chi Sqr probability surface",nSteps+1,-ctMIN/ca_fixed,-ctMAX/ca_fixed,nPSteps+1,-cptMIN/ca_fixed,-cptMAX/ca_fixed);
  hPChiSqr->GetXaxis()->SetTitle("C_T/C_A");
  hPChiSqr->GetYaxis()->SetTitle("C'_T/C_A");



  testChiSqr();

  for(int i = 0;i<=nPlusSteps+1;i++)
  {
    for(int j = 0; j <= nMinusSteps+1; j++)
    {
      CT_PLUS_CTp = CT_PLUS_CTpMIN + i*ctStep;
      CT_MINUS_CTp = CT_MINUS_CTpMIN + j*ctStep;
      
      //set ct and spr
      cT = (CT_PLUS_CTp + CT_MINUS_CTp)/2.;
      cpT = (CT_PLUS_CTp - CT_MINUS_CTp)/2.;
      //minimize chisquare with cs and cps floating
      NumericalMinimization("Minuit","Migrad",seconds);
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
