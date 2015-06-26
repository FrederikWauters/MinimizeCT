/**

This class performs the minimazation

Minimizer method taken from http://root.cern.ch/drupal/content/numerical-minimization 

**/


//Std includes
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cmath>


#include "Analyzer.h"
#include "Globals.h"
#include "DataManager.h"
#include "TData.h"

#include "Math/Factory.h"
#include "Math/Functor.h"

#include "TMath.h"
#include "TRandom2.h"

using std::cout;
using std::endl;
using std::vector;

//DECAY_CONSTANTS constants;
extern DECAY_CONSTANTS constants;

Analyzer::Analyzer( int dof , DataManager*  datamanager,const char* minName = "Minuit2", const char* algoName = "" )
{
  cout << "Construct analyzer ... " << endl;

  manager = datamanager;
  InitMinimizer(dof, minName,algoName);
  cv=constants.cv;
  ca_SM=constants.ca_fixed[0]; // only used to normalized ct+ctp

  data = manager->GetData();
}

Analyzer::~Analyzer()
{
  cout << "delete analyzer ... " << endl;
}

int Analyzer::InitMinimizer(int dof, const char * minName = "Minuit2", const char *algoName = "" )
{
    cout << "Initialize Minimizer" << endl;
    minimizer = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

    // set tolerance , etc...
    minimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
    minimizer->SetMaxIterations(10000);  // for GSL 
    minimizer->SetTolerance(0.0001);
    minimizer->SetPrintLevel(-1);
    
    ROOT::Math::Functor f2(this,&Analyzer::ChiSqrFunction2,2);
    ROOT::Math::Functor f2bis(this,&Analyzer::ChiSqrFunction2bis,2);
    ROOT::Math::Functor f4(this,&Analyzer::ChiSqrFunction4,4);
    ROOT::Math::Functor f4bis(this,&Analyzer::ChiSqrFunction4bis,4);  //for CS, C'S scan
    
    
    if(dof==2) functor = f2;
    else if(dof==-2) functor = f2bis;
    else if(dof==4) functor = f4;
    else if(dof==-4) functor = f4bis;
    minimizer->SetFunction(functor); 
   
    InitVariables(dof);
    
  return 1;
}

void Analyzer::InitVariables(int dof)
{
  int degree_of_freedom = abs(dof);
  double step[degree_of_freedom]; 
  double variable[degree_of_freedom];// = ca and Vud
    
  TRandom2 r(0);
  if(dof==2 || dof == -2)
  {
    step[0] = 0.0001; step[1] = 0.00001;
    variable[0] = r.Uniform(-1.26,-1.28);
    variable[1] = r.Uniform(constants.vUD[0]-0.001,constants.vUD[0]+0.001);

    minimizer->SetLimitedVariable(0,"Ca",variable[0], step[0],-1.5,-1.);
   minimizer->SetLimitedVariable(1,"Vud",variable[1],step[1],constants.vUD[0]-0.1,constants.vUD[0]+0.1);
    //minimizer->SetFixedVariable(1,"Vud",constants.vUD[0]);

  }
  if(dof==4)
  {
    double step[4]; step[0] =0.0001; step[1] = 0.00001; step[2] = 0.0001; step[3] = 0.0001;  

    variable[0] = r.Uniform(-1.26,-1.28);
    variable[1] = r.Uniform(constants.vUD[0]-0.001,constants.vUD[0]+0.001);
    variable[2] = r.Uniform(-0.02,0.02);
    variable[3] = r.Uniform(-0.02,0.02);

    minimizer->SetLimitedVariable(0,"Ca",variable[0], step[0],-1.5,-1.);
    minimizer->SetLimitedVariable(1,"Vud",variable[1],step[1],constants.vUD[0]-0.1,constants.vUD[0]+0.1);
    //minimizer->SetFixedVariable(1,"Vud",constants.vUD[0]);
    minimizer->SetLimitedVariable(2,"Cs",variable[2],step[2],-1,1);
    minimizer->SetLimitedVariable(3,"Csp",variable[3],step[3],-1,1);
  }
  if(dof==-4)
  {
    double step[4]; step[0] =0.0001; step[1] = 0.00001; step[2] = 0.0001; step[3] = 0.0001;  

    variable[0] = r.Uniform(-1.26,-1.28);
    variable[1] = r.Uniform(constants.vUD[0]-0.001,constants.vUD[0]+0.001);
    variable[2] = r.Uniform(-0.02,0.02);
    variable[3] = r.Uniform(-0.02,0.02);

    minimizer->SetLimitedVariable(0,"Ca",variable[0], step[0],-1.5,-1.);
    minimizer->SetLimitedVariable(1,"Vud",variable[1],step[1],constants.vUD[0]-0.1,constants.vUD[0]+0.1);
    //minimizer->SetFixedVariable(1,"Vud",constants.vUD[0]);
    minimizer->SetLimitedVariable(2,"Ct",variable[2],step[2],-1,1);
    minimizer->SetLimitedVariable(3,"Ctp",variable[3],step[3],-1,1);
  }
}


double Analyzer::ChiSqrFunction2(const double *xx )
{
  double ca = xx[0];
  double vud = xx[1];
  double cs = par2; double csp = par2;
  double ct = par1; double ctp = par1;
  //cout << "call chisqr " << endl;
  double value = GetChiSqr(ca,cs,csp,ct,ctp,vud);

  return value;
}

double Analyzer::ChiSqrFunction2bis(const double *xx )
{
  double ca = xx[0];
  double vud = xx[1];
  double cs = par2; double csp = -par2;
  double ct = par1; double ctp = -par1;
  //cout << "call chisqr " << endl;
  double value = GetChiSqr(ca,cs,csp,ct,ctp,vud);

  return value;
}

double Analyzer::ChiSqrFunction4(const double *xx )
{
  double ca = xx[0];
  double vud = xx[1];
  double cs = xx[2];
  double csp = xx[3];
  double ct = (par1 + par2)/2.;
  double ctp = (par1 - par2)/2.;
  double value = GetChiSqr(ca,cs,csp,ct,ctp,vud);
  return value;
}

double Analyzer::ChiSqrFunction4bis(const double *xx )
{
  double ca = xx[0];
  double vud = xx[1];
  double ct = xx[2];
  double ctp = xx[3];
  double cs = (par1 + par2)/2.;
  double csp = (par1 - par2)/2.;
  double value = GetChiSqr(ca,cs,csp,ct,ctp,vud);
  return value;
}


double Analyzer::GetChiSqr(double ca, double cs, double csp, double ct, double ctp, double vud)
{
  double value = 0.;
  for(int i = 0; i<data.size() ; i++ )
  {
    if(!data.at(i)->Use()) continue;
    double expValue=data.at(i)->ExpValue();
    double expError=data.at(i)->Error();
    std::vector<double> expectation = data.at(i)->GetExpectation(cv,ca,cs,csp,ct,ctp,vud);    
    double totErrorSq = expError*expError + expectation.at(1)*expectation.at(1);
    value+= pow(expValue-expectation.at(0),2)/totErrorSq;     
  }  
  return value;
}


int Analyzer::Run(double low1, double low2, double high1, double high2, double step, int dof)
{
  
  if(data.size()<1) return 0;

  int totalSteps = (int)((high1-low1)/step * (high2-low2)/step );
  cout << "Running minimizer with " << totalSteps << " steps ... " << endl;
  
  par1 = low1; 
  par2 = low2;  
  int nstep = 0;

  while(par1<high1)
  {  
    while(par2<high2)
    {
      nstep++;
      InitVariables(dof); //resetting the variables is imperative, otherwise they tend to drift away  
      minimizer->Minimize(); 
      const double *xs = minimizer->X();
      double minChiSqr = (double)minimizer->MinValue();
      manager->SetOutput(minChiSqr,xs,par1,par2);
      /*if(nstep%100==0) 
      {
        cout << "step nr " << nstep <<  "  ca = " << xs[0] << endl;
        cout << "par1 " << par1 << "   par2 " << par2 << " xs[0] " << xs[0] << "  xs[1] " << xs[1] << " xs[2] " << xs[2] << "  xs[3] " << xs[3] << endl; 
      }*/
      
      //cout << par1 << par2 << endl;
      par2+=step;
      int modulus = (int)(totalSteps/10.);
      if(nstep%modulus==0) { cout << "step nr " << nstep <<  "  ca = " << xs[0] << endl; }
    }
    par1+=step; par2 = low2;
  }

  return 1;
}

int Analyzer::TestRun(double x, double y, int dof)
{
  double ca = constants.ca_fixed[0];
  double vud = constants.vUD[0];
  double cs, csp, ct, ctp;

  cout << " data size " << data.size() << endl;

  if(dof==2)
  {
    cs = y; csp = y;
    ct = x; ctp = x;
  }
  if(dof==-2)
  {
    cs = y; csp = -y;
    ct = x; ctp = -x;
  }
  if(dof ==4)
  {
    ct = (x + y)/2.;
    ctp = (x - y)/2.;
    cs = 0.; csp = 0.;
  }
  if(dof == -4)
  {
    cs = (x + y)/2.;
    csp = (x - y)/2.;
    ct = 0.; ctp = 0.;
  }

  std::cout.precision(4);

  if(data.size()<1) return 0;

  cout << endl;
  cout << "***************************************" << endl;
  cout << "************* Test Run **************" << endl;
  cout << "***************************************" << endl << endl;
  cout << "Ca = " << ca << "	Vud = " << vud << "	cp = " << ct <<"	ctp = " << ctp << "	cs = " << cs << "	cs = " << csp << endl << endl;   

  cout << "Isotope	Par.	ExpVal.	  Error	Expect.	Error	DChiSqr" << endl << endl;

  for(int i = 0; i<data.size() ; i++ )
  {
    if( !data.at(i)->Use() ) continue;
    std::vector<double> expectation = data.at(i)->GetExpectation(cv,ca,cs,csp,ct,ctp,vud);  
    double expError = data.at(i)->Error();
    double expValue = data.at(i)->ExpValue();
    double totErrorSq = expError*expError + expectation.at(1)*expectation.at(1);
    double dChiSqr = pow(expValue-expectation.at(0),2)/totErrorSq;
    cout << data.at(i)->Isotope() << "	" << data.at(i)->Parameter() << "	" << data.at(i)->ExpValue() << "	  " << data.at(i)->Error()  << "	" << expectation.at(0) << "	" << expectation.at(1) << "	" << dChiSqr << endl;
    
  }

  return 1;

}



