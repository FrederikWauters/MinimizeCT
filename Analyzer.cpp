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

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TMath.h"
#include "TRandom2.h"

using std::cout;
using std::endl;
using std::vector;

//DECAY_CONSTANTS constants;
extern DECAY_CONSTANTS constants;

Analyzer::Analyzer( int dof , DataManager* datamanager,const char* minName = "Minuit2", const char* algoName = "" )
{
  cout << "Construct analyzer ... " << endl;

  manager = datamanager;
  InitMinimizer(dof,minName,algoName);
}

Analyzer::~Analyzer()
{
  cout << "delete analyzer ... " << endl;
}

int Analyzer::InitMinimizer( int dof,const char * minName = "Minuit2", const char *algoName = "" )
{
  cout << "Initialize Minimizer" << endl;
  minimizer = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

  // set tolerance , etc...
  minimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  minimizer->SetMaxIterations(10000);  // for GSL 
  minimizer->SetTolerance(0.0001);
  minimizer->SetPrintLevel(-1);

  if(dof==2)
  {
    
    ROOT::Math::Functor f(this,&Analyzer::ChiSqrFunction2,dof);
    minimizer->SetFunction(f); 
    double step[dof]; step[0] = 0.0001; step[1] = 0.00001;

    double variable[dof];// = ca and Vud
    TRandom2 r(0);
    variable[0] = r.Uniform(-1.26,-1.28);
    variable[1] = r.Uniform(0.95,1.05);
 
    // Set the free variables to be minimized!
    minimizer->SetLimitedVariable(0,"Ca",variable[0], step[0],-1.5,-1.);
    minimizer->SetLimitedVariable(1,"Vud",variable[1],step[1],0.7,1.3);
   
  }

  if(dof==4)
  {
    ROOT::Math::Functor f(this,&Analyzer::ChiSqrFunction4,dof);
    minimizer->SetFunction(f); 
    double step[dof]; step[0] =0.0001; step[1] = 0.00001; step[2] = 0.0001; step[3] = 0.0001;

    double variable[dof];// = ca, cs and csp and Vud
    TRandom2 r(0);
    variable[0] = r.Uniform(-1.26,-1.28);
    variable[1] = r.Uniform(0.95,1.05);
    variable[2] = r.Uniform(-0.05,0.05);
    variable[3] = r.Uniform(-0.05,0.05);
 
    // Set the free variables to be minimized!
    minimizer->SetLimitedVariable(0,"Ca",variable[0], step[0],-1.5,-1.);
    minimizer->SetLimitedVariable(1,"Vud",variable[1],step[1],0.7,1.3);
    minimizer->SetLimitedVariable(2,"Cs",variable[2],step[2],-1,1);
    minimizer->SetLimitedVariable(3,"Csp",variable[3],step[3],-1,1);
  }


  return 1;
}

double Analyzer::ChiSqrFunction2(const double *xx )
{
  double value;
  return value;
}

double Analyzer::ChiSqrFunction4(const double *xx )
{
  double value;
  return value;
}

int Analyzer::Run()
{
  data = manager->GetData();
  return 1;
}
