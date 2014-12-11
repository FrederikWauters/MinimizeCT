#ifndef Analyzer_h
#define Analyzer_h

#include <string>
#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>

#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "TH2.h"
#include "TGraph.h"

class DataManager;

class TData;

class Analyzer {

  public:

    Analyzer(int dof, DataManager* datamanager,const char* minName, const char* algoName );
    virtual ~Analyzer();

    int InitMinimizer(int dof, const char * minName, const char *algoName);
    int Run(double low1, double low2, double high1, double high2, double step);
    int TestRun(double x, double y, int dof); //Get predicted values

    double ChiSqrFunction2(const double *xx);
    double ChiSqrFunction4(const double *xx);
    double GetChiSqr(double ca, double cs, double csp, double ct, double ctp, double vud);



  private:
 

    ROOT::Math::Minimizer* minimizer;
    ROOT::Math::Functor functor;
    std::vector<TData*> data;
    DataManager* manager;

    double par1, par2; //The coupling constants to loop over. Either cs and csp, or (ct+ctp)/ca and (ct-ctp)/ca 
    double cv, ca_SM;

    


};
#endif
