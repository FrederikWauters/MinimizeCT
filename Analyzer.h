#ifndef Analyzer_h
#define Analyzer_h

#include <string>
#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>

#include "Math/Minimizer.h"

class DataManager;

class TData;

class Analyzer {

  public:

    Analyzer(int dof , DataManager* datamanager,const char* minName, const char* algoName );
    virtual ~Analyzer();

    int InitMinimizer( int dof,const char * minName, const char *algoName);
    int Run();

    double ChiSqrFunction2(const double *xx);
    double ChiSqrFunction4(const double *xx);

  private:
 

    ROOT::Math::Minimizer* minimizer;
    std::vector<TData*>* data;
    DataManager* manager;



};
#endif
