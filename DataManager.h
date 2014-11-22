#ifndef DataManager_h
#define DataManager_h


#include <string>
#include <iostream>
#include <cstring>
#include <fstream>

#include "TFile.h"



class TData;

class DataManager {


  private:
  
    int ProcessInputFile(std::ifstream* f);
    int SetOutputFile(char* outfile);

    TFile *fout;

  public:
    
    DataManager(char* fin, char* fout);
    virtual ~DataManager();

    std::vector<TData*> data;
    void PrintData();
    std::vector<TData*>* GetData() {return &data;};
    

    
};

#endif

