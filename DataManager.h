#ifndef DataManager_h
#define DataManager_h


#include <string>
#include <iostream>
#include <cstring>
#include <fstream>
#include <TGFrame.h>


#include "TFile.h"

#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMultiGraph.h"



class TData;

class DataManager {


  private:
  
    int ProcessInputFile(std::ifstream* f);
    int SetOutputFile(char* outfile);

    TFile *fout;
    int dof; //determines which mode we're running in

    TH2F* hCa;
    TH2F* hVud;
    TH2F* hChiSqr;
    TH2F* hPDF;
    
    TH2F* hC1;
    TH2F* hC2;
    
    TH1F* hPDF_X;
    TH1F* hPDF_Y;
    
    TH1F* hCL_X;
    TH1F* hCL_Y;
    
    TGraph* g68CL;
    TGraph* g90CL;
    TGraph* g95CL;
    TGraph* gScale;
    TMultiGraph* mga;

  public:
    
    DataManager(char* fin, char* fout, int d);
    virtual ~DataManager();

    std::vector<TData*> data;
    void PrintData();
    std::vector<TData*> GetData() {return data;}
    
    void SetOutput(double min, const double *xs,double par1, double par2);

    double PDFValue(double x,int n);
    void InitHistos(double low1, double low2, double high1, double high2, double step);
    void Plot();
    void WriteOutput();

    //post processing
    void ConstructContour(TH2* h1, TH2* h2, TGraph* g);
    void MakeCLContours(int nPoints, double dChiSqr);
    void Make1DContours(double CL, int nPoints);




    
};

#endif

