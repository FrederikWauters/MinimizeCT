//Frederik Wauters, Wept 2012

//Take the output from MinimizeCT.cpp and draw confidence contours


#include "TRandom2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TDatime.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraphSmooth.h"

//gROOT->Reset();

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

//input histograms
TH2F *hchiSqr;
TH2F *hPchiSqr;
TH2F *hPDFchiSqr;

//CL contour plots
TGraph *g68CL = new TGraph();
TGraph *g90CL;
TGraph *g95CL;
TGraph *gSmooth;

TGraph *gChiSqrPlus1;
TGraph *gChiSqrPlus2;
TGraph *gChiSqrPlus2p30;
TGraph *gChiSqrPlus2p71;
TGraph *gChiSqrPlus3;
TGraph *gChiSqrPlus4;
TGraph *gChiSqrPlus4p61;
TGraph *gChiSqrPlus5p99;

//pibeta limits
TGraph *gPionLow;
TGraph *gPionHigh;

TGraph *gScale; //stupid graph to fix the scale on the plots


//C.L. scan parameters
const double step = 0.05; // stepping through the chi sqr values
const int nPoints = 200; // minchiSqr +  nPoints*step will be scanned 

//Process input*************************
int ReadInput()
{
  gStyle->SetPalette(1);
  char fname[1024];

  cout << "give file name:" << endl;
  cin >> fname;

  TFile *f1 = new TFile(fname,"READ");
  if( f1->IsZombie() ) 
  {
    cerr << "***ERROR! Cannot open file [" << fname << "]" << endl;
    delete f1;
    return 0;
  }

  hchiSqr = (TH2F*)f1->FindObjectAny("ChiSqr");
  hPchiSqr = (TH2F*)f1->FindObjectAny("PChiSqr");
  hPDFchiSqr = (TH2F*)f1->FindObjectAny("PDFChiSqr");

  TTree* tree = (TTree *)f1->Get("tExpInputdata");
  
  char isotope[10];
  //const char* parameter;
  char parameter[10];
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
  
  tree->SetBranchAddress("isotope", &isotope);
  tree->SetBranchAddress("parameter", &parameter);
  tree->SetBranchAddress("jInit", &jInit);
  tree->SetBranchAddress("jFinal", &jFinal);
  tree->SetBranchAddress("zDaughter", &zDaughter);
  tree->SetBranchAddress("sign", &sign);
  tree->SetBranchAddress("MF", &MF);
  tree->SetBranchAddress("MGT", &MGT);
  tree->SetBranchAddress("mE", &mE);
  tree->SetBranchAddress("expValue", &expValue);
  tree->SetBranchAddress("error", &error);
  tree->SetBranchAddress("use", &use);

  Int_t nentries = tree->GetEntries();

  cout << "*******************************" << endl;
  cout << "Input data : " << endl;

  for(Int_t i = 0; i < nentries; i++)
  {
    tree->GetEvent(i);
    cout <<  isotope << "  "<< parameter  << "  "<< jInit  << "  "<< jFinal << "  " << zDaughter << "  " << sign << "  "<< MF << " " <<MGT << "  "<< mE << "  "<< expValue<< "  " << error << " " << use << endl;
  }
  cout << "*******************************" << endl;

  //delete f1;
  return 1;



}

void InitCLContours()
{
  g68CL = new TGraph();
  g68CL->SetLineColor(2);
  g68CL->SetLineWidth(2);

  g90CL = new TGraph();
  g90CL->SetLineColor(4);
  g90CL->SetLineWidth(2);

  g95CL = new TGraph();
  g95CL->SetLineColor(6);
  g95CL->SetLineWidth(2);

  gChiSqrPlus1 = new TGraph();
  gChiSqrPlus1->SetLineColor(2);
  gChiSqrPlus1->SetLineWidth(2);

  gChiSqrPlus2 = new TGraph();
  gChiSqrPlus2->SetLineColor(4);
  gChiSqrPlus2->SetLineWidth(2);

  gChiSqrPlus2p30 = new TGraph();
  gChiSqrPlus2p30->SetLineColor(2);
  gChiSqrPlus2p30->SetLineWidth(2);

  gChiSqrPlus2p71 = new TGraph();
  gChiSqrPlus2p71->SetLineColor(4);
  gChiSqrPlus2p71->SetLineWidth(2);

  gChiSqrPlus3 = new TGraph();
  gChiSqrPlus3->SetLineColor(6);
  gChiSqrPlus3->SetLineWidth(2);

  gChiSqrPlus4 = new TGraph();
  gChiSqrPlus4->SetLineColor(6);
  gChiSqrPlus4->SetLineWidth(2);

  gChiSqrPlus4p61 = new TGraph();
  gChiSqrPlus4p61->SetLineColor(4);
  gChiSqrPlus4p61->SetLineWidth(2);

  gChiSqrPlus5p99 = new TGraph();
  gChiSqrPlus5p99->SetLineColor(6);
  gChiSqrPlus5p99->SetLineWidth(2);


  Double_t y[2] = {-0.3,0.3};
  Double_t x1[2] = {-0.00856/2.,-0.00856/2.};
  Double_t x2[2] = {0.0069/2.,0.0069/2.};
  gPionLow = new TGraph(2,x1,y);
  gPionLow->SetLineColor(3);
  gPionLow->SetLineWidth(3);
  gPionHigh = new TGraph(2,x2,y);
  gPionHigh->SetLineColor(3);
  gPionHigh->SetLineWidth(3);
  

  double gx[2] = {-0.27,0.27};
  double gy[2] = {-0.27,0.27};
  gScale = new TGraph(2,gx,gy);

}



void ConstructContour(TH2* h1, TH2* h2, TGraph* g)
{
  bool inCLContour;
  Int_t np;
  double x,y;
  double firstx, firsty;
  bool filledFirst = false;
  
  for(Int_t k = 1; k <= h1->GetNbinsY(); k ++)
  {
    inCLContour = false;

    for(Int_t j = 1; j <= h1->GetNbinsX(); j++)
    {
      if(h1->GetXaxis()->GetBinCenter(j) > -0.3)
      {
        if(h2->GetBinContent(j,k) > 0 && h1->GetBinContent(j,k) <= 0 && !inCLContour) //Its not good to == floating points
        {
          x = h1->GetXaxis()->GetBinCenter(j);
          y = h1->GetYaxis()->GetBinCenter(k);
          if(!filledFirst){ firstx = x; firsty = y; filledFirst = true; }
          //cout << " Bin Content : " << h1->GetBinContent(j,k) << " x " << x << " y " << y << endl;
          np = g->GetN(); 
          g->SetPoint(np,x,y);
          inCLContour = true;
        }
      }
    }
  }

  for(Int_t k = h1->GetNbinsY(); k >= 1; k--)
  {
    inCLContour = false;

    for(Int_t j = h1->GetNbinsX(); j >= 1; j--)
    {

      if(h2->GetBinContent(j,k) > 0 && h1->GetBinContent(j,k) <= 0 && !inCLContour) //Its not good to == floating points
      {
        x = h1->GetXaxis()->GetBinCenter(j);
        y = h1->GetYaxis()->GetBinCenter(k);
        //cout << " Bin Content : " << h1->GetBinContent(j,k) << " x " << x << " y " << y << endl;
        np = g->GetN();
        g->SetPoint(np,x,y);
        inCLContour = true;
      }
    }
  }
  //close contour
  np = g->GetN();
  g->SetPoint(np,firstx,firsty);       

}




void CalculateCLContours()
{
  double maxPValue = hPchiSqr->GetBinContent(hPchiSqr->GetMaximumBin());
  cout << "maximum p value = " << maxPValue << endl;

  double minchiSqrValue = hchiSqr->GetBinContent(hPchiSqr->GetMaximumBin());
  cout << "min chi sqr value = " << minchiSqrValue<< endl;

  TH2F* hClone = (TH2F*)hPDFchiSqr->Clone();
  TH2F* hCloneB = (TH2F*)hPDFchiSqr->Clone(); //Will set negative values to detect C.L. in "ConstructContour". Bit a dirty trick bit if it works ...
  double totalInt = hClone->Integral();
  cout << "total Int : " << totalInt <<endl;

  double newInt;
  double chiSqr;

  bool filled68 = false;
  bool filled90 = false;
  bool filled95 = false;

  bool filledCS1 = false;
  bool filledCS2 = false;
  bool filledCS3 = false;
  bool filledCS4 = false;
  bool filledCS2p30 = false;
  bool filledCS2p71 = false;
  bool filledCS4p61 = false;
  bool filledCS5p99 = false;

  cout << "Step value : " << step << " and nPoints : " << nPoints << endl;

  for(Int_t i = 1; i <= nPoints ;i++)
  {
    chiSqr = minchiSqrValue + (double)i*step;

    for(Int_t j = 1; j <= hClone->GetNbinsX(); j++)
    {
      for(Int_t k = 1; k <= hClone->GetNbinsY(); k ++)
      {
        if(hchiSqr->GetBinContent(j,k) < chiSqr){  hClone->SetBinContent(j,k,0); }
      }
    }

    newInt = hClone->Integral();

    //cout << "chiSqr : " << chiSqr << " newInt : " << newInt << "  C.L. : " << 100.-newInt/totalInt*100. << endl;

    if(100.-newInt/totalInt*100. > 68. && !filled68)
    {
      cout << "Construct 68% C.L. contour, delta chi sqr =  " << chiSqr - minchiSqrValue<< endl;
      ConstructContour(hClone,hCloneB,g68CL);
      filled68=true;
    }

    if(100.-newInt/totalInt*100. > 90. && !filled90)
    {
      cout << "Construct 90% C.L. contour, delta chi sqr =  " << chiSqr - minchiSqrValue<< endl;
      ConstructContour(hClone,hCloneB,g90CL);
      filled90=true;
    } 

    if(100.-newInt/totalInt*100. > 95. && !filled95)
    {
      cout << "Construct 95% C.L. contour, delta chi sqr =  " << chiSqr - minchiSqrValue<< endl;
      ConstructContour(hClone,hCloneB,g95CL);
      filled95=true;
    } 

    if(chiSqr > (minchiSqrValue + 1.) && !filledCS1)
    {
      cout << "Construct Chi Sqr +1  C.L. contour " << endl;
      ConstructContour(hClone,hCloneB,gChiSqrPlus1);
      filledCS1=true;
    } 

    if(chiSqr > (minchiSqrValue + 2.) && !filledCS2)
    {
      cout << "Construct Chi Sqr +2  C.L. contour " << endl;
      ConstructContour(hClone,hCloneB,gChiSqrPlus2);
      filledCS2=true;
    } 

    if(chiSqr > (minchiSqrValue + 3.) && !filledCS3)
    {
      cout << "Construct Chi Sqr +3  C.L. contour " << endl;
      ConstructContour(hClone,hCloneB,gChiSqrPlus3);
      filledCS3=true;
    } 

    if(chiSqr > (minchiSqrValue + 4.) && !filledCS4)
    {
      cout << "Construct Chi Sqr +4  C.L. contour " << endl;
      ConstructContour(hClone,hCloneB,gChiSqrPlus4);
      filledCS4=true;
    }

    if(chiSqr > (minchiSqrValue + 2.30) && !filledCS2p30)
    {
      cout << "Construct Chi Sqr + 2.30  C.L. contour " << endl;
      ConstructContour(hClone,hCloneB,gChiSqrPlus2p30);
      filledCS2p30=true;
    }

    if(chiSqr > (minchiSqrValue + 2.71) && !filledCS2p71)
    {
      cout << "Construct Chi Sqr + 2.71  C.L. contour " << endl;
      ConstructContour(hClone,hCloneB,gChiSqrPlus2p71);
      filledCS2p71=true;
    }

    if(chiSqr > (minchiSqrValue + 4.61) && !filledCS4p61)
    {
      cout << "Construct Chi Sqr + 4.61  C.L. contour " << endl;
      ConstructContour(hClone,hCloneB,gChiSqrPlus4p61);
      filledCS4p61=true;
    }

    if(chiSqr > (minchiSqrValue + 5.99) && !filledCS5p99)
    {
      cout << "Construct Chi Sqr + 5.99  C.L. contour " << endl;
      ConstructContour(hClone,hCloneB,gChiSqrPlus5p99);
      filledCS5p99=true;
    }


    //cout << "i  " <<  i << " CL: " << xValues[i-1] << " value: " << yValues[i-1] << endl;
  }

}


void DrawContours()
{
  //g68CL->Draw("AP");
  TCanvas *c1 = new TCanvas("c1","68%, 90% and 95% C.L. contours",800,800);
  c1->cd();
  gPad->SetGridx();
  gPad->SetGridy();

  TMultiGraph *mga = new TMultiGraph();
  mga->Add(gScale,"");
  mga->Add(g95CL,"l");
  mga->Add(g90CL,"l");
  mga->Add(g68CL,"l");
  //mga->Add(gPionLow,"l");
  //mga->Add(gPionHigh,"l");
  mga->Draw("A");
  mga->GetXaxis()->SetTitle("C_{T}/C_{A}");
  mga->GetYaxis()->SetTitle("C_{S}/C_{V}");
  mga->GetXaxis()->CenterTitle();
  mga->GetYaxis()->CenterTitle();
  mga->GetYaxis()->SetRangeUser(-0.1,0.1);
  mga->GetXaxis()->SetRangeUser(-0.1,0.1);
  gPad->Modified();

   //TGraphSmooth *gs = new TGraphSmooth("normal");
   //gSmooth = gs->Approx(g68CL,"linear", 50, 0, 0, 0, 1, 0.5, "tied");


  /*TCanvas *c2 = new TCanvas("c2","Chi Sqr +2.3,3.6,6.0 C.L. contours",800,800);
  c2->cd();
  gPad->SetGridx();
  gPad->SetGridy();

  TMultiGraph *mgb = new TMultiGraph();
  mgb->Add(gScale,"");
  mgb->Add(gChiSqrPlus2p30,"l");
  mgb->Add(gChiSqrPlus4p61,"l");
  mgb->Add(gChiSqrPlus5p99,"l");
  //mgb->Add(gPionLow,"l");
  //mgb->Add(gPionHigh,"l");
  mgb->Draw("A");
  mgb->GetXaxis()->SetTitle("C_{T}/C_{A}");
  mgb->GetYaxis()->SetTitle("C_{S}/C_{V}");
  mgb->GetXaxis()->CenterTitle();
  mgb->GetYaxis()->CenterTitle();
  mgb->GetYaxis()->SetRangeUser(-0.1,0.1);
  mgb->GetXaxis()->SetRangeUser(-0.1,0.1);
  gPad->Modified();*/

 /* TCanvas *c3 = new TCanvas("c3","Chi Sqr +1,2.71,4.00 C.L. contours",800,800);
  c3->cd();
  gPad->SetGridx();
  gPad->SetGridy();

  TMultiGraph *mgc = new TMultiGraph();
  mgc->Add(gScale,"");
  mgc->Add(gChiSqrPlus1,"l");
  mgc->Add(gChiSqrPlus2p71,"l");
  mgc->Add(gChiSqrPlus4,"l");
  mgc->Add(gPionLow,"l");
  mgc->Add(gPionHigh,"l");
  mgc->Draw("A");
  mgc->GetXaxis()->SetTitle("#frac{(C_{T}+C'_{T})}{C_{A}}");
  mgc->GetYaxis()->SetTitle("#frac{(C_{T}-C'_{T})}{C_{A}}");
  mgc->GetXaxis()->CenterTitle();
  mgc->GetYaxis()->CenterTitle();
  mgc->GetYaxis()->SetRangeUser(-0.27,0.27);
  mgc->GetXaxis()->SetRangeUser(-0.03,0.03);
  gPad->Modified();

  TCanvas *c4 = new TCanvas("c4","Chi Sqr +2.30,4.61,6.17 C.L. contours",800,800);
  c4->cd();
  gPad->SetGridx();
  gPad->SetGridy();

  TMultiGraph *mgd = new TMultiGraph();
  mgd->Add(gScale,"");
  mgd->Add(gChiSqrPlus2p30,"l");
  mgd->Add(gChiSqrPlus4p61,"l");
  mgd->Add(gChiSqrPlus6p17,"l");
  mgd->Add(gPionLow,"l");
  mgd->Add(gPionHigh,"l");
  mgd->Draw("A");
  mgd->GetXaxis()->SetTitle("#frac{(C_{T}+C'_{T})}{C_{A}}");
  mgd->GetYaxis()->SetTitle("#frac{(C_{T}-C'_{T})}{C_{A}}");
  mgd->GetXaxis()->CenterTitle();
  mgd->GetYaxis()->CenterTitle();
  mgd->GetYaxis()->SetRangeUser(-0.27,0.27);
  mgd->GetXaxis()->SetRangeUser(-0.03,0.03);
  gPad->Modified();*/
}

//*********************************************
//*******************main program**************
//*********************************************

int MakeCLCont()
{

  if(!ReadInput()) return 0;

  InitCLContours();
  CalculateCLContours();
  DrawContours();
  


  return 1;
}
