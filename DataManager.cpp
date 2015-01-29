/**

Process input

**/

//Std includes
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cmath>


#include "TData.h"
#include "DataManager.h"
#include "TCanvas.h"
#include "Globals.h"
#include "TApplication.h"
#include <TGClient.h>
#include <TF1.h>
#include <TRandom.h>
#include "TMath.h"
#include "TTree.h"

//#include <TGButton.h>
//#inc/lude <TGFrame.h>
//#include <TRootEmbeddedCanvas.h>
//#include <RQ_OBJECT.h>
//#include "TGClient.h"

using std::cout;
using std::endl;

extern DECAY_CONSTANTS constants;


//constructor
DataManager::DataManager(char* fin, char* fout, int d)
{
  cout << "Construct data manager ... " << endl;

  data.clear();
  
  std::ifstream *inputfile = new std::ifstream(fin);
  if(inputfile->is_open()) cout << fin << " open ... " << endl;
    else cout << "failed to open inputfile ... " << endl; 
  if(!ProcessInputFile(inputfile)) cout << "input not properly processed" << endl; 
  inputfile->close();
 
  dof = d;

  if(SetOutputFile(fout)==1) cout << fout << " open ... " << endl;
  
  //InitHistograms();

}

//Destructor
DataManager::~DataManager()
{
  cout << "Delete data manager ... " << endl;
}

int DataManager::ProcessInputFile(std::ifstream* file)
{
  if(file==NULL) return 0;
  std::string line;
  int nInput = 0;
  
  while(std::getline(*file,line))
  {
    if(line[0]=='#' || line.size()<2) continue;
    
    else  
    {
      std::stringstream ss;
      ss << line;
      TData* newData = new TData();
      newData->SetData(ss);
      if(newData->Use()) data.push_back(newData);
  
      nInput++;
    }
  }

  cout<< "size " << data.size() << endl;
  return 1;
}

int DataManager::SetOutputFile(char* outfile)
{
  fout = new TFile ( outfile, "ReCreate");

  if ( fout->IsZombie() ) {
    std::cerr << "**ERROR! Cannot open file [" << outfile << "]" << endl;
    return 0;
  }
  fout->cd();
  return 1;
}

void DataManager::PrintData()
{
  cout << endl;
  cout << "***************************************" << endl;
  cout << "************* Input Data **************" << endl;
  cout << "***************************************" << endl << endl;
  if(data.size()<1) return;
  cout << "isotope" << "	" << "par" << "	" << "jInit"  << "	"<< "jFinal" << "	" << "zDaught" << "	" << "b +/-" << "	"<< "mF" << "	" << "mGT" << "	"<< "mE" << "	"<< "expVal" << "	" << "error" << "	" << "use" << endl<< endl;
  for(int i = 0; i<data.size(); i++) data.at(i)->Print();

  cout << endl;
   
}

double DataManager::PDFValue(double x,int n)
{
  return ( TMath::Exp(-x/2.) * pow(x,(n/2.-1)) ) / ( TMath::Gamma(n/2.) * pow(2,n/2.) );
}

void DataManager::SetOutput(double min, const double *xs,double par1, double par2)
{
  Int_t bin;
  if (dof==2 || dof == -2) bin = hCa->FindBin(par1/constants.ca_fixed[0],par2);
  else if (dof==4) bin = hCa->FindBin(par1/constants.ca_fixed[0],par2/constants.ca_fixed[0]);
  else if (dof==-4) bin = hCa->FindBin(par1/constants.cv,par2/constants.cv);
 
  double ca = xs[0];
  double vud = xs[1];

  int nDOF = data.size() - 2;

  hCa->SetBinContent(bin,ca);
  hVud->SetBinContent(bin,vud);

  //hChiSqr->SetBinContent(bin,min/nDOF);
  hChiSqr->SetBinContent(bin,min);
  hPDF->SetBinContent(bin,PDFValue(min,2));
  
  if(dof == 4) { hC1->SetBinContent(bin,xs[2]+xs[3]); hC2->SetBinContent(bin,xs[2]-xs[3]); }
  if(dof == -4) { hC1->SetBinContent(bin,(xs[2]+xs[3])/constants.ca_fixed[0]); hC2->SetBinContent(bin,(xs[2]-xs[3])/constants.ca_fixed[0]); }
}

void DataManager::InitHistos(double low1, double low2, double high1, double high2, double step)
{
  //Histograms
  int xBins = (int)((high1-low1)/step)+1;
  int yBins = (int)((high2-low2)/step)+1; cout << xBins << " " << yBins << " " <<step << endl;

  if(dof==2 || dof==-2)
  { 
    double ctca_min = high1/constants.ca_fixed[0];
    double ctca_max = low1/constants.ca_fixed[0];
    double cscv_min = low2/constants.cv;
    double cscv_max = high2/constants.cv;
    hCa = new TH2F("hCa","Minimum Ca; #frac{C_{T}}{C_{A}};#frac{C_{S}}{C_{V}}",xBins,ctca_min,ctca_max,yBins,cscv_min,cscv_max);
    hVud = new TH2F("hVud","Minimum Vud;#frac{C_{T}}{C_{A}};#frac{C_{S}}{C_{V}}",xBins,ctca_min,ctca_max,yBins,cscv_min,cscv_max);
    hChiSqr = new TH2F("hChiSqr","ChiSqr;#frac{C_{T}}{C_{A}};#frac{C_{S}}{C_{V}}",xBins,ctca_min,ctca_max,yBins,cscv_min,cscv_max);
    hPDF = new TH2F("hPDF","PDF surface;#frac{C_{T}}{C_{A}};#frac{C_{S}}{C_{V}}",xBins,ctca_min,ctca_max,yBins,cscv_min,cscv_max);
  }
  
  if(dof==4)
  { 
    double ct_plus_ctpca_min = high1/constants.ca_fixed[0];
    double ct_plus_ctpca_max = low1/constants.ca_fixed[0];
    double ct_minus_ctpca_min = high2/constants.ca_fixed[0];
    double ct_minus_ctpca_max = low2/constants.ca_fixed[0];
    hCa = new TH2F("hCa","Minimum Ca; #frac{C_{T}+C'_{T}}{C_{A}};#frac{C_{T}-C'_{T}}{C_{A}}",xBins,ct_plus_ctpca_min,ct_plus_ctpca_max,yBins,ct_minus_ctpca_min,ct_minus_ctpca_max);
    hVud = new TH2F("hVud","Minimum Vud;#frac{C_{T}+C'_{T}}{C_{A}};#frac{C_{T}-C'_{T}}{C_{A}}",xBins,ct_plus_ctpca_min,ct_plus_ctpca_max,yBins,ct_minus_ctpca_min,ct_minus_ctpca_max);
    hChiSqr = new TH2F("hChiSqr","ChiSqr;#frac{C_{T}+C'_{T}}{C_{A}};#frac{C_{T}-C'_{T}}{C_{A}}",xBins,ct_plus_ctpca_min,ct_plus_ctpca_max,yBins,ct_minus_ctpca_min,ct_minus_ctpca_max);
    hPDF = new TH2F("hPDF","PDF surface;#frac{C_{T}+C'_{T}}{C_{A}};#frac{C_{T}-C'_{T}}{C_{A}}",xBins,ct_plus_ctpca_min,ct_plus_ctpca_max,yBins,ct_minus_ctpca_min,ct_minus_ctpca_max);
    hC1 = new TH2F("hCS_plus","Minimum C_S + C'_S; #frac{C_{T}+C'_{T}}{C_{A}};#frac{C_{T}-C'_{T}}{C_{A}}",xBins,ct_plus_ctpca_min,ct_plus_ctpca_max,yBins,ct_minus_ctpca_min,ct_minus_ctpca_max);
    hC2 = new TH2F("hCS_minus","Minimum C_S - C'_{S}; #frac{C_{T}+C'_{T}}{C_{A}};#frac{C_{T}-C'_{T}}{C_{A}}",xBins,ct_plus_ctpca_min,ct_plus_ctpca_max,yBins,ct_minus_ctpca_min,ct_minus_ctpca_max);
  }
  
  if(dof==-4)
  {
    double cs_plus_cspcv_min = low1/constants.cv;
    double cs_plus_cspcv_max = high1/constants.cv;
    double cs_minus_cspcv_min = low2/constants.cv;
    double cs_minus_cspcv_max = high2/constants.cv;
    hCa = new TH2F("hCa","Minimum Ca; #frac{C_{S}+C'_{S}}{C_{S}};#frac{C_{S}-C'_{S}}{C_{V}}",xBins,cs_plus_cspcv_min,cs_plus_cspcv_max,yBins,cs_minus_cspcv_min,cs_minus_cspcv_max);
    hVud = new TH2F("hVud","Minimum Vud;#frac{C_{S}+C'_{S}}{C_{V}};#frac{C_{S}-C'_{S}}{C_{V}}",xBins,cs_plus_cspcv_min,cs_plus_cspcv_max,yBins,cs_minus_cspcv_min,cs_minus_cspcv_max);
    hChiSqr = new TH2F("hChiSqr","ChiSqr;#frac{C_{S}+C'_{S}}{C_{V}};#frac{C_{S}-C'_{S}}{C_{V}}",xBins,cs_plus_cspcv_min,cs_plus_cspcv_max,yBins,cs_minus_cspcv_min,cs_minus_cspcv_max);
    hPDF = new TH2F("hPDF","PDF surface;#frac{C_{S}+C'_{S}}{C_{V}};#frac{C_{S}-C'_{S}}{C_{V}}",xBins,cs_plus_cspcv_min,cs_plus_cspcv_max,yBins,cs_minus_cspcv_min,cs_minus_cspcv_max);
    hC1 = new TH2F("hCT_plus","Minimum C_S + C'_S; #frac{C_{S}+C'_{S}}{C_{V}};#frac{C_{S}-C'_{S}}{C_{V}}",xBins,cs_plus_cspcv_min,cs_plus_cspcv_max,yBins,cs_minus_cspcv_min,cs_minus_cspcv_max);
    hC2 = new TH2F("hCT_minus","Minimum C_S - C'_{S}; #frac{C_{S}+C'_{S}}{C_{V}};#frac{C_{S}-C'_{S}}{C_{V}}",xBins,cs_plus_cspcv_min,cs_plus_cspcv_max,yBins,cs_minus_cspcv_min,cs_minus_cspcv_max);
  }
  
   //Graphs (confidence contours)
   
  g68CL = new TGraph();
  g68CL->SetLineColor(13);
  g68CL->SetLineWidth(2);
  g68CL->SetEditable(kFALSE);
  g68CL->SetName("g68CL");

  g90CL = new TGraph();
  g90CL->SetLineColor(1);
  g90CL->SetLineWidth(2);
  g90CL->SetEditable(kFALSE);
  g90CL->SetName("g90CL");

  g95CL = new TGraph();
  g95CL->SetLineColor(27);
  g95CL->SetLineWidth(2);
  g95CL->SetLineStyle(7);
  g95CL->SetEditable(kFALSE);
  g95CL->SetName("g95CL");
  
  double gx[2] = {-0.27,0.27};
  double gy[2] = {-0.27,0.27};
  gScale = new TGraph(2,gx,gy);
  gScale->SetLineWidth(0);
  
  hPDF_X = new TH1F();
  hPDF_Y= new TH1F();
    
  hCL_X= new TH1F();
  hCL_Y= new TH1F();

  
}

void DataManager::Plot()
{
    //pop the ghui
    //MyMainFrame* myMainFrame = new MyMainFrame(gClient->GetRoot(),200,200); 
    //const TGWindow* main = gClient->GetRoot();
   // *fMain = new TGMainFrame(main ,200,200);
    //TRootEmbeddedCanvas *fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fMain,200,200);
    //fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
   
   //TCanvas *canvas = new TCanvas("MyCanvas","Test Canvas",10,10,800,800);
   //TGHorizontalFrame *hFrame = new TGHorizontalFrame(this);
   //hFrame->AddFrame(canvas, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX,10, 10, 10, 2));
   //TGTextButton *exitButton = new TGTextButton(hFrame,"Exit");
    
    TCanvas *canvas = new TCanvas("MyCanvas","Test Canvas",10,10,800,800);
    TApplication app("a", 0, 0);
    
    canvas->cd();
    hCa->Draw();
    app.Run(0);
    app.Terminate();
    cout << "Press ENTER to continue ... " << endl;
    std::cin.get();
    
}

void DataManager::WriteOutput()
{
  TCanvas *c1 = new TCanvas("c1","68%, 90% and 95% C.L. contours",800,800);
  c1->cd();
  gPad->SetGridx();
  gPad->SetGridy();
  gScale->SetLineColor(0);	
  mga = new TMultiGraph("mga","2D confidence contours");
  mga->Add(gScale,"");
  mga->Add(g95CL,"l");
  mga->Add(g90CL,"l");
  mga->Add(g68CL,"l");
  mga->Draw("A");
  
  if(dof==2 )
  {
    mga->GetXaxis()->SetTitle("#frac{C_{T}}{C_{A}}");
    mga->GetYaxis()->SetTitle("#frac{C_{S}}{C_{V}}");
    mga->GetYaxis()->SetRangeUser(-0.1,0.1);
    mga->GetXaxis()->SetRangeUser(-0.1,0.1);
  }
  if(dof==-2 )
  {
    mga->GetXaxis()->SetTitle("#frac{C_{T}}{C_{A}}");
    mga->GetYaxis()->SetTitle("#frac{C_{S}}{C_{V}}");
    mga->GetYaxis()->SetRangeUser(-0.01,0.01);
    mga->GetXaxis()->SetRangeUser(-0.01,0.01);
  }
  else if(dof==4)
  {
    mga->GetXaxis()->SetTitle("#frac{C_{T}+C'_{T}}{C_{A}}");
    mga->GetYaxis()->SetTitle("#frac{C_{T}-C'_{T}}{C_{A}}");
    mga->GetYaxis()->SetRangeUser(-0.3,0.3);
    mga->GetXaxis()->SetRangeUser(-0.1,0.1);
  }
  else if(dof==-4)
  {
    mga->GetXaxis()->SetTitle("#frac{C_{S}+C'_{S}}{C_{V}}");
    mga->GetYaxis()->SetTitle("#frac{C_{S}-C'_{S}}{C_{V}}");
    mga->GetYaxis()->SetRangeUser(-0.3,0.3);
    mga->GetXaxis()->SetRangeUser(-0.1,0.1);
  }
  mga->GetXaxis()->CenterTitle();
  mga->GetYaxis()->CenterTitle();
  
  gPad->Update();
  
  g68CL->Write();
  g90CL->Write();
  g95CL->Write();
  c1->Write();
  
  TData* datapoint = 0;
  TTree* tree = new TTree("InputData","Input data to this minimizer run");
  tree->Branch("datapoint","TData",&datapoint,0,64000);
  
  if(data.size()>0)
  {
    for(int i = 0; i < data.size(); i++)
    {
      datapoint = data.at(i);
      tree->Fill();
    }
  }
  data.at(0)->Write();
  //tree->Write();
  
  fout->Write();
  delete tree;
}

void DataManager::ConstructContour(TH2* h1, TH2* h2, TGraph* g)
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

void DataManager::MakeCLContours(int nPoints, double deltaChiSqr)
{

  int maxbin = hPDF->GetMaximumBin();

  double minchiSqrValue = hChiSqr->GetBinContent(maxbin);
  int nDOF = data.size() - 2;
  cout << endl << "min chi sqr value / DOF = " << minchiSqrValue/nDOF<< " DOF = " << nDOF << endl;
  cout << "Min Ca value = " << hCa->GetBinContent(maxbin) << endl;
  cout << "Min Vud value = " << hVud->GetBinContent(maxbin) << endl << endl;

  TH2F* hClone = (TH2F*)hPDF->Clone(); //
  TH2F* hCloneB = (TH2F*)hPDF->Clone(); // Will set negative values to detect C.L. in "ConstructContour". Bit a dirty trick bit if it works ...
  
  double totalInt = hClone->Integral();
  //cout << "total Int : " << totalInt <<endl;

  double newInt;
  double chiSqr;

  bool filled68 = false;
  bool filled90 = false;
  bool filled95 = false;

  double step = deltaChiSqr/(nPoints*1.);
  cout << "Step value : " << step << " and nPoints : " << nPoints << endl;

  for(Int_t i = 1; i <= nPoints ;i++)
  {
    chiSqr = minchiSqrValue + (double)i*step;

    for(Int_t j = 1; j <= hClone->GetNbinsX(); j++)
    {
      for(Int_t k = 1; k <= hClone->GetNbinsY(); k ++)
      {
        if(hChiSqr->GetBinContent(j,k) < chiSqr){  hClone->SetBinContent(j,k,0); }
      }
    }

    newInt = hClone->Integral();

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
  }
  delete hClone;
  delete hCloneB;
}

void DataManager::Make1DContours(double CL, int nPoints)
{
  hPDF_X = (TH1F*)hPDF->ProjectionX("hPDF_XProjection");
  hPDF_Y = (TH1F*)hPDF->ProjectionY("hPDF_YProjection");
  
  double totalInt1 = hPDF_X->Integral();
  double totalInt2 = hPDF_Y->Integral();


  //main loop   
  hCL_X = (TH1F*)hPDF_X->Clone();
  hCL_Y = (TH1F*)hPDF_Y->Clone();
  hCL_X->SetName("hCL_X");
  hCL_Y->SetName("hCL_Y");

  double maximum = hPDF_X->GetBinContent(hPDF_X->GetMaximumBin());
  double min = 0.0;
  double ratio, newInt;  

  for(Int_t i = 0; i < nPoints ;i++)
  {  
    min = 0 + i * maximum/(double)nPoints; 
    for(Int_t j = 1; j <= hCL_X->GetNbinsX(); j++)
    {
      if(hPDF_X->GetBinContent(j) < min){  hCL_X->SetBinContent(j,0); }
    }
    newInt = hCL_X->Integral();
    ratio = newInt/totalInt1;    
    if (ratio < CL) break;
  }
  
  maximum = hPDF_Y->GetBinContent(hPDF_Y->GetMaximumBin());
  min = 0.0;
  
  for(Int_t i = 0; i < nPoints ;i++)
  {  
    min = 0 + i * maximum/(double)nPoints; 
    for(Int_t j = 1; j <= hCL_Y->GetNbinsX(); j++)
    {
      if(hPDF_Y->GetBinContent(j) < min){  hCL_Y->SetBinContent(j,0); }
    }
    newInt = hCL_Y->Integral();
    ratio = newInt/totalInt2;    
    if (ratio < CL) break;
  }
  
  //Print CL
  
  double low = 0.; double high = 0.;
  
  int nBins_X = hCL_X->GetNbinsX();
  int nBins_Y = hCL_Y->GetNbinsX();
  
  for(int i = 1; i < nBins_X; i++)
  {
    if(hCL_X->GetBinContent(i) > 0) 
    {
      low = hCL_X->GetBinCenter(i);
      break;
    }
    else continue;    
  }
  
  for(int i = nBins_X; i > 0; i--)
  {
    if(hCL_X->GetBinContent(i) > 0) 
    {
      high = hCL_X->GetBinCenter(i);
      break;
    }
    else continue;
  }
  
  cout << CL << " confidence level (X) is [ " << low << " ; " << high << " ] " << endl;
  
  for(int i = 1; i < nBins_Y; i++)
  {
    if(hCL_Y->GetBinContent(i) > 0) 
    {
      low = hCL_Y->GetBinCenter(i);
      break;
    }
    else continue;    
  }
  
  for(int i = nBins_Y; i > 0; i--)
  {
    if(hCL_Y->GetBinContent(i) > 0) 
    {
      high = hCL_Y->GetBinCenter(i);
      break;
    }
    else continue;
  }
  
  
  cout << CL << " confidence level (Y) is [ " << low << " ; " << high << " ] " << endl;
  
}


