

//#include"iostream.h"
#include<TCanvas.h>
#include<TGraph.h>
#include<TGraphErrors.h>
#include<TF1.h>

int main()
{

  const int nPoints = 13;
  
  //isotopes            10C ,    14O ,	  22Mg ,   34Ar ,   26mAl,   34Cl,    38mK,    42Sc ,   46V ,    50Mn ,   54Co ,   62Ga ,   74Rb
  //gamma*mE/gamma  	    
  double x[nPoints] = { 0.61859, 0.43763, 0.30721, 0.21056, 0.29925, 0.23200, 0.21086, 0.19806, 0.18033, 0.16639, 0.15385, 0.13748, 0.12003 };
  //FT
  double y[nPoints] = { 3076.7   , 3071.4 , 3077.9 , 3065.6 , 3072.9 , 3070.7 , 3071.6 , 3072.4 , 3074.1 , 3071.2 , 3069.9 , 3071.5 , 3076 };
  //err FT
  double err[nPoints]={ 4.5    , 3.2    , 7.3    , 8.4    , 1      , 1.8    , 2      , 2.3    , 2      , 2.4    , 2.6    , 6.7    , 11   };
  
  //dummyxerr
  double errx[nPoints]={0.0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  TGraphErrors *graph = new TGraphErrors{nPoints, x, y,errx, err};
  
  TF1* f = new TF1("function","[0]/(1+[1]*x)");
  graph->Fit(f);
  
  TCanvas* c = new TCanvas("Fts","Fts",800,400);
  graph->GetXaxis()->SetRangeUser(0,1);
  graph->SetMarkerStyle(20);
  graph->Draw("AP");
  
  
  return 1;

}



