/**

This class is a data-type for experimental data + functionalities

**/

//Std includes
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cmath>


#include "TData.h"
#include "Globals.h"

using std::cout;
using std::endl;
using std::vector;

extern DECAY_CONSTANTS constants;
//CONSTANTS constants;

ClassImp(TData)

//constructor
TData::TData()
{
  Reset();

}

void TData::Reset()
{
  isotope ="";
  parameter = "";
  jInit = 0;
  jFinal = 0;
  zDaughter = 0;
  sign = 0;
  mF = 0;
  mGT = 0;
  mE = 0;
  expValue = 0;
  error = 0;
  use = 0;
}

//Destructor
TData::~TData()
{
}

void TData::SetData(std::iostream& ss)
{
  //cout << " Set data point " << endl;
  char iso[10];// = "Ca";
  char par[10];
  double jI = -1.0;
  double jF = -1.0;
  int zD = 0;
  std::string plmn;
  double mFermi;
  double mGamTel;
  double avME = 0.0;
  double expV = 0.0;
  double err = 0.0;
  int useIt = 0;
  
  ss >> iso >> par >> jI  >> jF  >> zD >> plmn >> mFermi >> mGamTel >> avME >> expV >>  err >> useIt;

  std::string plus("+");

  isotope = (char*)iso;
  parameter = (char*)par;
  jInit = jI;
  jFinal = jF;
  zDaughter = zD;
  if(plus.compare(plmn)==0) sign = -1; else sign = +1;
  mF = mFermi;
  mGT = mGamTel;
  mE = avME;
  expValue = expV;
  error = err;
  use = useIt;
}

void TData::Print()
{
  std::string betaType;
  if(sign==1) betaType = "-"; else betaType="+"; 
  cout <<  isotope << "	" << parameter  << "	"<< jInit  << "	"<< jFinal << "	" << zDaughter << "	" << betaType << "	"<< mF << "	" << mGT << "	"<< mE << "	"<< expValue<< "	" << error << "	" << use << endl;
}

//*************************************************
//******* Formulas for correlation coefficients ***
//*************************************************

double TData::Gamma()
{
  double alpha = constants.alpha;
  return pow(1-alpha*alpha*zDaughter*zDaughter,0.5);
}

double TData::Lambda()
{
  if( jFinal == jInit-1 )  return 1.; 
  if( jFinal == jInit)  return (1./(1.+(double)jInit));
  else return  (-1./(1.+(double)jInit));
}

int TData::Delta()
{
  if( jFinal == jInit) return 1;
  else return 0;
}

double TData::Xi(double cv, double ca,  double cs, double cps, double ct, double cpt)
{
  return pow(mF,2)*( cs*cs + cv*cv + cps*cps + cv*cv) + pow(mGT,2)*( ct*ct + ca*ca + cpt*cpt + ca*ca ); 
}

double TData::aXi(double cv, double ca, double cs, double cps, double ct, double cpt)
{
  return pow(mF,2)*( -cs*cs + cv*cv + -cps*cps + cv*cv ) + (pow(mGT,2)/3)*( ct*ct - ca*ca + cpt*cpt - ca*ca );
}

double TData::bXi(double cv, double ca, double cs, double cps, double ct, double cpt)
{
  return sign*2*Gamma()*(mF*mF*( cs*cv + cps*cv ) + mGT*mGT*( ct*ca + cpt*ca ));
}

double TData::AXi(double cv, double ca, double cs, double cps, double ct, double cpt)
{
  return mGT*mGT*Lambda()*( 2*sign*( ct*cpt - ca*ca) ) + Delta()*mF*mGT*sqrt((double)jInit/(1.+(double)jInit))*( 2*( cs*cpt + cps*ct - 2*ca*cv ));
  
}

double TData::b(double cv, double ca, double cs, double cps, double ct, double cpt)
{
  return bXi(cv,ca,cs,cps,ct,cpt)/Xi(cv,ca,cs,cps,ct,cpt);
}


double TData::a(double cv, double ca, double cs, double cps, double ct, double cpt)
{
  double value = aXi(cv,ca,cs,cps,ct,cpt)/Xi(cv,ca,cs,cps,ct,cpt);
  value = value / (1 + mE*b(cv,ca,cs, cps, ct,cpt));
  return value;
}


double TData::A(double cv, double ca, double cs, double cps, double ct, double cpt)
{
  double value = AXi(cv,ca,cs,cps,ct,cpt)/Xi(cv,ca,cs,cps,ct,cpt);
  value = value / (1 + mE*b(cv,ca,cs, cps, ct,cpt));
  return value;
}


vector<double> TData::tN_NF(double cv, double ca, double cs, double cps, double ct, double cpt, double vud) //neutron lifetime value, not fierz corrected yet
{
  vector<double> vOut; //contains value and error

  double strength = constants.neutronStrength[0];
  double errorstrength = constants.neutronStrength[1];

  double tInv = vud*vud*( mF*mF*( 1 + (cs*cs+cps*cps)/(2*cv*cv) ) +  mGT*mGT*( (ca*ca)/(cv*cv) + (ct*ct+cpt*cpt)/(2*cv*cv) )  );
  
  vOut.push_back(strength/tInv);
  vOut.push_back(errorstrength/tInv);

  return vOut;
}


vector<double> TData::tN(double cv, double ca, double cs, double cps, double ct, double cpt, double vud)
{
  vector<double> vOut;

  if(isotope.compare("n")!=0) { cout << " neutron lifetime requested for another nucleus! " << endl; vOut.clear(); return vOut;  }

  vector<double> noFierz = tN_NF(cv,ca,cs,cps,ct,cpt,vud);
  double value_NF = noFierz.at(0);
  double error_NF = noFierz.at(1);
  
  double value = value_NF / (1 + mE*b(cv,ca,cs, cps, ct,cpt) );
  double error = value*error_NF/value_NF;

  vOut.push_back(value);
  vOut.push_back(error);

  return vOut;
}

vector<double> TData::ftSuper_NF(double cv, double ca, double cs, double cps, double ct, double cpt, double vud)
{
  vector<double> vOut; //contains value and error

  double ftInv = vud*vud*mF*mF*( 1 + (cs*cs)/(cv*cv) );
  
  double strength = constants.supperStrength[0]; 
  double errorstrength = constants.supperStrength[1];

  vOut.push_back(strength/ftInv);
  vOut.push_back( errorstrength/ftInv);

  return vOut;
}

vector<double> TData::FTSupper(double cv, double ca, double cs, double cps, double ct, double cpt, double vud)
{
  vector<double> vOut;

  if(parameter.compare("Ft")!=0) { cout << " ft 0+->0+ requested for another nucleus! " << endl; vOut.clear(); return vOut; }

  vector<double> noFierz = ftSuper_NF(cv,ca,cs,cps,ct,cpt,vud); 
  double value_NF = noFierz.at(0);
  double error_NF = noFierz.at(1);
  
  double value = value_NF / (1 + mE*b(cv,ca,cs, cps, ct,cpt) );
  double error = value*error_NF/value_NF;

  vOut.push_back(value);
  vOut.push_back(error);

  return vOut;
}

double TData::RelPol(double cv, double ca, double cs, double cps, double ct, double cpt)
{
  if(parameter.compare("R")!=0) { cout << " Warning, called when not a relative polarization data point, will corrupt data" << endl; return 0.; }

  //do GT part
  mGT = 1; mF =0; sign = 1;
  double bGT = b(cv,ca,cs, cps, ct,cpt);
  //do MF part
  mGT = 0; mF =1; sign = -1;
  double bF = b(cv,ca,cs, cps, ct,cpt);

  return  mE*( bF +  bGT) + 1.;
}

std::vector<double> TData::GetExpectation(double cv, double ca, double cs, double cps, double ct, double cpt, double vud)
{
  double value=0.;
  double error=0.;

  if(parameter.compare("a")==0) value = a(cv,ca,cs,cps,ct,cpt);

  else if(parameter.compare("A")==0) value = A(cv,ca,cs,cps,ct,cpt);

  else if(parameter.compare("b")==0) value = b(cv,ca,cs,cps,ct,cpt);

  else if(parameter.compare("R")==0) value = RelPol(cv,ca,cs,cps,ct,cpt);

  else if(parameter.compare("tn")==0) {
    vector<double> values = tN(cv,ca,cs,cps,ct,cpt,vud);
    value = values.at(0);
    error = values.at(1);
  }

  else if(parameter.compare("Ft")==0) {
    vector<double> values = FTSupper(cv,ca,cs,cps,ct,cpt,vud);
    value = values.at(0); 
    error = values.at(1);
  }

  vector<double> vOut;
  vOut.push_back(value);
  vOut.push_back(error);

  return vOut;
  
}
