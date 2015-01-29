#ifndef TData_h
#define TData_h

#include <string>
#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include "TObject.h"


class TData : public TObject {

  private:
  
	std::string isotope;
	std::string parameter;
	double jInit;
	double jFinal;
	int zDaughter;
	int sign; //is + for beta- and - for beta+!
	double mF;
	double mGT;
	double mE;
	double expValue;
	double error;
	bool use;

	//outsiders have nu business with these intermediate steps
	double Gamma();  //Formulas from Rev. Mod. Phys. 78, 991, appendix B
	double Lambda();
	int Delta(); //checks if i-f spins are equal
	double Xi(double cv, double ca,  double cs, double cps, double ct, double cpt);
	double aXi(double cv, double ca, double cs, double cps, double ct, double cpt);
	double bXi(double cv, double ca, double cs, double cps, double ct, double cpt);
	double AXi(double cv, double ca, double cs, double cps, double ct, double cpt);
	std::vector<double> tN_NF(double cv, double ca, double cs, double cps, double ct, double cpt, double vud); //neutron lifetime value, not fierz corrected yet
	std::vector<double> ftSuper_NF(double cv, double ca, double cs, double cps, double ct, double cpt, double vud);


  public:
    
        TData();
	void Reset();
	virtual ~TData();

	std::string Isotope() { return isotope; }
	std::string Parameter() { return parameter; }
	double JInit() { return jInit; }
	double JFinal() { return jFinal; }
	int ZDaughter() { return zDaughter; }
	int Sign() { return sign; }
	double MF() { return mF; }
	double MGT() { return mGT; }
	double ME() { return mE; }
	double ExpValue() { return expValue; }
	double Error() { return error; }
	bool Use() { return use; }

	void SetData(std::iostream& ss);
	void Print();

	double a(double cv, double ca, double cs, double cps, double ct, double cpt);
	double b(double cv, double ca, double cs, double cps, double ct, double cpt);
	double A(double cv, double ca, double cs, double cps, double ct, double cpt);
	std::vector<double> tN(double cv, double ca, double cs, double cps, double ct, double cpt, double vud); //first entry is the value, second is it's error
	std::vector<double> FTSupper(double cv, double ca, double cs, double cps, double ct, double cpt, double vud);
	double RelPol(double cv, double ca, double cs, double cps, double ct, double cpt);

	std::vector<double> GetExpectation(double cv, double ca, double cs, double cps, double ct, double cpt, double vud); //probe for the expectation value for a specific set of couplings constants. Firts entry is the value, second is the error
	


 ClassDef(TData,3);
    
};

#endif

