#include "Globals.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;


DECAY_CONSTANTS::DECAY_CONSTANTS(){

  //Includes GF, m_e, f = 1.6887 from Nilkinson NPA 377 (1982) 474 and a factor (1+RC) = 1.03886 from Abele, Prog. P&N Phys. 60 (2008) 1-81  
  neutronStrength[0] = 4908.5;
  neutronStrength[1] = 1.9;
  
  //Towner and Hardy, ArXiv 2013, check
  supperStrength[0]=2915.64;
  supperStrength[0]=1.08;

  //constants
  cv = 1.;// C'V = CV
  ca_fixed[0] = -1.2723;//set C'A = CA , PDG2014 value
  ca_fixed[1] = 0.0023;//set C'A = CA , PDG2014 value
  alpha = 0.0072973525698;
  vUD[0] = 0.97425; //I S Towner and J C Hardy 2010 Rep. Prog. Phys. 73 046301, Now floating!
  vUD[1] = 0.00022;
}
DECAY_CONSTANTS constants;
