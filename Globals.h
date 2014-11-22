#ifndef Globals_h
#define Globals_h 1


class DECAY_CONSTANTS {

  public:
  DECAY_CONSTANTS();
  double neutronStrength[2]; //value and error
  double supperStrength[2];

  //constants
  double cv;
  double ca_fixed[2];
  double alpha;
  double vUD[2]; 


};
#endif
