#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>

#include "common.h"


struct Structure_Gamma_Spectrum{

  int iterations;
  int redshift;
  vector<double> Gamma_Energy;
  vector<double> Gamma_Spectrum;



};

struct Structure_Electron_Spectrum{

  int iterations;
  int redshift;
  vector<double> Electron_Energy;
  vector<double> Electron_Spectrum;



};


struct Structure_Particle_Physics_Model{

  double  M_x;
  double  E_0;
  double  Zeta_x;
  double  tau_x;
  double  T_x;
  double  z_x;

};
