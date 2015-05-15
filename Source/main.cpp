#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

#include "../Include/common.h"
#include "../Include/structures.h"
#include "../Include/EM_cascade.h"
#include "EM_cascade.cpp"
#include "../Include/Injected_spectrum.h"
#include "Injected_spectrum.cpp"
#include "../Include/tools.h"
#include "tools.cpp"
#include "../Include/BBN_constraints.h"
#include "BBN_constraints.cpp"

int main(void){

  struct Structure_Gamma_Spectrum Cascade_Spectrum;
  struct Structure_Particle_Physics_Model Particle_Physics_Model;
  double M_x = 50;
  double tau_x = pow(10,6);
  double Zeta_x = pow(10,-6);
  double z = 1000000;
  string nuclei ("4He");
  ofstream Spectrum("Output/Universal_spectrum.dat");

  Fill_Structure_Particle_Physics_Model(M_x, Zeta_x, tau_x, &Particle_Physics_Model);
  double z_initial = 5*Particle_Physics_Model.z_x;
  cout << " z initial " << z_initial << endl;
  double z_final = z_initial/2.;
  // print_spectrum_from_function(Spectrum, Universal_Spectrum, z, &Particle_Physics_Model);
  // Cascade_Spectrum_Calculation_From_Function(Dirac_Spectrum_After_One_Iteration, z, &Particle_Physics_Model, &Cascade_Spectrum, 1000, 10);
  Constraints_from_destruction_only(nuclei,  z_initial, z_final, 2, &Cascade_Spectrum, &Particle_Physics_Model, 1000);
  //





  return 0;
}
