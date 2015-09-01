#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

#include "../Include/common.h"
#include "../Include/structures.h"
#include "../Include/EM_cascade.h"
#include "../Include/injected_spectrum.h"
#include "../Include/tools.h"
#include "../Include/BBN_constraints.h"



int main(int argc, char** argv){
  if(argc>2){
    cout <<"Error : you have given too many files as an input.\n Please, restart cBBNFast with at most one '.ini' file." << endl;
    return 0;
  }

  // maps<string,unsigned> map_parameters;
  //
  // fill_default_parameters(map_parameters);

  /******* Precision parameters : control the number of step in integration scheme, ***************
  ********** as well as the number of iterations for computing the cascade spectrum **************/
  int number_iterations_photon = 5;      //Number of iterations for computing the photon spectrum
  int number_iterations_electron = 30;    //Number of iterations for computing the electron spectrum
  int z_step = 80;         //Number of redshift steps between injection time and the minimal redshift of integration
  int n_step = 200;        //Number of steps in the simpson algorithm for integrations, to be chosen small, it is adapted inside the code when needed.
  /***********************************************************************************************/

  /******* Some string useful in the program whatever calcutation is done*******/
  string nuclei ("4He");                       //This needs to be a nuclei in the list : "2H", "4He", "3He", "7Be", "7Li".
  string print_result ("yes");
  string calculation_mode ("triangular");      //This can be either "triangular" or "iterative".
  string photon_spectrum_choice ("Dirac");     //This can be either "Dirac", "universal" or user specified "from_file", "from_function".
  string electron_spectrum_choice ("none");    //This can be either "Dirac", "universal" or user specified "from_file", "from_function".
  string spectrum_mode ("writing");            //This can be either "writing", "reading" or "nothing".
  string inverse_compton_scattering ("no");    //In case you choose the iterative mode, it is possible to switch on and off the inverse_compton_scattering.
  string results_files ("automatic");
  string spectrum_files ("automatic");
  /*****************************************************************************/
  int task = 2; // 1 : Print spectrum in file, 2 : Compute_Constraints_from_destruction_only, 3 : Compute_constraints_from_destruction_and_production

  // /************To print cascade spectrum in a file*********/
  if(task==1){
  struct Structure_Spectrum Cascade_Spectrum;
  struct Structure_Particle_Physics_Model Particle_Physics_Model;
  struct Structure_Spectrum_and_Precision_Parameters Spectrum_and_Precision_Parameters;

  double M_x = 140;
  double T = 0.0001;
  // double tau_x = 3*pow(10,6);
  double tau_x = pow(T/T_0,-2)/(2*H_r);
  double Zeta_x = pow(10,-3);
  // double z = T/T_0-1;
  double z = 648069;
  fill_structure_particle_physics_model(M_x, Zeta_x, tau_x, &Particle_Physics_Model); // MANDATORY STEP
  fill_structure_spectrum_and_precision_parameters(number_iterations_photon,
                                                  number_iterations_electron,
                                                  z_step,
                                                  n_step,
                                                  calculation_mode,
                                                  photon_spectrum_choice,
                                                  electron_spectrum_choice,
                                                  spectrum_mode,
                                                  inverse_compton_scattering,
                                                  no_photons_injected,
                                                  no_electrons_injected,
                                                  &Spectrum_and_Precision_Parameters);

  // double z = 5*Particle_Physics_Model.z_x;
  mkdir("Output",01777);
  ofstream Spectrum("Output/universal_spectrum.dat");
  // print_spectrum_from_function(Spectrum, universal_spectrum, &Particle_Physics_Model);
  Cascade_Spectrum_Calculation(z,
                               &Particle_Physics_Model,
                               &Cascade_Spectrum,
                               &Spectrum_and_Precision_Parameters);
  }
  /*******************************************************/


  // /**************To compute constraints in the zeta-tau plane***********/
  if(task==2 || task==3){
  struct Structure_Particle_Physics_Model Particle_Physics_Model;
  struct Structure_Spectrum_and_Precision_Parameters Spectrum_and_Precision_Parameters;
  struct Structure_Scan_Parameters Scan_Parameters;
  struct Structure_Output_Options Output_Options;
  double tau_min = 1e4;
  double tau_max = 1e10;
  double tau_step = 100;
  double zeta_min = 1e-12;
  double zeta_max = 1e-3;
  double zeta_step = 100;
  double M_x =140;
  if(verbose>0){
  cout << "********************* Hello ! Thanks for using cBBNfast, ***********************" << endl;
  cout << "***** You're computing the constraints from BBN on EM-decaying particles *******" << endl;
  cout << "***************************** in the zeta-tau plane. ***************************" << endl;
  cout << "********** you have chosen the following range for your parameters *************" << endl;
  cout << "==> tau in ["<<tau_min<<","<<tau_max<<"]"<<" with " << tau_step << " steps.    " << endl;
  cout << "==> zeta in ["<<zeta_min<<","<<zeta_max<<"]"<<" with " << zeta_step << " steps." << endl;
  cout << "*************************** I now start computing ! ****************************" << endl;
  }
  mkdir("Output",01777);

  fill_structure_particle_physics_model(M_x, zeta_min, tau_min, &Particle_Physics_Model); // MANDATORY STEP
  fill_structure_spectrum_and_precision_parameters(number_iterations_photon,
                                                  number_iterations_electron,
                                                  z_step,
                                                  n_step,
                                                  calculation_mode,
                                                  photon_spectrum_choice,
                                                  electron_spectrum_choice,
                                                  spectrum_mode,
                                                  inverse_compton_scattering,
                                                  no_photons_injected,
                                                  no_electrons_injected,
                                                  &Spectrum_and_Precision_Parameters);
  fill_structure_scan_parameters(nuclei, tau_min, tau_max, tau_step, zeta_min, zeta_max, zeta_step, &Scan_Parameters);
  fill_output_options(print_result, results_files, spectrum_files, &Output_Options);
  if(task==2)Compute_Constraints_from_destruction_only(&Particle_Physics_Model,
                                            &Spectrum_and_Precision_Parameters,
                                            &Scan_Parameters,
                                            &Output_Options);
  if(task==3)Compute_constraints_from_destruction_and_production(&Particle_Physics_Model,
                                                      &Spectrum_and_Precision_Parameters,
                                                      &Scan_Parameters,
                                                      &Output_Options);
  }
  // // /********************************************************************************/



  cout << "I'm done ! Thanks for using cBBNFast" << endl;

  return 0;
}
