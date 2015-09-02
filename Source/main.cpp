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
    cout <<"Error : you have given too many files as an input! \nPlease, restart cBBNFast with only one '.ini' file." << endl;
    return 0;
  }

  mkdir("Output",777);
  mkdir("Output/Cascade_Spectrum_Folder",777);
  mkdir("Output/Result_Scan_Folder",777);

  map_parameters map_parameters;
  ifstream file_default("default_param.ini");
  ifstream file_input(argv[1]);
  fill_default_parameters(file_default,map_parameters);
  string task = "task";
  if(argc==2)get_parameter_from_file(file_input,task,task);
  else get_parameter_from_file(file_default,task,task);
  if(task=="task"){
    get_parameter_from_file(file_default,task,task);
  }
  cout << map_parameters["spectrum_files"]<<" " << map_parameters["z_step"] << endl;
  cout << task << endl;


  struct Structure_Particle_Physics_Model Particle_Physics_Model;
  struct Structure_Spectrum_and_Precision_Parameters Spectrum_and_Precision_Parameters;
  struct Structure_Scan_Parameters Scan_Parameters;
  struct Structure_Output_Options Output_Options;



  // /************To print cascade spectrum in a file*********/
  if(map_parameters["task"]=="compute_spectrum"){
  struct Structure_Spectrum Cascade_Spectrum;
  string redshift = "redshift";
  fill_structure_particle_physics_model(file_input, map_parameters, &Particle_Physics_Model); // MANDATORY STEP
  fill_structure_spectrum_and_precision_parameters(file_input, map_parameters, &Spectrum_and_Precision_Parameters);
  if(argc==2)get_parameter_from_file(file_input,redshift,redshift);
  else get_parameter_from_file(file_default,redshift,redshift);
  if(redshift=="redshift"){
    get_parameter_from_file(file_default,task,task);
  }
  ofstream Spectrum("Output/universal_spectrum.dat");
  Cascade_Spectrum_Calculation(atof(redshift.c_str()),
                               &Particle_Physics_Model,
                               &Cascade_Spectrum,
                               &Spectrum_and_Precision_Parameters);
  }
  /*******************************************************/


  // /**************To compute constraints in the zeta-tau plane***********/
  else if(map_parameters["task"]=="compute_constraints_from_destruction_only"||map_parameters["task"]=="compute_constraints_from_destruction_and_production"){

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

  fill_structure_particle_physics_model(file_input,map_parameters,&Particle_Physics_Model); // MANDATORY STEP
  fill_structure_spectrum_and_precision_parameters(file_input,map_parameters,&Spectrum_and_Precision_Parameters);
  fill_structure_scan_parameters(file_input,map_parameters,&Scan_Parameters);
  fill_output_options(file_input,map_parameters,&Output_Options);
  if(map_parameters["task"]=="compute_constraints_from_destruction_only")Compute_Constraints_from_destruction_only(&Particle_Physics_Model,
                                            &Spectrum_and_Precision_Parameters,
                                            &Scan_Parameters,
                                            &Output_Options);
  if(map_parameters["task"]=="compute_constraints_from_destruction_and_production")Compute_constraints_from_destruction_and_production(&Particle_Physics_Model,
                                                      &Spectrum_and_Precision_Parameters,
                                                      &Scan_Parameters,
                                                      &Output_Options);
  }
  // // /********************************************************************************/


  file_default.close();
  file_input.close();
  cout << "I'm done ! Thanks for using cBBNFast." << endl;

  return 0;
}
