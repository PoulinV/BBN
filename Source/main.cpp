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

#include "../include/common.h"
#include "../include/structures.h"
#include "../include/EM_cascade.h"
#include "../include/injected_spectrum.h"
#include "../include/tools.h"
#include "../include/BBN_constraints.h"


int main(int argc, char** argv){
  if(argc!=2){
    cout <<"Error : you have not given an input file! \nPlease, restart cBBNFast with one '.ini' file. If you want to use default parameters, just give 'default_param.ini'." << endl;
    return 0;
  }

  mkdir("Output",01777);
  mkdir("Output/Cascade_Spectrum_Folder",01777);
  mkdir("Output/Result_Scan_Folder",01777);

  map_parameters map_parameters;
  ifstream file_default("default_param.ini");
  ifstream file_input(argv[1]);
  fill_default_parameters(file_default,map_parameters);
  string task = "task";
  if(argc==2)get_parameter_from_file(file_input,task);
  if(task=="default" || argc==1) task=map_parameters["task"];


  struct Structure_Particle_Physics_Model Particle_Physics_Model;
  struct Structure_Spectrum_and_Precision_Parameters Spectrum_and_Precision_Parameters;
  struct Structure_Scan_Parameters_and_Results Scan_Parameters_and_Results;
  struct Structure_Output_Options Output_Options;



  // /************To print cascade spectrum in a file*********/
  if(task=="compute_cascade_spectrum"){

    cout << "************************** Hello ! Thanks for using cBBNfast, ***************************" << endl;
    cout << "********* You're computing the cascade spectrum from a decaying particles with **********" << endl;
    cout << "***************************** m_x = " << map_parameters["m_x"] << "MeV at z = " << map_parameters["redshift"] <<" *******************************" << endl;
    cout << "******************************** I now start computing ! ********************************" << endl;

  struct Structure_Spectrum Cascade_Spectrum;
  string redshift = "redshift";
  fill_structure_spectrum_and_precision_parameters(file_input, map_parameters, &Spectrum_and_Precision_Parameters);
  fill_structure_particle_physics_model(file_input, map_parameters, &Particle_Physics_Model); // MANDATORY STEP
  fill_structure_output_options(file_input, map_parameters, &Output_Options);
  if(argc==2)get_parameter_from_file(file_input,redshift);
  if(redshift=="default" || argc==1)redshift=map_parameters["redshift"];
  ofstream Spectrum("Output/universal_spectrum.dat");
  Cascade_Spectrum_Calculation(atof(redshift.c_str()),
                               &Output_Options,
                               &Particle_Physics_Model,
                               &Cascade_Spectrum,
                               &Spectrum_and_Precision_Parameters);
  }
  /*******************************************************/


  // /**************To compute constraints in the zeta-tau plane***********/
  else if(task=="compute_constraints_from_destruction_only"||task=="compute_constraints_from_destruction_and_production"){
  cout << "********************* Hello ! Thanks for using cBBNfast, ***********************" << endl;
  cout << "***** You're computing the constraints from BBN on EM-decaying particles *******" << endl;
  cout << "***************************** in the zeta-tau plane. ***************************" << endl;
  cout << "********** you have chosen the following range for your parameters *************" << endl;
  cout << "==> tau in ["<<map_parameters["tau_min"]<<","<<map_parameters["tau_max"]<<"]"<<" with " << map_parameters["tau_step"] << " steps.    " << endl;
  cout << "==> zeta in ["<<map_parameters["zeta_min"]<<","<<map_parameters["zeta_max"]<<"]"<<" with " << map_parameters["zeta_step"] << " steps." << endl;
  cout << "*************************** I now start computing ! ****************************" << endl;


  fill_structure_particle_physics_model(file_input,map_parameters,&Particle_Physics_Model); // MANDATORY STEP
  fill_structure_spectrum_and_precision_parameters(file_input,map_parameters,&Spectrum_and_Precision_Parameters);
  fill_structure_scan_parameters_and_results(file_input,map_parameters,&Scan_Parameters_and_Results);
  fill_structure_output_options(file_input,map_parameters,&Output_Options);
  cout << "here " <<endl;
  if(task=="compute_constraints_from_destruction_only")Compute_Constraints_from_destruction_only(&Particle_Physics_Model,
                                            &Spectrum_and_Precision_Parameters,
                                            &Scan_Parameters_and_Results,
                                            &Output_Options);
  if(task=="compute_constraints_from_destruction_and_production")Compute_constraints_from_destruction_and_production(&Particle_Physics_Model,
                                                      &Spectrum_and_Precision_Parameters,
                                                      &Scan_Parameters_and_Results,
                                                      &Output_Options);
  }
  // // /********************************************************************************/


  file_default.close();
  file_input.close();
  cout << "I'm done ! Thanks for using cBBNFast." << endl;

  return 0;
}
