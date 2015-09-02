#ifndef __tools__
#define __tools__


#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "common.h"
/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

/*
*
*This module contains all the numerical tools necesserary to the code.
*Most of the algorithm comes from Numerical Recipes.
*
*/

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/


void linearint(std::vector<double> &xa, std::vector<double> &ya, int n, double x, double &y);
void print_spectrum(std::ofstream &file, Structure_Spectrum * pt_Cascade_Spectrum, Structure_Particle_Physics_Model * pt_Particle_Model);
void print_spectrum_from_function(std::ifstream &file, double (*func)(double,double,double), Structure_Particle_Physics_Model * pt_Particle_Model);
void print_spectrum_automatic_names(int number_files, Structure_Spectrum * pt_Cascade_Spectrum, Structure_Particle_Physics_Model * pt_Particle_Model);
void fill_structure_particle_physics_model(std::ifstream &file, std::map<std::string,std::string> &map_parameters, Structure_Particle_Physics_Model * pt_Particle_Model);
void fill_structure_spectrum_and_precision_parameters(std::ifstream &file, std::map<std::string,std::string> &map_parameters, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
void fill_default_parameters(std::ifstream &file, std::map<std::string,std::string> &map_parameters);
void fill_structure_scan_parameters(std::ifstream &file, std::map<std::string,std::string> &map_parameters, Structure_Scan_Parameters * pt_Scan_Parameters);
void fill_output_options(std::ifstream &file, std::map<std::string,std::string> &map_parameters, Structure_Output_Options * pt_Structure_Output_Options);
void attribute_name_and_value(const std::string &line,std::string &name,std::string &value);
void check_value_and_name_error(const std::string &name,std::string &error_name, const std::string &value,std::string &error_value);
void check_parameter_error(const std::string &parameter, std::string &error_parameter);
void get_parameter_from_file(std::ifstream &file, std::string &parameter, std::string &value);
void check_energy_conservation(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                              Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                              Structure_Spectrum * pt_Gamma_Spectrum,
                              Structure_Spectrum * pt_Electron_Spectrum,
                              double &integrale);


#endif
