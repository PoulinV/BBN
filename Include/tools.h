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
void fill_structure_particle_physics_model(std::ifstream &file, const std::maps<string,string> &map_parameters, Structure_Particle_Physics_Model * pt_Particle_Model);
void fill_structure_spectrum_and_precision_parameters(std::ifstream &file,const std::maps<string,string> &map_parameters, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
void fill_default_parameters(std::ifstream &file, const std::maps<string,string> &map_parameters);
void fill_structure_scan_parameters(std::ifstream &file, const std::maps<string,string> &map_parameters, Structure_Scan_Parameters * pt_Scan_Parameters);
void fill_output_options(std::ifstream &file, const std::maps<string,string> &map_parameters, Structure_Output_Options * pt_Structure_Output_Options);
void check_energy_conservation(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                              Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                              Structure_Spectrum * pt_Gamma_Spectrum,
                              Structure_Spectrum * pt_Electron_Spectrum,
                              double &integrale);

#endif
