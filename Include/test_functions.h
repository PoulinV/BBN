#include <iostream>
#include <fstream>
#include <sstream>
#include "common.h"
#include "../include/structures.h"

double integrate_dsigma_phph(double E_MIN, double E_MAX, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,Structure_Output_Options * pt_Output_Options);
double integrate_dsigma_compton(double E_MIN, double E_MAX, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,Structure_Output_Options * pt_Output_Options);
double integrator_simpson_dsigma_pair_creation(double z,
																							 double E_ini,
																							 double E_max,
																							 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
																						 	 Structure_Output_Options * pt_Output_Options);
double Function_Integrand_Spectre_Compton_times_bb_spectrum(double E_e, double E_gamma,  double E_gamma_bar);
