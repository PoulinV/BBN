#ifndef __BBN_PHOTONS_H__
#define __BBN_PHOTONS_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include "bbn/common.h"
#include "bbn/structures.h"


double  dsigma_compton(double  x, double  z, double x_prime, Structure_Output_Options * pt_Output_Options);
double  dsigma_phph(double  x, double  z,  double x_prime, Structure_Output_Options * pt_Output_Options);
double dsigma_NPC(double E_gamma, double z, double E_e, Structure_Output_Options * pt_Output_Options);
double dsigma_pair_creation(double z,
                            double E_e,
                            double E_gamma,
                            Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                            Structure_Output_Options * pt_Output_Options);
double dsigma_pair_creation_v2(double z,
                               double E_e,
                               double E_gamma,
                               Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                               Structure_Output_Options * pt_Output_Options);
double rate_compton(double  x, double  z);
double rate_NPC(double  x, double  z);
double rate_gg_scattering(double  x, double  z);
double rate_pair_creation(double E, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
double rate_pair_creation_v2(double E_gamma, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
double function_integrand_pair_creation(double E_e, double E_gamma, double E_gamma_bar);
double function_integrand_pair_creation_v2(double x_gamma, double gamma_e, double x_bb);
double integrand_rate_pair_creation_v2(double E_gamma, double E_gamma_bar);
double integrand_rate_pair_creation_v3(double x_gamma, double x_gamma_bar, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
double integrand_rate_pair_creation(double s);



#endif // __BBN_PHOTONS_H__
