#ifndef __BBN_ELECTRONS_H__
#define __BBN_ELECTRONS_H__


#include <iostream>
#include <fstream>
#include <sstream>
#include "bbn/common.h"
#include "bbn/structures.h"



double Analytical_form_inverse_compton(double s, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
double rate_electron_inverse_compton(double z,
        double E_e,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options);
double dsigma_inverse_compton_electron_spectrum_v3(double z,
        double E_e,
        double E_e_prime,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options);
double dsigma_inverse_compton_electron_spectrum_v4(double z,
        double E_e,
        double E_e_prime,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options);
double gamma_inverse_compton_analytical(double gamma_e,
                                        double E_gamma,
                                        double z,
                                        int N,
                                        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                        Structure_Output_Options * pt_Output_Options);
double gamma_inverse_compton_analytical_v2(double gamma_e, double E_gamma, double z,Structure_Output_Options * pt_Output_Options);
double dsigma_inverse_compton_electron_spectrum(double z,
        double gamma_e,
        double gamma_prime,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options);
double dsigma_inverse_compton_electron_spectrum_v2(double z,
        double gamma_e,
        double gamma_prime,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options);

#endif // __BBN_ELECTRONS_H__
