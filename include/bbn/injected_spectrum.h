#ifndef __BBN_INJECTED_SPECTRUM_H__
#define __BBN_INJECTED_SPECTRUM_H__

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "bbn/common.h"

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

/*
*This module has for purposes to take the injected gamma and/or e^+e^- spectra from your the model of interest.
*
*Some typical (most common) spectrum are already implemented : the "Universal Spectra" for high-energy electrons or photons injection,
*and the "Dirac spectrum" for below pair-production threshold injection. See e.g.  arXiv:1503.04852 for more details.
*
*/

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/





double universal_spectrum(double  E, double  z, double E_0, Structure_Output_Options * pt_Output_Options);
double no_electrons_injected(double x, double z, double E_0, Structure_Output_Options * pt_Output_Options);
double no_photons_injected(double x, double z, double E_0, Structure_Output_Options * pt_Output_Options);
double Dirac_Spectrum_After_One_Iteration(double  x, double  z, double E_0, Structure_Output_Options * pt_Output_Options);
void attribute_map_spectrum(std::map<std::string,Spectrum> &map_spectrum,std::map<std::string,std::string> &map_parameters);
void check_name_spectrum(const std::string &value, std::string &error_value);

#endif // __BBN_INJECTED_SPECTRUM_H__
