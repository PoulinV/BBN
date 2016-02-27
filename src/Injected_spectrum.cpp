#include "../include/EM_cascade.h"
#include "../include/injected_spectrum.h"
#include "../include/structures.h"
#include "../include/BBN_constraints.h"
#include "../include/tools.h"
#include "../include/test_functions.h"

using namespace std;

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

/*
*This module has for purposes to take the injected gamma and/or e^+e^- spectra from your model of interest.
*
*Some typical (most common) spectrum are already implemented : the "Universal Spectra" for high-energy electrons or photons injection,
*and the "Dirac spectrum" for below pair-production threshold injection. See e.g.  arXiv:1503.04852 for more details.
*
*/

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/



double universal_spectrum(double  E, double  z, double E_0){

	double  E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
	double E_phph = m_e*m_e/(T_0*(1+z));
	double  K_0 = E_0/(pow(E_x,2)*(2+log(E_c/E_x)));
	double  f;
	double Rate_photons_E_g = rate_NPC(E,z)+rate_compton(E,z);
	if(E < E_phph) Rate_photons_E_g += rate_gg_scattering(E,z);
	if(E < E_x) f = K_0*pow(E_x/E,1.5)/(Rate_photons_E_g);
	else if(E > E_x && E < E_c) f =  K_0*pow(E_x/E,2)/(Rate_photons_E_g);
	else {f = 0;}
	// cout << " f = " << f << " E = "<< E <<" E_c = " << E_c << " z = " << z << endl;
	return f;
}

double no_electrons_injected(double x, double z, double E_0){
	return 0;
}

double no_photons_injected(double x, double z, double E_0){
	return 0;
}




/******************************************************************************************************************************************/
/*************************** A few functions that you need to slightly modify if you add an injection spectrum ****************************/
/***************************************** Those are the only functions you will have to modify ! *****************************************/
/******************************************************************************************************************************************/

void attribute_map_spectrum(map_spectrum &map_spectrum,map_parameters &map_parameters){
	map_spectrum["no_photons_injected"]=no_photons_injected;
	map_spectrum["no_electrons_injected"]=no_electrons_injected;
	/*You need to add here your own spectrum of the type "Spectrum" as it is defined at the beginning of "structures.h", it means your_function(double E, double z, double E_0).*/
}
void check_name_spectrum(const string &value, string &error_value){
	if(value == "Dirac" || value == "none" || value == "universal" || value == "from_file"){
		error_value = "no";
	}
	/* If you want to add a spectrum, you need to add also a "error checking" line. Here is a template you can use.
	else if{value == "name_of_your_function"}{
		error_value = "no";
	}
	*/
	else{
		error_value = "yes";
	}
}
