#include "../Include/EM_cascade.h"
#include "../Include/injected_spectrum.h"
#include "../Include/structures.h"
#include "../Include/BBN_constraints.h"
#include "../Include/tools.h"
using namespace std;

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



double universal_spectrum(double  E, double  z, double E_0){

	double  E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
	double  K_0 = E_0/(pow(E_x,2)*(2+log(E_c/E_x)));
	double  f;

	if(E < E_x) f = K_0*pow(E_x/E,1.5)/(rate_NPC(E,z)+rate_compton(E,z)+rate_gg_scattering(E,z));
	else if(E > E_x && E < E_c) f =  K_0*pow(E_x/E,2)/(rate_NPC(E,z)+rate_compton(E,z)+rate_gg_scattering(E,z));
	else {f = 0;}
	// cout << " f = " << f << " E = "<< E <<" E_c = " << E_c << " z = " << z << endl;
	return f;
}

double photon_dirac_spectrum_after_one_iteration(double  x, double  z, double E_0){


	double Gamma_tot = rate_NPC(E_0,z)+rate_compton(E_0,z)+rate_gg_scattering(E_0,z);

	double T = T_0*(1+z);
	double int_BB = 8./63.*pow(pi,4)*pow(T_0*(1+z),6);

	double spectre_gamma_gamma = 1112./(10125*pi)*pow(ALPHA*r_e,2)*pow(m_e,-6)*pow(E_0,2)*pow(1-x/E_0+pow(x/E_0,2),2)*int_BB;
	double spectre_compton = pi*r_e*r_e*m_e*pow(E_0,-2)*(x/E_0+E_0/x+pow(m_e/x-m_e/E_0,2)-2*m_e*(1/x-1/E_0))*n_e*pow(1+z,3);
	double f = (spectre_gamma_gamma+spectre_compton)/(Gamma_tot);
	if(x>E_0)f=0;
	return f;

}
double electron_dirac_spectrum_after_one_iteration(double  E, double  z, double E_0){
			 	double E_gamma_bb = 2.701*T_0*(1+z);
				double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	 			double resultat = 0, Gamma_electron_Dirac = 0, E_e;
				double dE = (E_0 - E_min)/ (double) (Electron_Table_Size-1);
				// cout << " E_gamma_bb = " << E_gamma_bb << " intbb = " << int_bb << endl;
				double F1, F2, F3, F4, F5, F6, F7;
				double q_min=0.0001, q1, q2, q3, q4, q5, q6, q7, dq, h;
				int n_step = 200;
				dq = (1-q_min)/ (double) (n_step-1);
				h = dq/6.;
				q1=q_min;
				for(int i=0;i<n_step-1;i++){
					if(i!=0){
						q1=q7;
					}
					q2=q1+h;
					q3=q2+h;
					q4=q3+h;
					q5=q4+h;
					q6=q5+h;
					q7=q6+h;

					F1 = Function_Integrand_Spectre_Compton_version_q(q1,E, E_gamma_bb);
					F2 = Function_Integrand_Spectre_Compton_version_q(q2,E, E_gamma_bb);
					F3 = Function_Integrand_Spectre_Compton_version_q(q3,E, E_gamma_bb);
					F4 = Function_Integrand_Spectre_Compton_version_q(q4,E, E_gamma_bb);
					F5 = Function_Integrand_Spectre_Compton_version_q(q5,E, E_gamma_bb);
					F6 = Function_Integrand_Spectre_Compton_version_q(q6,E, E_gamma_bb);
					F7 = Function_Integrand_Spectre_Compton_version_q(q7,E, E_gamma_bb);

					Gamma_electron_Dirac+= dq/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
				}
				Gamma_electron_Dirac*=2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb*E*E);
				resultat = 2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb*E_0*E_0*Gamma_electron_Dirac)*Function_Integrand_Spectre_Compton(E_0+E_gamma_bb-E,E_0,E_gamma_bb);

				return resultat;
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
	map_spectrum["electron_dirac_spectrum_after_one_iteration"]=electron_dirac_spectrum_after_one_iteration;
	map_spectrum["photon_dirac_spectrum_after_one_iteration"]=photon_dirac_spectrum_after_one_iteration;
	/*You need to add here your own spectrum of the type "Spectrum" as it is defined at the beginning of "structures.h", it means your_function(double E, double z, double E_0).*/
}
void check_name_spectrum(const string &value, string &error_value){
	if(value == "Dirac" || value == "none" || value == "universal"){
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
