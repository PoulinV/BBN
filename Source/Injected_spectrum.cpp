#include "../Include/Injected_spectrum.h"
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



double Universal_Spectrum(double  E, double  z, double E_0){

	double  E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
	double  K_0 = E_0/(pow(E_x,2)*(2+log(E_c/E_x)));
	double  f;

	if(E < E_x) f = K_0*pow(E_x/E,1.5);
	else if(E > E_x && E < E_c) f =  K_0*pow(E_x/E,2);
	else {f = 0;}
	// cout << " f = " << f << " E = "<< E <<" E_c = " << E_c << " z = " << z << endl;
	return f;
}

double Dirac_Spectrum_After_One_Iteration(double  x, double  z, double E_0){


	double Gamma_tot = gamma_NPC(E_0,z)+gamma_compton(E_0,z)+gamma_phph(E_0,z);

	double T = T_0*(1+z);
	double int_BB = 8./63.*pow(pi,4)*pow(T_0*(1+z),6);

	double spectre_gamma_gamma = 1112./(10125*pi)*pow(a*r_e,2)*pow(m_e,-6)*pow(E_0,2)*pow(1-x/E_0+pow(x/E_0,2),2)*int_BB;
	double spectre_compton = pi*r_e*r_e*m_e*pow(E_0,-2)*(x/E_0+E_0/x+pow(m_e/x-m_e/E_0,2)-2*m_e*(1/x-1/E_0))*n_e*pow(1+z,3);
	double f = (spectre_gamma_gamma+spectre_compton)/(Gamma_tot);
	if(x>E_0)f=0;
	return f;

}
double No_Electrons_Injected(double x, double z, double E_0){
	return 0;
}
double No_Photons_Injected(double x, double z, double E_0){
	return 0;
}
// void Initialize_Spectrum(double (*func)(double,double,struct),struct Structure_Particle_Model * pt_Particle_Model,struct Structure_Gamma_Spectrum * pt_Injected_Spectrum){
//
// 	pt_Cascade_Spectrum->Gamma_Energy.resize(Gamma_Table_Size);
//
//
// }
