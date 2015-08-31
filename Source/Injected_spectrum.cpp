#include "../Include/injected_spectrum.h"
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


	double Gamma_tot = rate_NPC(E_0,z)+rate_compton(E_0,z)+rate_gg_scattering(E_0,z);

	double T = T_0*(1+z);
	double int_BB = 8./63.*pow(pi,4)*pow(T_0*(1+z),6);

	double spectre_gamma_gamma = 1112./(10125*pi)*pow(a*r_e,2)*pow(m_e,-6)*pow(E_0,2)*pow(1-x/E_0+pow(x/E_0,2),2)*int_BB;
	double spectre_compton = pi*r_e*r_e*m_e*pow(E_0,-2)*(x/E_0+E_0/x+pow(m_e/x-m_e/E_0,2)-2*m_e*(1/x-1/E_0))*n_e*pow(1+z,3);
	double f = (spectre_gamma_gamma+spectre_compton)/(Gamma_tot);
	if(x>E_0)f=0;
	return f;

}
void Electron_dirac_spectrum_after_one_iteration(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
								                                   struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
								                                   struct Structure_Spectrum * pt_Electron_Spectrum){
			  double z = pt_Electron_Spectrum->redshift;
				double E_0 = pt_Particle_Physics_Model->E_0;
			 	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
			 	double E_gamma_bb = 2.701*T_0*(1+z);
				double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	 			double resultat = 0, Gamma_electron_Dirac = 0, E_e;
				double dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Electron_Table_Size-1);
				// cout << " E_gamma_bb = " << E_gamma_bb << " intbb = " << int_bb << endl;
				Gamma_electron_Dirac = Rate_Inverse_Compton(E_0,z,pt_Spectrum_and_Precision_Parameters);
				for(int i = 0 ; i < Electron_Table_Size ; i++){
					E_e = E_min + i*dE;
					pt_Electron_Spectrum->Energy[i] = E_e;
					resultat = 2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb*E_0*E_0*Gamma_electron_Dirac)*Function_Integrand_Spectre_Compton(E_0+E_gamma_bb-E_e,E_0,E_gamma_bb);
					pt_Electron_Spectrum->Spectrum[i] = resultat/Rate_Inverse_Compton(E_e,z,pt_Spectrum_and_Precision_Parameters);
					// cout << " E = " << E_e << " resultat = " << pt_Electron_Spectrum->Spectrum[i] << endl;
				}
}

void No_Electrons_Injected(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                           struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                           struct Structure_Spectrum * pt_Electron_Spectrum){

		 double dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Electron_Table_Size-1);
		 for(int i = 0 ; i < Electron_Table_Size ; i++){
				pt_Electron_Spectrum->Energy[i] = E_min + i*dE;
				pt_Electron_Spectrum->Spectrum[i] = 0;
			}
}

double No_Photons_Injected(double x, double z, double E_0){
	return 0;
}
