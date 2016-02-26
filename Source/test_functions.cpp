#include "../include/EM_cascade.h"
#include "../include/injected_spectrum.h"
#include "../include/structures.h"
#include "../include/BBN_constraints.h"
#include "../include/tools.h"
#include "../include/test_functions.h"


using namespace std;




double integrate_dsigma_phph(double E_MIN, double E_MAX, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,Structure_Output_Options * pt_Output_Options){

	double  q1, q2, q3, q4, q5, q6, q7,F1, F2, F3, F4, F5, F6, F7, dq, h, rate = 0;
	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
	dq = (E_MAX-pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
	h = dq/6.;
	q1=pt_Spectrum_and_Precision_Parameters->E_min_table;
	for(int i=0;i<pt_Spectrum_and_Precision_Parameters->n_step-1;i++){
		if(i!=0){
			q1=q7;
		}
		q2=q1+h;
		q3=q2+h;
		q4=q3+h;
		q5=q4+h;
		q6=q5+h;
		q7=q6+h;

		F1 = dsigma_phph(E_MAX,z, q1,pt_Output_Options);
		F2 = dsigma_phph(E_MAX,z, q2,pt_Output_Options);
		F3 = dsigma_phph(E_MAX,z, q3,pt_Output_Options);
		F4 = dsigma_phph(E_MAX,z, q4,pt_Output_Options);
		F5 = dsigma_phph(E_MAX,z, q5,pt_Output_Options);
		F6 = dsigma_phph(E_MAX,z, q6,pt_Output_Options);
		F7 = dsigma_phph(E_MAX,z, q7,pt_Output_Options);

		rate+= dq/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
		// cout << "q = " << q1 << " rate = " << rate << endl;
	}
	return rate;
}
double integrate_dsigma_compton(double E_MIN, double E_MAX, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,Structure_Output_Options * pt_Output_Options){

	double  q1, q2, q3, q4, q5, q6, q7,F1, F2, F3, F4, F5, F6, F7, dq, h, rate = 0;
	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
	dq = (E_MAX-pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (10000*pt_Spectrum_and_Precision_Parameters->n_step-1);
	h = dq/6.;
	q1=pt_Spectrum_and_Precision_Parameters->E_min_table;
	rate = 0.;

	for(int i=0;i<10000*pt_Spectrum_and_Precision_Parameters->n_step-1;i++){
		if(i!=0){
			q1=q7;
		}
		q2=q1+h;
		q3=q2+h;
		q4=q3+h;
		q5=q4+h;
		q6=q5+h;
		q7=q6+h;

		F1 = dsigma_compton(E_MAX,z, q1,pt_Output_Options);
		if(isnan(F1)==1)F1=0.;
		F2 = dsigma_compton(E_MAX,z, q2,pt_Output_Options);
		F3 = dsigma_compton(E_MAX,z, q3,pt_Output_Options);
		F4 = dsigma_compton(E_MAX,z, q4,pt_Output_Options);
		F5 = dsigma_compton(E_MAX,z, q5,pt_Output_Options);
		F6 = dsigma_compton(E_MAX,z, q6,pt_Output_Options);
		F7 = dsigma_compton(E_MAX,z, q7,pt_Output_Options);
		// cout << " F1 = " << F1 << " F2 = " << F2 <<" F3 = " << F3 <<" F4 = " << F4 <<" F5 = " << F5 <<" F6 = " << F6 <<" F7 = " << F7 << endl;
		// cout << " pt_Spectrum_and_Precision_Parameters->E_min_table = " << pt_Spectrum_and_Precision_Parameters->E_min_table << " dq = " << dq << endl;
		rate += dq/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
		// cout << "q = " << q1 << " rate = " <<  rate << endl;
	}
	return rate;
}
double integrator_simpson_dsigma_pair_creation(double z,
																							 double E_ini,
																							 double E_max,
																							 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
																						 	 Structure_Output_Options * pt_Output_Options){

	double T=T_0*(1+z);
	double result;
	double h, f[pt_Spectrum_and_Precision_Parameters->eval_max],E[pt_Spectrum_and_Precision_Parameters->eval_max];
	// double T=T_0*(1+z), E_max = 10*2.701*T;

	double dE = (E_max- E_ini)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
	// cout << "E_max "<< E_max << " E_ini " << E_ini <<" dE = " << dE << endl;
	int y=0;
	while(dE>E_ini){
		dE/=10.;
		y++;
	}
		h = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
		// cout << "(rate_pair_creation_v2 :) dE " << dE << " E_ini " << E_ini << " y " << y << endl;
		result=0;
		for(int i=0; i<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;i++){

					// cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
					for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
					{
						if(eval == 0){
						if(i==0)	E[eval]=E_ini;
						else E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
						}
						else{
							E[eval]=E[0]+eval*h;
						}
						f[eval]=dsigma_pair_creation(z,E[eval],E_max,pt_Spectrum_and_Precision_Parameters,pt_Output_Options);
						// f[eval]=integrand_rate_pair_creation_v2(E_gamma,E[eval])/(exp(E[eval]/T)-1);
						result += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
					// cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
				}
		}

	return result;
}
double Function_Integrand_Spectre_Compton_times_bb_spectrum(double E_e, double E_gamma,  double E_gamma_bar){
	double f;
 	double Gamma_e = 4*E_gamma_bar*E_e/(m_e*m_e);
	double q = E_gamma/(Gamma_e*(E_e-E_gamma));
	double n = E_gamma_bar*E_gamma/(m_e*m_e);

	if(q>=0. && q<=1.) {
			f = 2*q*log(q)
					+(1+2*q)*(1-q)
			  	+2*n*q*(1-q);
}
// 	if(q>=0. && q<=1.) {
// 			f = 2*q*log(q)
// 					+(1+2*q)*(1-q)
// 			  	+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q);
// }
	else f=0;
	if(f<0){
		// cout << "(Function_Integrand_Spectre_Compton : ) Eg= "<<E_gamma<<" Ee =" << E_e << " f = " << f << " q=" << q <<" Gamma_e = " << Gamma_e<< " Ecmb="<< E_gamma_bar<< endl;
		f=0;
	}
	cout << " q = " << q << " Gamma_e = " << Gamma_e << " f = " << f << " E_gamma_bar =" << E_gamma_bar<<endl;
	return f*E_gamma_bar/(exp(E_gamma_bar/0.0001)-1);
}
