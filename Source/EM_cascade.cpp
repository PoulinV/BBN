#include "../include/EM_cascade.h"
#include "../include/injected_spectrum.h"
#include "../include/structures.h"
#include "../include/BBN_constraints.h"
#include "../include/tools.h"

using namespace std;

double  dsigma_compton(double  x, double  z, double g){

	double  dsigma =0;
	dsigma = pi*pow(r_e,2)*m_e/pow(x,2)*(x/g+g/x+pow(m_e/g-m_e/x,2)-2*m_e*(1/g-1/x))*n_e*pow(1+z,3)*(1+Y/2)/(1+Y);
	// if(g<x && g>(x/(1+2*x)))dsigma = pi*pow(r_e,2)*m_e/pow(x,2)*(x/g+g/x+pow(m_e/g-m_e/x,2)-2*m_e*(1/g-1/x))*n_e*pow(1+z,3)/(1+Y);
 	return dsigma;

}
double  dsigma_phph(double  x, double  z,  double g){

	double  dsigma = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(ALPHA*r_e,2)*pow(m_e,-6)*pow(x,2)*pow(1-g/x+pow(g/x,2),2)/(63*10125*pi);
	return dsigma;

}

double dsigma_NPC(double E_gamma, double z, double E_e){

	double p = sqrt(E_e*E_e-m_e*m_e);
	double E_pos = E_gamma-E_e-m_e;
	double p_pos = E_pos*E_pos-m_e*m_e;
	if(p_pos>0){
		p_pos = sqrt(p_pos);
	}
	else return 0;
	double L = log((E_pos*E_e+p_pos*p+m_e*m_e)/(E_pos*E_e-p_pos*p+m_e*m_e));
	double l = log((E_e+p)/(E_e-p));
	double l_pos = log((E_pos+p_pos)/(E_pos-p_pos));

	double result =  ALPHA*pow(r_e,2)*pow(1+z,3)*eta*n_y_0*p*p_pos/pow(E_gamma,3)
				 *(-4./3-2*E_pos*E_e*(p_pos*p_pos+p*p)/(p*p*p_pos*p_pos)+m_e*m_e*(l*E_pos/pow(p,3)+l_pos*E_e/pow(p_pos,3)-l*l_pos/(p_pos*p))
				 +L*(-8*E_pos*E_e/(3*p_pos*p)+E_gamma*E_gamma/pow(p_pos*p,3)*(pow(E_pos*E_e,2)+pow(p_pos*p,2)-m_e*m_e*E_pos*E_e)-m_e*m_e*E_gamma/(2*p_pos*p)*(l_pos*(E_pos*E_e-pow(p_pos,2))/pow(p_pos,3)+l*(E_pos*E_e-p*p)/pow(p,3))));
	// result = 0;
	return 2*result;
}

// double dsigma_pair_creation(double E_gamma, double z, double s){
// s = pow(E_e/m_e,2);
// }

double  rate_compton(double  x, double  z){

	double X = 2*x/m_e;
	double sigma_cs = 2*pi*pow(r_e,2)/X*((1-4/X-8/pow(X,2))*log(1+X)+1/2+8/X+1/(2*pow(1+X,2)));
	double Gamma = sigma_cs*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3);
	return Gamma;

}
double  rate_NPC(double  x, double  z){

	double  k = x/m_e;
	double  rho = (2*k-4)/(k+2+2*pow(2*k,0.5));
	double  sigma_PCN;
	sigma_PCN = ALPHA*pow(r_e,2)*(28/9*log(2*k)-218/27) ;
	double  Gamma_2 = sigma_PCN*pow(1+z,3)*eta*n_y_0;
	return Gamma_2;

}

double  rate_gg_scattering(double  x, double  z){

	double  Gamma_3;
	Gamma_3 = 0.1513*pow(ALPHA,4)*m_e*pow(x/m_e,3)*pow(T_0*(1+z)/m_e,6);
	return Gamma_3;

}
// double rate_pair_creation(double E, double z){
// 	double E_gamma_bb = 2.701*T_0*(1+z);
// 	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
// 	return 1/(8*E)*int_bb/pow(E_gamma_bb,2)*integrator_Weddle_Hardy(function_integrand_pair_creation,4*m_e*m_e,4*E*E_gamma_bb);
// }




double Function_Integrand_Spectre_Compton(long double E_gamma, double E_e, long double E_gamma_bar){
	double f;
 	long double Gamma_e = 4*E_gamma_bar*E_e/(m_e*m_e);
	long double q = E_gamma/(Gamma_e*(E_e-E_gamma));

	if(q>=0. && q<=1.) {f = 2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q);
}
	else f=0;
	// cout << "Eg= "<<E_gamma<<"Ee =" << E_e << "f = " << f << "q=" << q <<" Gamma_e = " << Gamma_e<< "Ecmb="<< E_gamma_bar<< endl;
	return f;
}

double Function_Integrand_Spectre_Compton_version_q(double q, double E_e, double E_gamma_bar){
	double f;
	double Gamma_e = 4*E_gamma_bar*E_e/(m_e*m_e);
	double E_gamma=Gamma_e*E_e*q/(1+Gamma_e*q);

	if(q>=0. && q<=1.) {f = (2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q))*pow(E_e-E_gamma,2)*Gamma_e/E_e;
}
	else f=0;
	if(f<=0)f=0;
	// cout <<"q = "<< q << " Gamma_e = "<<Gamma_e<<" E_gamma = "<<E_gamma<<" E_e ="<<E_e<<" f = " << f << endl;

	return f;
}
void  Spectre_electron_compton(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
													 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
													 Structure_Spectrum * pt_Gamma_Spectrum,
													 Structure_Spectrum * pt_Electron_Spectrum){

 	double z = pt_Electron_Spectrum->redshift;
	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
	double E_0 = pt_Particle_Physics_Model->E_0;
	double E_gamma_bb = 2.701*T_0*(1+z);
	double E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	double resultat = 0, Gamma_electron = 0;
	double dE, dE_2, h;
	double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7;
	double E_gamma, E_e;
	double E_gamma_1, E_gamma_2, E_gamma_3;
	double resultat_monochromatique_1, resultat_monochromatique_2;
	dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Electron_Table_Size-1);
	for(int j =0; j<Electron_Table_Size;j++){
		resultat=0;
		Gamma_electron=0;
		E_e = E_min+j*dE;
		dE_2 = (pt_Particle_Physics_Model->E_0 - (E_e))/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
		h = dE_2/6.;

		for(int i=0;i<n_step-1;i++){


			if(i==0){
				E1=E_e;
				}
			else{
				E1=E7;
			}

			E2=E1 + h;
			E3=E1 + 2*h;
			E4=E1 + 3*h;
			E5=E1 + 4*h;
			E6=E1 + 5*h;
			E7=E1 + 6*h;
			if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E1, f1);
			else f1=0;
			if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E2, f2);
			else f2=0;
			if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E3, f3);
			else f3=0;
			if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E4, f4);
			else f4=0;
			if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E5, f5);
			else f5=0;
			if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E6, f6);
			else f6=0;
			if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E7, f7);
			else f7=0;

			f1*=(dsigma_compton(E1,z,(E1+m_e-E_e))+dsigma_NPC(E1+m_e,z,E_e));
			f2*=(dsigma_compton(E2,z,(E2+m_e-E_e))+dsigma_NPC(E2+m_e,z,E_e));
			f3*=(dsigma_compton(E3,z,(E3+m_e-E_e))+dsigma_NPC(E3+m_e,z,E_e));
			f4*=(dsigma_compton(E4,z,(E4+m_e-E_e))+dsigma_NPC(E4+m_e,z,E_e));
			f5*=(dsigma_compton(E5,z,(E5+m_e-E_e))+dsigma_NPC(E5+m_e,z,E_e));
			f6*=(dsigma_compton(E6,z,(E6+m_e-E_e))+dsigma_NPC(E6+m_e,z,E_e));
			f7*=(dsigma_compton(E7,z,(E7+m_e-E_e))+dsigma_NPC(E7+m_e,z,E_e));
			resultat += dE_2/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);

		}
		Gamma_electron=Rate_Inverse_Compton(E_e,z,pt_Spectrum_and_Precision_Parameters);
		if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice=="Dirac")resultat += (dsigma_compton(E_0,z,(E_0+m_e-E_e))+dsigma_NPC(E_0+m_e,z,E_0))/(rate_NPC(E_0,z)+rate_compton(E_0,z)+rate_gg_scattering(E_0,z)) ;
		pt_Electron_Spectrum->Energy[j]=E_e;
		pt_Electron_Spectrum->Spectrum[j]=resultat;
		pt_Electron_Spectrum->Spectrum[j]/=(Gamma_electron);
		// cout << " E = " << pt_Electron_Spectrum->Energy[j] << " resultat = " << pt_Electron_Spectrum->Spectrum[j] << endl;

	}

}

double Rate_Inverse_Compton(double E_e, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){
	double rate = 0;
 	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	double E_gamma_bb = 2.701*T_0*(1+z);
	double F1, F2, F3, F4, F5, F6, F7;
	double q_min=0.0001, q1, q2, q3, q4, q5, q6, q7, dq, h;
	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
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

		F1 = Function_Integrand_Spectre_Compton_version_q(q1,E_e, E_gamma_bb);
		F2 = Function_Integrand_Spectre_Compton_version_q(q2,E_e, E_gamma_bb);
		F3 = Function_Integrand_Spectre_Compton_version_q(q3,E_e, E_gamma_bb);
		F4 = Function_Integrand_Spectre_Compton_version_q(q4,E_e, E_gamma_bb);
		F5 = Function_Integrand_Spectre_Compton_version_q(q5,E_e, E_gamma_bb);
		F6 = Function_Integrand_Spectre_Compton_version_q(q6,E_e, E_gamma_bb);
		F7 = Function_Integrand_Spectre_Compton_version_q(q7,E_e, E_gamma_bb);

		rate+= dq/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
	}
	rate*=2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb*E_e*E_e);
	// cout << "Ee = " << E_e << " rate = " << rate << endl;
	return rate;
}

void Spectrum_electron_scattered(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
																 Structure_Spectrum * pt_Input_Electron_Spectrum,
															 	 Structure_Spectrum * pt_Output_Electron_Spectrum){

	  	double z = pt_Input_Electron_Spectrum->redshift;
		 	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
		 	double E_0 = pt_Particle_Physics_Model->E_0;
		 	double E_gamma_bb = 2.701*T_0*(1+z);
		 	double E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
		 	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
		 	double resultat = 0, Gamma_electron = 0, Gamma_electron_Dirac=0, Gamma_electron_test=0;
		 	double dE, dE_2, h;
		 	double E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7, g1, g2, g3, g4, g5, g6, g7;
		 	double E_gamma, E_e;
		 	double E_gamma_1, E_gamma_2, E_gamma_3;
		 	double resultat_monochromatique_1, resultat_monochromatique_2;
		 	dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Electron_Table_Size-1);



		 	for(int j =0; j<Electron_Table_Size;j++){
		 		resultat=0;
		 		Gamma_electron=0;
				E_e = E_min+j*dE;
				n_step = pt_Spectrum_and_Precision_Parameters->n_step;
				if(E_e<15){
					n_step *=50;
				}

				dE_2 = (pt_Particle_Physics_Model->E_0 - (E_e))/ (double) (n_step-1);

				h = dE_2/6.;

				for(int i=0;i<n_step-1;i++){

					if(i==0){
						E1=E_e;
						}
					else{
						E1=E7;
					}

					E2=E1 + h;
					E3=E1 + 2*h;
					E4=E1 + 3*h;
					E5=E1 + 4*h;
					E6=E1 + 5*h;
					E7=E1 + 6*h;

					if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E1, g1);
          else g1=0;
          if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E2, g2);
          else g2=0;
          if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E3, g3);
          else g3=0;
          if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E4, g4);
          else g4=0;
          if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E5, g5);
          else g5=0;
          if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E6, g6);
          else g6=0;
          if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E7, g7);
          else g7=0;

          F1=g1*Function_Integrand_Spectre_Compton(E1+E_gamma_bb-E_e,E1, E_gamma_bb)/((E1)*(E1));
          F2=g2*Function_Integrand_Spectre_Compton(E2+E_gamma_bb-E_e,E2, E_gamma_bb)/((E2)*(E2));
          F3=g3*Function_Integrand_Spectre_Compton(E3+E_gamma_bb-E_e,E3, E_gamma_bb)/((E3)*(E3));
          F4=g4*Function_Integrand_Spectre_Compton(E4+E_gamma_bb-E_e,E4, E_gamma_bb)/((E4)*(E4));
          F5=g5*Function_Integrand_Spectre_Compton(E5+E_gamma_bb-E_e,E5, E_gamma_bb)/((E5)*(E5));
          F6=g6*Function_Integrand_Spectre_Compton(E6+E_gamma_bb-E_e,E6, E_gamma_bb)/((E6)*(E6));
          F7=g7*Function_Integrand_Spectre_Compton(E7+E_gamma_bb-E_e,E7, E_gamma_bb)/((E7)*(E7));

          resultat+= dE_2/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);

					}

		 		pt_Output_Electron_Spectrum->Energy[j]=E_e;
				// cout << " E = " << pt_Output_Electron_Spectrum->Energy[j] << " resultat before constant = " << resultat << endl;

		 		Gamma_electron=Rate_Inverse_Compton(E_e,z,pt_Spectrum_and_Precision_Parameters);
				pt_Output_Electron_Spectrum->Spectrum[j]=resultat*2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb);
				// cout << " E = " << pt_Output_Electron_Spectrum->Energy[j] << " resultat spectre diffuse avant = " << pt_Output_Electron_Spectrum->Spectrum[j] << " Gamma_electron = " << Gamma_electron << " Gamma_test = " << Gamma_electron_test<<  endl;

				pt_Output_Electron_Spectrum->Spectrum[j]/=(Gamma_electron);
				// cout << " E = " << pt_Output_Electron_Spectrum->Energy[j] << " resultat spectre diffuse = " << pt_Output_Electron_Spectrum->Spectrum[j] << endl;
		}

}



void Spectre_gamma_compton(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
													 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
													 Structure_Spectrum * pt_Electron_Spectrum,
													 Structure_Spectrum * pt_Gamma_Spectrum){

	double z = pt_Gamma_Spectrum->redshift;

	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
	double E_gamma_bb = 2.701*T_0*(1+z);
	double E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	double resultat = 0;
	double dE, dE_2, h;
	double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7, F1, F2, F3, F4, F5, F6, F7;
	double E_gamma, E_e;
	double q_min=0.0001, q1, q2, q3, dq;
	double Gamma_electron=0;
	double E_0=pt_Particle_Physics_Model->E_0;
	dq = (1-q_min)/ (double) (n_step-1);
	dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Gamma_Table_Size-1);

		for(int j =0; j<Gamma_Table_Size;j++){
			resultat=0;
			E_gamma = E_min + j*dE;
			dE_2 = (pt_Particle_Physics_Model->E_0 - (E_gamma+m_e))/ (double) (n_step-1);
			h = dE_2/6.;
			for(int i=0;i<n_step-1;i++){
				if(i==0){
					E1=E_gamma+m_e;
				}
				else{
					E1=E7;
				}


				E2=E1 + h;
				E3=E1 + 2*h;
				E4=E1 + 3*h;
				E5=E1 + 4*h;
				E6=E1 + 5*h;
				E7=E1 + 6*h;
				if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E1, f1);
				else f1=0;
				if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E2, f2);
				else f2=0;
				if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E3, f3);
				else f3=0;
				if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E4, f4);
				else f4=0;
				if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E5, f5);
				else f5=0;
				if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E6, f6);
				else f6=0;
				if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E7, f7);
				else f7=0;
				F1 = Function_Integrand_Spectre_Compton(E_gamma,E1, E_gamma_bb);
				F2 = Function_Integrand_Spectre_Compton(E_gamma,E2, E_gamma_bb);
				F3 = Function_Integrand_Spectre_Compton(E_gamma,E3, E_gamma_bb);
				F4 = Function_Integrand_Spectre_Compton(E_gamma,E4, E_gamma_bb);
				F5 = Function_Integrand_Spectre_Compton(E_gamma,E5, E_gamma_bb);
				F6 = Function_Integrand_Spectre_Compton(E_gamma,E6, E_gamma_bb);
				F7 = Function_Integrand_Spectre_Compton(E_gamma,E7, E_gamma_bb);


				f1*=2*F1/(E1*E1);
				f2*=2*F2/(E2*E2);
				f3*=2*F3/(E3*E3);
				f4*=2*F4/(E4*E4);
				f5*=2*F5/(E5*E5);
				f6*=2*F6/(E6*E6);
				f7*=2*F7/(E7*E7);
				// if(f1!=0){cout << "E1 = " << E1<< " f1 = " << f1 << endl;
				// }
				// cout << "E = " << E_min+j*dE<<"F1 = "<< F1 << " f1 apres = " << f1 << endl;

				resultat+= dE_2/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);

			}
			pt_Gamma_Spectrum->Energy[j]=E_gamma ;
			pt_Gamma_Spectrum->Spectrum[j]=resultat*2*pi*r_e*r_e*m_e*m_e*int_bb/E_gamma_bb/(rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z));
			if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice=="Dirac"){
				Gamma_electron = Rate_Inverse_Compton(E_0,z,pt_Spectrum_and_Precision_Parameters);
				pt_Gamma_Spectrum->Spectrum[j]+=2*pi*r_e*r_e*m_e*m_e*2*int_bb*Function_Integrand_Spectre_Compton(E_gamma,E_0,E_gamma_bb)/(E_gamma_bb*Gamma_electron*E_0*E_0)/(rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z));
			}
			// cout << "E = "  << E_min+j*dE<< "resultat = " << 	pt_Gamma_Spectrum->Spectrum[j] << endl;

		}


}
void Spectrum_gamma_scattered(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
																 Structure_Spectrum * pt_Input_Gamma_Spectrum,
															 	 Structure_Spectrum * pt_Output_Gamma_Spectrum){


 double dE, dE_2, h, h2;
 double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7,E_gamma;
 double resultat ;
 double z = pt_Input_Gamma_Spectrum->redshift;
 dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Gamma_Table_Size-1);

	for(int i=0;i<Gamma_Table_Size;i++){
		E_gamma = E_min+i*dE;
		resultat = 0;
					for(int j=0; j<pt_Spectrum_and_Precision_Parameters->n_step-1;j++){
						dE_2 = (pt_Particle_Physics_Model->E_0 - (E_gamma))/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
						h2 = dE_2/6.;

						if(j==0){
							E1=E_gamma;
							}
						else{
							E1=E7;
						}

						E2=E1 + h2;
						E3=E1 + 2*h2;
						E4=E1 + 3*h2;
						E5=E1 + 4*h2;
						E6=E1 + 5*h2;
						E7=E1 + 6*h2;
						if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E1, f1);
						else f1=0;
						if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E2, f2);
						else f2=0;
						if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E3, f3);
						else f3=0;
						if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E4, f4);
						else f4=0;
						if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E5, f5);
						else f5=0;
						if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E6, f6);
						else f6=0;
						if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E7, f7);
						else f7=0;
						// f1*=(dsigma_compton(E1,z,E_gamma));
						// f2*=(dsigma_compton(E2,z,E_gamma));
						// f3*=(dsigma_compton(E3,z,E_gamma));
						// f4*=(dsigma_compton(E4,z,E_gamma));
						// f5*=(dsigma_compton(E5,z,E_gamma));
						// f6*=(dsigma_compton(E6,z,E_gamma));
						// f7*=(dsigma_compton(E7,z,E_gamma));
						f1*=(dsigma_phph(E1,z,E_gamma)+dsigma_compton(E1,z,E_gamma));
						f2*=(dsigma_phph(E2,z,E_gamma)+dsigma_compton(E2,z,E_gamma));
						f3*=(dsigma_phph(E3,z,E_gamma)+dsigma_compton(E3,z,E_gamma));
						f4*=(dsigma_phph(E4,z,E_gamma)+dsigma_compton(E4,z,E_gamma));
						f5*=(dsigma_phph(E5,z,E_gamma)+dsigma_compton(E5,z,E_gamma));
						f6*=(dsigma_phph(E6,z,E_gamma)+dsigma_compton(E6,z,E_gamma));
						f7*=(dsigma_phph(E7,z,E_gamma)+dsigma_compton(E7,z,E_gamma));

						resultat += dE_2/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);
						// cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << " i = " << i << endl;
					}
					pt_Output_Gamma_Spectrum->Energy[i]=E_gamma;
					pt_Output_Gamma_Spectrum->Spectrum[i]=resultat/(rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z));
				}
}
void Cascade_Spectrum_Reading_From_File(double z,
																				Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																				Structure_Spectrum * pt_Spectrum,
																				Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){
  ostringstream os;
  string name;
	double tmp_Energy, tmp_Spectrum;
  pt_Spectrum->redshift = z;

	if(pt_Spectrum->species == "photon"){
		if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="iterative")  os << "Output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << pt_Spectrum_and_Precision_Parameters->number_iterations_photon <<"iterations.dat";
		else if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="triangular")  os << "Output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << "triangular.dat";
	}
	else if(pt_Spectrum->species == "electron"){
		if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="iterative")  os << "Output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << pt_Spectrum_and_Precision_Parameters->number_iterations_electron <<"iterations.dat";
		else if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="triangular")  os << "Output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << "triangular.dat";
	}
	  name = os.str();
    ifstream file(name);
    if(file)cout << "Importing file " << name << " in structure Cascade_Spectrum." << endl;
		else{
			cout << "I couldn't recognize cascade spectrum file. Please check that it is present in the folder Cascade_Specrum_Folder with proper name : 'Spectrum_mXXX_zXXX_XXXiterations.dat' corresponding to the value of m, z and iterations you are using."<<endl;
			return;
		}
    while(file){
			string line;
	    getline(file, line);
	    // stringstream is(ligne);
	    if(line[0] == '#' or line[0] == '\0') continue;
      file >> tmp_Energy >> tmp_Spectrum ;
			pt_Spectrum->Energy.push_back(tmp_Energy);
			pt_Spectrum->Spectrum.push_back(tmp_Spectrum*(rate_NPC(tmp_Energy,z)+rate_compton(tmp_Energy,z)+rate_gg_scattering(tmp_Energy,z)));
			  // cout << z << "  " << tmp_Energy << "  " << tmp_Spectrum<<endl;

		}

		file.close();


}


void Triangular_Spectrum(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
																 Structure_Spectrum * pt_Cascade_Spectrum,
															 	 Structure_Spectrum * pt_Electron_Spectrum){
	double E_e_minus_1, E_e_plus_1, E_j, E_j_minus_1, E_j_plus_1, dE_j, dE;
	double resultat, E_e, E_g;
	double z = pt_Cascade_Spectrum->redshift;
	double E_gamma_bb = 2.701*T_0*(1+z);
	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	double E_0 = pt_Particle_Physics_Model->E_0;
			for(double i = (Electron_Table_Size-1); i>=0 ;i--){
			resultat=0;
			E_e = E_min*pow(E_0/E_min,(double) i/(Electron_Table_Size-1));
			E_g = E_e;
			E_e_plus_1 = E_min*pow(E_0/E_min,((double) i+1)/(Electron_Table_Size-1));
			E_e_minus_1 = E_min*pow(E_0/E_min,((double) i-1)/(Electron_Table_Size-1));
			dE = (E_e_plus_1 - E_e_minus_1)/2.;
			pt_Electron_Spectrum->Energy[i] = E_e;
			pt_Cascade_Spectrum->Energy[i] = E_g;
			if(i==(Electron_Table_Size-1)){
				if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "Dirac"){
					pt_Electron_Spectrum->Spectrum[i]=1/(dE*Rate_Inverse_Compton(E_e,z,pt_Spectrum_and_Precision_Parameters));
				}
				if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "Dirac"){
					pt_Cascade_Spectrum->Spectrum[i]=1/(dE*(rate_NPC(E_g,z)+rate_compton(E_g,z)+rate_gg_scattering(E_g,z)));
				}
				else{
					pt_Electron_Spectrum->Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum(E_e,z,pt_Particle_Physics_Model->E_0)/Rate_Inverse_Compton(E_e,z,pt_Spectrum_and_Precision_Parameters);
					pt_Cascade_Spectrum->Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E_g,z,pt_Particle_Physics_Model->E_0)/(rate_NPC(E_g,z)+rate_compton(E_g,z)+rate_gg_scattering(E_g,z));
				}
			}
			else{
				for (int j = i+1 ; j < Electron_Table_Size ; j ++){
					E_j = E_min*pow(E_0/E_min,(double) j/(Electron_Table_Size-1));
					E_j_plus_1 = E_min*pow(E_0/E_min,((double) j+1)/(Electron_Table_Size-1));
					E_j_minus_1 = E_min*pow(E_0/E_min,((double) j-1)/(Electron_Table_Size-1));
					dE_j = (E_j_plus_1 - E_j_minus_1)/2.;
					pt_Electron_Spectrum->Spectrum[i]+=pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum(E_j,z,pt_Particle_Physics_Model->E_0)/Rate_Inverse_Compton(E_j,z,pt_Spectrum_and_Precision_Parameters);
					pt_Cascade_Spectrum->Spectrum[i]+=pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E_j,z,pt_Particle_Physics_Model->E_0)/(rate_NPC(E_j,z)+rate_compton(E_j,z)+rate_gg_scattering(E_j,z));

					pt_Electron_Spectrum->Spectrum[i] += dE_j * pt_Cascade_Spectrum->Spectrum[j] * (dsigma_compton(E_j,z,(E_j+m_e-E_e))+dsigma_NPC(E_j+m_e,z,E_e))/Rate_Inverse_Compton(E_e,z,pt_Spectrum_and_Precision_Parameters);
					pt_Electron_Spectrum->Spectrum[i] += dE_j * pt_Electron_Spectrum->Spectrum[j] * 2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb) * (Function_Integrand_Spectre_Compton(E_j+E_gamma_bb-E_e,E_j, E_gamma_bb)/((E_j)*(E_j)))/Rate_Inverse_Compton(E_e,z,pt_Spectrum_and_Precision_Parameters);
					// cout << "Function_Integrand_Spectre_Compton(E_j+E_gamma_bb-E_e,E_j, E_gamma_bb) " << Function_Integrand_Spectre_Compton(E_j+E_gamma_bb-E_e,E_j, E_gamma_bb) << endl;
					pt_Cascade_Spectrum->Spectrum[i] += dE_j * pt_Cascade_Spectrum->Spectrum[j] * (dsigma_phph(E_j,z,E_g)+dsigma_compton(E_j,z,E_g))/(rate_NPC(E_g,z)+rate_compton(E_g,z)+rate_gg_scattering(E_g,z));
					pt_Cascade_Spectrum->Spectrum[i] += dE_j * pt_Electron_Spectrum->Spectrum[j] * (2*Function_Integrand_Spectre_Compton(E_g,E_j, E_gamma_bb)/(E_j*E_j))*(2*pi*r_e*r_e*m_e*m_e*int_bb/E_gamma_bb)/(rate_NPC(E_g,z)+rate_compton(E_g,z)+rate_gg_scattering(E_g,z));
				}
			}
			// cout << " E = "<< pt_Electron_Spectrum->Energy[i] << " pt_Electron_Spectrum->Spectrum[i] = " << pt_Electron_Spectrum->Spectrum[i] << " pt_Cascade_Spectrum->Spectrum[i] = " << pt_Cascade_Spectrum->Spectrum[i]  <<endl;
			}

}

void  Cascade_Spectrum_Calculation(double z,
																	 Structure_Output_Options * pt_Output_Options,
																	 Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																	 Structure_Spectrum * pt_Cascade_Spectrum,
																	 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){
	double dE, E1, E_e, E_gamma;
	double resultat, integrale;
	double E_c = E_c_0/(1+z);
	Structure_Spectrum Electron_Spectrum;
	Structure_Spectrum Inverse_Compton_Spectrum;
	Structure_Spectrum Diffused_Gamma_Spectrum;
	Structure_Spectrum Diffused_Electron_Spectrum;
	Structure_Spectrum Compton_Electron_Spectrum;
	Structure_Spectrum Tmp_Electron_Spectrum;
	double E_0 = pt_Particle_Physics_Model->E_0;

	/*****Initialization****/
	pt_Cascade_Spectrum->Energy.resize(Gamma_Table_Size);
	pt_Cascade_Spectrum->Spectrum.resize(Gamma_Table_Size);
	pt_Cascade_Spectrum->species = "photon";
	pt_Cascade_Spectrum->spectrum_name = "total_photon_";
	pt_Cascade_Spectrum->redshift = z;

	Inverse_Compton_Spectrum.Energy.resize(Gamma_Table_Size);
	Inverse_Compton_Spectrum.Spectrum.resize(Gamma_Table_Size);
	Inverse_Compton_Spectrum.species = "photon";
	Inverse_Compton_Spectrum.spectrum_name = "ICS_";
	Inverse_Compton_Spectrum.redshift=z;

	Diffused_Gamma_Spectrum.Energy.resize(Gamma_Table_Size);
	Diffused_Gamma_Spectrum.Spectrum.resize(Gamma_Table_Size);
	Diffused_Gamma_Spectrum.species = "photon";
	Diffused_Gamma_Spectrum.spectrum_name = "diffused_photon_";
	Diffused_Gamma_Spectrum.redshift=z;

	Electron_Spectrum.Energy.resize(Electron_Table_Size);
	Electron_Spectrum.Spectrum.resize(Electron_Table_Size);
	Electron_Spectrum.species="electron";
	Electron_Spectrum.spectrum_name = "total_electron_";
	Electron_Spectrum.redshift=z;

	Diffused_Electron_Spectrum.Energy.resize(Electron_Table_Size);
	Diffused_Electron_Spectrum.Spectrum.resize(Electron_Table_Size);
	Diffused_Electron_Spectrum.species = "electron";
	Diffused_Electron_Spectrum.spectrum_name = "diffused_electron_";
	Diffused_Electron_Spectrum.redshift=z;

	Compton_Electron_Spectrum.Energy.resize(Electron_Table_Size);
	Compton_Electron_Spectrum.Spectrum.resize(Electron_Table_Size);
	Compton_Electron_Spectrum.species = "electron";
	Compton_Electron_Spectrum.spectrum_name = "compton_electron_";
	Compton_Electron_Spectrum.redshift=z;

	Tmp_Electron_Spectrum.Energy.resize(Electron_Table_Size);
	Tmp_Electron_Spectrum.Spectrum.resize(Electron_Table_Size);
	Tmp_Electron_Spectrum.species = "electron";
	Tmp_Electron_Spectrum.spectrum_name = "tmp_electron_";
	Tmp_Electron_Spectrum.redshift=z;

	dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Gamma_Table_Size-1);


/*****************Start of the computation****************/

	if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "universal"){
					for(int i=0;i<Gamma_Table_Size;i++){
						E1=E_min+i*dE;
						pt_Cascade_Spectrum->Energy[i]=E1;
						pt_Cascade_Spectrum->Spectrum[i]=universal_spectrum(E1,z,pt_Particle_Physics_Model->E_0);
						Electron_Spectrum.Energy[i]=E1;
						Electron_Spectrum.Spectrum[i]=0;

					}
					check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);

	}
	else{
					if(E_c <= pt_Particle_Physics_Model->E_0){


								for(int i=0;i<Gamma_Table_Size;i++){
									E1=E_min+i*dE;
									pt_Cascade_Spectrum->Energy[i]=E1;
									pt_Cascade_Spectrum->Spectrum[i]=universal_spectrum(E1,z,pt_Particle_Physics_Model->E_0);
									Electron_Spectrum.Energy[i]=E1;
									Electron_Spectrum.Spectrum[i]=0;
								}
								// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);
								if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
										if(pt_Output_Options->verbose>1)cout <<" I will now print the spectrum in files." << endl;
										print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
								}
					}


					else{




						if(pt_Spectrum_and_Precision_Parameters->calculation_mode == "triangular"){
							for(int i=0;i<Electron_Table_Size;i++){
									E1=E_min+i*dE;
									Electron_Spectrum.Energy[i]=E1;
									Electron_Spectrum.Spectrum[i]=0;
									// Electron_Spectrum.Spectrum[i]=Tmp_Electron_Spectrum.Spectrum[i];
							}
							for(int i=0;i<Gamma_Table_Size;i++){
									E1=E_min+i*dE;
									pt_Cascade_Spectrum->Energy[i]=E1;
									pt_Cascade_Spectrum->Spectrum[i]=0;
									// Electron_Spectrum.Spectrum[i]=Tmp_Electron_Spectrum.Spectrum[i];
							}

							Electron_Spectrum.spectrum_name = "total_electron_triangular_";
							pt_Cascade_Spectrum->spectrum_name = "total_photon_triangular_";
							Triangular_Spectrum(pt_Particle_Physics_Model,
																					 pt_Spectrum_and_Precision_Parameters,
																					 pt_Cascade_Spectrum,
																					 &Electron_Spectrum);

						//  check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);
						print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Electron_Spectrum, pt_Particle_Physics_Model);
						print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
						}



						else if(pt_Spectrum_and_Precision_Parameters->calculation_mode == "iterative"){
							/********First step : compute initial ICS spectrum from the electon spectrum injected**********/

							for(int i=0;i<Electron_Table_Size;i++){
									E1=E_min+i*dE;
									Tmp_Electron_Spectrum.Energy[i]=E1;
									Tmp_Electron_Spectrum.Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum(E1,z,pt_Particle_Physics_Model->E_0);
									Electron_Spectrum.Energy[i]=E1;
									Electron_Spectrum.Spectrum[i]=Tmp_Electron_Spectrum.Spectrum[i];
							}

						for(int i=0;i<Gamma_Table_Size;i++){
								E1=E_min+i*dE;
								pt_Cascade_Spectrum->Energy[i]=E1;
								pt_Cascade_Spectrum->Spectrum[i]=0;
								// Electron_Spectrum.Spectrum[i]=Tmp_Electron_Spectrum.Spectrum[i];
						}

						// print_spectrum_automatic_names(0, &Electron_Spectrum, pt_Particle_Physics_Model);

						if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice != "none"){

									for(int i = 0; i<pt_Spectrum_and_Precision_Parameters->number_iterations_electron;i++){
											cout << " iteration electrons : " << i+1 << endl;
											Spectrum_electron_scattered(pt_Particle_Physics_Model,
																								pt_Spectrum_and_Precision_Parameters,
																								&Electron_Spectrum,
																								&Diffused_Electron_Spectrum);
											for(int j=0;j<Electron_Table_Size;j++){
											E1=E_min+j*dE;
											Electron_Spectrum.Energy[j]=E1;
											Electron_Spectrum.Spectrum[j]=Tmp_Electron_Spectrum.Spectrum[j]+Diffused_Electron_Spectrum.Spectrum[j];
											cout <<"E = " << Electron_Spectrum.Energy[j] << " Tmp_Electron_Spectrum = " << Tmp_Electron_Spectrum.Spectrum[j] << " Diffused_Electron_Spectrum = " << Diffused_Electron_Spectrum.Spectrum[j]<<endl;

											}
											print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Electron_Spectrum, pt_Particle_Physics_Model);

									}

									if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
											Diffused_Electron_Spectrum.spectrum_name = "electron_diffused_from_dirac_";
											print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Diffused_Electron_Spectrum, pt_Particle_Physics_Model);
									}
									Spectre_gamma_compton(pt_Particle_Physics_Model,
																				pt_Spectrum_and_Precision_Parameters,
																				&Electron_Spectrum,
																				&Inverse_Compton_Spectrum);
									if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
											Inverse_Compton_Spectrum.spectrum_name = "ICS_from_e_injection_";
											if(pt_Output_Options->verbose>1)cout <<" I will now print the spectrum in files." << endl;
											print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Inverse_Compton_Spectrum, pt_Particle_Physics_Model);
									}



						}
						/*********Second step : compute the gamma spectrum from gg->gg, ge->ge, gN->Nee for a certain number of iterations.*********/

							for(int i=0;i<Gamma_Table_Size;i++){
										pt_Cascade_Spectrum->spectrum_name="Initial_photon_spectrum_";
										E1=E_min+i*dE;
										pt_Cascade_Spectrum->Energy[i]=E1;
										if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice=="none")pt_Cascade_Spectrum->Spectrum[i]=0;
										else pt_Cascade_Spectrum->Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E1,z,pt_Particle_Physics_Model->E_0)/(rate_NPC(E1,z)+rate_compton(E1,z)+rate_gg_scattering(E1,z));
										pt_Cascade_Spectrum->Spectrum[i]+=Inverse_Compton_Spectrum.Spectrum[i];
								}
								// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);

								if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
										pt_Cascade_Spectrum->spectrum_name = "Cascade_";
										if(pt_Output_Options->verbose>1)cout <<" I will now print the spectrum in files." << endl;
										print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
								}



							for(int k = 0; k<pt_Spectrum_and_Precision_Parameters->number_iterations_photon;k++){
										if(pt_Output_Options->verbose>1)cout<<"iteration : " << k+1 << endl;


								Spectrum_gamma_scattered(pt_Particle_Physics_Model,
																				 pt_Spectrum_and_Precision_Parameters,
																				 pt_Cascade_Spectrum,
																				 &Diffused_Gamma_Spectrum);
								 for(int i=0;i<Gamma_Table_Size;i++){
										E_gamma = E_min+i*dE;
										pt_Cascade_Spectrum->Spectrum[i] = Diffused_Gamma_Spectrum.Spectrum[i] + Inverse_Compton_Spectrum.Spectrum[i];
										if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice!="none")pt_Cascade_Spectrum->Spectrum[i]+=pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E_gamma,z,pt_Particle_Physics_Model->E_0)/(rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z));
									}

									// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);

									if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
										pt_Cascade_Spectrum->spectrum_name = "Cascade_";
										print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
									}

						}
						// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);

						/**********Third Step : Compute the associated electron spectrum and the gamma spectrum from ICS.**********/
						if(pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering == "yes"){
							// pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice = "from_cascade";

							Spectre_electron_compton(pt_Particle_Physics_Model,
																			 pt_Spectrum_and_Precision_Parameters,
																			 pt_Cascade_Spectrum,
																			 &Compton_Electron_Spectrum);

						for(int j=0; j<Electron_Table_Size; j++){
							Tmp_Electron_Spectrum.Energy[j]=Compton_Electron_Spectrum.Energy[j];
							Tmp_Electron_Spectrum.Spectrum[j]=Compton_Electron_Spectrum.Spectrum[j];
						}
						for(int i = 0 ; i < pt_Spectrum_and_Precision_Parameters->number_iterations_electron; i++){
							cout << "iteration : " << i+1 << endl;
						 Spectrum_electron_scattered(pt_Particle_Physics_Model,
																				pt_Spectrum_and_Precision_Parameters,
																				&Tmp_Electron_Spectrum,
																				&Diffused_Electron_Spectrum);
						for(int j=0; j<Electron_Table_Size; j++){
							cout <<"E = " << Electron_Spectrum.Energy[j] << " Compton_Electron_Spectrum = " << Compton_Electron_Spectrum.Spectrum[j] << " Diffused_Electron_Spectrum = " << Diffused_Electron_Spectrum.Spectrum[j]<<endl;
							Tmp_Electron_Spectrum.Spectrum[j]= Compton_Electron_Spectrum.Spectrum[j]+Diffused_Electron_Spectrum.Spectrum[j];
							cout << " Electron_Spectrum = " << Electron_Spectrum.Spectrum[j] << endl;
						}
						print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Tmp_Electron_Spectrum, pt_Particle_Physics_Model);

					}
					for(int j=0; j<Electron_Table_Size; j++){
						Electron_Spectrum.Energy[j]=Tmp_Electron_Spectrum.Energy[j];
						Electron_Spectrum.Spectrum[j]+=Tmp_Electron_Spectrum.Spectrum[j];
					}
					print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Electron_Spectrum, pt_Particle_Physics_Model);

						if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "Dirac"){
							pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice = "from_cascade";
						}
							Spectre_gamma_compton(pt_Particle_Physics_Model,
																		pt_Spectrum_and_Precision_Parameters,
																		&Tmp_Electron_Spectrum,
																		&Inverse_Compton_Spectrum);
						if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "from_cascade"){
							pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice = "Dirac";
						}
							for(int l = 0; l<Gamma_Table_Size ; l++){
								cout <<"pt_Cascade_Spectrum->Spectrum[l] = " << pt_Cascade_Spectrum->Spectrum[l] << " Inverse_Compton_Spectrum.Spectrum[l] = " << Inverse_Compton_Spectrum.Spectrum[l] << endl;
								pt_Cascade_Spectrum->Spectrum[l]+=Inverse_Compton_Spectrum.Spectrum[l];
							}

							if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
								Compton_Electron_Spectrum.spectrum_name = "compton_electron_";
								print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Compton_Electron_Spectrum, pt_Particle_Physics_Model);
							}
							if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
								pt_Cascade_Spectrum->spectrum_name = "total_";
								print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
							}

							// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);

							if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
								Inverse_Compton_Spectrum.spectrum_name = "ICS_from_cascade_";
								print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Inverse_Compton_Spectrum, pt_Particle_Physics_Model);
							}
						}
						}

					}
	}



}
