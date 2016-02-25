#include "../include/EM_cascade.h"
#include "../include/injected_spectrum.h"
#include "../include/structures.h"
#include "../include/BBN_constraints.h"
#include "../include/tools.h"

using namespace std;

double  dsigma_compton(double  x, double  z, double g){

	double  dsigma =0;
	dsigma = pi*pow(r_e,2)*m_e/pow(x,2)*(x/g+g/x+pow(m_e/g-m_e/x,2)-2*m_e*(1/g-1/x))*n_e*pow(1+z,3);
	// if(g<x && g>(x/(1+2*x)))dsigma = pi*pow(r_e,2)*m_e/pow(x,2)*(x/g+g/x+pow(m_e/g-m_e/x,2)-2*m_e*(1/g-1/x))*n_e*pow(1+z,3)/(1+Y);
	// cout << "dsigma = " << dsigma << endl;
 	return dsigma;

}
double  dsigma_phph(double  x, double  z,  double g){

	// double  dsigma = pow(T_0*(1+z),6)*8*pow(pi,4)*10000*pow(ALPHA*r_e,2)*pow(m_e,-6)*pow(x,2)*pow(1-g/x+pow(g/x,2),2)/(63*10125*pi);
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
// n_H+(Z^2 = 4)*n_He = (1-Y)n_b+Y*n_He = n_b
	return 2*result;
}


double  rate_compton(double  x, double  z){

	double X = 2*x/m_e;
	double sigma_cs = 2*pi*pow(r_e,2)/X*((1-4./X-8./pow(X,2))*log(1+X)+0.5+8./X+1./(2*pow(1+X,2)));
	double Gamma = sigma_cs*n_e*pow(1+z,3);
	return Gamma;

}
double  rate_NPC(double  x, double  z){

	double	k = x/m_e;
	double  rho = (2*k-4)/(k+2+2*pow(2*k,0.5));
	double  sigma_PCN;
	if(k >= 4){
	sigma_PCN = ALPHA*pow(r_e,2)
								 *(28./9*log(2*k)-218./27
								 +pow(2/k,2)*(2./3*pow(log(2*k),3)-pow(log(2*k),2)+(6-pi*pi/3)*log(2*k)+2*1.20205+pi*pi/6-7./2.)
								 -pow(2/k,4)*(3./16.*log(2*k)+1./8.)
								 -pow(2/k,6)*(29./2304*log(2*k)-77./13824.));}
	else sigma_PCN = ALPHA*pow(r_e,2)*2*pi/3.*pow((k-2)/k,3)*(1+0.5*rho+23./40.*pow(rho,2)+11./60.*pow(rho,3)+29./960.*pow(rho,4));
	double Gamma = sigma_PCN*pow(1+z,3)*eta*n_y_0;			// n_H+(Z^2 = 4)*n_He = (1-Y)n_b+Y*n_He = n_b
	if(Gamma<0)Gamma = 0;
	return Gamma;
}

double  rate_gg_scattering(double  x, double  z){

	double  Gamma;
	// if(x < (m_e*m_e)/(T_0*(1+z)))
	Gamma = 1946./(50625*pi)*pow(ALPHA*r_e,2)*pow(x,3)*8*pow(pi,4)*pow(T_0*(1+z)/m_e,6)/63.;
	// else Gamma = 0;
	// Gamma = 0.1513*pow(ALPHA,4)*m_e*pow(x/m_e,3)*pow(T_0*(1+z)/m_e,6);
	return Gamma;

}
double integrate_dsigma_phph(double E_MIN, double E_MAX, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

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

		F1 = dsigma_phph(E_MAX,z, q1);
		F2 = dsigma_phph(E_MAX,z, q2);
		F3 = dsigma_phph(E_MAX,z, q3);
		F4 = dsigma_phph(E_MAX,z, q4);
		F5 = dsigma_phph(E_MAX,z, q5);
		F6 = dsigma_phph(E_MAX,z, q6);
		F7 = dsigma_phph(E_MAX,z, q7);

		rate+= dq/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
		// cout << "q = " << q1 << " rate = " << rate << endl;
	}
	return rate;
}
double integrate_dsigma_compton(double E_MIN, double E_MAX, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

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

		F1 = dsigma_compton(E_MAX,z, q1);
		if(isnan(F1)==1)F1=0.;
		F2 = dsigma_compton(E_MAX,z, q2);
		F3 = dsigma_compton(E_MAX,z, q3);
		F4 = dsigma_compton(E_MAX,z, q4);
		F5 = dsigma_compton(E_MAX,z, q5);
		F6 = dsigma_compton(E_MAX,z, q6);
		F7 = dsigma_compton(E_MAX,z, q7);
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
																							 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

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
						f[eval]=dsigma_pair_creation(z,E[eval],E_max,pt_Spectrum_and_Precision_Parameters);
						// f[eval]=integrand_rate_pair_creation_v2(E_gamma,E[eval])/(exp(E[eval]/T)-1);
						result += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
					// cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
				}
		}
	return result;
}

double integrand_rate_pair_creation(double s){
	double beta = pow(1-(4*m_e*m_e/s),0.5);
	double result;
	if(beta>=0){
		beta = sqrt(beta);
		result =  s*0.5*pi*r_e*r_e*(1-(beta*beta))*((3-pow(beta,4))*log((1+beta)/(1-beta))-2*beta*(2-beta*beta));
	}
	else {
		beta = 0;
		result = 0;
	}
	// cout << "(integrand_rate_pair_creation : ) result = " << result << " s = " << s << " beta = " << beta << endl;
	return result;
}
double rate_pair_creation(double E_gamma, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){
	double E_gamma_bb = 2.701*T_0*(1+z);
	// if(E_gamma > 10000)E_gamma = 10000.;
	// cout << "E_gamma_bb = " << E_gamma_bb << endl;
	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	double T = T_0*(1+z);
	double resultat = 0;
	double s[pt_Spectrum_and_Precision_Parameters->eval_max-1],f[pt_Spectrum_and_Precision_Parameters->eval_max-1],E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
	double h, ds;
	double f_gamma_bb;
	vector<double> cmb_spectrum_convoluted_with_cross_section_energy;
	vector<double> cmb_spectrum_convoluted_with_cross_section;
	double E_cmb_min = m_e*m_e/E_gamma;
	double E_cmb_max = 10*E_gamma_bb;
	int y;
	for(int j = 0 ; j < pt_Spectrum_and_Precision_Parameters->z_step ; j ++){

				E_gamma_bb = (E_cmb_min)*pow(E_cmb_max/E_cmb_min,(double) j/(pt_Spectrum_and_Precision_Parameters->z_step-1));
				// cout << " E_gamma_bb " << E_gamma_bb << " E_cmb_min " << E_cmb_min << " E_cmb_max " << E_cmb_max << endl;
		 		f_gamma_bb= E_gamma_bb*E_gamma_bb/(pi*pi)/(exp(E_gamma_bb/T)-1);
				ds = (4*E_gamma_bb*E_gamma - 4*m_e*m_e)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
				y = 0;
				while(ds>4*m_e*m_e){
					ds/=10.;
					y++;
				}
				h = ds/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

				// cout << " ds = " << ds << "s_min = " << 4*m_e*m_e << " y " << y << endl;
				// ds = (E*E_gamma_bb/(m_e*m_e) - 1)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);

				for(int i=0; i<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step;i++){

					// cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
					for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
					{
						if(eval == 0){
						if(i==0)	s[eval]=4*m_e*m_e;
						else s[eval]=s[pt_Spectrum_and_Precision_Parameters->eval_max-1];
						}
						else{
							s[eval]=s[0]+eval*h;
						}

					f[eval]=integrand_rate_pair_creation(s[eval]);
					resultat += ds/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
					// cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
				}
			}
				// cout << "(function rate pair creation : ) resultat intermediaire = " << resultat << endl;

				resultat*=f_gamma_bb/(E_gamma_bb*E_gamma_bb);
				cmb_spectrum_convoluted_with_cross_section_energy.push_back(E_gamma_bb);
				cmb_spectrum_convoluted_with_cross_section.push_back(resultat);
	}
  ds = (E_cmb_max-E_cmb_min)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
	y = 0;
	while(ds>E_cmb_min){
		ds/=10.;
		y++;
	}
	h = ds/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
	// cout << " dS " << ds << " E_cmb_min " << E_cmb_min << endl;
	resultat=0;
	for(int i=0; i<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;i++){

				// cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
				for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
				{
					if(eval == 0){
					if(i==0)	s[eval]=E_cmb_min;
					else s[eval]=s[pt_Spectrum_and_Precision_Parameters->eval_max-1];
					}
					else{
						s[eval]=s[0]+eval*h;
					}

					linearint(cmb_spectrum_convoluted_with_cross_section_energy, cmb_spectrum_convoluted_with_cross_section, cmb_spectrum_convoluted_with_cross_section_energy.size(), s[eval], f[eval]);
					resultat += ds/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
				// cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
			}
	}

	// cout << "(function rate pair creation : ) resultat = " << 1/(8*E*E)*resultat<< " energy = " << E << endl;
	return 1/(8*E_gamma*E_gamma)*resultat;
	// return 1/(8*E*E)*int_bb/pow(E_gamma_bb,2)*resultat;
}

double integrand_rate_pair_creation_v2(double E_gamma, double E_gamma_bar){
	double s_0;
	double beta_0, beta_0_squared;
	double w_0;
	double L_0, L_prime_0;
	double result;
	int N = 100;
	s_0 = E_gamma*E_gamma_bar/(m_e*m_e);
	beta_0_squared = 1-1/s_0;
if(beta_0_squared>=0){
	beta_0 = pow(beta_0_squared,0.5);
	w_0 = (1+beta_0)/(1-beta_0);
	L_prime_0 = pi*pi/12;
	for(int i = 1; i <= N ; i++){
		L_prime_0 -= pow(-1,i-1)*pow(i,-2)*pow(w_0,-i);
	}
	L_0 = 0.5*pow(log(w_0),2)+L_prime_0;
	result = (1+beta_0_squared)/(1-beta_0_squared)*log(w_0)-beta_0_squared*log(w_0)-pow(log(w_0),2)-4*beta_0/(1-beta_0_squared)
	+2*beta_0+4*log(w_0)*log(w_0+1)-L_0;
}
else result = 0;

	return result;

}
double integrand_rate_pair_creation_v3(double x_gamma, double x_gamma_bar, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

	double v;
	double w;
	double result;

	v = x_gamma*x_gamma_bar-1;

if(v>=0){
	w = (pow(1+v,0.5)+pow(v,0.5))/(pow(1+v,0.5)-pow(v,0.5));
	result =(1+2*v+2*v*v)/(1+v)*log(w)-(2*pow(v,0.5)*(1+2*v))/pow(1+v,0.5)-pow(log(w),2)+2*pow(log(1+w),2)+4*polylog_2(1/(1+w),pt_Spectrum_and_Precision_Parameters)-pi*pi/3;
}
else result = 0;
	// cout << " v " << v << " w " << w << " result " << result << endl;
	return result;

}
double rate_pair_creation_v2(double E_gamma, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){
	double h, dE, E[pt_Spectrum_and_Precision_Parameters->eval_max],f[pt_Spectrum_and_Precision_Parameters->eval_max], T=T_0*(1+z), result;
	double E_gamma_bb = 2.701*T_0*(1+z);
	double E_cmb_min = m_e*m_e/E_gamma;
	double E_cmb_max = 10*E_gamma_bb;
	double A = 8*pi*pow(m_e/(2*pi),3);
	int y;

	dE = (E_cmb_max-E_cmb_min)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
	y = 0;
	while(dE>E_cmb_min){
		dE/=10.;
		y++;
	}
	h = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
	// cout << "(rate_pair_creation_v2 :) dE " << dE << " E_cmb_min " << E_cmb_min << " y " << y << endl;
	result=0;
	for(int i=0; i<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;i++){

				// cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
				for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
				{
					if(eval == 0){
					if(i==0)	E[eval]=E_cmb_min;
					else E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
					}
					else{
						E[eval]=E[0]+eval*h;
					}
					f[eval]=integrand_rate_pair_creation_v3(E_gamma/m_e,E[eval]/m_e,pt_Spectrum_and_Precision_Parameters)/(exp(E[eval]/T)-1);
					// f[eval]=integrand_rate_pair_creation_v2(E_gamma,E[eval])/(exp(E[eval]/T)-1);
					result += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
				// cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
			}
	}
	if(result < 0 ) result = 0;
	// cout << "(rate_pair_creation_v2 : ) resultat = " << pow(r_e,2)*pow(m_e,4)/(pi*E_gamma*E_gamma)*result<< " energy = " << E_gamma << endl;
	// return pow(r_e,2)*pow(m_e,4)/(pi*E_gamma*E_gamma)*result;
	return 3*sigma_T*A*m_e/(8*E_gamma*E_gamma)*result;
	// return 3*sigma_T/(E_gamma)*result;
}

double Function_Integrand_Spectre_Compton(double E_e, double E_gamma,  double E_gamma_bar){
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
	// cout << " q = " << q << " Gamma_e = " << Gamma_e << " f = " << f << " E_gamma_bar =" << E_gamma_bar<<endl;
	return f/E_gamma_bar;
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
double Function_Integrand_Spectre_Compton_version_q(double q, double E_e, double E_gamma_bar){
	double f;
	double Gamma_e = 4*E_gamma_bar*E_e/(m_e*m_e);
	double E_gamma=Gamma_e*E_e*q/(1+Gamma_e*q);

	if(q>=0. && q<=1.) {
		f = (2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q))*(E_e*Gamma_e/pow(1+Gamma_e*q,2));
		// f = (2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q))*pow(E_e-E_gamma,2)*Gamma_e/E_e;
}
	else f=0;
	if(f<=0)f=0;
	// cout <<"q = "<< q << " Gamma_e = "<<Gamma_e<<" E_gamma = "<<E_gamma<<" E_e ="<<E_e<<" f = " << f << endl;

	return f;
}

double Analytical_form_inverse_compton(double E_e, double E_gamma_bar, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){
	double gamma_e = E_e/m_e;
	E_gamma_bar /= m_e;
	double s = 4*gamma_e*E_gamma_bar;
	double resultat;

	resultat = (s+9+8./s)*log(1+s)-(16+18*s+s*s)/(2*(1+s))+4*polylog_2(-s,pt_Spectrum_and_Precision_Parameters);

		return resultat/(s*s);
}
/* void  Spectre_electron_pair_creation(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
													 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
													 Structure_Spectrum * pt_Gamma_Spectrum,
													 Structure_Spectrum * pt_Electron_Spectrum){
			double E_e;
			double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7;
			double dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);
		 	for(int j =0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;j++){
		 		resultat=0;
		 		Gamma_electron=0;
		 		E_e = pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE;
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

		 			f1*=(function_integrand_electrom_pair_creation(E_e,E1,E_gamma_bb))/pow(E1,3);
		 			f2*=(function_integrand_electrom_pair_creation(E_e,E2,E_gamma_bb))/pow(E2,3);
		 			f3*=(function_integrand_electrom_pair_creation(E_e,E3,E_gamma_bb))/pow(E3,3);
		 			f4*=(function_integrand_electrom_pair_creation(E_e,E4,E_gamma_bb))/pow(E4,3);
		 			f5*=(function_integrand_electrom_pair_creation(E_e,E5,E_gamma_bb))/pow(E5,3);
		 			f6*=(function_integrand_electrom_pair_creation(E_e,E6,E_gamma_bb))/pow(E6,3);
		 			f7*=(function_integrand_electrom_pair_creation(E_e,E7,E_gamma_bb))/pow(E7,3);
		 			resultat += dE_2/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);

		 		}
				pt_Electron_Spectrum->s
			}

}*/

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
	dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);
	for(int j =0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;j++){
		resultat=0;
		Gamma_electron=0;
		E_e = pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE;
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
	double gamma = E_e/m_e;
	double q_min=(1/(4*gamma*gamma)), q1, q2, q3, q4, q5, q6, q7, dq, h;
	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
	dq = (1-q_min)/ (double) (10*n_step-1);
	h = dq/6.;
	q1=q_min;
	for(int i=0;i<10*n_step-1;i++){
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
		 	dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);



		 	for(int j =0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;j++){
		 		resultat=0;
		 		Gamma_electron=0;
				E_e = pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE;
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
	dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);

		for(int j =0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;j++){
			resultat=0;
			E_gamma = pt_Spectrum_and_Precision_Parameters->E_min_table + j*dE;
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
				// cout << "E = " << pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE<<"F1 = "<< F1 << " f1 apres = " << f1 << endl;

				resultat+= dE_2/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);

			}
			pt_Gamma_Spectrum->Energy[j]=E_gamma ;
			pt_Gamma_Spectrum->Spectrum[j]=resultat*2*pi*r_e*r_e*m_e*m_e*int_bb/E_gamma_bb/(rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z));
			if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice=="Dirac"){
				Gamma_electron = Rate_Inverse_Compton(E_0,z,pt_Spectrum_and_Precision_Parameters);
				pt_Gamma_Spectrum->Spectrum[j]+=2*pi*r_e*r_e*m_e*m_e*2*int_bb*Function_Integrand_Spectre_Compton(E_gamma,E_0,E_gamma_bb)/(E_gamma_bb*Gamma_electron*E_0*E_0)/(rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z));
			}
			// cout << "E = "  << pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE<< "resultat = " << 	pt_Gamma_Spectrum->Spectrum[j] << endl;

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
 dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);

	for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
		E_gamma = pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
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
		// if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="iterative")  os << "Output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << pt_Spectrum_and_Precision_Parameters->number_iterations_electron <<"iterations.dat";
		// else if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="triangular")  os << "Output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << "triangular.dat";
		if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_file_name=="automatic")  os << "Output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << "triangular.dat";
		else os << pt_Spectrum_and_Precision_Parameters->electron_spectrum_file_name;
		// cout  << pt_Spectrum_and_Precision_Parameters->electron_spectrum_file_name;
	}
	  name = os.str();
    ifstream file(name);
    if(file)cout << "Importing file " << name << " in structure " << pt_Spectrum->spectrum_name << "." <<endl;
		else{
			cout << "I couldn't recognize cascade spectrum file. Please check that it is present in the folder Cascade_Specrum_Folder with proper name : 'Spectrum_mXXX_zXXX_XXXiterations.dat' corresponding to the value of m, z and iterations you are using."<<endl;
			return;
		}
    while(file){
			string line;
	    getline(file, line);
	    // stringstream is(ligne);
	    // if(line[0] == '#' or line[0] == '\0') continue;
      file >> tmp_Energy >> tmp_Spectrum ;
			cout << z << "  " << tmp_Energy << "  " << tmp_Spectrum<<endl;

			pt_Spectrum->Energy.push_back(tmp_Energy);
			pt_Spectrum->Spectrum.push_back(tmp_Spectrum);
			// pt_Spectrum->Spectrum.push_back(tmp_Spectrum*(rate_NPC(tmp_Energy,z)+rate_compton(tmp_Energy,z)+rate_gg_scattering(tmp_Energy,z)));

		}

		file.close();


}
double integrator_simpson_gamma_inverse_compton(double (*func)(double,double,double),
																													 double z,
																													 double E_ini,
																													 double E_max,
																													 double E_e,
																													 double E_gamma,
																													 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

	double T=T_0*(1+z);
	// double T=T_0*(1+z), E_max = 10*2.701*T;

	double dE = (E_max- E_ini)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
	int y=0;
	while(dE>E_ini){
		dE/=10.;
		y++;
	}
	// cout << "E_max "<< E_max << " E_ini " << E_ini <<" dE = " << dE << " y " << y <<endl;
	double E[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h2;
	double resultat = 0, res_initial=0;
	double f_initial=0, precision;
	precision = 1e-4;
								for(int j=0; j<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;j++){

									h2 = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
									// cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
									for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
									{
										if(eval == 0){
										if(j==0)	E[eval]=E_ini;
										else E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
										}
										else{
											E[eval]=E[0]+eval*h2;
										}

									if((exp(E[eval]/T)-1)!=0.)f[eval]=(*func)(E_e,E_gamma,E[eval])*E[eval]*E[eval]/(exp(E[eval]/T)-1);
									else f[eval] = 0;
									resultat += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
									// cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
								}

								// if(res_initial==0 && resultat !=0)res_initial = resultat;
								// if(res_initial!=0 && resultat/res_initial<precision)break;
								// cout << E1 << " f1 " << f1 << " (exp(E1/T)-1) "  << (exp(E1/T)-1) << " f7 " << f7 << endl;


								// 	cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << endl;
							}

							// if(resultat != 0.)cout << "resultat = " << resultat << endl;
 							return resultat/(pi*pi);
}
double gamma_inverse_compton_analytical_v2(double gamma_e, double E_gamma, double z){
	double T = T_0*(1+z);
	double eta_c = T*E_gamma/(m_e*m_e);
	double eta_0 = E_gamma*E_gamma/(4*gamma_e*m_e*(gamma_e*m_e-E_gamma));
	// cout << "eta_c " << eta_c << "eta_0 " << eta_0 << endl;

	double I = pi*pi/6*eta_c*(exp(-5./4.*pow(eta_0/eta_c,0.5))+2*eta_0*exp(-5./7*pow(eta_0/eta_c,0.7)))*exp(-2./3.*eta_0/eta_c);
	// return sigma_T*T*T*E_gamma/(8*gamma_e*gamma_e);
	return 3*sigma_T*T*m_e*m_e/(4*pi*pi*gamma_e*gamma_e)*I;
}

double gamma_inverse_compton_analytical(double gamma_e,
																				double E_gamma,
																				double z,
																				int N,
																				Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){
	double A = 8*pi*pow(m_e/(2*pi),3);
	double T = T_0*(1+z);
	E_gamma /= m_e;
	double E_n = E_gamma/(4*gamma_e*(gamma_e-E_gamma));
	double theta = T/m_e;
	double x = E_n/theta;
	double q = 2.257;
	double f_0,f_1,f_minus_1,f_ln;
	double B[22];
	double result;
	B[0] = 1, B[1]= -0.5, B[2] = 1./6., B[4] = -1./30, B[6] = 1./42, B[8] = -1./30, B[10] = 5./66, B[12] = -691./2730, B[14] = 7./6, B[16] = -3617/510;
	f_0 = -log(1-exp(-x));
	f_1 = polylog_2(exp(-x),pt_Spectrum_and_Precision_Parameters)+x*f_0;

  if(x <= q)	f_minus_1 = 0.0366377 + pow(x,-1) + log(x)/2 -pow(q,-1)-log(q)/2;
	else f_minus_1 = 0;
	if(x<= q)	f_ln =0.5*(q*(1-log(q))+pow(log(q),2)) - 0.5*(x*(1-log(x))+pow(log(x),2)) + 0.125505;
	else f_ln = 0;
	for(int i = 1; i <= N ; i++){
		if(x <= q) f_minus_1 +=pow(q,2*i-1)*B[2*i]/((2*i-1)*factorial(2*i)) - pow(x,2*i-1)*B[2*i]/((2*i-1)*factorial(2*i));
		else f_minus_1 += expint(1,i*x);
		// else f_minus_1 += exponential_integral(i*x,pt_Spectrum_and_Precision_Parameters);
		if(x <= q) f_ln = pow(q,2*i)*(2*i*log(q)-1)*B[2*i]/(4*i*i*factorial(2*i))-pow(x,2*i)*(2*i*log(x)-1)*B[2*i]/(4*i*i*factorial(2*i));
		else f_ln = (expint(1,i*x)+exp(-i*x)*log(x))/i;
		// else f_ln = (exponential_integral(i*x,pt_Spectrum_and_Precision_Parameters)+exp(-i*x)*log(x))/i;
	}
	// double I = pi*pi/6*eta_c*(exp(-5./4.*pow(eta_0/eta_c,0.5))+2*eta_0*exp(-5./7*pow(eta_0/eta_c,0.7)))*exp(-2./3.*eta_0/eta_c);
	// return sigma_T*T*T*E_gamma/(8*gamma_e*gamma_e);
	// cout << " x = " << x << endl;
	// cout << "f1 = " << f_1 << " f_0 = " << f_0 << " f_ln = " << f_ln << " f_minus_1 = " << f_minus_1 << endl;
	result = 3*sigma_T*A/(4*gamma_e*gamma_e)*((1-2*E_gamma*E_n+2*log(E_n/theta))*E_n*theta*f_0+(1+2*E_gamma*E_n)*theta*theta*f_1-2*E_n*E_n*f_minus_1-2*E_n*E_n*f_minus_1-2*E_n*theta*f_ln);
	if(result < 0) result = 0;
	return result;
}



double integrator_simpson_scattered_electron_inverse_compton(double z,
																													 double E_e_prime,
																													 double E_e,
																													 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

	double T=T_0*(1+z);
	double gamma_e, E_cmb_min, E_cmb_max, E_gamma_bb;
	E_gamma_bb = 2.701*T;
	gamma_e = E_e/m_e;
	E_cmb_min = (E_e/E_e_prime-1)/(4*gamma_e)*m_e;
	E_cmb_max = 10*E_gamma_bb;

	double dE = (E_cmb_max- E_cmb_min)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
	int y=0;
	while(dE>E_cmb_min){
		dE/=10.;
		y++;
	}
	// cout << "E_max "<< E_cmb_max << " E_ini " << E_cmb_min <<" dE = " << dE << " y = " <<y << endl;
	double E[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h2;
	double resultat = 0, res_initial=0;
	double f_initial = 0, precision;
	precision = 1e-4;

								for(int j=0; j<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;j++){

									h2 = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
									// cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
									for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
									{
										if(eval == 0){
										if(j==0)	E[eval]=E_cmb_min;
										else E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
										}
										else{
											E[eval]=E[0]+eval*h2;
										}
									if((exp(E[eval]/T)-1)!=0.)f[eval]=Function_Integrand_Spectre_Compton(E_e,E_e+E[eval]-E_e_prime,E[eval])*E[eval]*E[eval]/(exp(E[eval]/T)-1);
									else f[eval] = 0;
									resultat += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
									// cout <<"E_e "<<E_e<< " E_e+E[eval]-E_e_prime " << E_e+E[eval]-E_e_prime << " eval " << eval << " E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
								}

								// if(res_initial==0 && resultat !=0)res_initial = resultat;
								// if(res_initial!=0 && resultat/res_initial<precision)break;
								// cout << E1 << " f1 " << f1 << " (exp(E1/T)-1) "  << (exp(E1/T)-1) << " f7 " << f7 << endl;


								// 	cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << endl;
		 					}
							// if(resultat != 0.)cout << "resultat = " << resultat << endl;
 							return resultat/(pi*pi);
}
double Analytical_form_scattered_electron_from_inverse_compton(double z, double gamma_e, double gamma_prime, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){
	double T = T_0*(1+z);
	double A = 8*pi*pow(m_e/(2*pi),3);
	double a, b;
	double r = 0.5*(gamma_e/gamma_prime + gamma_prime/gamma_e);
	int N = 3;
	double theta = T/m_e;
	double k = (gamma_e/gamma_prime-1)/(4*gamma_e*theta);
	double Eps_0 = 0, Eps_1 = 0;
	double result, LI2 = 0;
	for(int i = 1; i <= N ; i++){
		// Eps_0 += exponential_integral(i*k,pt_Spectrum_and_Precision_Parameters);
		// Eps_1 += exponential_integral(i*k,pt_Spectrum_and_Precision_Parameters)/i;
		Eps_0 += expint(1,i*k);
		Eps_1 += expint(1,i*k)/i;
	}
	// for(int i = 1 ; i < 200*N; i++){
	// 	LI2 += pow(exp(-k),i)/(i*i);
	// }
	//
	// a = r*LI2;
	a = r*polylog_2(exp(-k),pt_Spectrum_and_Precision_Parameters);
	b = 2*k*(log(1-exp(-k))+k*Eps_0+Eps_1);
	result = 3*sigma_T*A*theta*theta/(4*gamma_e*gamma_e)*(a-b);
	if(result < 0) result = 0;
	// cout << " exp(-k) "<< exp(-k) << " Eps_0 " << Eps_0 << " Eps_1 " << Eps_1  << " result " << result << " a/r = " << a/r <<" LI2 = " << LI2<< " b = " << b <<endl;
	return result;
}

double dsigma_pair_creation(double z,
													 double E_e,
													 double E_gamma,
													 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

	double gamma_e, gamma_prime, x_j;
  double E_gamma_bb,E_cmb_min, E_cmb_max, E_s;
	double E[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h2;
	double T=T_0*(1+z);
	double dE;
	double resultat = 0;
	double res_initial =0 , precision;
	int y=0;

	E_gamma_bb = 2.701*T_0*(1+z);
	gamma_e = E_gamma/m_e;
	gamma_prime = E_e/m_e;
  E_cmb_max = 10*E_gamma_bb;
	x_j = E_gamma/m_e;
	gamma_e = E_e/m_e;
	gamma_prime = x_j - gamma_e;
	E_s = x_j*x_j/(4*gamma_e*gamma_prime);
	E_cmb_min = E_s/x_j*m_e;
	while(E_cmb_min>E_cmb_max)E_cmb_max*=10.;
	dE = (E_cmb_max- E_cmb_min)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
	// cout << "E_cmb_max "<< E_cmb_max << " E_cmb_min " << E_cmb_min <<" dE = " << dE << endl;

	if(E_s > 1 ){

			while(dE>E_cmb_min){
				dE/=10.;
				y++;
			}
			// cout << "(dsigma_pair_creation :) after y = "<<y <<" dE = " << dE << endl;
			h2 = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

			// cout << "E_max "<< E_cmb_max << " E_ini " << E_cmb_min <<" dE = " << dE << endl;
			// double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7, h2;

			precision = 1e-4;
				 					for(int j=0; j<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;j++){
										// cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
										for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
										{
											if(eval == 0){
											if(j==0)	E[eval]=E_cmb_min;
											else E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
											}
											else{
												E[eval]=E[0]+eval*h2;
											}

										if((exp(E[eval]/T)-1)!=0.)f[eval]=function_integrand_pair_creation(E_e,E_gamma,E[eval])/(exp(E[eval]/T)-1);
										else f[eval] = 0;
										resultat += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
										// cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
									}


										if(res_initial==0 && resultat !=0)res_initial = resultat;
										if(res_initial!=0 && resultat/res_initial<precision)break;

				 					}
	}
	else resultat = 0;

							// cout << "resultat = " << resultat << endl;
 							return resultat/(pi*pi);
}



double integrator_simpson_blackbody_spectrum_over_scattered_inverse_compton(double (*func)(double,double,double),
																													 double z,
																													 double E_ini,
																													 double E_max,
																													 double E_e,
																													 double E_e_prime,
																													 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

	double T=T_0*(1+z);
	// double T=T_0*(1+z), E_max = 10*2.701*T;

	double dE = (E_max- E_ini)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
	// cout << "E_max "<< E_max << " E_ini " << E_ini <<" dE = " << dE << endl;
	int y=0;
	while(dE>E_ini){
		dE/=10.;
		y++;
	}
	// cout << "after y = "<<y <<" dE = " << dE << endl;
	double E[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h2;
	double resultat = 0;
		 					for(int j=0; j<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;j++){

								h2 = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
								// cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
								for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
								{
									if(eval == 0){
									if(j==0)	E[eval]=E_ini;
									else E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
									}
									else{
										E[eval]=E[0]+eval*h2;
									}

								if((exp(E[eval]/T)-1)!=0.)f[eval]=(*func)(E_e,E_e+E[eval]-E_e_prime,E[eval])*E[eval]*E[eval]/(exp(E[eval]/T)-1);
								else f[eval] = 0;
								resultat += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
								// cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
							}
								// 	cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << endl;
		 					}
							// if(resultat != 0)cout << "resultat = " << resultat << endl;
 							return 2*pi*r_e*r_e*m_e*m_e*resultat/(pi*pi);
}

double integrator_simpson_rate_inverse_compton(double (*func)(double,double),
																																			 double z,
																																			 double E_ini,
																																			 double E_max,
																																			 double E_e,
																																			 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

	double T=T_0*(1+z);
	double dE = (E_max- E_ini)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
	double E[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h2;

	int y=0;
	while(dE>E_ini){
		dE/=10.;
		y++;
	}
	// cout << "E_max "<< E_max << " E_ini " << E_ini <<" dE = " << dE << endl;
	double resultat = 0;
	// cout << "here "<< dE << endl;

								for(int j=0; j<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;j++){

									h2 = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
									// cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
									for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
									{
										if(eval == 0){
										if(j==0)	E[eval]=E_ini;
										else E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
										}
										else{
											E[eval]=E[0]+eval*h2;
										}

									if((exp(E[eval]/T)-1)!=0.)f[eval]=(*func)(E_e,E[eval])*E[eval]/(exp(E[eval]/T)-1);
									else f[eval] = 0;
									resultat += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
									// cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
								}

								// if(res_initial==0 && resultat !=0)res_initial = resultat;
								// if(res_initial!=0 && resultat/res_initial<precision)break;
								// cout << E1 << " f1 " << f1 << " (exp(E1/T)-1) "  << (exp(E1/T)-1) << " f7 " << f7 << endl;


								// 	cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << endl;
							}

							resultat*=2*pi*r_e*r_e*m_e*m_e/(E_e*E_e*pi*pi);
							// if(resultat != 0)cout << "(integrator_simpson_blackbody_spectrum_over_rate_inverse_compton : ) resultat = " << resultat << endl;
 							return resultat;
}
double integrator_simpson_rate_inverse_compton(double z,
																									 double E_ini,
																									 double E_max,
																									 double E_e,
																									 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

	double T=T_0*(1+z);
	double A = 8*pi*pow(m_e/(2*pi),3);
	double dE = (E_max- E_ini)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1)/m_e;
	double theta = T/m_e;
	double E[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h2;

	int y=0;
	while(dE>E_ini){
		dE/=10.;
		y++;
	}
	// cout << "(integrator_simpson_rate_inverse_compton :) E_max "<< E_max << " E_ini " << E_ini <<" dE = " << dE << " y = " << y << endl;
	h2 = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
	double resultat = 0;
	// cout << "here "<< dE << endl;

								for(int j=0; j<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;j++){

									// cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
									for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
									{
										if(eval == 0){
										if(j==0)	E[eval]=E_ini/m_e;
										else E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
										}
										else{
											E[eval]=E[0]+eval*h2;
										}

									if((exp(E[eval]/theta)-1)!=0.)f[eval]=Analytical_form_inverse_compton(E_e,E[eval],pt_Spectrum_and_Precision_Parameters)*E[eval]*E[eval]/(exp(E[eval]/theta)-1);
									else f[eval] = 0;
									resultat += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
									// cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
								}

								// if(res_initial==0 && resultat !=0)res_initial = resultat;
								// if(res_initial!=0 && resultat/res_initial<precision)break;
								// cout << E1 << " f1 " << f1 << " (exp(E1/T)-1) "  << (exp(E1/T)-1) << " f7 " << f7 << endl;


								// 	cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << endl;
							}

							resultat*=3*sigma_T/2.*A;
							// resultat*=2*pi*r_e*r_e*m_e*m_e/(E_e*E_e*pi*pi);
							// if(resultat != 0)cout << "(integrator_simpson_blackbody_spectrum_over_rate_inverse_compton : ) resultat = " << resultat << endl;
 							return resultat;
}
double integrator_simpson_rate_inverse_compton_full(double z,
																										double E_e,
																										double E_max,
																										Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

	double T=T_0*(1+z);
	double dE, E_ini, E_gamma;
	double E_gamma_bb = 2.701*T;
	double E[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],q_tab[pt_Spectrum_and_Precision_Parameters->eval_max],dq,h;
	vector<double> cmb_spectrum_convoluted_with_cross_section_energy;
	vector<double> cmb_spectrum_convoluted_with_cross_section;
	double E_cmb_min = E_gamma_bb/100.;
	double E_cmb_max = 1e2*E_gamma_bb;
	double q, q_min = m_e*m_e/(4*E_e*E_e);
	int y;
	// cout << "E_max "<< E_max << " E_ini " << E_ini <<" dE = " << dE << endl;
	double tmp_res = 0, resultat = 0;
	dE = (E_cmb_max - E_cmb_min)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
	y = 0;
	while(dE>E_cmb_min){
		dE/=10.;
		y++;
		// cout << "dE " << dE << " E_cmb_min " << E_cmb_min << endl;
	}
	h = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

							for(int j = 0 ; j < pt_Spectrum_and_Precision_Parameters->z_step ; j ++){

										q = (q_min)*pow(1/q_min,(double) j/(pt_Spectrum_and_Precision_Parameters->z_step-1));
										// cout << " E_gamma_bb " << E_gamma_bb << " E_cmb_min " << E_cmb_min << " E_cmb_max " << E_cmb_max << endl;

										// cout << " q = " << q << "q_min = " << q_min << " y " << y << endl;

										// dE = (E*E_gamma_bb/(m_e*m_e) - 1)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);

										for(int i=0; i<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step;i++){

											// cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
											for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
											{
												if(eval == 0){
												if(i==0)	E[eval]=E_cmb_min;
												else E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
												}
												else{
													E[eval]=E[0]+eval*h;
												}

											f[eval]=Function_Integrand_Spectre_Compton_version_q(q, E_e,E[eval])*E[eval]*E[eval]/(pi*pi)/(exp(E[eval]/T)-1);
											tmp_res = dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
											resultat += tmp_res;
											// cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat <<" (exp(E[eval]/T)-1) " << (exp(E[eval]/T)-1) << endl;
										}
										if(tmp_res == 0) break;
									}
										// cout << "(function rate electron : ) resultat intermediaire = " << tmp_res << endl;

										cmb_spectrum_convoluted_with_cross_section_energy.push_back(q);
										cmb_spectrum_convoluted_with_cross_section.push_back(resultat);
							}


							dq = (1 - q_min)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
							y = 0;
							while(dq>1000*q_min){
								dq/=10.;
								y++;
								cout << "dq " << dq << " q_min " << q_min << endl;

							}
							// h = dq/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

								for(int j=0; j<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;j++){


									// cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
									for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
									{
										if(eval == 0){
										if(j==0)	q_tab[eval]=q_min;
										else q_tab[eval]=q_tab[pt_Spectrum_and_Precision_Parameters->eval_max-1];
										}
										else{
											q_tab[eval]=q_tab[0]+eval*h;
										}


									linearint(cmb_spectrum_convoluted_with_cross_section_energy, cmb_spectrum_convoluted_with_cross_section, cmb_spectrum_convoluted_with_cross_section_energy.size(), q_tab[eval], f[eval]);
									resultat += dq/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
									// cout << "eval " << eval << "q = " << q_tab[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
								}

								// if(res_initial==0 && resultat !=0)res_initial = resultat;
								// if(res_initial!=0 && resultat/res_initial<precision)break;
								// cout << E1 << " f1 " << f1 << " (exp(E1/T)-1) "  << (exp(E1/T)-1) << " f7 " << f7 << endl;


								// 	cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << endl;
							}

							resultat*=2*pi*r_e*r_e*m_e*m_e/(E_e*E_e);
							// if(resultat != 0)cout << "(integrator_simpson_blackbody_spectrum_over_rate_inverse_compton : ) resultat = " << resultat << endl;
 							return resultat;
}
double print_interaction_rate(double z,
															double E_MIN,
															double E_MAX,
															Structure_Output_Options * pt_Output_Options,
													  	Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

																ostringstream os;
																string name;
																if(pt_Output_Options->interaction_rate_files=="automatic"){
																	os << "Output/Interaction_Rate_Folder/Interaction_Rate_at_z"<< z <<".dat";
																}
																else os << pt_Output_Options->interaction_rate_files << ".dat";
																name = os.str();

																ofstream file(name);
																cout << "Printing in file " << name <<"."<<  endl;


																	double E_e, E_g;
																	double E_c = E_c_0/(1+z);
																	double E_phph = m_e*m_e/(T_0*(1+z));
																	double E_gamma_bb = 2.701*T_0*(1+z);
																	double Rate_photons_E_g = 0, Rate_electrons_E_e = 0;
																	double rate_PP= 0, rate_COM= 0, rate_DP= 0, rate_DP_2= 0, rate_NP = 0;
																	double E_cmb_min = E_gamma_bb/100., E_cmb_max = E_gamma_bb*100.;
																	double PP = 0, CS= 0, ICS_g = 0, ICS_e= 0, NPC= 0, COM= 0, DP= 0;
																			for(double i = (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1); i>=0 ;i--){
																			Rate_photons_E_g = 0;
																			E_g = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_MAX/pt_Spectrum_and_Precision_Parameters->E_min_table,(double) i/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
																			if(pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering == "yes")Rate_electrons_E_e = integrator_simpson_rate_inverse_compton(z,E_cmb_min,E_cmb_max,E_g,pt_Spectrum_and_Precision_Parameters);
																			else Rate_electrons_E_e = 0.;
																			if(pt_Spectrum_and_Precision_Parameters->pair_creation_in_nuclei =="yes")rate_NP= rate_NPC(E_g,z);
																			else rate_NP = 0.;
																			if(pt_Spectrum_and_Precision_Parameters->compton_scattering =="yes")rate_COM = rate_compton(E_g,z);
																			else rate_COM = 0.;
																			if(pt_Spectrum_and_Precision_Parameters->photon_photon_diffusion =="yes" && E_g < E_phph) rate_PP=rate_gg_scattering(E_g,z);
																			else rate_PP = 0.;
																			if(pt_Spectrum_and_Precision_Parameters->double_photon_pair_creation=="yes" && E_g > E_c/10.) rate_DP=rate_pair_creation_v2(E_g,z,pt_Spectrum_and_Precision_Parameters);
																			else rate_DP = 0.;

																			Rate_photons_E_g += rate_NP;
																			Rate_photons_E_g += rate_COM;
																			Rate_photons_E_g+=rate_PP;
																			Rate_photons_E_g+=rate_DP;

																			file << E_g << " " << rate_NP << "  " << rate_COM << " " << rate_PP << " " << rate_DP << "  " << Rate_photons_E_g << " " << Rate_electrons_E_e << endl;
																			if(i == (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1))cout << "E_g" << " " << "rate_NP" << "  " << "rate_COM" << " " << "rate_PP" << " " << "rate_DP" << "  " << "Rate_photons_E_g" << " " << "Rate_electrons_E_e" << endl;
																			cout << E_g << " " << rate_NP << "  " << rate_COM << " " << rate_PP << " " << rate_DP << "  " <<  Rate_photons_E_g << " " << Rate_electrons_E_e << endl;
																			// file << 100/E_g << " " << 2*dsigma_pair_creation(z,100,E_g,pt_Spectrum_and_Precision_Parameters)/rate_DP*E_g/m_e  << endl;
																			// cout << 100/E_g << " "  << endl;
																			// cout << "check : " << rate_PP/integrate_dsigma_phph(0,E_g,z,pt_Spectrum_and_Precision_Parameters) << "  " << rate_PP/integrator_simpson_dsigma_pair_creation(z,100,E_g,pt_Spectrum_and_Precision_Parameters) << endl;
													}
													return 0;
}

double integrator_simpson_blackbody_spectrum(double z,
																						 double E_ini,
																						 double E_max,
																						 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

	double T=T_0*(1+z);
	double dE = (E_max- E_ini)/ (double) (10*pt_Spectrum_and_Precision_Parameters->n_step-1);
	// cout << "E_max "<< E_max << " E_ini " << E_ini <<" dE = " << dE << endl;
	double E[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h2;
	double resultat = 0;
	// cout << "here "<< dE << endl;
	int y = 1;
	for(int j=0; j<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;j++){
										h2 = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
										// cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
										for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
										{
											if(eval == 0){
											if(j==0)	E[eval]=E_ini;
											else E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
											}
											else{
												E[eval]=E[0]+eval*h2;
											}

										if((exp(E[eval]/T)-1)!=0.)f[eval]=E[eval]*E[eval]/(exp(E[eval]/T)-1);
										else f[eval] = 0;
										resultat += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
										// cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
									}

									// if(res_initial==0 && resultat !=0)res_initial = resultat;
									// if(res_initial!=0 && resultat/res_initial<precision)break;
									// cout << E1 << " f1 " << f1 << " (exp(E1/T)-1) "  << (exp(E1/T)-1) << " f7 " << f7 << endl;


									// 	cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << endl;
								}
							resultat*=1./(pi*pi);
							// if(resultat != 0)cout << "(integrator_simpson_blackbody_spectrum_over_rate_inverse_compton : ) resultat = " << resultat << endl;
 							return resultat;
}
double function_integrand_pair_creation(double E_e, double E_gamma, double E_gamma_bar){
	double E_e_prime = E_gamma+E_gamma_bar-E_e, A, B, C, D;

	  A = 4*pow(E_e+E_e_prime,2)/(E_e*E_e_prime)*log(4*E_gamma_bar*E_e*E_e_prime/(m_e*m_e*(E_e+E_e_prime)));
		B =	-(1-m_e*m_e/(E_gamma_bar*(E_e+E_e_prime)))*pow(E_e+E_e_prime,4)/pow(E_e*E_e_prime,2);
		C =	2*(2*E_gamma_bar*(E_e+E_e_prime)-m_e*m_e)*pow(E_e+E_e_prime,2)/(m_e*m_e*E_e*E_e_prime);
		D = -8*E_gamma_bar*E_gamma/(m_e*m_e);
		// if(A < 0) A = 0;
		// if(B < 0) B = 0;
		// if(C < 0) C = 0;
		// if(D < 0) D = 0;
	double result = A + B +C +D;
	if(result < 0)result = 0;
	if(result < 0)cout << "(function_integrand_pair_creation : ) A =" << A <<" B = " << B << " C = " << C << " D = " << D << " result = " << result << " E_e = " << E_e << " E_gamma = " << E_gamma << " E_e_prime = " << E_e_prime << endl;
	return 2*result;			// We treat electron and positron on an equal footing, hence we simply multiply here by 2.
}
void Triangular_Spectrum(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
												 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
												 Structure_Spectrum * pt_Cascade_Spectrum,
											 	 Structure_Spectrum * pt_Electron_Spectrum,
											 	 Structure_Output_Options * pt_Output_Options){


	double E_e_minus_1, E_e_plus_1, E_j, E_j_minus_1, E_j_plus_1, dE_j, dE;
	double resultat, E_e, E_g, f_e, check;
	double z = pt_Cascade_Spectrum->redshift;

	Structure_Spectrum Tmp_Electron_Spectrum;
	Structure_Spectrum Tmp_Photon_Spectrum;
	Tmp_Electron_Spectrum.spectrum_name = "cascade_electrons";
	Tmp_Electron_Spectrum.redshift = z;
	Tmp_Electron_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Tmp_Electron_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);


	Tmp_Photon_Spectrum.spectrum_name = "cascade_photons";
	Tmp_Photon_Spectrum.redshift = z;
	Tmp_Photon_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Tmp_Photon_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);



	double E_c = E_c_0/(1+z), E_x = E_x_0/(1+z);
	double E_phph = m_e*m_e/(T_0*(1+z));
	double E_gamma_bb = 2.701*T_0*(1+z), E_e_ICS, E_g_lim;
	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	double E_0 = pt_Particle_Physics_Model->E_0;
	double integrale=E_0;
	double tau_x = pt_Particle_Physics_Model->tau_x;
	double Rate_photons_E_j = 0, Rate_photons_E_g = 0, Rate_electrons_E_j = 0, Rate_electrons_E_e = 0, Rate_electrons_E_e_2;
	double rate_PP= 0, rate_COM= 0, rate_DP= 0, rate_NP = 0;
	double n_x = n_y_0*pow(1+z,3)*exp(-1/(2*H_r*pow(1+z,2)*tau_x))/E_0;
	double E_max = E_0, E_cmb_min, E_cmb_max, gamma_e;
	double E_s, x_j, gamma_prime;
	double PP= 0, CS= 0, ICS_g = 0, ICS_e= 0, NPC= 0, COM= 0, DP= 0;
	bool approximation =0;

	if(pt_Output_Options->EM_cascade_verbose > 3)cout << "E_c = " << E_c << "E_max = " << E_max << " E_phph = " << E_phph << " int_bb = " << int_bb << endl;

			for(int i = (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1); i>=0 ;i--){
			resultat=0;
			Rate_photons_E_g = 0;
			E_cmb_max = 10*E_gamma_bb;
			E_cmb_min = E_gamma_bb/10.;
			E_e = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,(double) i/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
			E_g = E_e;
			E_e_plus_1 = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,((double) i+1)/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
			E_e_minus_1 = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,((double) i-1)/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
			dE = (E_e_plus_1 - E_e_minus_1)/2.;
			pt_Electron_Spectrum->Energy[i] = E_e;
			pt_Cascade_Spectrum->Energy[i] = E_g;
			Tmp_Electron_Spectrum.Energy[i] = E_e;
			Tmp_Photon_Spectrum.Energy[i] = E_g;

			// if(approximation == 0)Rate_electrons_E_e = integrator_simpson_rate_inverse_compton(z,E_cmb_min,E_cmb_max,E_g,pt_Spectrum_and_Precision_Parameters);
			// if(Rate_Inverse_Compton(E_e,z,pt_Spectrum_and_Precision_Parameters)<Rate_electrons_E_e || approximation == 1){
			// 	Rate_electrons_E_e =Rate_Inverse_Compton(E_e,z,pt_Spectrum_and_Precision_Parameters);
			// 	approximation = 1;
			// }
			Rate_electrons_E_e = integrator_simpson_rate_inverse_compton(z,E_cmb_min,E_cmb_max,E_g,pt_Spectrum_and_Precision_Parameters);

		if(pt_Output_Options->EM_cascade_verbose > 3)cout << "(Rate electrons : ) at E = " << E_e << " tot = " << Rate_electrons_E_e << endl;


			if(pt_Spectrum_and_Precision_Parameters->pair_creation_in_nuclei == "yes") rate_NP= rate_NPC(E_g,z);
			else rate_NP = 0;
			if(pt_Spectrum_and_Precision_Parameters->compton_scattering == "yes") rate_COM = rate_compton(E_g,z);
			else rate_COM = 0;

			if(pt_Spectrum_and_Precision_Parameters->double_photon_pair_creation=="yes" && E_g >= E_c)rate_DP=rate_pair_creation_v2(E_g,z,pt_Spectrum_and_Precision_Parameters);
			else rate_DP = 0.;
			if(pt_Spectrum_and_Precision_Parameters->photon_photon_diffusion == "yes" && E_g < E_phph )rate_PP=rate_gg_scattering(E_g,z);
			else rate_PP = 0.;
			Rate_photons_E_g += rate_PP;
			Rate_photons_E_g += rate_NP;
			Rate_photons_E_g += rate_COM;
			Rate_photons_E_g += rate_DP;

		if(pt_Output_Options->EM_cascade_verbose > 3)cout << "(Rate photons : ) at E = " << E_g << " rate_NP = " << rate_NP << " rate_COM = " << rate_COM << " rate_PP = " << rate_PP << " rate_DP = " << rate_DP << " tot = " << Rate_photons_E_g << endl;


			if(i==(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1)){
				if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "Dirac"){
					pt_Electron_Spectrum->Spectrum[i]=1./(dE*Rate_electrons_E_e);
				}
				if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "Dirac"){
					// pt_Cascade_Spectrum->Spectrum[i]=n_x/(dE*(Rate_photons_E_g*tau_x));
					pt_Cascade_Spectrum->Spectrum[i]=1./(dE*Rate_photons_E_g);
				}
				// else{
				// 	pt_Electron_Spectrum->Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum(E_e,z,pt_Particle_Physics_Model->E_0)/Rate_electrons_E_e;
				// 	pt_Cascade_Spectrum->Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E_g,z,pt_Particle_Physics_Model->E_0)/Rate_photons_E_g;
				// }
			}
			else{

				for (int j = i+1 ; j < pt_Spectrum_and_Precision_Parameters->Energy_Table_Size; j ++){
					E_j = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,(double) j/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
					E_j_plus_1 = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,((double) j+1)/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
					E_j_minus_1 = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,((double) j-1)/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));

					dE_j = (E_j_plus_1 - E_j_minus_1)/2.;

					// if(pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E_j,z,pt_Particle_Physics_Model->E_0)!=0){
					// 	Rate_photons_E_j = rate_NPC(E_j,z)+rate_compton(E_j,z)+rate_gg_scattering(E_j,z);
					// 	if(pt_Spectrum_and_Precision_Parameters->double_photon_pair_creation=="yes" && E_j>E_c/2.){
					// 		Rate_photons_E_j+=rate_pair_creation(E_j,z,pt_Spectrum_and_Precision_Parameters);
					// 	}
					// 	pt_Cascade_Spectrum->Spectrum[i] +=pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E_j,z,pt_Particle_Physics_Model->E_0)/Rate_photons_E_j;
					// 	cout << "(source :) pt_Cascade_Spectrum->Spectrum[i] = " << pt_Cascade_Spectrum->Spectrum[i] << " E = " << E_g << " E_x = " << E_x << endl;
					// }
					// if(pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum(E_j,z,pt_Particle_Physics_Model->E_0)!=0){
					// 	Rate_electrons_E_j =	Rate_Inverse_Compton(E_j,z,pt_Spectrum_and_Precision_Parameters);
					// 	pt_Electron_Spectrum->Spectrum[i]+=pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum(E_j,z,pt_Particle_Physics_Model->E_0)/Rate_electrons_E_j;
					// }

					if(Rate_electrons_E_e!=0){
						if(pt_Spectrum_and_Precision_Parameters->compton_scattering == "yes")COM = dE_j * pt_Cascade_Spectrum->Spectrum[j] * dsigma_compton(E_j,z,(E_j+m_e-E_e));
						else COM = 0;
						if(pt_Spectrum_and_Precision_Parameters->pair_creation_in_nuclei == "yes")NPC = dE_j * pt_Cascade_Spectrum->Spectrum[j] * dsigma_NPC(E_j+m_e,z,E_e);
						else NPC = 0;
						if(pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering == "yes"){
						gamma_e = E_j/m_e;
						gamma_prime = E_e/m_e;

						while(E_cmb_max<E_cmb_min)E_cmb_max*=10.;
						// ICS_e =  dE_j * pt_Electron_Spectrum->Spectrum[j]  * 2*pi*r_e*r_e*m_e*m_e*integrator_simpson_scattered_electron_inverse_compton(z,E_e,E_j,pt_Spectrum_and_Precision_Parameters)/(E_j*E_j);
						ICS_e =  dE_j * pt_Electron_Spectrum->Spectrum[j]  * Analytical_form_scattered_electron_from_inverse_compton(z, gamma_e,  gamma_prime,  pt_Spectrum_and_Precision_Parameters);
						}
						else ICS_e = 0;
						if(pt_Spectrum_and_Precision_Parameters->double_photon_pair_creation=="yes" && E_j>=E_c){
							DP = dE_j * pt_Cascade_Spectrum->Spectrum[j] * 1/4.*pi*r_e*r_e*pow(m_e,4)*dsigma_pair_creation(z,E_e,E_j,pt_Spectrum_and_Precision_Parameters)/(E_j*E_j*E_j);
							if(DP<0) DP = 0.;
						}
						else DP = 0.;
						pt_Electron_Spectrum->Spectrum[i] += COM/(Rate_electrons_E_e);
						pt_Electron_Spectrum->Spectrum[i] += NPC/(Rate_electrons_E_e);
						pt_Electron_Spectrum->Spectrum[i] += ICS_e/(Rate_electrons_E_e);
						pt_Electron_Spectrum->Spectrum[i] += DP/(Rate_electrons_E_e);
						Tmp_Electron_Spectrum.Spectrum[i] = pt_Electron_Spectrum->Spectrum[i]*Rate_electrons_E_e;
					}
					else pt_Electron_Spectrum->Spectrum[i] += 0;

					if(pt_Output_Options->EM_cascade_verbose > 3)cout <<"(Scattering electrons : ) at E = " << E_e << " E_j = " <<  E_j <<  " COM = " << COM << " NPC = " << NPC << " ICS_e = " << ICS_e << " DP = " << DP << endl;
					if(Rate_photons_E_g!=0){
					if(	pt_Spectrum_and_Precision_Parameters->photon_photon_diffusion == "yes" &&  E_g < E_phph)	PP = dE_j * pt_Cascade_Spectrum->Spectrum[j] * (dsigma_phph(E_j,z,E_g));
					else PP = 0;
					if(pt_Spectrum_and_Precision_Parameters->compton_scattering == "yes") CS = dE_j * pt_Cascade_Spectrum->Spectrum[j] * (dsigma_compton(E_j,z,E_g));
					else CS = 0;




					if(pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering == "yes"){
						E_e_ICS = E_j;
						gamma_e = E_e_ICS/m_e;
						if(E_e_ICS <= E_max){
							linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E_e_ICS, f_e);
						}
						else f_e = 0;

						ICS_g = dE_j*f_e*gamma_inverse_compton_analytical_v2(gamma_e,E_g,z)/(E_g);
						// ICS_g = dE_j*f_e*gamma_inverse_compton_analytical(gamma_e,E_g,z,3,pt_Spectrum_and_Precision_Parameters);


					}
					else ICS_g = 0;
					pt_Cascade_Spectrum->Spectrum[i] += CS/(Rate_photons_E_g);
					pt_Cascade_Spectrum->Spectrum[i] += PP/(Rate_photons_E_g);
					pt_Cascade_Spectrum->Spectrum[i] += ICS_g/(Rate_photons_E_g);
					Tmp_Photon_Spectrum.Spectrum[i] = pt_Cascade_Spectrum->Spectrum[i]*(Rate_photons_E_g);
				if(pt_Output_Options->EM_cascade_verbose > 3)cout <<"(Scattering photons : ) at E = " << E_g << " E_j = " <<  E_j <<  " PP = " << PP << " CS = " << CS << " ICS_g = " << ICS_g << "f_e " << f_e << endl;
				}
				else pt_Cascade_Spectrum->Spectrum[i] += 0;
			}
			}
			if(pt_Output_Options->EM_cascade_verbose > 1)cout << " ************************************** E = "<< pt_Electron_Spectrum->Energy[i] << " pt_Electron_Spectrum->Spectrum[i] = " << pt_Electron_Spectrum->Spectrum[i] << " pt_Cascade_Spectrum->Spectrum[i] = " << pt_Cascade_Spectrum->Spectrum[i]  <<endl;
			}
			if(pt_Spectrum_and_Precision_Parameters->check_energy_conservation == "yes"){
				check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,&Tmp_Photon_Spectrum,&Tmp_Electron_Spectrum,integrale);

				for(int i = (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1); i>=0 ;i--){

					pt_Cascade_Spectrum->Spectrum[i]*=E_0/integrale;

				}
			}

			else cout << "*** No energy check energy conservation check requested ***" << endl;


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
	pt_Cascade_Spectrum->Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	pt_Cascade_Spectrum->Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	pt_Cascade_Spectrum->species = "photon";
	pt_Cascade_Spectrum->spectrum_name = "total_photon_";
	pt_Cascade_Spectrum->redshift = z;

	Inverse_Compton_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Inverse_Compton_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Inverse_Compton_Spectrum.species = "photon";
	Inverse_Compton_Spectrum.spectrum_name = "ICS_";
	Inverse_Compton_Spectrum.redshift=z;

	Diffused_Gamma_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Diffused_Gamma_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Diffused_Gamma_Spectrum.species = "photon";
	Diffused_Gamma_Spectrum.spectrum_name = "diffused_photon_";
	Diffused_Gamma_Spectrum.redshift=z;

	Electron_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Electron_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Electron_Spectrum.species="electron";
	Electron_Spectrum.spectrum_name = "total_electrons_";
	Electron_Spectrum.redshift=z;

	Diffused_Electron_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Diffused_Electron_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Diffused_Electron_Spectrum.species = "electron";
	Diffused_Electron_Spectrum.spectrum_name = "diffused_electron_";
	Diffused_Electron_Spectrum.redshift=z;

	Compton_Electron_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Compton_Electron_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Compton_Electron_Spectrum.species = "electron";
	Compton_Electron_Spectrum.spectrum_name = "compton_electron_";
	Compton_Electron_Spectrum.redshift=z;

	Tmp_Electron_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Tmp_Electron_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
	Tmp_Electron_Spectrum.species = "electron";
	Tmp_Electron_Spectrum.spectrum_name = "tmp_electron_";
	Tmp_Electron_Spectrum.redshift=z;
	double E_max = E_0;
	if(E_0>1.5*E_c){
		E_max = 1.5*E_c;
	}

	dE = (E_max - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);


/*****************Start of the computation****************/

	if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "universal"){
					for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
						E1=pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
						pt_Cascade_Spectrum->Energy[i]=E1;
						pt_Cascade_Spectrum->Spectrum[i]=universal_spectrum(E1,z,pt_Particle_Physics_Model->E_0);
						Electron_Spectrum.Energy[i]=E1;
						Electron_Spectrum.Spectrum[i]=0;
					}
					pt_Cascade_Spectrum->spectrum_name="universal_photon_";
					if(pt_Spectrum_and_Precision_Parameters->check_energy_conservation == "yes")check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);
					else cout << "*** No energy check energy conservation check requested ***" << endl;
					print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);

	}
	else{
					if(E_c <= pt_Particle_Physics_Model->E_0 && pt_Spectrum_and_Precision_Parameters->calculation_mode == "simplified"){


								for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
									E1=pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
									pt_Cascade_Spectrum->Energy[i]=E1;
									pt_Cascade_Spectrum->Spectrum[i]=universal_spectrum(E1,z,pt_Particle_Physics_Model->E_0);
									Electron_Spectrum.Energy[i]=E1;
									Electron_Spectrum.Spectrum[i]=0;
								}
								if(pt_Spectrum_and_Precision_Parameters->check_energy_conservation == "yes")check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);
								else cout << "*** No energy check energy conservation check requested ***" << endl;
								if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
										if(pt_Output_Options->EM_cascade_verbose>1)cout <<" I will now print the spectrum in files." << endl;
										print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
								}
					}


					else{




						if(pt_Spectrum_and_Precision_Parameters->calculation_mode == "triangular"){

							if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "from_file"){
								Electron_Spectrum.Energy.resize(0);
								Electron_Spectrum.Spectrum.resize(0);
								Cascade_Spectrum_Reading_From_File(z,pt_Particle_Physics_Model,&Electron_Spectrum,pt_Spectrum_and_Precision_Parameters);
							}
							else if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "from_function"){
								// to be filled
							}
							else{
								for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
										E1=pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
										Electron_Spectrum.Energy[i]=E1;
										Electron_Spectrum.Spectrum[i]=0;
										// Electron_Spectrum.Spectrum[i]=Tmp_Electron_Spectrum.Spectrum[i];
								}
							}

							if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "from_file"){

							}
							else if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "from_function"){
								// to be filled
							}
							else{
								for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
										E1=pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
										pt_Cascade_Spectrum->Energy[i]=E1;
										pt_Cascade_Spectrum->Spectrum[i]=0;
										// Electron_Spectrum.Spectrum[i]=Tmp_Electron_Spectrum.Spectrum[i];
								}
							}
							Electron_Spectrum.spectrum_name = "total_electrons_";
							pt_Cascade_Spectrum->spectrum_name = "total_photons_";
							Triangular_Spectrum(pt_Particle_Physics_Model,
																 pt_Spectrum_and_Precision_Parameters,
																 pt_Cascade_Spectrum,
																 &Electron_Spectrum,
															 	 pt_Output_Options);
						print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Electron_Spectrum, pt_Particle_Physics_Model);
						print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
						}



						else if(pt_Spectrum_and_Precision_Parameters->calculation_mode == "iterative"){
							/********First step : compute initial ICS spectrum from the electon spectrum injected**********/

							for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
									E1=pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
									Tmp_Electron_Spectrum.Energy[i]=E1;
									Tmp_Electron_Spectrum.Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum(E1,z,pt_Particle_Physics_Model->E_0);
									Electron_Spectrum.Energy[i]=E1;
									Electron_Spectrum.Spectrum[i]=Tmp_Electron_Spectrum.Spectrum[i];
							}

						for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
								E1=pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
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
											for(int j=0;j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;j++){
											E1=pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE;
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
											if(pt_Output_Options->EM_cascade_verbose>1)cout <<" I will now print the spectrum in files." << endl;
											print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Inverse_Compton_Spectrum, pt_Particle_Physics_Model);
									}



						}
						/*********Second step : compute the gamma spectrum from gg->gg, ge->ge, gN->Nee for a certain number of iterations.*********/

							for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
										pt_Cascade_Spectrum->spectrum_name="Initial_photon_spectrum_";
										E1=pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
										pt_Cascade_Spectrum->Energy[i]=E1;
										if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice=="none")pt_Cascade_Spectrum->Spectrum[i]=0;
										else pt_Cascade_Spectrum->Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E1,z,pt_Particle_Physics_Model->E_0)/(rate_NPC(E1,z)+rate_compton(E1,z)+rate_gg_scattering(E1,z));
										pt_Cascade_Spectrum->Spectrum[i]+=Inverse_Compton_Spectrum.Spectrum[i];
								}
								// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);

								if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
										pt_Cascade_Spectrum->spectrum_name = "Cascade_";
										if(pt_Output_Options->EM_cascade_verbose>1)cout <<" I will now print the spectrum in files." << endl;
										print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
								}



							for(int k = 0; k<pt_Spectrum_and_Precision_Parameters->number_iterations_photon;k++){
										if(pt_Output_Options->EM_cascade_verbose>1)cout<<"iteration : " << k+1 << endl;


								Spectrum_gamma_scattered(pt_Particle_Physics_Model,
																				 pt_Spectrum_and_Precision_Parameters,
																				 pt_Cascade_Spectrum,
																				 &Diffused_Gamma_Spectrum);
								 for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
										E_gamma = pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
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

						for(int j=0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size; j++){
							Tmp_Electron_Spectrum.Energy[j]=Compton_Electron_Spectrum.Energy[j];
							Tmp_Electron_Spectrum.Spectrum[j]=Compton_Electron_Spectrum.Spectrum[j];
						}
						for(int i = 0 ; i < pt_Spectrum_and_Precision_Parameters->number_iterations_electron; i++){
							cout << "iteration : " << i+1 << endl;
						 Spectrum_electron_scattered(pt_Particle_Physics_Model,
																				pt_Spectrum_and_Precision_Parameters,
																				&Tmp_Electron_Spectrum,
																				&Diffused_Electron_Spectrum);
						for(int j=0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size; j++){
							cout <<"E = " << Electron_Spectrum.Energy[j] << " Compton_Electron_Spectrum = " << Compton_Electron_Spectrum.Spectrum[j] << " Diffused_Electron_Spectrum = " << Diffused_Electron_Spectrum.Spectrum[j]<<endl;
							Tmp_Electron_Spectrum.Spectrum[j]= Compton_Electron_Spectrum.Spectrum[j]+Diffused_Electron_Spectrum.Spectrum[j];
							cout << " Electron_Spectrum = " << Electron_Spectrum.Spectrum[j] << endl;
						}
						print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Tmp_Electron_Spectrum, pt_Particle_Physics_Model);

					}
					for(int j=0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size; j++){
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
							for(int l = 0; l<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size ; l++){
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
