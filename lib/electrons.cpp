#include "bbn/EM_cascade.h"
#include "bbn/injected_spectrum.h"
#include "bbn/structures.h"
#include "bbn/BBN_constraints.h"
#include "bbn/tools.h"
#include "bbn/test_functions.h"
#include "bbn/photons.h"
#include "bbn/electrons.h"

using namespace std;

/**
* @brief  the electron module :
* Contains all processes undergone by the electron population.
* @file electrons.cpp
*/

double dsigma_inverse_compton_electron_spectrum(double z,
        double gamma_e,
        double gamma_prime,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options)
{
  /// Parameterization introduced in Zdziarski 1988 ApJ 335:786-802.
  /// Can be used to compute the gamma spectrum by setting gamma_prime -> gamma_e - (E_photon/m_e). It has been checked to work.
    double T = T_0*(1+z);
    double A = 8*pi*pow(m_e/(2*pi),3);
    double a, b;
    double r = 0.5*(gamma_e/gamma_prime + gamma_prime/gamma_e);
    int N = 100;
    double theta = T/m_e;
    double k = (gamma_e/gamma_prime-1)/(4*gamma_e*theta);
    double Eps_0 = 0, Eps_1 = 0;
    double result, LI2 = 0;

    for(int i = 1; i <= N ; i++) {
        // Eps_0 += exponential_integral(i*k,pt_Spectrum_and_Precision_Parameters);
        // Eps_1 += exponential_integral(i*k,pt_Spectrum_and_Precision_Parameters)/i;
        Eps_0 += expint(1,i*k);
        Eps_1 += expint(1,i*k)/i;
    }
    // for(int i = 1 ; i <= 10*N; i++){
    // 	LI2 += pow(exp(-k),i)/(i*i);
    // }
    //
    // a = r*LI2;
    a = r*polylog_2(exp(-k),pt_Spectrum_and_Precision_Parameters);
    b = 2*k*(log(1-exp(-k))+k*Eps_0+Eps_1);
    result = 3*sigma_T*A*theta*theta/(4*gamma_e*gamma_e)*(a-b);
    // if(result < 0) {
    //     result = 0;
    // }
    if(pt_Output_Options->EM_cascade_verbose > 2) {
        cout << "(dsigma_inverse_compton_electron_spectrum : )" <<" z " << z << " gamma_e "<< gamma_e << "gamma_prime "<< gamma_prime << " result "<< result << endl;
    }
    if(gamma_e == gamma_prime || k < 1) result = 0;
    // cout << " exp(-k) "<< exp(-k) << " Eps_0 " << Eps_0 << " Eps_1 " << Eps_1  << " result " << result << " a/r = " << a/r <<" LI2 = " << LI2<< " b = " << b <<endl;
    result/=(m_e*m_e); //Conversion factor from m_e unit to MeV.
    return result*4*pi;
}




double dsigma_inverse_compton_electron_spectrum_v2(double z,
        double gamma_e,
        double gamma_prime,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options){

  /// Taken from Zdziarski 1988 ApJ 335:786-802. Compute numerically the integral, only for checking i) the integration algorithm ii) the analytical formula.
  double T = T_0*(1+z);
  double A = 8*pi*pow(m_e/(2*pi),3);
  double r = 0.5*(gamma_e/gamma_prime + gamma_prime/gamma_e);
  double theta = T/m_e;
  double x[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h;
  double E_s = (gamma_e/gamma_prime-1)/4, E;
  double x_cmb_min = E_s / gamma_e, x_cmb_max = 1000*x_cmb_min;
  double result = 0;
  double dlogx = (log(x_cmb_max) - log(x_cmb_min))/(pt_Spectrum_and_Precision_Parameters->n_step);
  h = dlogx/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
    for(int j=0; j<pt_Spectrum_and_Precision_Parameters->n_step; j++) {
        // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
        for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {

            x[eval]=exp(log(x_cmb_min)+eval*h+j*dlogx);
            E = x[eval]*gamma_e;
            if(E>E_s) {
                f[eval]=(r+(2-r)*E_s/E-2*pow(E_s/E,2)-2*E_s/E*log(E/E_s))*x[eval]*x[eval]/(exp(x[eval]/theta)-1);
            } else {
                f[eval] = 0;
            }
            result += dlogx*x[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
            // cout << "eval " << eval << "E = " << x[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" result = " << result << endl;
        }

    }
result *=3*A*sigma_T/(4*gamma_e*gamma_e*m_e*m_e*m_e);
// cout << "result = " << result << endl;
return result ;
// return 0;

}
double dsigma_inverse_compton_electron_spectrum_v3(double z,
        double E_e,
        double E_e_prime,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options){

  /// Taken from Zdziarski 1988 ApJ 335:786-802. Compute numerically the integral, only for checking i) the integration algorithm ii) the analytical formula.
  double T = T_0*(1+z);
  double E[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h;
  double Gamma, E_gamma, q;
  double E_s = (E_e/E_e_prime-1)/4;
  double E_cmb_min = (E_s / E_e)*pow(m_e,2), E_cmb_max = 1000*E_cmb_min;
  double result = 0;
  double dlogE = (log(E_cmb_max) - log(E_cmb_min))/(pt_Spectrum_and_Precision_Parameters->n_step);
  h = dlogE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
    for(int j=0; j<pt_Spectrum_and_Precision_Parameters->n_step; j++) {
        // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
        for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {

            E[eval]=exp(log(E_cmb_min)+eval*h+j*dlogE);
            Gamma = 4*E[eval]*E_e/(m_e*m_e);
            E_gamma = E_e + E[eval] - E_e_prime;
            q = E_gamma/(Gamma*(E_e-E_gamma));
            if(q<=1 && q >= 0) {
                f[eval]=(2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma*q,2)*(1-q)/(2*(1-Gamma*q)))*E[eval]/(exp(E[eval]/T)-1);
            } else {
                f[eval] = 0;
            }
            if(f[eval]<0)f[eval]*=-1;
            result += dlogE*E[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
            // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" result = " << result << endl;
        }

    }
result *= 2*r_e*r_e*m_e*m_e/(E_e*E_e*pi);
cout << "result = " << result << endl;
return result ;
// return 0;

}

double gamma_inverse_compton_analytical_v2(double gamma_e, double E_gamma, double z,Structure_Output_Options * pt_Output_Options)
{
  /// Parameterization introduced in Petruk 2009 arXiv:0807.1969v2 [astro-ph].

    double T = T_0*(1+z);
    double eta_c = T*E_gamma/(m_e*m_e);
    double eta_0 = E_gamma*E_gamma/(4*gamma_e*m_e*(gamma_e*m_e-E_gamma));
    // cout << "eta_c " << eta_c << "eta_0 " << eta_0 << endl;

    double I = pi*pi/6*eta_c*(exp(-5./4.*pow(eta_0/eta_c,0.5))+2*eta_0*exp(-5./7*pow(eta_0/eta_c,0.7)))*exp(-2./3.*eta_0/eta_c);
    // return sigma_T*T*T*E_gamma/(8*gamma_e*gamma_e);
    if(pt_Output_Options->EM_cascade_verbose > 2) {
        cout << "(gamma_inverse_compton_analytical_v2 : )" <<"z " << z << " gamma_e "<< gamma_e << "E_gamma "<< E_gamma << " result " << 3*sigma_T*T*m_e*m_e/(4*pi*pi*gamma_e*gamma_e)*I << endl;
    }

    return 3*sigma_T*T*m_e*m_e/(4*pi*pi*gamma_e*gamma_e)*I/m_e;
}

double gamma_inverse_compton_analytical(double gamma_e,
                                        double E_gamma,
                                        double z,
                                        int N,
                                        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                        Structure_Output_Options * pt_Output_Options)
{
    ///Parameterization introduced in Zdziarski and Pjanka 201, arXiv:1307.6732v2 [astro-ph.HE].
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

    if(x <= q)	{
        f_minus_1 = 0.0366377 + pow(x,-1) + log(x)/2 -pow(q,-1)-log(q)/2;
    } else {
        f_minus_1 = 0;
    }
    if(x<= q)	{
        f_ln =0.5*(q*(1-log(q))+pow(log(q),2)) - 0.5*(x*(1-log(x))+pow(log(x),2)) + 0.125505;
    } else {
        f_ln = 0;
    }
    for(int i = 1; i <= N ; i++) {
        if(x <= q) {
            f_minus_1 +=pow(q,2*i-1)*B[2*i]/((2*i-1)*factorial(2*i)) - pow(x,2*i-1)*B[2*i]/((2*i-1)*factorial(2*i));
        } else {
            f_minus_1 += expint(1,i*x);
        }
        // else f_minus_1 += exponential_integral(i*x,pt_Spectrum_and_Precision_Parameters);
        if(x <= q) {
            f_ln = pow(q,2*i)*(2*i*log(q)-1)*B[2*i]/(4*i*i*factorial(2*i))-pow(x,2*i)*(2*i*log(x)-1)*B[2*i]/(4*i*i*factorial(2*i));
        } else {
            f_ln = (expint(1,i*x)+exp(-i*x)*log(x))/i;
        }
        // else f_ln = (exponential_integral(i*x,pt_Spectrum_and_Precision_Parameters)+exp(-i*x)*log(x))/i;
    }
    // double I = pi*pi/6*eta_c*(exp(-5./4.*pow(eta_0/eta_c,0.5))+2*eta_0*exp(-5./7*pow(eta_0/eta_c,0.7)))*exp(-2./3.*eta_0/eta_c);
    // return sigma_T*T*T*E_gamma/(8*gamma_e*gamma_e);
    // cout << " x = " << x << endl;
    // cout << "f1 = " << f_1 << " f_0 = " << f_0 << " f_ln = " << f_ln << " f_minus_1 = " << f_minus_1 << endl;
    result = 3*sigma_T*A/(4*gamma_e*gamma_e)*((1-2*E_gamma*E_n+2*log(E_n/theta))*E_n*theta*f_0+(1+2*E_gamma*E_n)*theta*theta*f_1-2*E_n*E_n*f_minus_1-2*E_n*E_n*f_minus_1-2*E_n*theta*f_ln);
    if(result < 0) {
        result = 0;
    }
    if(pt_Output_Options->EM_cascade_verbose > 2) {
        cout << "(gamma_inverse_compton_analytical : )" <<" z " << z << " gamma_e "<< gamma_e << "E_gamma "<< E_gamma << " result " << result << endl;
    }

    return result/(m_e*m_e);
}







double Analytical_form_inverse_compton(double s, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters)
{
  /// Taken from Zdziarski 1988 ApJ 335:786-802. To be used with rate_electron_inverse_compton.


    double result = 0;
    if(s<1){
      for(int n =2 ; n < 100 ; n++){
        result += 3./4.*(pow(n,4)-7*n*n+18*n-8)/((n-1)*n*n*(n+1))*pow(-s,n-2);
      }
    }
    // result = (s+9+8./s)*log(1+s)-(16+18*s+s*s)/(2*(1+s))+4*polylog_2(-s,pt_Spectrum_and_Precision_Parameters);
    else {
      result = (s+9+8./s)*log(1+s)-8-(2*s+s*s)/(2+2*s)+4*polylog_2(-s,pt_Spectrum_and_Precision_Parameters);
      result *= 3/(2*s*s);
    }
    if(-s > 1) {
        cout << "result " << result << " s " << s << endl;
    }
    return result;
}
double rate_electron_inverse_compton(double z,
        double E_e,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options)
{
  /// Taken from Zdziarski 1988 ApJ 335:786-802. To be used with Analytical_form_inverse_compton.

    double T=T_0*(1+z);
    double A = 8*pi*pow(m_e/(2*pi),3);
    double E_gamma_bb = 2.701*T;
    double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
    double x_ini = E_gamma_bb/10./m_e, x_max = E_gamma_bb*10./m_e;
    double dx = (x_max- x_ini)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1), dlogx;
    double theta = T/m_e, s;
    double x[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h2;

    int y=0;
    // while(dE>E_ini) {
    //     dE/=10.;
    //     y++;
    // }
    // cout << "(rate_electron_inverse_compton :) E_max "<< E_max << " E_ini " << E_ini <<" dE = " << dE << " y = " << y << endl;
    // h2 = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
    dlogx = (log(x_max)-log(x_ini))/(pt_Spectrum_and_Precision_Parameters->n_step-1);
    h2 = dlogx/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

    // x[eval]=exp(log(x_cmb_min)+eval*h2+j*dlogx);

    double result = 0;
    // cout << "here "<< dE << endl;

    for(int j=0; j<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1; j++) {

        // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
        for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {

            // E[eval]=E_ini+eval*h2+j*dE;
            x[eval]=exp(log(x_ini)+eval*h2+j*dlogx);
            s = 4*E_e*x[eval]/m_e;
            if((exp(x[eval]/theta)-1)!=0.) {
                f[eval]=Analytical_form_inverse_compton(s,pt_Spectrum_and_Precision_Parameters)*x[eval]*x[eval]/(exp(x[eval]/theta)-1);
            } else {
                f[eval] = 0;
            }
            result += dlogx*x[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
            // result += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
            // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" result = " << result << endl;
        }

    }

    result*=sigma_T*A;
    result/=(m_e);      //Conversion factor from unit m_e to unit MeV.
    if(pt_Output_Options->EM_cascade_verbose > 2) {
        cout << "(rate_electron_inverse_compton : )" <<" z " << z << " E_e "<< E_e << " result "<<  result/(pi*pi) << endl;
    }

    return result;
}
