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
* Processus undergone by electron population
*/

double dsigma_inverse_compton_electron_spectrum(double z,
        double gamma_e,
        double gamma_prime,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options)
{
    double T = T_0*(1+z);
    double A = 8*pi*pow(m_e/(2*pi),3);
    double a, b;
    double r = 0.5*(gamma_e/gamma_prime + gamma_prime/gamma_e);
    int N = 5;
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
    for(int i = 1 ; i < 200*N; i++){
    	LI2 += pow(exp(-k),i)/(i*i);
    }

    a = r*LI2;
    // a = r*polylog_2(exp(-k),pt_Spectrum_and_Precision_Parameters);
    b = 2*k*(log(1-exp(-k))+k*Eps_0+Eps_1);
    result = 3*sigma_T*A*theta*theta/(4*gamma_e*gamma_e)*(a-b);
    // if(result < 0) {
    //     result = 0;
    // }
    if(pt_Output_Options->EM_cascade_verbose > 2) {
        cout << "(dsigma_inverse_compton_electron_spectrum : )" <<" z " << z << " gamma_e "<< gamma_e << "gamma_prime "<< gamma_prime << " result "<< result << endl;
    }
    if(gamma_e == gamma_prime) result = 0;
    // cout << " exp(-k) "<< exp(-k) << " Eps_0 " << Eps_0 << " Eps_1 " << Eps_1  << " result " << result << " a/r = " << a/r <<" LI2 = " << LI2<< " b = " << b <<endl;
    return result;
}




double dsigma_inverse_compton_electron_spectrum_v2(double z,
        double gamma_e,
        double gamma_prime,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options){

  double T = T_0*(1+z);
  double A = 8*pi*pow(m_e/(2*pi),3);
  double r = 0.5*(gamma_e/gamma_prime + gamma_prime/gamma_e);
  double theta = T/m_e;
  double x[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h;
  double E_s = (gamma_e/gamma_prime-1)/4, E;
  double x_cmb_min = E_s / gamma_e, x_cmb_max = 10*(2.701*T)/m_e;
  double result = 0;
  double dlogx = (log(x_cmb_max) - log(x_cmb_min))/(pt_Spectrum_and_Precision_Parameters->n_step);
  h = dlogx/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
  if(E>E_s){
    for(int j=0; j<pt_Spectrum_and_Precision_Parameters->n_step; j++) {
        // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
        for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {

            x[eval]=exp(log(x_cmb_min)+eval*h+j*dlogx);
            E = x[eval]*gamma_e;
            if(exp(x[eval]/theta)-1!=0.) {
                f[eval]=(r+(2-r)*E_s/E-2*pow(E_s/E,2)-2*E_s/E*log(E_s/E))*x[eval]*x[eval]/(exp(x[eval]/theta)-1)/(gamma_e*m_e);
            } else {
                f[eval] = 0;
            }
            result += dlogx*x[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
            // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" result = " << result << endl;
        }

    }
  }
  else result = 0;
return 3*A*sigma_T*result/(4);

}

double gamma_inverse_compton_analytical_v2(double gamma_e, double E_gamma, double z,Structure_Output_Options * pt_Output_Options)
{
    double T = T_0*(1+z);
    double eta_c = T*E_gamma/(m_e*m_e);
    double eta_0 = E_gamma*E_gamma/(4*gamma_e*m_e*(gamma_e*m_e-E_gamma));
    // cout << "eta_c " << eta_c << "eta_0 " << eta_0 << endl;

    double I = pi*pi/6*eta_c*(exp(-5./4.*pow(eta_0/eta_c,0.5))+2*eta_0*exp(-5./7*pow(eta_0/eta_c,0.7)))*exp(-2./3.*eta_0/eta_c);
    // return sigma_T*T*T*E_gamma/(8*gamma_e*gamma_e);
    if(pt_Output_Options->EM_cascade_verbose > 2) {
        cout << "(gamma_inverse_compton_analytical_v2 : )" <<"z " << z << " gamma_e "<< gamma_e << "E_gamma "<< E_gamma << " result " << 3*sigma_T*T*m_e*m_e/(4*pi*pi*gamma_e*gamma_e)*I << endl;
    }

    return 3*sigma_T*T*m_e*m_e/(4*pi*pi*gamma_e*gamma_e)*I;
}

double gamma_inverse_compton_analytical(double gamma_e,
                                        double E_gamma,
                                        double z,
                                        int N,
                                        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                        Structure_Output_Options * pt_Output_Options)
{
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

    return result;
}







double Analytical_form_inverse_compton(double E_e, double E_gamma_bar, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters)
{
    double gamma_e = E_e/m_e;
    E_gamma_bar =  E_gamma_bar/m_e;
    double s = 4*gamma_e*E_gamma_bar;
    double result;

    // result = (s+9+8./s)*log(1+s)-(16+18*s+s*s)/(2*(1+s))+4*polylog_2(-s,pt_Spectrum_and_Precision_Parameters);
    result = (s+9+8./s)*log(1+s)-8-(2*s+s*s)/(2+2*s)+4*polylog_2(-s,pt_Spectrum_and_Precision_Parameters);
    if(-s > 1) {
        cout << "result " << result << " s " << s << endl;
    }
    return 3*result/(2*s*s);
}

double integrator_simpson_rate_inverse_compton(double z,
        double E_e,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options)
{

    double T=T_0*(1+z);
    double A = 8*pi*pow(m_e/(2*pi),3);
    double E_gamma_bb = 2.701*T;
    double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
    double E_ini = E_gamma_bb/10., E_max = E_gamma_bb*10.;
    double dE = (E_max- E_ini)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1), dlogE;
    double theta = T/m_e;
    double E[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h2;

    int y=0;
    // while(dE>E_ini) {
    //     dE/=10.;
    //     y++;
    // }
    // cout << "(integrator_simpson_rate_inverse_compton :) E_max "<< E_max << " E_ini " << E_ini <<" dE = " << dE << " y = " << y << endl;
    // h2 = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
    dlogE = (log(E_max)-log(E_ini))/(pt_Spectrum_and_Precision_Parameters->n_step-1);
    h2 = dlogE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

    // x[eval]=exp(log(x_cmb_min)+eval*h2+j*dlogx);

    double result = 0;
    // cout << "here "<< dE << endl;

    for(int j=0; j<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1; j++) {

        // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
        for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {

            // E[eval]=E_ini+eval*h2+j*dE;
            E[eval]=exp(log(E_ini)+eval*h2+j*dlogE);

            if((exp(E[eval]/T)-1)!=0.) {
                f[eval]=Analytical_form_inverse_compton(E_e,E[eval],pt_Spectrum_and_Precision_Parameters)*E[eval]*E[eval]/(exp(E[eval]/T)-1);
            } else {
                f[eval] = 0;
            }
            result += dlogE*E[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
            // result += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
            // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" result = " << result << endl;
        }

        // if(res_initial==0 && result !=0)res_initial = result;
        // if(res_initial!=0 && result/res_initial<precision)break;
        // cout << E1 << " f1 " << f1 << " (exp(E1/T)-1) "  << (exp(E1/T)-1) << " f7 " << f7 << endl;


        // 	cout << "Egamma = " << E_gamma << " E7 = " << E7 << " result = " << result<< " j = " << j << endl;
    }

    result*=sigma_T*A/(pow(m_e,3));
    // result*=2*pi*r_e*r_e*m_e*m_e/(E_e*E_e*pi*pi);
    // if(result != 0)cout << "(integrator_simpson_blackbody_spectrum_over_rate_inverse_compton : ) result = " << result << endl;
    if(pt_Output_Options->EM_cascade_verbose > 2) {
        cout << "(integrator_simpson_rate_inverse_compton : )" <<" z " << z << " E_e "<< E_e << " result "<<  result/(pi*pi) << endl;
    }

    return result;
}
double integrator_simpson_rate_inverse_compton_v2(double z,
        double gamma_e,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options)
{

    double T=T_0*(1+z);
    double A = 8*pi*pow(m_e/(2*pi),3);
    double E_gamma_bb = 2.701*T;
    double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
    double gamma_ini = 1, gamma_max = gamma_e;
    double dlogGamma = (log(gamma_max)- log(gamma_ini))/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
    double theta = T/m_e;
    double gamma[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h2;

    int y=0;
    // while(dE>E_ini) {
    //     dE/=10.;
    //     y++;
    // }
    // cout << "(integrator_simpson_rate_inverse_compton :) E_max "<< E_max << " E_ini " << E_ini <<" dE = " << dE << " y = " << y << endl;
    // h2 = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
    dlogGamma = (log(gamma_max)-log(gamma_ini))/(10*pt_Spectrum_and_Precision_Parameters->n_step-1);
    h2 = dlogGamma/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

    // x[eval]=exp(log(x_cmb_min)+eval*h2+j*dlogx);

    double result = 0;
    // cout << "here "<< dE << endl;

    for(int j=0; j<10*pt_Spectrum_and_Precision_Parameters->n_step-1; j++) {

        // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
        for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {

            // E[eval]=E_ini+eval*h2+j*dE;
            gamma[eval]=exp(log(gamma_ini)+eval*h2+j*dlogGamma);
            if(j == 10*pt_Spectrum_and_Precision_Parameters->n_step-2 && eval == pt_Spectrum_and_Precision_Parameters->eval_max-1 ) f[eval]=0;
            else f[eval]=dsigma_inverse_compton_electron_spectrum(z, gamma_e,  gamma[eval],  pt_Spectrum_and_Precision_Parameters, pt_Output_Options);

            result += dlogGamma*gamma[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
            // result += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
            // cout << "eval " << eval << " gamma_e " << gamma_e << "gamma = " << gamma[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" result = " << result << endl;
        }

        // if(res_initial==0 && result !=0)res_initial = result;
        // if(res_initial!=0 && result/res_initial<precision)break;
        // cout << E1 << " f1 " << f1 << " (exp(E1/T)-1) "  << (exp(E1/T)-1) << " f7 " << f7 << endl;


        // 	cout << "Egamma = " << E_gamma << " E7 = " << E7 << " result = " << result<< " j = " << j << endl;
    }

    // result*=2*pi*r_e*r_e*m_e*m_e/(E_e*E_e*pi*pi);
    // if(result != 0)cout << "(integrator_simpson_blackbody_spectrum_over_rate_inverse_compton : ) result = " << result << endl;
    if(pt_Output_Options->EM_cascade_verbose > 1) {
        cout << "(integrator_simpson_rate_inverse_compton v2: )" <<" z " << z << " gamma_e "<< gamma_e << " result "<<  result << endl;
    }

    return result;
}
