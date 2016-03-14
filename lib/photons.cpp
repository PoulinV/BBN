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
* @brief The photons module :
* Contains all processes udergone by the photon population.
* @file photons.cpp
*/

double  dsigma_compton(double  x, double  z, double x_prime, Structure_Output_Options * pt_Output_Options)
{
  ///Text book result taken from Kawasaki and Moroi 1995 ApJ 452:506-514.
    double  dsigma = pi*pow(r_e,2)*m_e/pow(x,2)*(x/x_prime+x_prime/x+pow(m_e/x_prime-m_e/x,2)-2*m_e*(1/x_prime-1/x))*n_e*pow(1+z,3);
    return dsigma;

}
double  dsigma_phph(double  x, double  z,  double x_prime, Structure_Output_Options * pt_Output_Options)
{
  /// Taken from Svensson and Zdziarski 1990, ApJ, 349:415-428.

    double  dsigma = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(ALPHA*r_e,2)*pow(m_e,-6)*pow(x,2)*pow(1-x_prime/x+pow(x_prime/x,2),2)/(63*10125*pi);
    return dsigma;

}

double dsigma_NPC(double E_gamma, double z, double E_e, Structure_Output_Options * pt_Output_Options)
{
  ///Taken from Kawasaki and Moroi 1995 ApJ 452:506-514. It has been derived in Berestetskii, Lifschitz and Pitaevskii 1971, Relativistic Quantum Theory (Oxford: Pergamon Press).

    double p = sqrt(E_e*E_e-m_e*m_e);
    double E_pos = E_gamma-E_e-m_e;
    double p_pos = E_pos*E_pos-m_e*m_e;
    if(p_pos>0) {
        p_pos = sqrt(p_pos);
    } else {
        return 0;
    }
    double L = log((E_pos*E_e+p_pos*p+m_e*m_e)/(E_pos*E_e-p_pos*p+m_e*m_e));
    double l = log((E_e+p)/(E_e-p));
    double l_pos = log((E_pos+p_pos)/(E_pos-p_pos));

    double result =  ALPHA*pow(r_e,2)*pow(1+z,3)*eta*n_y_0*p*p_pos/pow(E_gamma,3)
                     *(-4./3-2*E_pos*E_e*(p_pos*p_pos+p*p)/(p*p*p_pos*p_pos)+m_e*m_e*(l*E_pos/pow(p,3)+l_pos*E_e/pow(p_pos,3)-l*l_pos/(p_pos*p))
                       +L*(-8*E_pos*E_e/(3*p_pos*p)+E_gamma*E_gamma/pow(p_pos*p,3)*(pow(E_pos*E_e,2)+pow(p_pos*p,2)-m_e*m_e*E_pos*E_e)-m_e*m_e*E_gamma/(2*p_pos*p)*(l_pos*(E_pos*E_e-pow(p_pos,2))/pow(p_pos,3)+l*(E_pos*E_e-p*p)/pow(p,3))));
// n_H+(Z^2 = 4)*n_He = (1-Y)n_b+Y*n_He = n_b
    return 2*result;
}

double  rate_compton(double  x, double  z)
{
  ///Text book result taken from Kawasaki and Moroi 1995 ApJ 452:506-514.
    double X = 2*x/m_e;
    double sigma_cs = 2*pi*pow(r_e,2)/X*((1-4./X-8./pow(X,2))*log(1+X)+0.5+8./X-1./(2*pow(1+X,2)));
    double Gamma = sigma_cs*n_e*pow(1+z,3);
    return Gamma;

}
double  rate_NPC(double  x, double  z)
{
  ///Taken from Kawasaki and Moroi 1995 ApJ 452:506-514. It has been derived in Maximon 1968, J. Res. NBS, 72(B), 79.
    double	k = x/m_e;
    double  rho = (2*k-4)/(k+2+2*pow(2*k,0.5));
    double  sigma_PCN;
    if(k >= 4) {
        sigma_PCN = ALPHA*pow(r_e,2)
                    *(28./9*log(2*k)-218./27
                      +pow(2/k,2)*(2./3*pow(log(2*k),3)-pow(log(2*k),2)+(6-pi*pi/3)*log(2*k)+2*1.20205+pi*pi/6-7./2.)
                      -pow(2/k,4)*(3./16.*log(2*k)+1./8.)
                      -pow(2/k,6)*(29./2304*log(2*k)-77./13824.));
    } else {
        sigma_PCN = ALPHA*pow(r_e,2)*2*pi/3.*pow((k-2)/k,3)*(1+0.5*rho+23./40.*pow(rho,2)+11./60.*pow(rho,3)+29./960.*pow(rho,4));
    }
    double Gamma = sigma_PCN*pow(1+z,3)*eta*n_y_0;			// n_H+(Z^2 = 4)*n_He = (1-Y)n_b+Y*n_b = n_b
    if(Gamma<0) {
        Gamma = 0;
    }
    return Gamma;
}

double  rate_gg_scattering(double  x, double  z)
{
    /// Taken from Svensson and Zdziarski 1990, ApJ, 349:415-428.
    double  Gamma;
    Gamma = 1946./(50625*pi)*pow(ALPHA*r_e,2)*pow(x,3)*8*pow(pi,4)*pow(T_0*(1+z)/m_e,6)/63.;

    return Gamma;
}



double integrand_dsigma_pair_creation_v2(double x_gamma, double gamma_e, double x_bb)
{
  /// It is taken from Zdziarski 1988 ApJ 335:786-802. To be used with dsigma_pair_creation_v2.

    double gamma_prime = x_gamma - gamma_e;
    double r = 1./2*(gamma_e/gamma_prime+gamma_prime/gamma_e);
    double E_s = x_gamma*x_gamma/(4*gamma_e*gamma_prime);
    double E = x_gamma*x_bb;
    double A = 8*pi*pow(m_e/(2*pi),3);
    double result;

    result = r-(2+r)*E_s/E+2*pow(E_s/E,2)+2*E_s/E*log(E/E_s);
    result *= A*3*sigma_T/(4*E*x_gamma);

    return result;

}

double dsigma_pair_creation_v2(double z,
                               double E_e,
                               double E_gamma,
                               Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                               Structure_Output_Options * pt_Output_Options)
{
  /// It is taken from Zdziarski 1988 ApJ 335:786-802. To be used with integrand_dsigma_pair_creation_v2.

    double gamma_e, gamma_prime, x_gamma;
    double E_gamma_bb,x_cmb_min, x_cmb_max, E_s;
    double x[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h2;
    double T=T_0*(1+z);
    double theta = T/m_e;
    double dx, dlogx;
    double result = 0;
    double res_initial =0 , precision;
    int y=0;

    E_gamma_bb = 2.701*T_0*(1+z);
    x_cmb_max = 1000*E_gamma_bb/m_e;
    x_gamma = E_gamma/m_e;
    gamma_e = E_e/m_e;
    gamma_prime = x_gamma - gamma_e;
    E_s = x_gamma*x_gamma/(4*gamma_e*gamma_prime);
    x_cmb_min = E_s/x_gamma;
    while(x_cmb_min>x_cmb_max) {
        x_cmb_max*=10.;
    }
    // dx = (x_cmb_max- x_cmb_min)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
    // cout << "x_cmb_max "<< x_cmb_max << " x_cmb_min " << x_cmb_min <<" dx = " << dx << endl;

    // if(E_s > 1 ){

    // while(dx>x_cmb_min) {
    //     dx/=10.;
    //     y++;
    // }
    // h2 = dx/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

    dlogx = (log(x_cmb_max)-log(x_cmb_min))/(pt_Spectrum_and_Precision_Parameters->n_step-1);
    h2 = dlogx/(pt_Spectrum_and_Precision_Parameters->eval_max-1);


    // precision = 1e-4;
    for(int j=0; j<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1; j++) {
        // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
        for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {

            // x[eval]=x_cmb_min+eval*h2+j*dx;
            x[eval]=exp(log(x_cmb_min)+eval*h2+j*dlogx);

            if((exp(x[eval]/theta)-1)!=0.) {
                f[eval]=integrand_dsigma_pair_creation_v2(x_gamma,gamma_e,x[eval])*x[eval]*x[eval]/(exp(x[eval]/theta)-1);
            } else {
                f[eval] = 0;
            }
            // result += dx/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
            if(isnan(f[eval])==1)f[eval]=0;
            result += dlogx*x[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
            // cout << "eval " << eval << "E = " << x[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" result = " << result << endl;
        }


        // if(res_initial==0 && result !=0) {
        //     res_initial = result;
        // }
        // if(res_initial!=0 && result/res_initial<precision) {
        //     break;
        // }

    }
    // }
    // else result = 0;

    if(pt_Output_Options->EM_cascade_verbose > 2) {
        cout << "(dsigma_pair_creation_v2 : )" <<" z " << z << " E_e "<< E_e << "E_gamma "<< E_gamma << " result "<<  result << endl;
    }
    return 2*result/(m_e*m_e); //We treat positron and electron on an equal footing, hence the factor 2.
}

double integrand_dsigma_pair_creation(double E_e, double E_gamma, double E_gamma_bar)
{
  /// Taken from Kawasaki and Moroi 1995 ApJ 452:506-514. To be used with dsigma_pair_creation.

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
    if(result < 0) {
        result = 0;
    }
    if(result < 0) {
        cout << "(function_integrand_pair_creation : ) A =" << A <<" B = " << B << " C = " << C << " D = " << D << " result = " << result << " E_e = " << E_e << " E_gamma = " << E_gamma << " E_e_prime = " << E_e_prime << endl;
    }
    return 2*result;			// We treat electron and positron on an equal footing, hence we simply multiply here by 2.
}


double dsigma_pair_creation(double z,
                            double E_e,
                            double E_gamma,
                            Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                            Structure_Output_Options * pt_Output_Options)
{
    /// Diffential pair creation cross section as given in Kawasaki and Moroi 1995 ApJ 452:506-514. To be used with integrand_dsigma_pair_creation.
    double gamma_e, gamma_prime, x_j;
    double E_gamma_bb,E_cmb_min, E_cmb_max, E_s;
    double E[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max],h2;
    double T=T_0*(1+z);
    double dE;
    double result = 0;
    double res_initial =0 , precision;
    int y=0;

    E_gamma_bb = 2.701*T_0*(1+z);
    E_cmb_max = 10*E_gamma_bb;
    x_j = E_gamma/m_e;
    gamma_e = E_e/m_e;
    gamma_prime = x_j - gamma_e;
    E_s = x_j*x_j/(4*gamma_e*gamma_prime);
    E_cmb_min = E_s/x_j*m_e;
    while(E_cmb_min>E_cmb_max) {
        E_cmb_max*=10.;
    }
    dE = (E_cmb_max- E_cmb_min)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
    // cout << "E_cmb_max "<< E_cmb_max << " E_cmb_min " << E_cmb_min <<" dE = " << dE << endl;

    if(E_s > 1 ) {

        while(dE>E_cmb_min) {
            dE/=10.;
            y++;
        }
        // cout << "(dsigma_pair_creation :) after y = "<<y <<" dE = " << dE << endl;
        h2 = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

        // cout << "E_max "<< E_cmb_max << " E_ini " << E_cmb_min <<" dE = " << dE << endl;
        // double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7, h2;

        precision = 1e-4;
        for(int j=0; j<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1; j++) {
            // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
            for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {

                E[eval]=E_cmb_min+eval*h2+j*dE;

                if((exp(E[eval]/T)-1)!=0.) {
                    f[eval]=integrand_dsigma_pair_creation(E_e,E_gamma,E[eval])/(exp(E[eval]/T)-1);
                } else {
                    f[eval] = 0;
                }
                result += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
                // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" result = " << result << endl;
            }


            if(res_initial==0 && result !=0) {
                res_initial = result;
            }
            if(res_initial!=0 && result/res_initial<precision) {
                break;
            }

        }
    } else {
        result = 0;
    }

    if(pt_Output_Options->EM_cascade_verbose > 2) {
        cout << "(dsigma_pair_creation : )" <<" z " << z << " E_e "<< E_e << "E_gamma "<< E_gamma << " result "<<  result/(pi*pi) << endl;
    }
    return  1/4.*r_e*r_e*pow(m_e,4)*result/(E_gamma*E_gamma*E_gamma*pi);
}








double integrand_rate_pair_creation(double s)
{
  /// Taken from Kawasaki and Moroi 1995 ApJ 452:506-514. To be used with rate_pair_creation.

    double beta = pow(1-(4*m_e*m_e/s),0.5);
    double result;
    if(beta>=0) {
        beta = sqrt(beta);
        result =  s*0.5*pi*r_e*r_e*(1-(beta*beta))*((3-pow(beta,4))*log((1+beta)/(1-beta))-2*beta*(2-beta*beta));
    } else {
        beta = 0;
        result = 0;
    }
    return result;
}
double rate_pair_creation(double E_gamma, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters)
{
  /// Pair creation rate as given in Kawasaki and Moroi 1995 ApJ 452:506-514. To be used with integrand_rate_pair_creation.
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
    for(int j = 0 ; j < pt_Spectrum_and_Precision_Parameters->z_step ; j ++) {

        E_gamma_bb = (E_cmb_min)*pow(E_cmb_max/E_cmb_min,(double) j/(pt_Spectrum_and_Precision_Parameters->z_step-1));
        // cout << " E_gamma_bb " << E_gamma_bb << " E_cmb_min " << E_cmb_min << " E_cmb_max " << E_cmb_max << endl;
        f_gamma_bb= E_gamma_bb*E_gamma_bb/(pi*pi)/(exp(E_gamma_bb/T)-1);
        ds = (4*E_gamma_bb*E_gamma - 4*m_e*m_e)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
        y = 0;
        while(ds>4*m_e*m_e) {
            ds/=10.;
            y++;
        }
        h = ds/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

        // cout << " ds = " << ds << "s_min = " << 4*m_e*m_e << " y " << y << endl;
        // ds = (E*E_gamma_bb/(m_e*m_e) - 1)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);

        for(int i=0; i<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step; i++) {

            // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
            for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {
                if(eval == 0) {
                    if(i==0)	{
                        s[eval]=4*m_e*m_e;
                    } else {
                        s[eval]=s[pt_Spectrum_and_Precision_Parameters->eval_max-1];
                    }
                } else {
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
    while(ds>E_cmb_min) {
        ds/=10.;
        y++;
    }
    h = ds/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
    // cout << " dS " << ds << " E_cmb_min " << E_cmb_min << endl;
    resultat=0;
    for(int i=0; i<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1; i++) {

        // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
        for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {
            if(eval == 0) {
                if(i==0)	{
                    s[eval]=E_cmb_min;
                } else {
                    s[eval]=s[pt_Spectrum_and_Precision_Parameters->eval_max-1];
                }
            } else {
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

double integrand_rate_pair_creation_v2(double x_gamma, double x_gamma_bar, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters)
{
   /// Parameterization introduced in Zdziarski 1988 ApJ 335:786-802.

    double v;
    double w;
    double result;

    v = x_gamma*x_gamma_bar-1;

    if(v>=0) {
        w = (pow(1+v,0.5)+pow(v,0.5))/(pow(1+v,0.5)-pow(v,0.5));
        result =(1+2*v+2*v*v)/(1+v)*log(w)-(2*pow(v,0.5)*(1+2*v))/pow(1+v,0.5)-pow(log(w),2)+2*pow(log(1+w),2)+4*polylog_2(1/(1+w),pt_Spectrum_and_Precision_Parameters)-pi*pi/3;
    } else {
        result = 0;
    }
    // cout << " v " << v << " w " << w << " result " << result << endl;
    return result;

}
double rate_pair_creation_v2(double E_gamma, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters)
{
    /// It is taken from Zdziarski 1988 ApJ 335:786-802. To be used with integrand_rate_pair_creation_v2.
    double h, dE, dlogE, E[pt_Spectrum_and_Precision_Parameters->eval_max],f[pt_Spectrum_and_Precision_Parameters->eval_max], T=T_0*(1+z), result;
    double E_gamma_bb = 2.701*T_0*(1+z);
    double E_cmb_min = m_e*m_e/E_gamma;
    double E_cmb_max = 10*E_gamma_bb;
    double A = 8*pi*pow(m_e/(2*pi),3);
    int y;

    dE = (E_cmb_max-E_cmb_min)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
    y = 0;
    // while(dE>E_cmb_min) {
    //     dE/=10.;
    //     y++;
    // }
    // h = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
    dlogE = (log(E_cmb_max)-log(E_cmb_min))/(pt_Spectrum_and_Precision_Parameters->n_step-1);
    h = dlogE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

            // x[eval]=exp(log(x_cmb_min)+eval*h2+j*dlogx);

    // cout << "(rate_pair_creation_v2 :) dE " << dE << " E_cmb_min " << E_cmb_min << " y " << y << endl;
    result=0;
    for(int i=0; i<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1; i++) {

        // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
        for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {
            E[eval]=exp(log(E_cmb_min)+i*dlogE+eval*h);

            f[eval]=integrand_rate_pair_creation_v2(E_gamma/m_e,E[eval]/m_e,pt_Spectrum_and_Precision_Parameters)/(exp(E[eval]/T)-1);
            result += dlogE*E[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
            // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
        }
    }
    if(result < 0 ) {
        result = 0;
    }

    result *=3*sigma_T*A*m_e/(8*E_gamma*E_gamma) ;
    result /= m_e; // Conversion factor
    // cout << "(rate_pair_creation_v2 : ) resultat = " << result<< " energy = " << E_gamma << endl;

    return result;
    // return 3*sigma_T/(E_gamma)*result;
}
