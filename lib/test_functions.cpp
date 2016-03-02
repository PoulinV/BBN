#include "bbn/EM_cascade.h"
#include "bbn/injected_spectrum.h"
#include "bbn/structures.h"
#include "bbn/BBN_constraints.h"
#include "bbn/tools.h"
#include "bbn/test_functions.h"


using namespace std;




double integrate_dsigma_phph(double E_MIN, double E_MAX, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,Structure_Output_Options * pt_Output_Options)
{

    double  q1, q2, q3, q4, q5, q6, q7,F1, F2, F3, F4, F5, F6, F7, dq, h, rate = 0;
    int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
    dq = (E_MAX-pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
    h = dq/6.;
    q1=pt_Spectrum_and_Precision_Parameters->E_min_table;
    for(int i=0; i<pt_Spectrum_and_Precision_Parameters->n_step-1; i++) {
        if(i!=0) {
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
double integrate_dsigma_compton(double E_MIN, double E_MAX, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,Structure_Output_Options * pt_Output_Options)
{

    double  q1, q2, q3, q4, q5, q6, q7,F1, F2, F3, F4, F5, F6, F7, dq, h, rate = 0;
    int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
    dq = (E_MAX-pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (10000*pt_Spectrum_and_Precision_Parameters->n_step-1);
    h = dq/6.;
    q1=pt_Spectrum_and_Precision_Parameters->E_min_table;
    rate = 0.;

    for(int i=0; i<10000*pt_Spectrum_and_Precision_Parameters->n_step-1; i++) {
        if(i!=0) {
            q1=q7;
        }
        q2=q1+h;
        q3=q2+h;
        q4=q3+h;
        q5=q4+h;
        q6=q5+h;
        q7=q6+h;

        F1 = dsigma_compton(E_MAX,z, q1,pt_Output_Options);
        if(isnan(F1)==1) {
            F1=0.;
        }
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
        Structure_Output_Options * pt_Output_Options)
{

    double T=T_0*(1+z);
    double result;
    double h, f[pt_Spectrum_and_Precision_Parameters->eval_max],E[pt_Spectrum_and_Precision_Parameters->eval_max];
    // double T=T_0*(1+z), E_max = 10*2.701*T;

    double dE = (E_max- E_ini)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
    // cout << "E_max "<< E_max << " E_ini " << E_ini <<" dE = " << dE << endl;
    int y=0;
    while(dE>E_ini) {
        dE/=10.;
        y++;
    }
    h = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
    // cout << "(rate_pair_creation_v2 :) dE " << dE << " E_ini " << E_ini << " y " << y << endl;
    result=0;
    for(int i=0; i<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1; i++) {

        // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
        for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {
            if(eval == 0) {
                if(i==0)	{
                    E[eval]=E_ini;
                } else {
                    E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
                }
            } else {
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
double Function_Integrand_Spectre_Compton_times_bb_spectrum(double E_e, double E_gamma,  double E_gamma_bar)
{
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
    else {
        f=0;
    }
    if(f<0) {
        // cout << "(Function_Integrand_Spectre_Compton : ) Eg= "<<E_gamma<<" Ee =" << E_e << " f = " << f << " q=" << q <<" Gamma_e = " << Gamma_e<< " Ecmb="<< E_gamma_bar<< endl;
        f=0;
    }
    cout << " q = " << q << " Gamma_e = " << Gamma_e << " f = " << f << " E_gamma_bar =" << E_gamma_bar<<endl;
    return f*E_gamma_bar/(exp(E_gamma_bar/0.0001)-1);
}


double  print_func_kawmor(double  z, double  E_0, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters)
{
    ofstream file("spectre_kawmor.dat");
    double  a_pp[3],a_low[3],N_pp[3],N_low[3];
    int h;
    if(T_0*(1+z) == 1e-6) {
        h =0;
    } else if(T_0*(1+z) == 1e-5) {
        h =1;
    } else if(T_0*(1+z) == 1e-4) {
        h =2;
    } else {
        cout << "(print_func_kawmor :) Please choose T = 1e-6, 1e-5, 1e-4 MeV!" << endl;
        exit(0);
    }
    if(E_0 == 1e7) {
        a_pp[0]=-5.10 ;
        a_pp[1]=-5.20 ;
        a_pp[2]=-4.84 ;
        a_low[0]=-1.57 ;
        a_low[1]=-1.34 ;
        a_low[2]=-1.22 ;
        N_pp[0]=6.9*pow(10.,-18) ;
        N_pp[1]=6.0*pow(10.,-18) ;
        N_pp[2]=1.1*pow(10.,-17) ;
        N_low[0]=1.6*pow(10.,8) ;
        N_low[1]=5.4*pow(10.,8) ;
        N_low[2]=1.7*pow(10.,9) ;
    }
    if(E_0 == 1e6) {
        a_pp[0]=-5.07 ;
        a_pp[1]=-5.17 ;
        a_pp[2]=-4.79 ;
        a_low[0]=-1.56 ;
        a_low[1]=-1.34 ;
        a_low[2]=-1.22 ;
        N_pp[0]=6.2*pow(10.,-18) ;
        N_pp[1]=5.5*pow(10.,-18) ;
        N_pp[2]=1.0*pow(10.,-17) ;
        N_low[0]=1.4*pow(10.,8) ;
        N_low[1]=4.9*pow(10.,8) ;
        N_low[2]=1.4*pow(10.,9) ;
    } else if(E_0 == 1e5) {
        a_pp[0]=-5.01 ;
        a_pp[1]=-5.15 ;
        a_pp[2]=-4.74 ;
        a_low[0]=-1.56 ;
        a_low[1]=-1.33 ;
        a_low[2]=-1.22 ;
        N_pp[0]=5.7*pow(10.,-18) ;
        N_pp[1]=5.3*pow(10.,-18) ;
        N_pp[2]=1.1*pow(10.,-17) ;
        N_low[0]=1.4*pow(10.,8) ;
        N_low[1]=4.7*pow(10.,8) ;
        N_low[2]=1.3*pow(10.,9) ;
    } else if(E_0 == 1e4) {
        a_pp[0]=-5.01 ;
        a_pp[1]=-5.12 ;
        a_pp[2]=-4.77 ;
        a_low[0]=-1.56 ;
        a_low[1]=-1.33 ;
        a_low[2]=-1.22 ;
        N_pp[0]=5.7*pow(10.,-18) ;
        N_pp[1]=5.5*pow(10.,-18) ;
        N_pp[2]=9.6*pow(10.,-18) ;
        N_low[0]=1.4*pow(10.,8) ;
        N_low[1]=4.5*pow(10.,8) ;
        N_low[2]=1.3*pow(10.,9) ;
    } else {
        cout << "(print_func_kawmor :) Please choose E_0 = 1e7, 1e6, 1e5 or 1e4 MeV!" << endl;
        exit(0);
    }

    double  f, x ;
    double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));

    double  T = T_0*(1+z);
    for(int i = 0; i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size; i++) {
        x = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_0/pt_Spectrum_and_Precision_Parameters->E_min_table,(double) i/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
        if(x < E_x) {
            f = E_0*N_low[h]*pow(T*pow(10.,-3),-3)*pow(x*pow(10.,-3),a_low[h])*pow(10.,-3);
        } else if(x > E_x && x < E_c) {
            f = E_0*N_pp[h]*pow(T*pow(10.,-3),-6)*pow(x*pow(10.,-3),a_pp[h])*pow(10.,-3);
        } else {
            f = 0;
        }

        file << x << "  " << f*1e-6 << endl;
        cout << x << "  " << f*1e-6 << endl;
    }
    file.close();
    return 0;
}

void check_energy_conservation(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                               Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                               Structure_Output_Options * pt_Output_Options,
                               Structure_Spectrum * pt_Gamma_Spectrum,
                               Structure_Spectrum * pt_Electron_Spectrum,
                               double &integrale)
{
  if(pt_Output_Options->Test_functions_verbose > 1){
      #pragma omp critical(print)
      {
        cout << " *** Currently checking energy conservation *** " << endl;
      }
  }
    integrale = 0.;
    double integrale_electrons = 0;
    double dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table) / (double) (pt_Spectrum_and_Precision_Parameters->n_step - 1);
    double h;
    double dE_2, h2;
    double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7, g1, g2, g3, g4, g5, g6, g7, E_gamma, E_e;
    double rate_E1,rate_E2,rate_E3,rate_E4,rate_E5,rate_E6,rate_E7;
    double resultat_1 = 0, resultat_2 = 0, resultat_3 = 0, resultat_4 = 0, resultat_5 = 0;
    double F1, F2, F3, F4, F5, F6, F7;
    double z = pt_Gamma_Spectrum->redshift;
    double E_c = E_c_0/(1+z);
    int y = 0;
    double E_gamma_bb = 2.701*T_0*(1+z);
    double E_cmb_max = 10*E_gamma_bb;
    double E_cmb_min = E_gamma_bb/100.;
    vector<double> Gamma_Spectrum_Integrated_Over_Kernel;
    vector<double> Gamma_Spectrum_Integrated_Over_Kernel_energy;
    vector<double> Electron_Spectrum_Integrated_Over_Kernel;
    vector<double> Electron_Spectrum_Integrated_Over_Kernel_energy;
    double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
    double dE_3 = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table) / (double) (pt_Spectrum_and_Precision_Parameters->z_step - 1);
    while(dE > pt_Spectrum_and_Precision_Parameters->E_min_table) {
        dE/=10.;
        y++;
    }
    h = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

    if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice=="universal") {
      {
      int end = pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;
      int end2 = pt_Spectrum_and_Precision_Parameters->eval_max;
      #pragma omp parallel for ordered schedule(dynamic)

        for(int j=0; j<end; j++) {
          if(pt_Output_Options->Test_functions_verbose > 1){
            #pragma omp critical(print)
            {
              cout << " step j : " << j << " still " << pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1-j << " to go."<< endl;
            }
          }
          double E[end2],f[end2], rate_E;
          double resultat_photons = 0;
            // cout << "end2 = " << end2 << " h2 " << h2 << endl;
            for(int eval=0; eval < end2; eval++) {
                // if(eval == 0) {
                //     if(j==0)	{
                //         E[eval]=pt_Spectrum_and_Precision_Parameters->E_min_table;
                //     } else {
                //         E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
                //     }
                // } else {
                //     E[eval]=E[0]+eval*h;
                // }
                E[eval]=pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE+eval*h;

                if(E[eval]<pt_Particle_Physics_Model->E_0) {
                    linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E[eval], f[eval]);
                } else {
                    f[eval]=0;
                }
                rate_E = rate_NPC(E[eval],z)+rate_compton(E[eval],z)+rate_gg_scattering(E[eval],z);
                // rate_E = 1;
                f[eval]*=rate_E*E[eval];
                // f[eval]=universal_spectrum(E[eval],   z, pt_Particle_Physics_Model->E_0)*E[eval];
                resultat_photons += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];

            }

            // if(res_initial==0 && resultat !=0)res_initial = resultat;
            // if(res_initial!=0 && resultat/res_initial<precision)break;
            // cout << E1 << " f1 " << f1 << " (exp(E1/T)-1) "  << (exp(E1/T)-1) << " f7 " << f7 << endl;

            #pragma omp critical(dataupdate)
            {
              integrale += resultat_photons;
            }
            // 	cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << endl;
        }

      }
    } else if(pt_Spectrum_and_Precision_Parameters->calculation_mode == "triangular") {
      {
      int end = pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;
      int end2 = pt_Spectrum_and_Precision_Parameters->eval_max;
      #pragma omp parallel for ordered schedule(dynamic)

        for(int j=0; j<end; j++) {

          double E[end2],f[end2],g[end2],rate_E;
          double resultat_photons = 0, resultat_electrons = 0 ;
          if(pt_Output_Options->Test_functions_verbose > 1){
            #pragma omp critical(print)
            {
              cout << " step j : " << j << " still " << pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1-j << " to go."<< endl;
            }
          }
            // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
            for(int eval=0; eval < end2; eval++) {
                // if(eval == 0) {
                //     if(j==0)	{
                //         E[eval]=1;
                //     } else {
                //         E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
                //     }
                // } else {
                //     E[eval]=E[0]+eval*h;
                // }
                E[eval]=pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE+eval*h;

                if(E[eval]<pt_Particle_Physics_Model->E_0) {

                    linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E[eval], f[eval]);

                } else {
                    f[eval]=0;
                }
                // rate_E = rate_NPC(E[eval],z)+rate_compton(E[eval],z)+rate_gg_scattering(E[eval],z);
                // if(pt_Spectrum_and_Precision_Parameters->double_photon_pair_creation=="yes" && E[eval] >= E_c) {
                //     rate_E+=rate_pair_creation_v2(E[eval],z,pt_Spectrum_and_Precision_Parameters);
                // }
                rate_E = 1;
                f[eval]*=rate_E*E[eval];
                resultat_photons += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];

                if(E[eval]<pt_Particle_Physics_Model->E_0) {

                    linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E[eval], g[eval]);

                } else {
                    g[eval]=0;
                }
                g[eval]*=E[eval];
                // g[eval]*=E[eval]*(integrator_simpson_rate_inverse_compton_v2(z,E_cmb_min,E_cmb_max,E[eval],pt_Spectrum_and_Precision_Parameters));
                resultat_electrons += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*g[eval];
                // #pragma omp critical(print)
                // {
                //   cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat_photons << endl;
                // }
            }

            // if(res_initial==0 && resultat !=0)res_initial = resultat;
            // if(res_initial!=0 && resultat/res_initial<precision)break;
            // cout << E1 << " f1 " << f1 << " (exp(E1/T)-1) "  << (exp(E1/T)-1) << " f7 " << f7 << endl;
            #pragma omp critical(dataupdate)
            {
              integrale += resultat_photons;
              integrale_electrons += resultat_electrons;
            }
            // 	cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << endl;
        }
      }



    }

    if(pt_Output_Options->Test_functions_verbose > 0){
      #pragma omp critical(print)
      {
          cout << "At z = "<< z  << ", the total energy contained in " <<   pt_Gamma_Spectrum->spectrum_name  << "spectrum is " << integrale << " MeV";
          cout << " and in " << pt_Electron_Spectrum->spectrum_name << "spectrum is " << integrale_electrons << " MeV ";
          cout << "for a total of "<< integrale+integrale_electrons <<" MeV, you had injected " << pt_Particle_Physics_Model->E_0 << " MeV." << endl;
      }
    }

}
