#include "bbn/EM_cascade.h"
#include "bbn/injected_spectrum.h"
#include "bbn/structures.h"
#include "bbn/BBN_constraints.h"
#include "bbn/tools.h"
#include "bbn/test_functions.h"
#include "bbn/photons.h"
#include "bbn/electrons.h"

using namespace std;




double integrate_dsigma_phph(double E_MIN,
                             double E_MAX,
                             double z,
                             Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                             Structure_Output_Options * pt_Output_Options)
{
  double h, dE, f[pt_Spectrum_and_Precision_Parameters->eval_max],E[pt_Spectrum_and_Precision_Parameters->eval_max], result, result_eval;
  int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
  // double int_BB = 8*pow(pi,4)*pow(T_0*(1+z),6)/63.;
  int y=0;

  dE = (E_MAX-E_MIN)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
  while(dE>E_MIN) {
      dE/=10.;
      y++;
  }
  h = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
  result = 0.;
  {
      int end = pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;
      #pragma omp parallel for ordered schedule(static) private(result_eval,f,E) shared(result)
      for(int i=0; i<end; i++) {
          result_eval = 0.;
          // double  f[pt_Spectrum_and_Precision_Parameters->eval_max],E[pt_Spectrum_and_Precision_Parameters->eval_max];
          for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {

              E[eval]=E_MIN+i*dE+eval*h;

              f[eval]=dsigma_phph(E_MAX,z,E[eval],pt_Output_Options);
              // f[eval]=integrand_rate_pair_creation_v2(E_gamma,E[eval])/(exp(E[eval]/T)-1);
              result_eval += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
              // #pragma  omp critical(print)
              // {
              //   cout << "i " << i <<  " eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" result = " << result_eval << endl;
              // }
          }
          #pragma omp atomic
          result +=result_eval;
      }
  }
  cout << "(integrate_dsigma_phph: ) E_max = " << E_MAX << " result " << result << " ratio = " << rate_gg_scattering(E_MAX,z)/result << endl;

  return result;
}
double integrate_dsigma_compton(double E_MIN,
                                double E_MAX,
                                double z,
                                Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                Structure_Output_Options * pt_Output_Options)
{

    double h, dlogE, f[pt_Spectrum_and_Precision_Parameters->eval_max],E[pt_Spectrum_and_Precision_Parameters->eval_max], result;
    double  result_eval = 0.;
    //  E_MIN = E_MAX*pow(10,-3);
    int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
    int y=3;

    dlogE = (log(E_MAX)-log(E_MIN))/ (double) (pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step);
    // while(dE>E_MIN) {
    //     dE/=10.;
    //     y++;
    // }
    h = dlogE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
    result = 0.;
    {
        int end = pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step;
        #pragma omp parallel for ordered schedule(static) private(result_eval,f,E) shared(result)
        for(int i=0; i<end; i++) {
          result_eval = 0;
            for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {

                E[eval]=exp(log(E_MIN)+i*dlogE+eval*h);

                f[eval]=dsigma_compton(E_MAX,z,E[eval],pt_Output_Options);
                result_eval += dlogE*E[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
                // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << result_eval << endl;
            }
            #pragma omp atomic
            result +=result_eval;
        }
    }
    cout << "(integrate_dsigma_compton: ) E_max = " << E_MAX << " result " << result << " ratio = " << rate_compton(E_MAX,z)/result << endl;

    return result;
}
double integrate_dsigma_NPC(double E_MIN,
                                double E_MAX,
                                double z,
                                Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                Structure_Output_Options * pt_Output_Options)
{

    double h, dlogE, f[pt_Spectrum_and_Precision_Parameters->eval_max],E[pt_Spectrum_and_Precision_Parameters->eval_max], result;
    double  result_eval = 0.;
    //  E_MIN = E_MAX*pow(10,-3);
    int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
    int y=2;

    dlogE = (log(E_MAX)-log(E_MIN))/ (double) (pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step);
    // while(dE>E_MIN) {
    //     dE/=10.;
    //     y++;
    // }
    h = dlogE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
    result = 0.;
    {
        int end = pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step;
        #pragma omp parallel for ordered schedule(static) private(result_eval,f,E) shared(result)
        for(int i=0; i<end; i++) {
          result_eval = 0;
            for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {

                E[eval]=exp(log(E_MIN)+i*dlogE+eval*h);

                f[eval]=dsigma_NPC(E_MAX,z,E[eval],pt_Output_Options);
                result_eval += dlogE*E[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
                // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << result_eval << endl;
            }
            #pragma omp atomic
            result +=result_eval;
        }
    }
    cout << "(integrate_dsigma_NPC: ) E_max = " << E_MAX << " result " << result << " ratio = " << rate_NPC(E_MAX,z)/result << endl;

    return result;
}
double integrate_dsigma_pair_creation(double E_MIN,
                                      double E_MAX,
                                      double z,
                                      Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                      Structure_Output_Options * pt_Output_Options)
{
  // bool flag = false;
  double h, dE, dlogE, f[pt_Spectrum_and_Precision_Parameters->eval_max],E[pt_Spectrum_and_Precision_Parameters->eval_max], result, result_eval;
  int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
  int y=0, count = 0;
  double result_old = 0.0, precision = 1e-6;
  // dE = (E_MAX-E_MIN)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
  dlogE = (log(E_MAX)-log(E_MIN))/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
  // while(dE>E_MIN) {
  //     dE/=10.;
  //     y++;
  // }
  // h = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
  h = dlogE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
  result = 0.;
  {
      int end = pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;
      #pragma omp parallel for ordered schedule(dynamic) private(result_eval) shared(result)
      for(int i=0; i<end; i++){


        double result_eval = 0.;
        double f[pt_Spectrum_and_Precision_Parameters->eval_max],E[pt_Spectrum_and_Precision_Parameters->eval_max];
        // if(flag)continue;
          result_eval = 0.;
          for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {

              E[eval]=exp(log(E_MIN)+i*dlogE+eval*h);

              f[eval]=dsigma_pair_creation_v2(z,E[eval],E_MAX,pt_Spectrum_and_Precision_Parameters,pt_Output_Options);
              // f[eval]=integrand_rate_pair_creation_v2(E_gamma,E[eval])/(exp(E[eval]/T)-1);
              result_eval += dlogE*E[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
              // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
          }
          #pragma omp critical(dataupdate)
          result +=result_eval;



      }
  }
  cout << "(integrate_dsigma_pair_creation: ) E_max = " << E_MAX << " result " << result << " ratio = " << rate_pair_creation_v2(E_MAX,z,pt_Spectrum_and_Precision_Parameters)/result << endl;

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
    if(pt_Output_Options->Test_functions_verbose > 1) {
        #pragma omp critical(print)
        {
            cout << " *** Currently checking energy conservation *** " << endl;
        }
    }
    integrale = 0.;
    double integrale_electrons = 0;
    // double dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table) / (double) (pt_Spectrum_and_Precision_Parameters->n_step - 1);
    double h;
    double dE_2, h2;
    double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7, g1, g2, g3, g4, g5, g6, g7, E_gamma, E_e;
    double rate_E1,rate_E2,rate_E3,rate_E4,rate_E5,rate_E6,rate_E7;
    double resultat_1 = 0, resultat_2 = 0, resultat_3 = 0, resultat_4 = 0, resultat_5 = 0;
    double F1, F2, F3, F4, F5, F6, F7;
    double z = pt_Gamma_Spectrum->redshift;
    double E_c = E_c_0/(1+z);
    double dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table) / (double) (pt_Spectrum_and_Precision_Parameters->n_step - 1);

    int y = 0;
    double E_gamma_bb = 2.701*T_0*(1+z);
    double E_cmb_max = 10*E_gamma_bb;
    double E_cmb_min = E_gamma_bb/100.;
    vector<double> Gamma_Spectrum_Integrated_Over_Kernel;
    vector<double> Gamma_Spectrum_Integrated_Over_Kernel_energy;
    vector<double> Electron_Spectrum_Integrated_Over_Kernel;
    vector<double> Electron_Spectrum_Integrated_Over_Kernel_energy;
    double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
    double dlogE;
    // while(dE > pt_Spectrum_and_Precision_Parameters->E_min_table) {
    //     dE/=10.;
    //     y++;
    // }
    // h = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
    dlogE = (log(pt_Particle_Physics_Model->E_0)-log(pt_Spectrum_and_Precision_Parameters->E_min_table))/(pt_Spectrum_and_Precision_Parameters->n_step-1);
    h = dlogE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
    if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice=="universal") {
        {
            int end = pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;
            int end2 = pt_Spectrum_and_Precision_Parameters->eval_max;
            #pragma omp parallel for ordered schedule(dynamic)

            for(int j=0; j<end; j++) {
                if(pt_Output_Options->Test_functions_verbose > 1) {
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
                    // E[eval]=pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE+eval*h;
                    E[eval]=exp(log(pt_Spectrum_and_Precision_Parameters->E_min_table)+j*dlogE+eval*h);

                    if(E[eval]<pt_Particle_Physics_Model->E_0) {
                        linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E[eval], f[eval]);
                    } else {
                        f[eval]=0;
                    }
                    rate_E = rate_NPC(E[eval],z)+rate_compton(E[eval],z)+rate_gg_scattering(E[eval],z);
                    // rate_E = 1;
                    f[eval]*=rate_E*E[eval];
                    // f[eval]=universal_spectrum(E[eval],   z, pt_Particle_Physics_Model->E_0)*E[eval];
                    resultat_photons += dlogE*E[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
                    // resultat_photons += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];

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
                if(pt_Output_Options->Test_functions_verbose > 1) {
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
                    // E[eval]=pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE+eval*h;
                    E[eval]=exp(log(pt_Spectrum_and_Precision_Parameters->E_min_table)+j*dlogE+eval*h);

                    if(E[eval]<pt_Particle_Physics_Model->E_0) {

                        linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E[eval], f[eval]);

                    } else {
                        f[eval]=0;
                    }
                    // rate_E = rate_NPC(E[eval],z)+rate_compton(E[eval],z)+rate_gg_scattering(E[eval],z);
                    // if(pt_Spectrum_and_Precision_Parameters->double_photon_pair_creation=="yes" && E[eval] >= E_c) {
                    //     rate_E+=rate_pair_creation_v2(E[eval],z,pt_Spectrum_and_Precision_Parameters);
                    // }
                    // rate_E = 1;
                    f[eval]*=E[eval];
                    resultat_photons += dlogE*E[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
                    // resultat_photons += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];

                    if(E[eval]<pt_Particle_Physics_Model->E_0) {

                        linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E[eval], g[eval]);

                    } else {
                        g[eval]=0;
                    }
                    g[eval]*=E[eval];
                    // g[eval]*=E[eval]*(integrator_simpson_rate_inverse_compton_v2(z,E_cmb_min,E_cmb_max,E[eval],pt_Spectrum_and_Precision_Parameters));
                    resultat_electrons += dlogE*E[eval]/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*g[eval];
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

    if(pt_Output_Options->Test_functions_verbose > 0) {
        #pragma omp critical(print)
        {
            cout << "At z = "<< z  << ", the total energy contained in " <<   pt_Gamma_Spectrum->spectrum_name  << "spectrum is " << integrale << " MeV";
            cout << " and in " << pt_Electron_Spectrum->spectrum_name << "spectrum is " << integrale_electrons << " MeV ";
            cout << "for a total of "<< integrale+integrale_electrons <<" MeV, you had injected " << pt_Particle_Physics_Model->E_0 << " MeV." << endl;
        }
    }

}


double print_interaction_rate(double z,
                              double E_MIN,
                              double E_MAX,
                              Structure_Output_Options * pt_Output_Options,
                              Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters)
{

    ostringstream os;
    string name;
    if(pt_Output_Options->interaction_rate_files=="automatic") {
        os << "output/Interaction_Rate_Folder/Interaction_Rate_at_z"<< z <<".dat";
    } else {
        os << pt_Output_Options->interaction_rate_files << ".dat";
    }
    name = os.str();

    ofstream file(name);
    if(file) {
        if(pt_Output_Options->BBN_constraints_verbose > 0) {
            cout << "Printing in file " << name <<  endl;
        }
    } else {
        cout << "I cannot open file " << name << ", please check that the folder exist." << endl;
        exit(0);
    }



    double E_e, E_g;
    double E_c = E_c_0/(1+z);
    double E_phph = m_e*m_e/(T_0*(1+z));
    double E_gamma_bb = 2.701*T_0*(1+z);
    double Rate_photons_E_g = 0, Rate_electrons_E_e = 0;
    double rate_PP= 0, rate_COM= 0, rate_DP= 0, rate_DP_2= 0, rate_NP = 0;
    double PP = 0, CS= 0, ICS_g = 0, ICS_e= 0, NPC= 0, COM= 0, DP= 0;
    for(double i = (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1); i>=0 ; i--) {
        Rate_photons_E_g = 0;
        E_g = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_MAX/pt_Spectrum_and_Precision_Parameters->E_min_table,(double) i/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
        #pragma omp parallel sections // starts a new team
        {
            #pragma omp section

            {
                Rate_electrons_E_e = integrator_simpson_rate_inverse_compton(z,E_g,pt_Spectrum_and_Precision_Parameters,pt_Output_Options);

                if(pt_Output_Options->Test_functions_verbose > 0)
                {
                    #pragma omp critical(print)
                    {
                        cout << "(Rate electrons : ) at E = " << E_e << " tot = " << Rate_electrons_E_e << endl;
                    }
                }

            }



            #pragma omp section
            {
                if(pt_Spectrum_and_Precision_Parameters->pair_creation_in_nuclei == "yes")
                {
                    rate_NP= rate_NPC(E_g,z);
                } else {
                    rate_NP = 0;
                }
                if(pt_Spectrum_and_Precision_Parameters->compton_scattering == "yes")
                {
                    rate_COM = rate_compton(E_g,z);
                } else {
                    rate_COM = 0;
                }

                if(pt_Spectrum_and_Precision_Parameters->double_photon_pair_creation=="yes" && E_g >= E_c)
                {
                    rate_DP=rate_pair_creation_v2(E_g,z,pt_Spectrum_and_Precision_Parameters);
                } else {
                    rate_DP = 0.;
                }
                if(pt_Spectrum_and_Precision_Parameters->photon_photon_diffusion == "yes" && E_g < E_phph )
                {
                    rate_PP=rate_gg_scattering(E_g,z);
                }
                // else rate_PP=rate_gg_scattering(E_phph,z);
                else {
                    rate_PP = 0.;
                }
                Rate_photons_E_g += rate_PP;
                Rate_photons_E_g += rate_NP;
                Rate_photons_E_g += rate_COM;
                Rate_photons_E_g += rate_DP;

                if(pt_Output_Options->Test_functions_verbose > 0)
                {
                    #pragma omp critical(print)
                    {
                        cout << "(Rate photons : ) at E = " << E_g << " rate_NP = " << rate_NP << " rate_COM = " << rate_COM << " rate_PP = " << rate_PP << " rate_DP = " << rate_DP << " tot = " << Rate_photons_E_g << endl;
                    }
                }

            }

        }

        file << E_g << " " << rate_NP << "  " << rate_COM << " " << rate_PP << " " << rate_DP << "  " << Rate_photons_E_g << " " << Rate_electrons_E_e << endl;
        // file << 100/E_g << " " << 2*dsigma_pair_creation(z,100,E_g,pt_Spectrum_and_Precision_Parameters)/rate_DP*E_g/m_e  << endl;
        // cout << 100/E_g << " "  << endl;
        // cout << "check : " << rate_PP/integrate_dsigma_phph(0,E_g,z,pt_Spectrum_and_Precision_Parameters) << "  " << rate_PP/integrator_simpson_dsigma_pair_creation(z,100,E_g,pt_Spectrum_and_Precision_Parameters,pt_Output_Options) << endl;
    }
    return 0;
}
