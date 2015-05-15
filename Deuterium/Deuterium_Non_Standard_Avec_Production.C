#define FUNC(x,a,i,g,E_O,Z_x,z_0) ((*func)(x,a,i,g,E_O,Z_x,z_0))
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <string>
using namespace std;
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#define EPS 1.0e-1
#define JMAX 10
#include <sstream>
#include "../Include/Abondancev5_2.h"
/******Programme*******/
/******Calcul graphe tau-zeta pour Deutérium AVEC PRODUCTION via Helium, spectre non standard******/
/*****ATTENTION - vérifier Abondancev5_2.h pour savoir quel spectre est utilisé*****/
/****** Pour les intervalles des abondances, voir arxiv:1307.6955*****/
int main()
{
    double Y_De_max_cyburt =1*pow(10,-4), Y_De_max_standard = 2.61*pow(10,-5),Y_De_min_cyburt=1*pow(10,-5),Y_De_min_standard = 2.45*pow(10,-5) ;
    // double Y_De_max_cyburt =1*pow(10,-4), Y_De_max_standard = 3.48*pow(10,-5),Y_De_min_cyburt=1*pow(10,-5),Y_De_min_standard = 2.56*pow(10,-5) ;
    float Emin[19];
    Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
    //~ for(int i=0;i<14;i++)Emin[i]+=0.;
    float T;
    double T_0 = 2.7255*0.862*pow(10,-10);
    double H_0 = 2.187*pow(10,-18);
    double H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5), Z_x = pow(10,-6) , n_y_0=3.154*pow(10,-30)  ;
    // double omega_b = 0.02273;
    double omega_b = 0.02225;
    double tau_n = 880.3, tau_n_0 = 880.3;
    double z,count1,count2,count3,Y_0 = 2.53*pow(10,-5)*pow(omega_b/0.02273,-1.62)*pow(tau_n/tau_n_0,0.41),Y_De,A,B;
    cout << " Y_2H_0 = " << Y_0 << endl;

    float tau, Ex , Ec, s,s2;
    double pi=3.14159;

    // ofstream output_De("OutputDe/results_Graphe_De_avec_production_avec_reinjection_30MeV_4Mai2015.dat");
    ofstream output_De("OutputDe/results_Graphe_De_avec_production_standard_5Mai2015.dat");
    int x =0;
    Attribution_avec_correction(x);
    cout << " avec_correction " << avec_correction << endl;
    ifstream file0("../TableIntegraleSpectre/results_standard_Destruc4He.dat");

    ifstream file("../TableIntegraleSpectre/results_Reinj_Monochromatique_Destruc4He_30MeV.dat");
    ifstream file2("../TableIntegraleSpectre/results_deux_iterations_Destruc4He_30MeV.dat");
    ifstream file3("../TableIntegraleSpectre/results_trois_iterations_Destruc4He_30MeV.dat");

    // ifstream file4("../TableIntegraleSpectre/results_Reinj_Monochromatique_Destruc2H_30MeV.dat");
    ifstream file5("../TableIntegraleSpectre/results_deux_iterations_Destruc2H_30MeV.dat");
    ifstream file6("../TableIntegraleSpectre/results_trois_iterations_Destruc2H_30MeV_v2.dat");
    // ifstream file7("../TableIntegraleSpectre/results_Reinj_Monochromatique_Produc2H_30MeV.dat");
    ifstream file8("../TableIntegraleSpectre/results_deux_iterations_Produc2H_30MeV.dat");
    ifstream file9("../TableIntegraleSpectre/results_trois_iterations_Produc2H_30MeV_v2.dat");
    ifstream file10("../TableIntegraleSpectre/results_quatre_iterations_Produc2H_30MeV_corrected.dat");
    ifstream file11("../TableIntegraleSpectre/results_quatre_iterations_Destruc2H_30MeV_corrected.dat");
    ifstream file12("../TableIntegraleSpectre/results_quatre_iterations_Destruc4He_30MeV_corrected.dat");
    double M_x,E_0;
    // M_x=60;
    // E_0=M_x/2;
    E_0 = 70.;
    //~ T=1*pow(10,-1);
    fill_table_from_file(file0,vector_z_s4He_destruc_standard,vector_s4He_destruc_standard);
    // for(int i=0;i<vector_z_s4He_destruc_standard.size();i++)cout << vector_z_s4He_destruc_standard[i]<< " "<< vector_s4He_destruc_standard[i] << endl;

      fill_table_from_file(file,vector_z_s4He_destruc,vector_s4He_destruc);
      // for(int i=0;i<vector_z.size();i++)output_Check << vector_z[i]<< " "<< vector_s4He_destruc[i] << endl;
      fill_table_from_file(file2,vector_z_s4He_destruc_deux_iterations,vector_s4He_destruc_deux_iterations);
      fill_table_from_file(file3,vector_z_s4He_destruc_trois_iterations,vector_s4He_destruc_trois_iterations);
      //
        // fill_table_from_file(file4,vector_z_s2H_destruc,vector_s2H_destruc);
        // for(int i=0;i<vector_z.size();i++)output_Check << vector_z[i]<< " "<< vector_s2H_destruc[i] << endl;
        fill_table_from_file(file5,vector_z_s2H_destruc_deux_iterations,vector_s2H_destruc_deux_iterations);
        fill_table_from_file(file6,vector_z_s2H_destruc_trois_iterations,vector_s2H_destruc_trois_iterations);
          // fill_table_from_file(file7,vector_z_s2H_produc,vector_s2H_produc);
          // for(int i=0;i<vector_z.size();i++)output_Check << vector_z[i]<< " "<< vector_s2H_produc[i] << endl;
          fill_table_from_file(file8,vector_z_s2H_produc_deux_iterations,vector_s2H_produc_deux_iterations);
          fill_table_from_file(file9,vector_z_s2H_produc_trois_iterations,vector_s2H_produc_trois_iterations);
          fill_table_from_file(file10,vector_z_s2H_produc_quatre_iterations,vector_s2H_produc_quatre_iterations);
          fill_table_from_file(file11,vector_z_s2H_destruc_quatre_iterations,vector_s2H_destruc_quatre_iterations);
          fill_table_from_file(file12,vector_z_s4He_destruc_quatre_iterations,vector_s4He_destruc_quatre_iterations);
          // for(int i=0;i<vector_z_s2H_destruc_quatre_iterations.size();i++)cout << vector_z_s2H_destruc_quatre_iterations[i]<< " "<< vector_s2H_destruc_quatre_iterations[i] << endl;


    tau=pow(10,4);double tau_moy=pow(10,6);double tau_moy2=pow(10,8);double tau_max=pow(10,10);
    //~ double T_min = 0.511*0.511/(1.6*22);
    double T_max = pow(tau_max*(2*H_r),-0.5)*T_0;
    //~ double T_max = 0.511*0.511/(2.2*22);
    //~ tau=pow(T_min/T_0,-2)/(2*H_r);tau_max=pow(T_max/T_0,-2)/(2*H_r);
    cout << " tau min = " << tau << " tau max = "<< tau_max << endl;
    double z_max = 35000;
    double h_moy = 100;
    double h_moy2 = 2*h_moy;
    double h_max = 3*h_moy;

    double Delta = (tau_moy-tau)/h_moy;
    double Delta2 = (tau_moy2-tau_moy)/h_moy;
    double Delta3 = (tau_max-tau_moy2)/h_moy;
    cout << "Delta 1 = " << Delta << " Delta 2 = " << Delta2 << " Delta 3 = " << Delta3 << endl;
    //~ tau=tau_moy;
    //~ double test;
    double  y, dy;

    for(int h=0;h<h_max;h++)
    {
        if(h==h_moy){Delta=Delta2;tau=tau_moy;}
        if(h==h_moy2){Delta=Delta3;tau=tau_moy2;}
        if(h>1)tau+=Delta;
        T=pow(0.2*tau*(2*H_r),-0.5)*T_0;
        z = T/T_0 - 1;
        count1 = 0;count2=0;count3=0;
        Ex = 0.261121/(80*T_0*(1+z)), Ec = 0.261121/(22*T_0*(1+z));

        cout << "****************h = "<<h<<"  tau = "<<tau<<"***************"<<endl;
        cout << "Ex = " << Ex << " Ec = " << Ec << " T= " << T << " tau = " << tau <<" z = "<< z << endl;
        //~ for(int n = 1 ; n<10;n++)cout << trapzd(func_x,z_max,z,z,n,14,E_0,Z_x)<<endl;
        //~ cout << qsimp(func_x,z_max,z,z,14,E_0,Z_x)<<endl;
        for(int l=11;l>1;l--)
          {
            for(double k=1;k<10;k+=0.1)
              {
                Z_x=k*pow(10,-l);
                A = qsimp_Z(K_Perte_sans_reinj,z_max,z,z,14,E_0,Z_x,z);
                B = qsimp_Z(S_Gain_sans_reinj,z_max,z,z,14,E_0,Z_x,z);
                // A = qsimp_Z(K_Perte_non_standard,z_max,z,z,14,E_0,Z_x,z);
                // B = qsimp_Z(S_Gain_non_standard,z_max,z,z,14,E_0,Z_x,z);
                // cout << " B = " << B << endl;
                // A = qsimp(S_Gain_Check_De,z_max,z,z,14,E_0,Z_x);
                // cout << " A = " << A << endl;
                //~ for(int n=1;n<2;n++)B = trapzd(S_Gain,z_max,z,z,2,14,E_0,Z_x);
                //~ test = trapzd(func_x,z_max,z,z,n,14,E_0,Z_x);

                //~ cout << " A = " << A << " B = " << B << endl;
                //~ cout << "Gain = " << S_Gain(z,z,14,z,E_0,Z_x)<< endl;
                Y_De=exp(-Z_x*n_y_0/(E_0*H_r*tau)*A)*(Y_0 + Z_x*n_y_0/(E_0*H_r*tau)*B);
                // Y_De=exp(-1*Z_x*n_y_0/(E_0*H_r*tau)*A)*Y_0;
                //~ cout << "Y_De = " << Y_De << endl;
                if(Y_De>Y_De_max_standard && count1==0)
                  {
                    // output_De << Z_x <<"  "<< tau <<"  "<<Y_De << "  " << 1 << endl;
                    // cout << " Y_De = " << Y_De << " A = " << A << " B = " << B << "  " <<  1  << endl;
                    output_De <<tau <<"  "<<  Z_x << "  " << endl;
                    count1++;
                    count3++;
                    cout << "count 1 = " << count1<<endl;
                  }
                  if(Y_De<Y_De_max_standard && count1==1)
                    {
                      // output_De << Z_x <<"  "<< tau <<"  "<<Y_De << "  " << 1 << endl;
                      // cout << " Y_De = " << Y_De << " A = " << A << " B = " << B << "  " <<  1  << endl;
                      output_De <<tau <<"  "<<  Z_x << "  " << endl;
                      count1++;
                      cout << "count 1 = " << count1<<endl;
                      count3++;
                    }
                    if(Y_De<Y_De_min_standard&& count2 == 0)
                      {
                        // output_De << Z_x <<"  "<< tau <<"  "<<Y_De << "  " << 2 << endl;
                        // cout << " Y_De = " << Y_De << " A = " << A << " B = " << B << "  " <<  2 <<  endl;
                        output_De <<tau <<"  "<<  Z_x << "  " << endl;
                        count2++;
                        cout << "count 2 = " << count2<<endl;
                        count3++;
                      }
                              if(count3==1){break;}
                      // if(count1==2 && count2==1)break;
                      // else if(count2==1&&count1==0)break;

                  }
                          if(count3 == 1)break;
                    // if(count1==2 && count2==1)break;
                    // else if(count2==1&&count1==0)break;

                  }
          }
    output_De.close();
    //~ output_He.close();
    //~ output_Be.close();
    return 0;
}
