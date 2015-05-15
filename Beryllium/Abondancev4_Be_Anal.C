#define FUNC(x,a,i,g,E_O) ((*func)(x,a,i,g,E_O))
#include <iostream>
using namespace std;
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#define EPS 1.0e-2
#define JMAX 20
#include <sstream>
#include "Abondancev4.h"
						/******Programme*******/
/******Calcul graphe tau-zeta pour Beryllium, spectre gaussien******/
/*****ATTENTION - vérifier Abondancev4.h pour savoir quel spectre est utilisé*****/
int main()
{

	float Emin[19];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	//~ for(int i=0;i<14;i++)Emin[i]+=0.;
	float T;
	double T_0 = 2.7255*0.862*pow(10,-10);
	double H_0 = 2.187*pow(10,-18);
	double H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5), Z_x = pow(10,-6) , n_y_0=3.154*pow(10,-30) , E_0 ;
	double B ;
	double z,s_Li,s_Be,s_He,s_De,A,dY_Li,dY_Be,dY_De,dY_He,A_theo;

	float tau, Ex , Ec, s,s2;
	double pi=3.14159;
	int writing_variable;
	double Z_x_old;
	//~ ofstream output_De("results_Graphe_De.dat");
	//~ ofstream output_He("results_Graphe_He.dat");
	ofstream output_Be("results_Graphe_Be_non_standard_2.2MeV.dat");
	// ofstream output_Be("results_Graphe_Be_Anal_02.dat");
	double M_x,E_max;
	M_x=4.4;
	E_0=M_x/2;
	//~ T=1*pow(10,-1);

	tau=pow(10,4);double tau_moy=pow(10,6);double tau_moy2=pow(10,8);double tau_max=pow(10,10);
	//~ double T_min = 0.511*0.511/(1.6*22);
	//~ double T_max = pow(tau_max*(2*H_r),-0.5)*T_0;
	//~ double T_max = 0.511*0.511/(2.2*22);
	//~ tau=pow(T_min/T_0,-2)/(2*H_r);
	//~ tau_max=pow(T_max/T_0,-2)/(2*H_r);
	cout << " tau min = " << tau << " tau max = "<< tau_max << endl;
	double z_max = 35000;
	double h_moy = 500;
	double h_moy2 = 2*h_moy;
	double h_max = 3*h_moy;
	//~ double Delta = (tau_max-tau)/h_moy;
	double Delta = (tau_moy-tau)/h_moy;
	double Delta2 = (tau_moy2-tau_moy)/h_moy;
	double Delta3 = (tau_max-tau_moy2)/h_moy;
	cout << "Delta 1 = " << Delta << " Delta 2 = " << Delta2 << " Delta 3 = " << Delta3 << endl;
	for(int h=1;h<h_max;h++)
	{
		writing_variable=0;
		if(h==h_moy)Delta=Delta2;
		if(h==h_moy2)Delta=Delta3;
		if(h>1)tau+=Delta;
		T=pow(0.2*tau*(2*H_r),-0.5)*T_0;
		z = T/T_0 - 1;


		cout << "****************h = "<<h<<"  tau = "<<tau<<"***************"<<endl;
		//~ cout << "Ex = " << Ex << " Ec = " << Ec << " T= " << T << " tau = " << tau <<" z = "<< z << endl;




		s_Li = 0; s_Be = 0;
		s_De = 0; s_He = 0;
		for(int i=7;i<8;i++)
		{
			s=0;
			//~ cout << "*******processus "<< i <<"*******"<< endl;
			//~ if(Ec>3)Ec=3;

				s=qsimp(func_Anal,z_max,z,z,i,E_0);

			if(i<7)
			{
				s_Li+=s;
			}
			if(i>=7 && i <14)
			{
				s_Be+=s;
			}
			if(i==14)
			{
				s_De+=s;
			}
			if(i>=15)
			{
				s_He+=s;
			}

			//~ cout << "s[" << i <<"] = " << s << " A = " << A << endl;

		}

		for(int l=10;l>1;l--)
		{

			for(double k=1;k<10;k+=0.01)
			{
				Z_x=k*pow(10,-l);
            //  cout << Z_x << endl;
				B=Z_x*n_y_0/(E_0*H_r*tau);
				dY_Be=s_Be*B;
				if(exp(-1*dY_Be)<0.6 && exp(-1*dY_Be)>0.3 && writing_variable == 0 )
				{
					output_Be <<tau <<"  "<<  Z_x << "  " ;
					writing_variable = 1;
				}
				if(exp(-1*dY_Be)<0.6 && exp(-1*dY_Be)>0.3 && writing_variable == 1 )
				{
					Z_x_old = Z_x;
				}
				else if(exp(-1*dY_Be)>0.6 || exp(-1*dY_Be)<0.3)
				{
					if(writing_variable == 1)
					{
						output_Be << Z_x_old << endl;
						writing_variable=2;
						break;
					}
				}
			}
			if(writing_variable==2)
			{

				break;
			}
		}
	}
	output_Be.close();
	return 0;
}
