#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include "nrutil.h"
#include "nrutil.c"
#include <iostream>
#include <string>
#include <vector>
double  gamma_tot_compton(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0);
double  gamma_tot_phph(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0);
double  gamma_tot_phph_2(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0);
double  gamma_tot_phph_3(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0);
double spectre_dirac(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0);
double integrale_quatre_interactions(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0);
double integrale_sept_interactions_table(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0);


double v2=1;
int avec_correction;

// g correspond au stockage de x_0 pour la prise en compte des gains dans le taux d'interaction relatif au photon
/*******Routine pour integrer, issu de Numerical Recipes*******/

vector<double>  vector_z_s2H_destruc, vector_z_s2H_produc,vector_z_s3He_destruc,vector_z_s3He_produc,vector_z_s4He_destruc,  vector_s2H_destruc, vector_s2H_produc, vector_s3He_destruc, vector_s3He_produc, vector_s4He_destruc ;
vector<double>  vector_z_s2H_destruc_deux_iterations, vector_z_s2H_produc_deux_iterations,vector_z_s3He_destruc_deux_iterations,vector_z_s3He_produc_deux_iterations,vector_z_s4He_destruc_deux_iterations, vector_s2H_destruc_deux_iterations, vector_s2H_produc_deux_iterations, vector_s3He_destruc_deux_iterations, vector_s3He_produc_deux_iterations, vector_s4He_destruc_deux_iterations ;
vector<double>  vector_z_s2H_destruc_trois_iterations, vector_z_s2H_produc_trois_iterations,vector_z_s3He_destruc_trois_iterations,vector_z_s3He_produc_trois_iterations,vector_z_s4He_destruc_trois_iterations, vector_s2H_destruc_trois_iterations, vector_s2H_produc_trois_iterations, vector_s3He_destruc_trois_iterations, vector_s3He_produc_trois_iterations, vector_s4He_destruc_trois_iterations ;
vector<double>  vector_z_s2H_destruc_quatre_iterations, vector_z_s2H_produc_quatre_iterations,vector_z_s3He_destruc_quatre_iterations,vector_z_s3He_produc_quatre_iterations,vector_z_s4He_destruc_quatre_iterations, vector_s2H_destruc_quatre_iterations, vector_s2H_produc_quatre_iterations, vector_s3He_destruc_quatre_iterations, vector_s3He_produc_quatre_iterations, vector_s4He_destruc_quatre_iterations ;
vector<double>  vector_z_s2H_destruc_standard, vector_z_s2H_produc_standard,vector_z_s3He_destruc_standard,vector_z_s3He_produc_standard,vector_z_s4He_destruc_standard, vector_s2H_destruc_standard, vector_s2H_produc_standard, vector_s3He_destruc_standard, vector_s3He_produc_standard, vector_s4He_destruc_standard ;


vector<double> vector_E, vector_Spectre;
vector<double> vector_E_cinq_iterations, vector_Spectre_cinq_iterations;
vector<double> vector_E_sept_iterations, vector_Spectre_sept_iterations;
vector<double> vector_E_correction, vector_correction;

void Attribution_avec_correction(int &x){
	avec_correction = x;
}
// float * vector_z, * vector_s_4He;
// int lignes=0;
void fill_table_from_file(ifstream &file, vector<double> &vector_z, vector<double> &vector_y){
	double tmp_z, tmp_y;
	// if(file==NULL)cout<<"error cannot open file"<<endl;
	while(file) //Tant qu'on n'est pas à la fin, on lit
	{
	   string ligne;
	   getline(file, ligne);
	   if(ligne[0] == '#' or ligne[0] == '\0') continue;
		file >> tmp_z >>  tmp_y;
		vector_z.push_back(tmp_z);
		vector_y.push_back(tmp_y);
	}
}

void linearint(vector<double> &xa, vector<double> &ya, int n, double x, double &y, double &dy){
				int i,m,ns=1;
				float x0,x1,y0,y1;
				double dift,dif;
				// cout << "xa[1]" << xa[1] << endl;
				dif=fabs(x-xa[1]);
				// cout << "dif = " << dif << endl;
				for (i=1;i<=n;i++) {
						if ( (dift=fabs(x-xa[i])) < dif) {
						ns=i;
						dif=dift;
						y0=ya[i];
						x0=xa[i];
						y=ya[i];
						}
						// if(ns!=i)break;
				}

	if(pow(10,y)!=0){
				if(fabs(x-xa[ns+1])<fabs(x-xa[ns-1])&&pow(10,ya[ns+1])!=0){
					x1=xa[ns+1];
					y1=ya[ns+1];
					if(x1-x0!=0.0)y=y0+(x-x0)/(x1-x0)*fabs(y1-y0);
				}
				else{
						if(pow(10,ya[ns-1])!=0){
						x1=xa[ns-1];
						y1=ya[ns-1];
						if(x1-x0!=0.0)y=y0-(x-x0)/(x1-x0)*fabs(y1-y0);
					}
					else if(pow(10,ya[ns+1])!=0){
						x1=xa[ns+1];
						y1=ya[ns+1];
						if(x1-x0!=0.0)y=y0+(x-x0)/(x1-x0)*fabs(y1-y0);
					}
				}
			}
			// else{cout<<" vivian z " << vector_z[ns] << " y " << y << endl; }

// cout << " y " << y << " (x1-x0) "<< (x1-x0) << endl;
}
void polint(vector<double> &xa, vector<double> &ya, int n, double x, double &y, double &dy){
				int i,m,ns=1;
				float den,dif,dift,ho,hp,w;
				float * d,*e;
				// cout << "xa[1]" << xa[1] << endl;
				dif=fabs(x-xa[1]); e=vector_num_rec(1,n); d=vector_num_rec(1,n);
				// cout << "dif = " << dif << endl;
				for (i=1;i<=n;i++) {
						if ( (dift=fabs(x-xa[i])) < dif) {
						ns=i;
						dif=dift;


						}
						e[i]=ya[i];
				    d[i]=ya[i];

				}
				y=ya[ns--];
				dy=0;
				if(pow(10,y)!=0){
				for (m=1;m<n;m++) {
							for (i=1;i<=n-m;i++) {
								ho=xa[i]-x;
								hp=xa[i+m]-x;
								w=e[i+1]-d[i];
								if ( (den=ho-hp) == 0.0){
									// nrerror("Error in routine polint");
							break;}
								den=w/den;
								d[i]=hp*den;
								e[i]=ho*den;
							}
							if(den==0.0)break;
							y += (dy=(2*ns < (n-m) ? e[ns+1] : d[ns--]));
				}
				}
				free_vector(d,1,n);
				free_vector(e,1,n);
}
double  spec_protheroe(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
double m_e=0.511;
double k;
k=1-4/(3*(log(2*x/(m_e))+1/2));
return k;
}
double  trapzdb(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int n, int i, double  E_0, double  Z_x, double z_0)
{
	double  x,tnma,sum,del;
	static double  s;
	int it,j;
	if (n == 1) {
	return (s = 0.5*(b-a)*(FUNC(a,d,i,a,E_0,Z_x,z_0)+FUNC(b,d,i,a,E_0,Z_x,z_0)));
	}

	else
	{
		for (it = 1,j = 1;j<n-1;j++) it <<=  1;
		tnma = it;
		del = (b-a)/tnma;
		// cout << " del = " << del<< endl;
		x = a+0.5*del;
		for (sum = 0.0,j = 1;j<= it;j++,x+= del) sum +=  FUNC(x,d,i,a,E_0,Z_x,z_0);
		s = 0.5*(s+(b-a)*sum/tnma);
		// cout << " s = " << s <<endl;
		return s;
	}

}
double  trapzd_3(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int n, int i, double  E_0, double  Z_x, double z_0)
{
	double  x,tnm,suma,dela;
	static double  s;
	int it,j;
	if (n == 1) {
		return (s = 0.5*(b-a)*(FUNC(a,d,i,a,E_0,Z_x,z_0)+FUNC(b,d,i,a,E_0,Z_x,z_0)));
	}

	else
		{
			for (it = 1,j = 1;j<n-1;j++) it <<=  1;
				tnm = it;
				dela = (b-a)/tnm;
				// cout << " del = " << del<< endl;
				x = a+0.5*dela;
				for (suma = 0.0,j = 1;j<= it;j++,x+= dela) suma +=  FUNC(x,d,i,a,E_0,Z_x,z_0);
					s = 0.5*(s+(b-a)*suma/tnm);
					// cout << " s = " << s <<endl;
					return s;
				}

			}
double  qsimp_3(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int i, double  E_0, double  Z_x, double z_0)
{

	double  trapzd_3(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int n, int i, double  E_0, double  Z_x, double z_0);
	int j;
	double  s,st,ost = 0.0,os = 0.0;
	for (j = 1;j<= JMAX;j++)
	{
		//~ if(j%10==0) cout << j*100/JMAX << "%" << endl;
		// cout << " j = " << j << endl;

		st = trapzd_3(func,a,b,d,j,i,E_0,Z_x,z_0);
		s = (4.0*st-ost)/3.0;
		// cout << "os = " << os << "s=" << s<<endl;
		if (j > 5)if (fabs(s-os) < EPS*fabs(os) ||(s == 0.0 && os == 0.0)) return s;
		os = s;
		ost = st;

	}
	return s;
}
double  trapzd_5(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int n, int i, double  E_0, double  Z_x, double z_0)
{
	double  x,tnm,suma,dela;
	static double  s;
	int it,j;
	if (n == 1) {
		return (s = 0.5*(b-a)*(FUNC(a,d,i,a,E_0,Z_x,z_0)+FUNC(b,d,i,a,E_0,Z_x,z_0)));
	}

	else
		{
			for (it = 1,j = 1;j<n-1;j++) it <<=  1;
				tnm = it;
				dela = (b-a)/tnm;
				// cout << " del = " << del<< endl;
				x = a+0.5*dela;
				for (suma = 0.0,j = 1;j<= it;j++,x+= dela) suma +=  FUNC(x,d,i,a,E_0,Z_x,z_0);
					s = 0.5*(s+(b-a)*suma/tnm);
					// cout << " s = " << s <<endl;
					return s;
				}

			}
double  qsimp_5(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int i, double  E_0, double  Z_x, double z_0)
{

	double  trapzd_5(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int n, int i, double  E_0, double  Z_x, double z_0);
	int j;
	double  s,st,ost = 0.0,os = 0.0;
	for (j = 1;j<= JMAX;j++)
	{
		//~ if(j%10==0) cout << j*100/JMAX << "%" << endl;
		// cout << " j = " << j << endl;

		st = trapzd_5(func,a,b,d,j,i,E_0,Z_x,z_0);
		s = (4.0*st-ost)/3.0;
		// cout << "os = " << os << "s=" << s<<endl;
		if (j > 5)if (fabs(s-os) < EPS*fabs(os) ||(s == 0.0 && os == 0.0)) return s;
		os = s;
		ost = st;

	}
	return s;
}
double  trapzd_Z(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int n, int i, double  E_0, double  Z_x, double z_0)
{
	double  x,tnm,suma,dela;
	static double  s;
	int it,j;
	if (n == 1) {
		return (s = 0.5*(b-a)*(FUNC(a,d,i,a,E_0,Z_x,z_0)+FUNC(b,d,i,a,E_0,Z_x,z_0)));
	}

	else
		{
			for (it = 1,j = 1;j<n-1;j++) it <<=  1;
				tnm = it;
				dela = (b-a)/tnm;
				// cout << " del = " << del<< endl;
				x = a+0.5*dela;
				for (suma = 0.0,j = 1;j<= it;j++,x+= dela) suma +=  FUNC(x,d,i,a,E_0,Z_x,z_0);
					s = 0.5*(s+(b-a)*suma/tnm);
					// cout << " s = " << s <<endl;
					return s;
				}

			}
double  qsimp_Z(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int i, double  E_0, double  Z_x, double z_0)
{

	double  trapzd_Z(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int n, int i, double  E_0, double  Z_x, double z_0);
	int j;
	double  s,st,ost = 0.0,os = 0.0;
	for (j = 1;j<= JMAX;j++)
	{
		//~ if(j%10==0) cout << j*100/JMAX << "%" << endl;
		// cout << " j = " << j << endl;

		st = trapzd_Z(func,a,b,d,j,i,E_0,Z_x,z_0);
		s = (4.0*st-ost)/3.0;
		// cout << "os = " << os << "s=" << s<<endl;
		if (j > 5)if (fabs(s-os) < EPS*fabs(os) ||(s == 0.0 && os == 0.0)) return s;
		os = s;
		ost = st;

	}
	return s;
}
double  trapzd(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int n, int i, double  E_0, double  Z_x, double z_0)
{
	double  x,tnm,suma,dela;
	static double  s;
	int it,j;
	if (n == 1) {
		return (s = 0.5*(b-a)*(FUNC(a,d,i,a,E_0,Z_x,z_0)+FUNC(b,d,i,a,E_0,Z_x,z_0)));
	}

	else
		{
			for (it = 1,j = 1;j<n-1;j++) it <<=  1;
				tnm = it;
				dela = (b-a)/tnm;
				// cout << " del = " << del<< endl;
				x = a+0.5*dela;
				for (suma = 0.0,j = 1;j<= it;j++,x+= dela) suma +=  FUNC(x,d,i,a,E_0,Z_x,z_0);
					s = 0.5*(s+(b-a)*suma/tnm);
					// cout << " s = " << s <<endl;
					return s;
				}

			}
double  qsimp(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int i, double  E_0, double  Z_x, double z_0)
{

	double  trapzd(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int n, int i, double  E_0, double  Z_x, double z_0);
	int j;
	double  s,st,ost = 0.0,os = 0.0;
	for (j = 1;j<= JMAX;j++)
	{
		//~ if(j%10==0) cout << j*100/JMAX << "%" << endl;
		// cout << " j = " << j << endl;

		st = trapzd(func,a,b,d,j,i,E_0,Z_x,z_0);
		s = (4.0*st-ost)/3.0;
		// cout << "os = " << os << "s=" << s<<endl;
		if (j > 5)if (fabs(s-os) < EPS*fabs(os) ||(s == 0.0 && os == 0.0)) return s;
		os = s;
		ost = st;

	}
	return s;
}
double  qsimp_E(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int i, double  E_0, double  Z_x, double z_0)
{
	double  trapzdb(double  (*func)(double,double,int,double,double,double,double), double  a, double  b, double  d, int n, int i, double  E_0, double  Z_x, double z_0);
	int j;
	double  t,st,ost = 0.0,os = 0.0;
	for (j = 1;j<= JMAX;j++)
		{
			//~ if(j%10==0) cout << j*100/JMAX << "%" << endl;
			// cout << " j = " << j << " JMAX = " << JMAX << endl;

			st = trapzdb(func,a,b,d,j,i,E_0,Z_x,z_0);
			t = (4.0*st-ost)/3.0;
			// cout << "os = " << os << "t=" << t<<endl;
			if (j > 5)if (fabs(t-os) < EPS*fabs(os) ||(t == 0.0 && os == 0.0)) return t;
			os = t;
			ost = st;

		}
		return t;
}
double  dsigma_compton(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  pi = 3.14159;
	double  m_e = 0.511;
	double  r_e = 1.42*pow(10.,-2);
	double  dsigma = pi*pow(r_e,2)*m_e/pow(x,2)*(x/g+g/x+pow(m_e/g-m_e/x,2)-2*m_e*(1/g-1/x));
	return dsigma;
}
double  dsigma_phph(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  dsigma = pow(x,2)*pow(1-g/x+pow(g/x,2),2);
	return dsigma;

}

double  func_gain_cs(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  n_y_0 = 3.154*pow(10.,-30) ;
	double  eta = 6.05*pow(10.,-10) ;
	double  Y = 0.25;
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double  m_e = 0.511;
	double  r_e = 1.42*pow(10.,-2);
	double  a = 0.0073;
	double  K_0 = E_0/(pow(E_x,2)*(2+log(E_c/E_x)));
	double  pi = 3.15159;
		double  f_gain_cs;
	//~ double  f_gain_cs =K_0*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3)*(pow(E_x,2)/(60*x*x)*(((15*x*E_x*(x*x+2*x*m_e-2*m_e*m_e)+12*x*x*m_e*m_e+20*m_e*E_x*E_x*(m_e-2*x)+30*x*pow(E_x,3))/pow(E_x,5) - 15*x*E_c*(x*x+2*x*m_e-2*m_e*m_e)+12*x*x*m_e*m_e+20*m_e*E_c*E_c*(m_e-2*x)+30*x*pow(E_c,3))/pow(E_c,5)) + 2*pow(E_x,3/2)/(315*pow(x,9/2))*((-pow(x,5/2)*(45*x*E_x*(x*x+2*x*m_e-2*m_e*m_e)+35*x*x*m_e*m_e+63*m_e*E_x*E_x*(m_e-2*x)-105*x*pow(E_x,3)))/pow(E_x,9/2)+150*pow(x,2)-36*x*m_e+8*pow(m_e,2)));
	//~ double  f_gain_cs = K_0*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3)*pi*pow(r_e,2)*m_e/(1260*pow(x,5)*pow(E_x,3)*pow(E_c,7))*(-30*pow(x,6)*(7*pow(E_x,7)*E_c+5*E_x*pow(E_c,7))-20*pow(x,5)*m_e*(3*pow(E_x,7)*(7*E_c+3*m_e)+15*E_x*pow(E_c,7)+5*pow(E_c,7)*m_e)+3*pow(x,4)*E_x*E_c*(20*pow(m_e,2)*(7*pow(E_x,6)+5*pow(E_c,6))+168*E_x*E_c*m_e*(pow(E_x,5)+pow(E_c,5))-35*pow(E_x,2)*pow(E_c,2)*(3*pow(E_x,4)+5*pow(E_c,4)))-252*pow(x,3)*pow(E_x,2)*pow(E_c,2)*pow(m_e,2)*(pow(E_x,5)+pow(E_c,5))+1200*pow(E_x,7)*pow(E_c,7)*pow(a/E_x,5/2)-288*pow(E_x,6)*pow(E_c,7)*m_e*pow(a/E_x,3/2)+64*pow(E_x,5)*pow(E_c,7)*pow(m_e,2)*pow(x/E_x,0.5));

	//~ f_gain_cs = K_0*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3)*pi*pow(r_e,2)*m_e/(1260*pow(x,5)*pow(E_x,3)*pow(E_c,7))*(-30*pow(x,6)*(7*pow(E_x,7)*E_c+5*E_x*pow(E_c,7))-20*pow(x,5)*m_e*(3*pow(E_x,7)*(7*E_c+3*m_e)+15*E_x*pow(E_c,7)+5*pow(E_c,7)*m_e)+3*pow(x,4)*E_x*E_c*(20*pow(m_e,2)*(7*pow(E_x,6)+5*pow(E_c,6))+168*E_x*E_c*m_e*(pow(E_x,5)+pow(E_c,5))-35*pow(E_x,2)*pow(E_c,2)*(3*pow(E_x,4)+5*pow(E_c,4)))-252*pow(x,3)*pow(E_x,2)*pow(E_c,2)*pow(m_e,2)*(pow(E_x,5)+pow(E_c,5))+1200*pow(E_x,7)*pow(E_c,7)*pow(x/E_x,5/2)-288*pow(E_x,6)*pow(E_c,7)*m_e*pow(x/E_x,3/2)+64*pow(E_x,5)*pow(E_c,7)*pow(m_e,2)*pow(x/E_x,0.5));

	//~ f_gain_cs = K_0*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3)*pi*pow(r_e,2)*m_e/(1260*pow(x,3))*(-45*pow(x,4)*(7*pow(E_x,4)+pow(E_c,4))/(pow(E_x,2)*pow(E_c,4))-2*pow(x,3)*m_e*(63*pow(E_x,5)*(5*E_c+2*m_e)+45*E_x*pow(E_c,5)+14*pow(E_c,5)*m_e)/(pow(E_x,3)*pow(E_c,5))+6*pow(x,2)*(-35*pow(E_x,2)*(3*pow(E_c,2)-4*E_c*m_e-3*pow(m_e,2))/pow(E_c,4)+15*m_e*m_e/pow(E_c,2)+28*m_e/E_x-35)+64*m_e*m_e*pow(E_x/x,3/2)-84*x*m_e*m_e*(5*pow(E_x,3)+pow(E_c,3))/(E_x*pow(E_c,3))-288*E_x*m_e*pow(E_x/x,0.5)+1200*x*E_x*pow(E_x/x,0.5));

	return f_gain_cs;
}

double  func_gain_phph(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  pi = 3.14159;
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double  m_e = 0.511;
	double  r_e = 1.42*pow(10.,-2);
	double  a = 0.0073;
	double  K_0 = E_0/(pow(E_x,2)*(2+log(E_c/E_x)));
	double  f_gain_phph ;
	//~ f_gain_phph = K_0*35584/(10125*pi)*pow(a*r_e,2)*pow(m_e,-6)*8*pow(pi,4)*pow(T_0*(1+z),6)/63*(pow(x,2)/3*(-pow(x,4)/pow(E_c,3)-3*pow(x,3)/pow(E_c,2)+3*pow(x,2)/E_c+(pow(x,4)+3*pow(x,3)*E_x-3*pow(x,2)*pow(E_x,2)-3*pow(E_x,4))/pow(E_x,3)+6*x*log(E_x/E_c)+3*E_c)+2/15*(23*pow(x,2)+(-3*pow(x,4)-10*pow(x,3)*E_x+15*pow(x,2)*pow(E_x,2)-30*x*pow(E_x,3)+5*pow(E_x,4))*(x/pow(E_x,3/2))/E_x));
	//~ f_gain_phph = 1/(30*E_x*pow(E_c,5))*(-6*pow(x,4)*(pow(E_x,5)+pow(E_c,5))+5*pow(x,3)*(3*pow(E_x,5)*E_c+5*E_x*pow(E_c,5))-30*pow(x,2)*pow(E_x,2)*pow(E_c,2)*(pow(E_x,3)+5*pow(E_c,3))+6*x*pow(E_x,3)*pow(E_c,3)*(pow(E_c,2)*(42*pow(x/E_x,0.5)-25)+5*pow(E_x,2))+10*pow(E_x,4)*pow(E_c,4)*(5*E_c-3*E_x));
	//~ f_gain_phph=K_0*35584/(10125*pi)*pow(a*r_e,2)*pow(m_e,-6)*8*pow(pi,4)*pow(T_0*(1+z),6)/63*f_gain_phph;}
	//~ f_gain_phph=-pow(x,4)*pow(E_x,2)/(3*pow(E_c,3))-pow(x,4)/(15*E_x)+pow(x,3)*pow(E_x,2)/pow(E_c,2)+42*pow(x,3)*pow(E_x/x,3/2)/5+pow(x,3)/3-3*pow(x,2)*pow(E_x,2)/E_c-3*pow(x,2)*E_x-2*x*pow(E_x,2)*log(E_c)-4*x*pow(E_x,2)+2*x*pow(E_x,2)*log(E_x)-pow(E_x,3)+pow(E_x,2)*E_c;
	//~ f_gain_phph=K_0*35584/(10125*pi)*pow(a*r_e,2)*pow(m_e,-6)*8*pow(pi,4)*pow(T_0*(1+z),6)/63*f_gain_phph;
	f_gain_phph=0;
	if(f_gain_phph<0) {f_gain_phph=0;}
	return f_gain_phph;
}

double  gamma_compton(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double Y = 0.25;
	double m_e = 0.511;
	double X = 2*x/m_e;
	double pi = 3.14159;
	double a = 0.0073;
	double r_e = 1.42*pow(10.,-2);
	double sigma_cs = 2*pi*pow(r_e,2)/X*((1-4/X-8/pow(X,2))*log(1+X)+1/2+8/X+1/(2*pow(1+X,2)));
	//~ double  sigma_cs = 2*pi*pow(r_e,2)/X*log(X);
	double Gamma = sigma_cs*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3);
	return Gamma;
}
double  gamma_NPC(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  n_y_0 = 3.154*pow(10.,-30);
	double  eta = 6.05*pow(10.,-10);
	double  Y = 0.25;
	double  m_e = 0.511;
	double  k = x/m_e;
	//~ cout << "x = " << x << endl;
	double  r_e = 1.42*pow(10.,-2);
	double  rho = (2*k-4)/(k+2+2*pow(2*k,0.5));
	double  sigma_PCN;
	double  a = 0.0073;
	double  pi = 3.14159;
	//~ if(k<4) sigma_PCN = a*pow(r_e,2)*2*pi/3*pow((k-2)/k,3)*(1+rho/2+23*pow(rho,2)/40+11*pow(rho,3)/60+20*pow(rho,4)/960);

	//~ else sigma_PCN = a*pow(r_e,2)*(28/9*log(2*k)-218/27+pow(2/k,2)*(2/3*pow(log(2*k),3)-pow(log(2*k),2)+(6-pi*pi/3)*log(2*k)+2*1.2021+pi*pi/6-7/2)-pow(2/k,4)*(3/16*log(2*k)+1/8)-pow(2/k,6)*(29/2304*log(2*k)-77/13824)) ;
	sigma_PCN = a*pow(r_e,2)*(28/9*log(2*k)-218/27) ;
	double  Gamma_2 = sigma_PCN*pow(1+z,3)*eta*n_y_0;
	//~ cout << "sigma_PCN = " << sigma_PCN << endl;
//~ if(Gamma_2!=0) cout << "Gamma NPC = "<< Gamma_2<<endl;
	return Gamma_2;
}

double  gamma_phph(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  m_e = 0.511;
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  a = 0.0073;
double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	//~ double  Gamma = 0.1513*pow(a,4)*m_e*pow(x/m_e,3)*pow(T_0/m_e,6)*pow(1+z,6);
	double  Gamma_3;
	Gamma_3 = 0.1513*pow(a,4)*m_e*pow(x/m_e,3)*pow(T_0*(1+z)/m_e,6);

	return Gamma_3;
}

double  func_spectre(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{

		double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double  K_0 = E_0/(pow(E_x,2)*(2+log(E_c/E_x)));
	double  f;

	if(x < E_x) f = K_0*pow(E_x/x,1.5);
	else if(x > E_x && x < E_c) f =  K_0*pow(E_x/x,2);
	else {f = 0;}
	return f;
}
double  func_spectre2(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{

	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double  K_0 = E_0/(pow(E_x,2)*(2+log(E_c/E_x)));
	double  K_1 = pow(E_0,1/2)/(2*pow(E_x,3/2));
	double  K_2 = E_0/(pow(E_x,2)*(2+log(E_0/E_x)));
	double  f;
	if(E_0<E_x)
	{
		if(x < E_0) f = K_1*pow(E_x/x,1.5);
		else {f = 0;}
		//~ cout << "here 0" <<endl;
	}
	else if(E_0>E_x && E_0<E_c)
	{
		if(x < E_x) f = K_2*pow(E_x/x,1.5);
		else if(x > E_x && x < E_0) f =  K_2*pow(E_x/x,2);
		else {f = 0;}
		//~ cout << "here 1" <<endl;
	}
	else
	{
		if(x < E_x) f = K_0*pow(E_x/x,1.5);
		else if(x > E_x && x < E_c) f =  K_0*pow(E_x/x,2);
		else {f = 0;}
		//~ cout << "here 2" <<endl;
	}
	return f;
}
double  func_spectre_gauss(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
double  T_0 = 2.7255*0.862*pow(10.,-10);
double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
double pi = 3.14159;
double sigma = 0.01;
double E_mean = E_0;
//~ cout << E_mean << endl;
//~ double  K = E_0/(exp(-E_mean*E_mean/(2*pow(sigma,2)))/pow(2*pi,0.5) + E_mean/2*(1+erf(E_mean/(pow(2,0.5)*sigma))));
//~ double  K = ((-exp(-(pow(E_c - E_mean,2)/(2*pow(sigma,2)))) + exp(-(pow(E_mean,2)/(2*pow(sigma,2)))))*sigma)/pow(2*pi,0.5) + 1/2*E_mean*(erf(E_mean/(pow(2,0.5)*sigma)) - erf((-E_c + E_mean)/(pow(2,0.5)*sigma)));
double  K = 2*E_0/(exp(-pow(E_0,2)/(2*pow(sigma,2)))*pow(2/pi,0.5)*sigma+E_0*(1+erf(E_0/(pow(2,0.5)*sigma))));
//~ cout << " K " << K << endl;
return exp(-pow(x-E_mean,2)/(2*pow(sigma,2)))/(K*sigma*sqrt(2*pi));

}

double Spectre_une_interaction_monochromatique(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
// cout << " x " << x << endl;
double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
// cout << " H_r " <<  H_r << endl;

	double Gamma_tot = gamma_NPC(E_0,z,i,g,E_0,Z_x, z_0)+gamma_compton(E_0,z,i,g,E_0,Z_x, z_0)+gamma_phph(E_0,z,i,g,E_0,Z_x, z_0);
	double conversion_s_en_MeV = 1.52*pow(10,15)*pow(10,6);
	double t_inj =  pow((1+z_0),-2)/(2*H_r);
	double tau = 5*t_inj;
	// cout << "1/Gamma_tot = " << 1/(Gamma_tot*tau*conversion_s_en_MeV)<< endl;
	double proba_ph = gamma_phph(E_0,z,i,g,E_0,Z_x, z_0)/Gamma_tot;
	double proba_compton = gamma_compton(E_0,z,i,g,E_0,Z_x, z_0)/Gamma_tot;
	// cout << " proba ph = " << proba_ph << " proba compton = " << proba_compton << endl;
	double E = x/0.511;
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double T = T_0*(1+z);
	// cout << " T " << T << endl;
	double m_e=0.511;
	double conversion = 3.2*pow(10,-20)*pow(1.52*pow(10,15),-1)*pow(10,-6);
	double alpha = 1./137;
	double r_e = alpha/m_e;
	double E_gamma = E_0/m_e;
	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double Y = 0.25;
	double pi = 3.14159;
	double n_e = eta*n_y_0*(1+Y/2);
	double int_BB = 8./63.*pow(pi,4)*pow(T_0*(1+z),6);
	// cout << " BB " << int_BB << endl;
	// double R = 1.83*pow(10,-27)*50*conversion*pow(1+z,6)*pow(E_gamma,3);
	// double spectre_gamma_gamma = R*(20/7)/E_gamma*pow(1-E/E_gamma+pow(E/E_gamma,2),2);
	double spectre_gamma_gamma = 1112./(10125*pi)*pow(alpha*r_e,2)*pow(m_e,-6)*pow(E_0,2)*pow(1-x/E_0+pow(x/E_0,2),2)*int_BB;
	// double spectre_gamma_gamma = 0;

	double spectre_compton = pi*r_e*r_e*m_e*pow(E_0,-2)*(x/E_0+E_0/x+pow(m_e/x-m_e/E_0,2)-2*m_e*(1/x-1/E_0))*n_e*pow(1+z,3);
	// double f = spectre_gamma_gamma/gamma_phph(E_0,z,i,g,E_0,Z_x, z_0);
	double f = (spectre_gamma_gamma+spectre_compton)/Gamma_tot;
	if(x>E_0)f=0;
	return f;
}
double  func_int_phph(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{

double  integrand ;
integrand = Spectre_une_interaction_monochromatique(x,z,k,g,E_0,Z_x, z_0)*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));

//~ else integrand = func_spectre(x,z,k,g,E_0,Z_x, z_0)*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_tot_compton(x,z,k-1,g,E_0,Z_x, z_0)+gamma_NPC(x,z,k-1,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,k-1,g,E_0,Z_x, z_0));
return integrand;
}
double  func_int_compton(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{

double  integrand ;
integrand = Spectre_une_interaction_monochromatique(x,z,k,g,E_0,Z_x, z_0)*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));

//~ else integrand = func_spectre(x,z,k,g,E_0,Z_x, z_0)*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_tot_compton(x,z,k-1,g,E_0,Z_x, z_0)+gamma_NPC(x,z,k-1,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,k-1,g,E_0,Z_x, z_0));
return integrand;

}
double  Spectre_deux_interactions(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{	double  a = 0.0073;
	double  pi = 3.14159;
	double  r_e = 1.42*pow(10.,-2);
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double Y=0.25;
	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double n_e = eta*n_y_0*(1+Y/2);
	double  m_e = 0.511;
	double  Gamma_phph, Gamma_compton;

	//~ double  GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0)/(1+pow(T_0*(1+z),6)*8*pow(pi,4)*35584*pow(a*r_e,2)*pow(m_e,-6)*qsimp(func_int_phph,x,E_c,z,k)/(63*10125*pi*func_spectre(x,z,k,g,E_0,Z_x, z_0)));
	Gamma_phph = qsimp_3(func_int_phph,x,E_0,z,k,E_0,Z_x, z_0);
	Gamma_compton =  qsimp_3(func_int_compton,x,E_0,z,k,E_0,Z_x, z_0);
	double  gain_phph = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(a*r_e,2)*pow(m_e,-6)*Gamma_phph/(63*10125*pi);
	double gain_compton = Gamma_compton*n_e*pow(1+z,3);


	//~ double  gain = 0;
		//~ cout << " gain num = " << gain << endl;
	// double  GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0)-gain;
	//  if(GammaTot<0){cout << "GammaTot = " << GammaTot << " Gamma phph = " <<gamma_phph(x,z,k,g,E_0,Z_x, z_0) << endl;GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0);}
	// cout <<" gain_compton " << gain_compton <<  "gain_phph " << gain_phph << endl;
	return gain_compton+gain_phph;
}
double  func_int_phph_2(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{

double  integrand ;
integrand = Spectre_deux_interactions(x,z,k,g,E_0,Z_x, z_0)*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));

//~ else integrand = func_spectre(x,z,k,g,E_0,Z_x, z_0)*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_tot_compton(x,z,k-1,g,E_0,Z_x, z_0)+gamma_NPC(x,z,k-1,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,k-1,g,E_0,Z_x, z_0));
return integrand;
}
double  func_int_compton_2(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{

double  integrand ;
integrand = Spectre_deux_interactions(x,z,k,g,E_0,Z_x, z_0)*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));

//~ else integrand = func_spectre(x,z,k,g,E_0,Z_x, z_0)*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_tot_compton(x,z,k-1,g,E_0,Z_x, z_0)+gamma_NPC(x,z,k-1,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,k-1,g,E_0,Z_x, z_0));
return integrand;

}
double  Spectre_trois_interactions(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{	double  a = 0.0073;
	double  pi = 3.14159;
	double  r_e = 1.42*pow(10.,-2);
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double Y=0.25;
	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double n_e = eta*n_y_0*(1+Y/2);
	double  m_e = 0.511;
	double  Gamma_phph, Gamma_compton;

	//~ double  GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0)/(1+pow(T_0*(1+z),6)*8*pow(pi,4)*35584*pow(a*r_e,2)*pow(m_e,-6)*qsimp(func_int_phph,x,E_c,z,k)/(63*10125*pi*func_spectre(x,z,k,g,E_0,Z_x, z_0)));
	Gamma_phph = qsimp_E(func_int_phph_2,x,E_0,z,k,E_0,Z_x, z_0);
	Gamma_compton =  qsimp_E(func_int_compton_2,x,E_0,z,k,E_0,Z_x, z_0);
	double  gain_phph = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(a*r_e,2)*pow(m_e,-6)*Gamma_phph/(63*10125*pi);
	double gain_compton = Gamma_compton*n_e*pow(1+z,3);


	//~ double  gain = 0;
		//~ cout << " gain num = " << gain << endl;
	// double  GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0)-gain;
	//  if(GammaTot<0){cout << "GammaTot = " << GammaTot << " Gamma phph = " <<gamma_phph(x,z,k,g,E_0,Z_x, z_0) << endl;GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0);}
	// cout <<" gain_compton " << gain_compton <<  "gain_phph " << gain_phph << endl;
	return (gain_compton+gain_phph);
}
double  func_int_phph_3(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{

double  integrand ;
integrand = Spectre_trois_interactions(x,z,k,g,E_0,Z_x, z_0)*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));

//~ else integrand = func_spectre(x,z,k,g,E_0,Z_x, z_0)*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_tot_compton(x,z,k-1,g,E_0,Z_x, z_0)+gamma_NPC(x,z,k-1,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,k-1,g,E_0,Z_x, z_0));
return integrand;
}
double  func_int_compton_3(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{

double  integrand ;
integrand = Spectre_trois_interactions(x,z,k,g,E_0,Z_x, z_0)*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));

//~ else integrand = func_spectre(x,z,k,g,E_0,Z_x, z_0)*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_tot_compton(x,z,k-1,g,E_0,Z_x, z_0)+gamma_NPC(x,z,k-1,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,k-1,g,E_0,Z_x, z_0));
return integrand;

}
double  Spectre_quatre_interactions(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0){

	double  a = 0.0073;
	double  pi = 3.14159;
	double  r_e = 1.42*pow(10.,-2);
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double Y=0.25;
	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double n_e = eta*n_y_0*(1+Y/2);
	double  m_e = 0.511;
	double  Gamma_phph, Gamma_compton;

	//~ double  GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0)/(1+pow(T_0*(1+z),6)*8*pow(pi,4)*35584*pow(a*r_e,2)*pow(m_e,-6)*qsimp(func_int_phph,x,E_c,z,k)/(63*10125*pi*func_spectre(x,z,k,g,E_0,Z_x, z_0)));
	Gamma_phph = qsimp(func_int_phph_3,x,E_0,z,k,E_0,Z_x, z_0);
	Gamma_compton =  qsimp(func_int_compton_3,x,E_0,z,k,E_0,Z_x, z_0);
	double  gain_phph = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(a*r_e,2)*pow(m_e,-6)*Gamma_phph/(63*10125*pi);
	double gain_compton = Gamma_compton*n_e*pow(1+z,3);


	//~ double  gain = 0;
		//~ cout << " gain num = " << gain << endl;
	// double  GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0)-gain;
	//  if(GammaTot<0){cout << "GammaTot = " << GammaTot << " Gamma phph = " <<gamma_phph(x,z,k,g,E_0,Z_x, z_0) << endl;GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0);}
	// cout <<" gain_compton " << gain_compton <<  "gain_phph " << gain_phph << endl;
	return (gain_compton+gain_phph);

}
double  func_int_phph_4(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0){

double  integrand ;
double y, dy;
polint(vector_E,vector_Spectre,vector_E.size(),x,y,dy);
double Spectre_quatre_interactions=y;
// integrand = Spectre_quatre_interactions*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));
integrand = Spectre_quatre_interactions*dsigma_phph(x,z,k,g,E_0,Z_x, z_0);
// integrand = Spectre_quatre_interactions(x,z,k,g,E_0,Z_x, z_0)*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));

//~ else integrand = func_spectre(x,z,k,g,E_0,Z_x, z_0)*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_tot_compton(x,z,k-1,g,E_0,Z_x, z_0)+gamma_NPC(x,z,k-1,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,k-1,g,E_0,Z_x, z_0));
return integrand;
}
double  func_int_compton_4(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{

double  integrand ;
double y, dy;
polint(vector_E,vector_Spectre,vector_E.size(),x,y,dy);
double Spectre_quatre_interactions=y;
// integrand = Spectre_quatre_interactions*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));
integrand = Spectre_quatre_interactions*dsigma_compton(x,z,k,g,E_0,Z_x, z_0);
// integrand = Spectre_quatre_interactions(x,z,k,g,E_0,Z_x, z_0)*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));

//~ else integrand = func_spectre(x,z,k,g,E_0,Z_x, z_0)*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_tot_compton(x,z,k-1,g,E_0,Z_x, z_0)+gamma_NPC(x,z,k-1,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,k-1,g,E_0,Z_x, z_0));
return integrand;

}
double  Spectre_cinq_interactions(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{	double  a = 0.0073;
	double  pi = 3.14159;
	double  r_e = 1.42*pow(10.,-2);
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double Y=0.25;
	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double n_e = eta*n_y_0*(1+Y/2);
	double  m_e = 0.511;
	double  Gamma_phph, Gamma_compton;

	//~ double  GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0)/(1+pow(T_0*(1+z),6)*8*pow(pi,4)*35584*pow(a*r_e,2)*pow(m_e,-6)*qsimp(func_int_phph,x,E_c,z,k)/(63*10125*pi*func_spectre(x,z,k,g,E_0,Z_x, z_0)));
	Gamma_phph = qsimp_5(func_int_phph_4,x,E_0,z,k,E_0,Z_x, z_0);
	Gamma_compton =  qsimp_5(func_int_compton_4,x,E_0,z,k,E_0,Z_x, z_0);
	double  gain_phph = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(a*r_e,2)*pow(m_e,-6)*Gamma_phph/(63*10125*pi);
	double gain_compton = Gamma_compton*n_e*pow(1+z,3);


	//~ double  gain = 0;
		//~ cout << " gain num = " << gain << endl;
	// double  GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0)-gain;
	//  if(GammaTot<0){cout << "GammaTot = " << GammaTot << " Gamma phph = " <<gamma_phph(x,z,k,g,E_0,Z_x, z_0) << endl;GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0);}
	// cout <<" gain_compton " << gain_compton <<  "gain_phph " << gain_phph << endl;
	return (gain_compton+gain_phph);
}
double  func_int_phph_5(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{

double  integrand ;
double y, dy;
linearint(vector_E_cinq_iterations,vector_Spectre_cinq_iterations,vector_E.size(),x,y,dy);
double Spectre_cinq_interactions=y;
integrand = Spectre_cinq_interactions*dsigma_phph(x,z,k,g,E_0,Z_x, z_0);
// integrand = Spectre_cinq_interactions*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));
// integrand = Spectre_cinq_interactions(x,z,k,g,E_0,Z_x, z_0)*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));

//~ else integrand = func_spectre(x,z,k,g,E_0,Z_x, z_0)*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_tot_compton(x,z,k-1,g,E_0,Z_x, z_0)+gamma_NPC(x,z,k-1,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,k-1,g,E_0,Z_x, z_0));
return integrand;
}
double  func_int_compton_5(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{

double  integrand ;
double y, dy;
linearint(vector_E_cinq_iterations,vector_Spectre_cinq_iterations,vector_E.size(),x,y,dy);
double Spectre_cinq_interactions=y;
integrand = Spectre_cinq_interactions*dsigma_compton(x,z,k,g,E_0,Z_x, z_0);
// integrand = Spectre_cinq_interactions*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));
// integrand = Spectre_cinq_interactions(x,z,k,g,E_0,Z_x, z_0)*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));

//~ else integrand = func_spectre(x,z,k,g,E_0,Z_x, z_0)*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_tot_compton(x,z,k-1,g,E_0,Z_x, z_0)+gamma_NPC(x,z,k-1,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,k-1,g,E_0,Z_x, z_0));
return integrand;

}
double  Spectre_six_interactions(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{	double  a = 0.0073;
	double  pi = 3.14159;
	double  r_e = 1.42*pow(10.,-2);
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double Y=0.25;
	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double n_e = eta*n_y_0*(1+Y/2);
	double  m_e = 0.511;
	double  Gamma_phph, Gamma_compton;

	//~ double  GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0)/(1+pow(T_0*(1+z),6)*8*pow(pi,4)*35584*pow(a*r_e,2)*pow(m_e,-6)*qsimp(func_int_phph,x,E_c,z,k)/(63*10125*pi*func_spectre(x,z,k,g,E_0,Z_x, z_0)));
	Gamma_phph = qsimp_Z(func_int_phph_5,x,E_0,z,k,E_0,Z_x, z_0);
	Gamma_compton =  qsimp_Z(func_int_compton_5,x,E_0,z,k,E_0,Z_x, z_0);
	double  gain_phph = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(a*r_e,2)*pow(m_e,-6)*Gamma_phph/(63*10125*pi);
	double gain_compton = Gamma_compton*n_e*pow(1+z,3);


	//~ double  gain = 0;
		//~ cout << " gain num = " << gain << endl;
	// double  GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0)-gain;
	//  if(GammaTot<0){cout << "GammaTot = " << GammaTot << " Gamma phph = " <<gamma_phph(x,z,k,g,E_0,Z_x, z_0) << endl;GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0);}
	// cout <<" gain_compton " << gain_compton <<  "gain_phph " << gain_phph << endl;
	return (gain_compton+gain_phph);
}
double  func_int_phph_6(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{

double  integrand ;
double y, dy;
// linearint(vector_E_cinq_iterations,vector_Spectre_cinq_iterations,vector_E.size(),x,y,dy);
// double Spectre_cinq_interactions=y;
// integrand = Spectre_cinq_interactions*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));
integrand = Spectre_six_interactions(x,z,k,g,E_0,Z_x, z_0)*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));

//~ else integrand = func_spectre(x,z,k,g,E_0,Z_x, z_0)*dsigma_phph(x,z,k,g,E_0,Z_x, z_0)/(gamma_tot_compton(x,z,k-1,g,E_0,Z_x, z_0)+gamma_NPC(x,z,k-1,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,k-1,g,E_0,Z_x, z_0));
return integrand;
}
double  func_int_compton_6(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{

double  integrand ;
double y, dy;
// linearint(vector_E_cinq_iterations,vector_Spectre_cinq_iterations,vector_E.size(),x,y,dy);
// double Spectre_cinq_interactions=y;
// integrand = Spectre_cinq_interactions*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));
integrand = Spectre_six_interactions(x,z,k,g,E_0,Z_x, z_0)*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,k,g,E_0,Z_x, z_0)+gamma_compton(x,z,k,g,E_0,Z_x, z_0)+gamma_phph(x,z,k,g,E_0,Z_x, z_0));

//~ else integrand = func_spectre(x,z,k,g,E_0,Z_x, z_0)*dsigma_compton(x,z,k,g,E_0,Z_x, z_0)/(gamma_tot_compton(x,z,k-1,g,E_0,Z_x, z_0)+gamma_NPC(x,z,k-1,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,k-1,g,E_0,Z_x, z_0));
return integrand;

}
double  Spectre_sept_interactions(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{	double  a = 0.0073;
	double  pi = 3.14159;
	double  r_e = 1.42*pow(10.,-2);
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double Y=0.25;
	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double n_e = eta*n_y_0*(1+Y/2);
	double  m_e = 0.511;
	double  Gamma_phph, Gamma_compton;

	//~ double  GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0)/(1+pow(T_0*(1+z),6)*8*pow(pi,4)*35584*pow(a*r_e,2)*pow(m_e,-6)*qsimp(func_int_phph,x,E_c,z,k)/(63*10125*pi*func_spectre(x,z,k,g,E_0,Z_x, z_0)));
	Gamma_phph = qsimp_E(func_int_phph_6,x,E_0,z,k,E_0,Z_x, z_0);
	Gamma_compton =  qsimp_E(func_int_compton_6,x,E_0,z,k,E_0,Z_x, z_0);
	double  gain_phph = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(a*r_e,2)*pow(m_e,-6)*Gamma_phph/(63*10125*pi);
	double gain_compton = Gamma_compton*n_e*pow(1+z,3);


	//~ double  gain = 0;
		//~ cout << " gain num = " << gain << endl;
	// double  GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0)-gain;
	//  if(GammaTot<0){cout << "GammaTot = " << GammaTot << " Gamma phph = " <<gamma_phph(x,z,k,g,E_0,Z_x, z_0) << endl;GammaTot = gamma_phph(x,z,k,g,E_0,Z_x, z_0);}
	// cout <<" gain_compton " << gain_compton <<  "gain_phph " << gain_phph << endl;
	return (gain_compton+gain_phph);
}
double  Integrand_gamma_e(double  x, double  z, int i, double  g, double  E_0, double  x_0, double z_0){
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double m_e=0.511;
	// cout << " x_0 " << x_0 << endl;
	double E_gamma_b = 2.701*T_0*(1+z);
	double Gamma_e = 4*E_gamma_b*x_0/(m_e*m_e);
	// double q = x/(Gamma_e*(x_0-x));
	double e_gamma=Gamma_e*x_0*x/(1+Gamma_e*x);
	double F = 0;
	// if(q<=1 && q >=0) F = 2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q);
	F = (2*x*log(x)+(1+2*x)*(1-x)+pow(Gamma_e*x,2)/(2*(1-Gamma_e*x))*(1-x))*pow(x_0-e_gamma,2)*Gamma_e/x_0;
	// cout << " F = " << F << endl;
	return F;
}
double  Integrand_spectre_electron_compton(double  x, double  z, int i, double  g, double  E_0, double  x_0, double z_0){
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	// cout << " H_r " <<  H_r << endl;
// 	double t_inj =  pow((1+z_0),-2)/(2*H_r);
// double t = pow((1+z),-2)/(2*H_r);

		double Gamma_tot = gamma_NPC(E_0,z,i,g,E_0,x_0, z_0)+gamma_compton(E_0,z,i,g,E_0,x_0, z_0)+gamma_phph(E_0,z,i,g,E_0,x_0, z_0);
		double conversion_s_en_MeV = 1.52*pow(10,15)*pow(10,6);
		// double Delta_t = (t-t_inj)*conversion_s_en_MeV;
		double T_0 = 2.7255*0.862*pow(10.,-10);
		double T = T_0*(1+z);
		// cout << " T " << T << endl;
		double m_e=0.511;
		double conversion = 3.2*pow(10,-20)*pow(1.52*pow(10,15),-1)*pow(10,-6);
		double alpha = 1./137;
		double r_e = alpha/m_e;
		double n_y_0 = 3.154*pow(10.,-30) ;
		double eta = 6.05*pow(10.,-10) ;
		double Y = 0.25;
		double pi = 3.14159;
		double n_e = eta*n_y_0*(1+Y/2);
		double spectre_compton = 0;
		double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
		double Gamma_electron;
		double E_gamma_b = 2.701*T_0*(1+z);
		if(E_0+m_e-x>=0 ) {
			Gamma_electron = 2*pi*r_e*r_e*m_e*m_e/(x*x)*qsimp_E(Integrand_gamma_e,0.0001,1,z,i,E_0,x, z_0)*int_bb/E_gamma_b;
			if(Gamma_electron>0)spectre_compton =  pi*r_e*r_e*m_e*pow(E_0,-2)*((E_0+m_e-x)/E_0+E_0/(E_0+m_e-x)+pow(m_e/(E_0+m_e-x)-m_e/E_0,2)-2*m_e*(1/(E_0+m_e-x)-1/E_0))*n_e*pow(1+z,3)/(Gamma_tot*Gamma_electron);
			else spectre_compton =0;
			// cout << "spectre_compton =  " << spectre_compton << endl;
		}
		// cout << "(E_0+m_e-x) =  " << (E_0+m_e-x) << endl;
		// double f = spectre_gamma_gamma/gamma_phph(E_0,z,i,g,E_0,Z_x, z_0);
		double Gamma_e = 4*E_gamma_b*x/(m_e*m_e);
		double q = x_0/(Gamma_e*(x-x_0));
		double F = 0;
		if(q<=1 && q >=0) F = 2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q);
		// cout << " F = " << F << endl;
		double f = 2*spectre_compton*F/(x*x);
		return f;
}
double  Spectre_electron_compton(double  x, double  z, int i, double  g, double  E_0, double  x_0, double z_0){


		double Gamma_tot = gamma_NPC(E_0,z,i,g,E_0,x_0, z_0)+gamma_compton(E_0,z,i,g,E_0,x_0, z_0)+gamma_phph(E_0,z,i,g,E_0,x_0, z_0);
		// double Delta_t = (t-t_inj)*conversion_s_en_MeV;
		double T_0 = 2.7255*0.862*pow(10.,-10);
		double T = T_0*(1+z);
		// cout << " T " << T << endl;
		double m_e=0.511;
		double alpha = 1./137;
		double r_e = alpha/m_e;
		double n_y_0 = 3.154*pow(10.,-30) ;
		double eta = 6.05*pow(10.,-10) ;
		double Y = 0.25;
		double pi = 3.14159;
		double n_e = eta*n_y_0*(1+Y/2);
		double spectre_compton = 0;
		double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
		double Gamma_electron;
		double E_gamma_b = 2.701*T_0*(1+z);
		if(E_0+m_e-x>=0 ) {
			Gamma_electron = 2*pi*r_e*r_e*m_e*m_e/(x*x)*qsimp_E(Integrand_gamma_e,0.0001,1,z,i,E_0,x, z_0)*int_bb/E_gamma_b;
			if(Gamma_electron>0)spectre_compton =  pi*r_e*r_e*m_e*pow(E_0,-2)*((E_0+m_e-x)/E_0+E_0/(E_0+m_e-x)+pow(m_e/(E_0+m_e-x)-m_e/E_0,2)-2*m_e*(1/(E_0+m_e-x)-1/E_0))*n_e*pow(1+z,3)/(Gamma_tot*Gamma_electron);
			else spectre_compton =0;
			// cout << "spectre_compton =  " << spectre_compton << endl;
		}
		// cout << "(E_0+m_e-x) =  " << (E_0+m_e-x) << endl;
		// double f = spectre_gamma_gamma/gamma_phph(E_0,z,i,g,E_0,Z_x, z_0);

		// cout << " F = " << F << endl;
		double f = spectre_compton;
		return f;
}
double Spectre_Gamma_compton(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0){

	double pi = 3.14159;
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double m_e=0.511;
	double alpha = 1./137;
	double r_e = alpha/m_e;
	double E_gamma_b = 2.701*T_0*(1+z);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	// cout << " Ec = " << E_c << endl;
	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	double result_int = 0;

	if(x+m_e<E_0) result_int = qsimp_3(Integrand_spectre_electron_compton,x+m_e,E_0,z,k,E_0,x, z_0);
	double f = 2*pi*r_e*r_e*m_e*m_e*int_bb*result_int/E_gamma_b;
	// cout << " f compton = " << result_int << endl;
	return f;
}
double  Integrand_gamma_e_avec_iterations(double  x, double  z, int i, double  g, double  E_0, double  x_0, double z_0){
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double m_e=0.511;
	double E_gamma_b = 2.701*T_0*(1+z);
	double Gamma_e = 4*E_gamma_b*x_0/(m_e*m_e);
	// double q = x/(Gamma_e*(x_0-x));
	double e_gamma=Gamma_e*x_0*x/(1+Gamma_e*x);
	double F = 0;
	// if(q<=1 && q >=0) F = 2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q);
	F = (2*x*log(x)+(1+2*x)*(1-x)+pow(Gamma_e*x,2)/(2*(1-Gamma_e*x))*(1-x))*pow(x_0-e_gamma,2)*Gamma_e/x_0;
	// cout << " F = " << F << endl;
	return F;
}
double  Integrand_Spectre_avec_iterations(double  x, double  z, int i, double  g, double  E_0, double  x_0, double z_0){
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double m_e=0.511;
	// cout << " x_0 " << x_0 << endl;

	// cout << " x_0 " << x_0 << endl;
	double E_gamma_b = 2.701*T_0*(1+z);
	double E_gamma = x+m_e-x_0;
	double F = 0, f_gamma=0;
	double y,dy;
	// cout << " x " << x << endl;
	double Gamma_tot = gamma_NPC(E_0,z,i,g,E_0,x_0, z_0)+gamma_compton(E_0,z,i,g,E_0,x_0, z_0)+gamma_phph(E_0,z,i,g,E_0,x_0, z_0);

 // 	f_gamma = spectre_dirac(x, z, 1, 1., E_0, 1., z_0);

		f_gamma = integrale_sept_interactions_table(x, z, 1, 1., E_0, 1., z_0);
	// if(q<=1 && q >=0) F = 2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q);
	F = f_gamma*1/(x*x)*(E_gamma/x+x/E_gamma+pow(m_e/E_gamma - m_e/x,2)-2*m_e*(1/E_gamma-1/x));
	if(F!=0){cout << " F = " << F << endl;cout << " x_0 " << x_0 << endl;
}
if(F<0)F=0;
return F;

	// if(x==E_0)F=pow(E_0,-2)*((E_0+m_e-x)/E_0+E_0/(E_0+m_e-x)+pow(m_e/(E_0+m_e-x)-m_e/E_0,2)-2*m_e*(1/(E_0+m_e-x)-1/E_0))/Gamma_tot;
}
double  Integrand_spectre_electron_compton_avec_iterations(double  x, double  z, int i, double  g, double  E_0, double  x_0, double z_0){
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	// cout << " H_r " <<  H_r << endl;
// 	double t_inj =  pow((1+z_0),-2)/(2*H_r);
// double t = pow((1+z),-2)/(2*H_r);

		double Gamma_tot = gamma_NPC(E_0,z,i,g,E_0,x_0, z_0)+gamma_compton(E_0,z,i,g,E_0,x_0, z_0)+gamma_phph(E_0,z,i,g,E_0,x_0, z_0);
		double conversion_s_en_MeV = 1.52*pow(10,15)*pow(10,6);
		// double Delta_t = (t-t_inj)*conversion_s_en_MeV;
		double T_0 = 2.7255*0.862*pow(10.,-10);
		double T = T_0*(1+z);
		// cout << " T " << T << endl;
		double m_e=0.511;
		double conversion = 3.2*pow(10,-20)*pow(1.52*pow(10,15),-1)*pow(10,-6);
		double alpha = 1./137;
		double r_e = alpha/m_e;
		double n_y_0 = 3.154*pow(10.,-30) ;
		double eta = 6.05*pow(10.,-10) ;
		double Y = 0.25;
		double pi = 3.14159;
		double n_e = eta*n_y_0*(1+Y/2);
		double spectre_compton = 0;
		double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
		double Gamma_electron;
		double E_gamma_b = 2.701*T_0*(1+z);
		double Integrale_spectre=0;
		double good_integral =0;
		// cout <<"E_0+m_e-x" << E_0+m_e-x << endl;
		// cout << " E_0 " << E_0 << endl;
		if(E_0+m_e-x>0 ) {
			// cout << " in " << endl;
			Integrale_spectre = qsimp_5(Integrand_Spectre_avec_iterations,x,E_0-1,z,i,E_0,x, z_0);
			Integrale_spectre += pow(E_0,-2)*((E_0+m_e-x)/E_0+E_0/(E_0+m_e-x)+pow(m_e/(E_0+m_e-x)-m_e/E_0,2)-2*m_e*(1/(E_0+m_e-x)-1/E_0))/Gamma_tot;
			good_integral = pow(E_0,-2)*((E_0+m_e-x)/E_0+E_0/(E_0+m_e-x)+pow(m_e/(E_0+m_e-x)-m_e/E_0,2)-2*m_e*(1/(E_0+m_e-x)-1/E_0))/Gamma_tot;
			if(Integrale_spectre!=good_integral)cout << " Integrale_spectre " << Integrale_spectre << " good integral << " << good_integral << "good x_0 << " << x <<endl;
			// Integrale_spectre = qsimp_Z(Integrand_Spectre_avec_iterations,x,E_0,z,i,E_0,x, z_0);
			Gamma_electron = qsimp_Z(Integrand_gamma_e_avec_iterations,0.0001,1,z,i,E_0,x, z_0);
			if(Gamma_electron>0)spectre_compton =  n_e*pow(1+z,3)*E_gamma_b*(x*x)/(2*m_e*int_bb)*Integrale_spectre/Gamma_electron;
			else spectre_compton =0;
			// cout << "spectre_compton =  " << spectre_compton << endl;
		}
		// cout << "(E_0+m_e-x) =  " << (E_0+m_e-x) << endl;
		// double f = spectre_gamma_gamma/gamma_phph(E_0,z,i,g,E_0,Z_x, z_0);
		double Gamma_e = 4*E_gamma_b*x/(m_e*m_e);
		double q = x_0/(Gamma_e*(x-x_0));
		double F = 0;
		if(q<=1 && q >=0) F = 2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q);
		// cout << " F = " << F << endl;
		double f = 2*spectre_compton*F/(x*x);
		return f;
}
double  Spectre_electron_compton_avec_iterations(double  x, double  z, int i, double  g, double  E_0, double  x_0, double z_0){


		double Gamma_tot = gamma_NPC(E_0,z,i,g,E_0,x_0, z_0)+gamma_compton(E_0,z,i,g,E_0,x_0, z_0)+gamma_phph(E_0,z,i,g,E_0,x_0, z_0);
		// double Delta_t = (t-t_inj)*conversion_s_en_MeV;
		double T_0 = 2.7255*0.862*pow(10.,-10);
		double T = T_0*(1+z);
		// cout << " T " << T << endl;
		double m_e=0.511;
		double alpha = 1./137;
		double r_e = alpha/m_e;
		double n_y_0 = 3.154*pow(10.,-30) ;
		double eta = 6.05*pow(10.,-10) ;
		double Y = 0.25;
		double pi = 3.14159;
		double n_e = eta*n_y_0*(1+Y/2);
		double spectre_compton = 0;
		double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
		double Gamma_electron;
		double E_gamma_b = 2.701*T_0*(1+z);
		if(E_0+m_e-x>=0 ) {
			Gamma_electron = 2*pi*r_e*r_e*m_e*m_e/(x*x)*qsimp_E(Integrand_gamma_e_avec_iterations,0.0001,1,z,i,E_0,x, z_0)*int_bb/E_gamma_b;
			if(Gamma_electron>0)spectre_compton =  pi*r_e*r_e*m_e*pow(E_0,-2)*((E_0+m_e-x)/E_0+E_0/(E_0+m_e-x)+pow(m_e/(E_0+m_e-x)-m_e/E_0,2)-2*m_e*(1/(E_0+m_e-x)-1/E_0))*n_e*pow(1+z,3)/(Gamma_tot*Gamma_electron);
			else spectre_compton =0;
			// cout << "spectre_compton =  " << spectre_compton << endl;
		}
		// cout << "(E_0+m_e-x) =  " << (E_0+m_e-x) << endl;
		// double f = spectre_gamma_gamma/gamma_phph(E_0,z,i,g,E_0,Z_x, z_0);

		// cout << " F = " << F << endl;
		double f = spectre_compton;
		return f;
}
double Spectre_Gamma_compton_avec_iterations(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{
	double pi = 3.14159;
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double m_e=0.511;
	double alpha = 1./137;
	double r_e = alpha/m_e;
	double E_gamma_b = 2.701*T_0*(1+z);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	// cout << " Ec = " << E_c << endl;
	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	double result_int = 0;

	if(x+m_e<E_0) result_int = qsimp_3(Integrand_spectre_electron_compton_avec_iterations,x+m_e,E_0,z,k,E_0,x, z_0);

	double f = 2*pi*r_e*r_e*m_e*m_e*int_bb*result_int/E_gamma_b;
	// cout << " f compton = " << result_int << endl;
	return f;
}
double  gamma_tot_compton(double  x, double  z, int k, double  g, double  E_0, double  Z_x, double z_0)
{
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));

	double  n_y_0 = 3.154*pow(10.,-30) ;
	double  eta = 6.05*pow(10.,-10) ;
	double  Y = 0.25;
	double  Gamma;

	//~ double  GammaTot = gamma_compton(x,z,k,g,E_0,Z_x, z_0)/(1+(eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3)*qsimp(func_int_compton,x,E_c,z,k))/func_spectre(x,z,k,g,E_0,Z_x, z_0));

	double  gain = eta*n_y_0*(1+Y/2)*pow(1+z,3)*qsimp(func_int_compton,x,E_c,k,z,E_0,Z_x, z_0)/(func_spectre_gauss(x,z,k,g,E_0,Z_x, z_0)*(1+Y));
	double  GammaTot = gamma_compton(x,z,k,g,E_0,Z_x, z_0)/(1+gain);
	//~ cout << "Gamma tot = " << GammaTot<< " Gamma compton = " << gamma_compton(x,z,k,g,E_0,Z_x, z_0)<<endl;
	return GammaTot;

}

double  func_kawmor(double  x, double  z, int h, double  E_0, double  Z_x, double z_0)
{	//Attention résultat en GEV
	double  a_pp[3],a_low[3],N_pp[3],N_low[3];

	//E0 = 10 TeV
	//~ a_pp[0]=-5.10 ; a_pp[1]=-5.20 ; a_pp[2]=-4.84 ;
	//~ a_low[0]=-1.57 ; a_low[1]=-1.34 ; a_low[2]=-1.22 ;
	//~ N_pp[0]=6.9*pow(10.,-18) ; N_pp[1]=6.0*pow(10.,-18) ; N_pp[2]=1.1*pow(10.,-17) ;
	//~ N_low[0]=1.6*pow(10.,8) ; N_low[1]=5.4*pow(10.,8) ; N_low[2]=1.7*pow(10.,9) ;
	//E0 = 10 GeV
	a_pp[0]=-5.01 ; a_pp[1]=-5.12 ; a_pp[2]=-4.77 ;
	a_low[0]=-1.56 ; a_low[1]=-1.33 ; a_low[2]=-1.22 ;
	N_pp[0]=5.7*pow(10.,-18) ; N_pp[1]=5.5*pow(10.,-18) ; N_pp[2]=9.6*pow(10.,-18) ;
	N_low[0]=1.4*pow(10.,8) ; N_low[1]=4.5*pow(10.,8) ; N_low[2]=1.3*pow(10.,9) ;

	double  f;
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));

	double  T = T_0*(1+z);


	if(x < E_x) f = E_0*N_low[h]*pow(T*pow(10.,-3),-3)*pow(x*pow(10.,-3),a_low[h])*pow(10.,-3);
	else if(x > E_x && x < E_c) f = E_0*N_pp[h]*pow(T*pow(10.,-3),-6)*pow(x*pow(10.,-3),a_pp[h])*pow(10.,-3);
	else f = 0;

	return f;
}



double  func_sigma(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{



/*******fonctions a integrer*******/

/****Processus 0 = 7Li(y,t)4He****/
/****Processus 1 = 7Li(y,n)6Li****/
/****Processus 2 = 7Li(y,2np)4He****/
/****Processus 3 = 7Li(y,2H)4He+n****/
/****Processus 4 = 7Li(y,p)6He, 6He(p,n)6Li***/
/****Processus 5 = 7Li(y,3H)3H+p****/
/****Processus 6 = 7Li(y,3H)3He+n****/
/****Processus 7 = 7Be(y,3He)4He****/
/****Processus 8 = 7Be(y,p)6Li****/
/****Processus 9 = 7Be(y,2pn)4He****/
/****Processus 10 = 7Be(y,2H)4He+p****/
/****Processus 11 = 7Be(y,n)6Be, 6Be(n,4He)2p+n****/
/****Processus 12 = 7Be(y,3He)3He+n****/
/****Processus 13 = 7Be(y,3H)3He+p****/
/****Processus 14 = d(y,n)p****/
/****Processus 15 = 4He(y,p)t****/
/****Processus 16 = 4He(y,n)3He****/
/****Processus 17 = 4He(y,d)d****/
/****Processus 18 = 4He(y,np)d****/
/****Processus 19 = 3He(y,p)d****/
/****Processus 20 = 3He(y,np)p****/

	double  y=0;

	if(i == 0)
	{
	double  Q = 2.467032;
	//y = (0.105*2371*pow(x,-2)*exp(-2.5954*pow(x-Q,-0.5))*exp(-2.056*(x-Q))*(1+2.2875*pow(x-Q,2)-1.1798*pow(x-Q,3)+2.5279*pow(x-Q,4)))*pow(x,d);
	y = 0.057*931.434*pow(x,-2)*exp(-2.59*pow(x-Q,-0.5));
	}
	else if(i == 1)
	{
	double  Q = 7.249962;
	y = (0.176*pow(Q,1.51)*pow(x-Q,0.49)*pow(x,-2)+1205*pow(Q,5.5)*pow(x-Q,5)*pow(x,-10.5)+0.06/(1+pow((x-7.46)/0.188,2)));
	}
	else if(i == 2)
	{
	double  Q = 10.948850;
	y = 122*pow(Q,4)*pow(x-Q,3)*pow(x,-7);
	}
	else if(i == 3)
	{
	double  Q0 = 8.725;
	double  Q1 = 23;
	y = 3.8 * pow(Q0,2.3)*(x-Q0)/pow(x,3.3);
	if(x>=Q1) y += 2.1*pow(Q1,1.5)*(x-Q1)/pow(x,2.5);

	}
	else if(i == 4)
	{
	double  Q = 9.98;
	y = 10.8 *pow(Q,2)*pow(x-Q,1.2)*(x-Q)/pow(x,3.2);
	}
	else if(i == 5)
	{
	double  Q = 22.28;
	y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
	}
	else if(i == 6)
	{
	double  Q = 23.05;
	y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
	}
	else if(i == 7)
	{
		double  Q = 1.586627;
		//y =  (0.504*2371*pow(x,-2)*exp(-5.1909*pow(x-Q,-0.5))*exp(-0.548*(x-Q))*(1-0.428*pow(x-Q,2)+0.543*pow(x-Q,3)-0.115*pow(x-Q,4)));
		if(x>1.59)y = 0.26*931.434*pow(x,-2)*exp(-5.19*pow(x-Q,-0.5));
		else y=0;
	}
	else if(i == 8)
	{
		double  Q = 5.605794;
		y = (32.6*pow(Q,10)*pow((x-Q),2)*pow(x,-12)+2.27*pow(10.,6)*pow(Q,8.8335)*pow((x-Q),13)*pow(x,-21.8335));
	}
	else if(i == 9)
	{
		double  Q = 9.30468;
		y = 133*pow(Q,4)*pow(x-Q,3)*pow(x,-7);
	}
	else if(i == 10)
	{
	double  Q0 = 7.08;
	double  Q1 = 23;
	y = 3.8 * pow(Q0,2.3)*(x-Q0)/pow(x,3.3);
	if(x>=Q1) y += 2.1*pow(Q1,1.5)*(x-Q1)/pow(x,2.5);
	}
	else if(i == 11)
	{
	double  Q = 10.68;
	y = 10.8 *pow(Q,2)*pow(x-Q,1.2)*(x-Q)/pow(x,3.2);
	}
	else if(i == 12)
	{
	double  Q = 22.17;
	y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
	}
	else if(i == 13)
	{
	double  Q = 21.4;
	y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
	}
	else if(i == 14)
	{
	double  Q = 2.224573;
	if(x>=Q)y = 18.75*(pow(pow(Q*(x-Q),0.5)/x,3)+0.007947*pow(pow(Q*(x-Q),0.5)/x,2)*pow(pow(Q,0.5)-pow(0.037,0.5),2)/(x-Q+0.037));
	}
	else if(i == 15)
	{
	double  Q = 19.813852;
	//~ y = 128.9*pow(Q,4.524)*pow(x-Q,2.512)/pow(x,4.524+2.512);
	if(x>=Q)y = 19.5*pow(Q,3.5)*pow(x-Q,1)/pow(x,4.5);
	}
	else if(i == 16)
	{
	double  Q = 20.577615;
	//~ y = 31.68*pow(Q,3.663)*pow(x-Q,1.580)/pow(x,3.663+1.580);
	if(x>=Q)y = 17.1*pow(Q,3.5)*pow(x-Q,1.)/pow(x,4.5);
	}
	else if(i == 17)
	{
	double  Q = 23.846527;
	if(x>=Q)y = 10.7*pow(Q,10.2)*pow(x-Q,3.4)/pow(x,13.6);
	}
	else if(i == 18)
	{
	double  Q = 26.0711;
	if(x>=Q)y = 21.7*pow(Q,4.0)*pow(x-Q,3.0)/pow(x,7.0);
	}
	else if(i == 19)
	{
	double  Q = 5.483485;
	if(x>=Q)y = 8.88*pow(Q,1.75)*pow(x-Q,1.65)/pow(x,3.4);
	}
	else if(i == 20)
	{
	double  Q = 7.718058;
	if(x>=Q)y = 16.7*pow(Q,1.95)*pow(x-Q,2.3)/pow(x,4.25);
	}
	y=y*2.569*pow(10.,-6);//Conversion mb en Mev^-2
	return y;

}
/***************EXTRA FONCTION, pour autres etudes******************/
/*******************************************************************/
double  func_z(double  z,double  T, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  tau =  pow(T/T_0,-2)/(2*H_r);
	double  A;
	A = exp(-1/(2*H_r*tau*(z+1)*(z+1)));
	return A;
}
double func_anal(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double t_inj =  pow((1+z_0),-2)/(2*H_r);
	double tau = 5*t_inj;
	double f=0;

	// if(E_c>2.2) {
	// 	f=exp(-1/(2*H_r*tau*(z+1)*(z+1)))*func_sigma(E_0,z,i,g,E_0)/(gamma_NPC(E_0,z,i,g,E_0)+gamma_compton(E_0,z,i,g,E_0)+gamma_phph(E_0,z,i,g,E_0));
	// 	//~ cout << " je suis la " << endl;
	//  }
	// else
	// {
	// 	f = exp(-1/(2*H_r*tau*(z+1)*(z+1)))*qsimp2(func_standard,1.6,E_c,z,i,E_0);
	// //~ cout << " je suis ici " << endl;
	// //~ cout << f << endl;
	//  }

	if(E_c>2.2)f=exp(-1/(2*H_r*tau*(z+1)*(z+1)))*func_sigma(E_0,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(E_0,z,i,g,E_0,Z_x, z_0)+gamma_compton(E_0,z,i,g,E_0,Z_x, z_0)+gamma_phph(E_0,z,i,g,E_0,Z_x, z_0));


	return f;
}


/*****************************************************************************************/
/************************************Fin Extra Fonctions**********************************/
/*****************************************************************************************/






/*****************************************************************************************/
/***********************Debut Fonction Deuterium Standard AVEC reinj**********************/
/*****************************************************************************************/
double  func_standard(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  f;
	f=func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));

	return f;
}
double  K_Perte(double z, double x, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double E_min = 2.23;
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  n_y_0 = 3.154*pow(10.,-30) ;

	double B= exp(-1/(2*H_r*tau*(z+1)*(z+1)));
	double K, s=0, s_2H = 0, s_2H_diffuse = 0;
	//~ cout << " tau = " << tau << " z = " << z <<endl;
	//~ if(E_c>E_min)s_2H=qsimp(func,E_min,E_c,z,14,E_0,Z_x, z_0);

	/****Fit realise avec E_0 = 5000 MeV soit M_X = 10000 MeV****/
/****	Fit SANS REINJ, ajout s_2H_diffusee pour tenir compte de la reinj après une iteration****/
	if(E_c>E_min){if(z<pow(10,5.91))s_2H = pow(z,-2.50882)*pow(10,36.8106);
	else if(z>pow(10,5.91) && z<pow(10,6.6)) s_2H = pow(z,-2.95593)*pow(10,39.4551);
	else if(z>pow(10,6.53) && z<pow(10,6.8)) s_2H = pow(z,-3.89315)*pow(10,45.5827);
	else if(z>pow(10,6.8) && z<pow(10,7.05)) s_2H = pow(z,-4.68502)*pow(10,50.966);
	else if(z>pow(10,7.05) && z<pow(10,7.20)) s_2H = pow(z,-7.04065)*pow(10,67.5667);
	else if(z>pow(10,7.2) && z<pow(10,7.31)) s_2H = pow(z,-11.7845)*pow(10,101.75);
	else s_2H=0;


	s_2H = E_0*s_2H/5000;
	cout << "s_2H " << s_2H << endl;
}
//
// 	if(E_c>E_min)s_2H_diffuse = qsimp_E(func_standard_deuterium,E_min,E_c,z,i,E_0,Z_x, z_0);
// if(s_2H_diffuse<0)	cout << " s_2H_diffuse " << s_2H_diffuse << endl;
// 	s_2H = s_2H_diffuse+s_2H;



	//~ else if(z>pow(10,7.17) && z<pow(10,7.24)) s_2H = pow(z,-10.2433)*pow(10,90.7182);
	//~ else {if(E_c>E_min)s_2H=qsimp(func,E_min,E_c,z,14,E_0,Z_x, z_0);}
	K = B*s_2H;
	//~ cout << " s = " << s << " s_int = " << s_int << endl;
	//~ cout << " s = " << s << " B = " << B <<  " Y = " << Y << endl;
	return K;
}

double  func_z_v2(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
{

	double  Emin[19];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  A,s_4He,s_2H;
	double y, dy;

	s_4He=0,s_2H=0;
	//~ cout << "E_c = " << E_c <<endl;
	//~ if(E_0 < 3 && E_c>3)E_c == 3;
	//~ for(int k=14;k<19;k++)
	//~ {
	//~ if(E_c>Emin[k])
	//~ {
	//~ s_int+=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);
	//~ cout << " s["<<k<<"] = " << s_int << endl;
	//~ }
	//~ }
	// if(E_c>Emin[15]){if(z<=pow(10,5.26))s_4He+=pow(z,-2.5297)*pow(10,36.9046);
	// 	else if(z>pow(10,5.26) && z<=pow(10,5.62))s_4He+=pow(z,-2.84119)*pow(10,38.5537);
	// 	else if(z>pow(10,5.62) && z<=pow(10,5.81))s_4He+=pow(z,-3.8207)*pow(10,44.0669);
	// 	else if(z>pow(10,5.81) && z<=pow(10,6.))s_4He+=pow(z,-4.78069 )*pow(10,49.6438);
	// 	else if(z>pow(10,6) && z<=pow(10,6.19))s_4He+=pow(z,-6.05892)*pow(10,57.309);
	// 	else if(z>pow(10,6.19) && z<=pow(10,6.39))s_4He+=pow(z,-7.50809)*pow(10,66.2643);
	// 	//~ else if(z>pow(10,6.19) && z<=pow(10,6.25))s_4He=pow(z,-7.50809)*pow(10,66.2643);
	// 	//~ else {for(int k = 15; k<19;k++)if(E_c>Emin[k])s_4He=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);}
	// 	else s_4He+=0;
	// 	s_4He = E_0*s_4He/5000;}
	if(E_c>Emin[15]){
		// for(int k=15;k<19;k++)s_4He_diffuse+=qsimp(func_trois_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
	linearint(vector_z_s4He_destruc_standard,vector_s4He_destruc_standard,vector_z_s4He_destruc_standard.size(),log10(z),y,dy);
		// output_Check << log10(z) << " "<< y << endl;
		s_4He=E_0*pow(10,y)/5000;
		// s_2H_compton = qsimp(func_spectre_compton,E_min,E_0,z,14,E_0,Z_x, z_0);
		// if(s_2H_diffuse <s_2H_compton) cout << "s_2H_compton =  " << s_2H_compton << endl;
		// cout << "s_4He_diffuse quatre iterations=  " << s_4He_diffuse << endl;
	}
		if(E_c>Emin[14]){if(z<pow(10,6.20))s_2H += pow(z,-2.54461)*pow(10,37.2207);
			else if(z>pow(10,6.2) && z<pow(10,6.6)) s_2H += pow(z,-3.14249)*pow(10,40.9313);
			else if(z>pow(10,6.6) && z<pow(10,6.8)) s_2H += pow(z,-4.1445)*pow(10,47.5602);
			else if(z>pow(10,6.8) && z<pow(10,7.)) s_2H += pow(z,-4.83139)*pow(10,52.2337);
			else if(z>pow(10,7) && z<pow(10,7.17)) s_2H += pow(z,-6.55167)*pow(10,64.2476);
			else if(z>pow(10,7.17) && z<pow(10,7.34)) s_2H += pow(z,-10.2433)*pow(10,90.7182);
			else s_2H+=0;
			s_2H = E_0*s_2H/5000;}

			//~ else if(z>pow(10,7.17) && z<pow(10,7.24)) s_2H = pow(z,-10.2433)*pow(10,90.7182);
			//~ else {if(E_c>Emin[14])s_2H=qsimp(func,Emin[14],E_c,z,14,E_0,Z_x, z_0);}
			//~ {
			//~ for(int k=14;k<19;k++)
			//~ {
			//~ if(E_c>Emin[k])
			//~ {
			//~ s+=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);
			//~
			//~ }
			//~ }
			//~ }

			//~ cout << "s int = " << s_int << " s = " << s << endl;
			A = exp(-1/(2*H_r*tau*(z+1)*(z+1)))*(s_4He-s_2H);
			//~ cout << " A = " << A << " z = " << z <<  endl;
			return A;
		}
double  Y_He_K_Perte(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double E_min = 2.23;
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  n_y_0 = 3.154*pow(10.,-30) ;
	double B=(Z_x*n_y_0)/(E_0*H_r*tau);
	double Y_He,s;
	s=qsimp_E(func_z_v2,z,z_0,z_0,14,E_0,Z_x, z_0);
	//~ cout << "z_0 " << z_0 << " z = " << z << endl;
	Y_He = exp(-B*s);
	//cout << " Y = " << Y_He << " B = " << B << " s = " << s << endl;
	return Y_He;
}

double  S_Gain(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
{	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  Emin[19];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  n_y_0 = 3.154*pow(10.,-30) ;
	//~ cout << " tau = " << tau << endl;
	double a=z;
	double B=(Z_x*n_y_0)/(E_0*H_r*tau);
	double Y=0;
	double Y_He_0 = 0.25,s=0;

	//~ for(int k=17;k<19;k++)
	//~ {
	//~ if(E_c>Emin[k])s_int=qsimp_E(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);
	//~ }
	//
	if(E_c>Emin[17]){if(z<=pow(10,5.3))s=pow(z,-2.60454)*pow(10,36.0398);
		else if(z>pow(10,5.3) && z<=pow(10,5.6))s=pow(z,-3.27568)*pow(10,39.5678);
		else if(z>pow(10,5.6) && z<=pow(10,6.15))s=pow(z,-5.83282)*pow(10,53.9883);
		else if(z>pow(10,6.2) && z<=pow(10,6.315))s=pow(z,-33.4332)*pow(10,224.703);
		else s = 0;}
		//~ {
		// for(int k=17;k<19;k++)
		// 	{
		// 		if(E_c>Emin[k])
		// 			{
		// 				if(k==17)s+=2*qsimp(func,Emin[k],E_c,a,k,E_0,Z_x, z_0);
		// 				if(k==18)s+=qsimp(func,Emin[k],E_c,a,k,E_0,Z_x, z_0);
		// 				// cout << " s = " << s << endl;
		// 			}
		// 		}
		//~ }
		//~
		//~ cout << " s = " << s << " s_int = " << s_int << endl;
		// Y =exp(-1/(2*H_r*pow(1+z,2)*tau))*s*Y_He_0;
		Y =exp(-1/(2*H_r*pow(1+z,2)*tau))*Y_He_K_Perte(a,z_0,i,z_0,E_0,Z_x, z_0)*s*Y_He_0;
		//~ Y =-B*exp(-1/(2*H_r*pow(1+z,2)*tau))*s*Y_He_0;
		//~ Y = 1;
		//~ cout << " Y = " << Y << "z = " << z << endl;
		//~ cout << " Y_He_K_Perte(z,z_0,i,z_0,E_0,Z_x, z_0) = " << Y_He_K_Perte(z,z_0,i,z_0,E_0,Z_x, z_0) << endl;
		return Y;

	}


/*****************************************************************************************/
/***********************Fin Fonction Deuterium Standard AVEC reinj************************/
/*****************************************************************************************/














/*****************************************************************************************/
/***********************Debut Fonction Deuterium Standard SANS reinj**********************/
/*****************************************************************************************/

double  func_sans_reinj(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  f;
	f=func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
	// f=func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+spec_protheroe(x,z,i,g,E_0,Z_x, z_0)*gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
	return f;
}
double  K_Perte_sans_reinj(double z, double x, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double E_min = 2.23;
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  n_y_0 = 3.154*pow(10.,-30) ;



	double B= exp(-1/(2*H_r*tau*(z+1)*(z+1)));
	double K, s_int=0, s_2H = 0;
	double s_2H_table = 0, tmp_z, tmp_sd_p, count_sd_p = 0;
	// if(E_c>E_min){s_2H=qsimp_E(func_sans_reinj,E_min,E_c,z,14,E_0,Z_x, z_0);}
	/****Fit realise avec E_0 = 5000 MeV soit M_X = 10000 MeV****/
	if(E_c>E_min){if(z<pow(10,5.91))s_2H = pow(z,-2.50882)*pow(10,36.8106);
		else if(z>pow(10,5.91) && z<pow(10,6.6)) s_2H = pow(z,-2.95593)*pow(10,39.4551);
		else if(z>pow(10,6.53) && z<pow(10,6.8)) s_2H = pow(z,-3.89315)*pow(10,45.5827);
		else if(z>pow(10,6.8) && z<pow(10,7.05)) s_2H = pow(z,-4.68502)*pow(10,50.966);
		else if(z>pow(10,7.05) && z<pow(10,7.20)) s_2H = pow(z,-7.04065)*pow(10,67.5667);
		else if(z>pow(10,7.2) && z<pow(10,7.31)) s_2H = pow(z,-11.7845)*pow(10,101.75);
		else s_2H=0;
		s_2H = E_0*s_2H/5000;
	}

	K = B*s_2H;
	// cout << " s = " << s_2H << " s_int = " << s_int << endl;
	//~ cout << " s = " << s << " B = " << B <<  " Y = " << Y << endl;
	return K;
}
double  func_z_v2_sans_reinj(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
{

	double  Emin[19];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  A;
	double y, dy;


	double s_4He=0,s_4He_int=0,s_2H=0,s_2H_int=0;
	double s_4He_table=0,s_2H_table=0;
	//    cout << " z = " << z << endl;
	// for(int k = 15; k<19;k++)
	// 	{
	// 		if(E_c>Emin[k])s_4He_int+=qsimp_E(func_sans_reinj,Emin[k],E_c,z,k,E_0,Z_x, z_0);
	// 	}
	// if(E_c>Emin[14])s_2H_int=qsimp_E(func_sans_reinj,Emin[14],E_c,z,14,E_0,Z_x, z_0);
	/****Fit realise avec E_0 = 5000 MeV soit M_X = 10000 MeV****/

	// if(E_c>Emin[15]){
	// 	if(z<=pow(10,5.33))s_4He+=pow(z,-2.51233)*pow(10,36.7075);
	// 	else if(z>pow(10,5.33) && z<=pow(10,5.78))s_4He+=pow(z,-3.2278)*pow(10,40.5342);
	// 	else if(z>pow(10,5.78) && z<=pow(10,6.1))s_4He+=pow(z,-5.14742)*pow(10,51.6394);
	// 	else if(z>pow(10,6.1) && z<=pow(10,6.2))s_4He+=pow(z,-6.04727)*pow(10,57.132);
	// 	else if(z>pow(10,6.2) && z<=pow(10,6.3))s_4He+=pow(z,-7.78083)*pow(10,67.8923);
	// 	else if(z>pow(10,6.3) && z<=pow(10,6.37))s_4He+=pow(z,-13.083)*pow(10,101.287);
	// 	//~ else if(z>pow(10,6.19) && z<=pow(10,6.25))s_4He=pow(z,-7.50809)*pow(10,66.2643);
	// 	//~ else {for(int k = 15; k<19;k++)if(E_c>Emin[k])s_4He=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);}
	// 	else s_4He+=0;
	// 	s_4He = E_0*s_4He/5000;}
	if(E_c>Emin[15]){
		// for(int k=15;k<19;k++)s_4He_diffuse+=qsimp(func_trois_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
	linearint(vector_z_s4He_destruc_standard,vector_s4He_destruc_standard,vector_z_s4He_destruc_standard.size(),log10(z),y,dy);
		// output_Check << log10(z) << " "<< y << endl;
		s_4He=E_0*pow(10,y)/5000;
		// s_2H_compton = qsimp(func_spectre_compton,E_min,E_0,z,14,E_0,Z_x, z_0);
		// if(s_2H_diffuse <s_2H_compton) cout << "s_2H_compton =  " << s_2H_compton << endl;
		// cout << "s_4He_diffuse quatre iterations=  " << s_4He_diffuse << endl;
	}
		/****Fit realise avec E_0 = 5000 MeV soit M_X = 10000 MeV****/
		if(E_c>Emin[14]){if(z<pow(10,5.91))s_2H = pow(z,-2.50882)*pow(10,36.8106);
			else if(z>pow(10,5.91) && z<pow(10,6.6)) s_2H = pow(z,-2.95593)*pow(10,39.4551);
			else if(z>pow(10,6.53) && z<pow(10,6.8)) s_2H = pow(z,-3.89315)*pow(10,45.5827);
			else if(z>pow(10,6.8) && z<pow(10,7.05)) s_2H = pow(z,-4.68502)*pow(10,50.966);
			else if(z>pow(10,7.05) && z<pow(10,7.20)) s_2H = pow(z,-7.04065)*pow(10,67.5667);
			else if(z>pow(10,7.2) && z<pow(10,7.31)) s_2H = pow(z,-11.7845)*pow(10,101.75);
			else s_2H=0;
			s_2H = E_0*s_2H/5000;}

			// cout << "s He int = " << s_4He_int << " s He = " << s_4He << endl;
			// cout << "s 2H int = " << s_2H_int << " s 2H = " << s_2H << endl;

			A = exp(-1/(2*H_r*tau*(z+1)*(z+1)))*(s_4He-s_2H);
			//~ cout << " A = " << A << " z = " << z <<  endl;
			return A;
		}



	double  Y_He_K_Perte_sans_reinj(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double E_min = 2.23;
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  n_y_0 = 3.154*pow(10.,-30) ;
	double B=(Z_x*n_y_0)/(E_0*H_r*tau);
	double Y_He,s;
	s=qsimp_E(func_z_v2_sans_reinj,z,z_0,z_0,14,E_0,Z_x, z_0);
	//~ cout << "z_0 " << z_0 << " z = " << z << endl;
	Y_He = exp(-B*s);
	cout << " Y = " << Y_He << " B = " << B << " s = " << s << endl;
	return Y_He;
}
double  S_Gain_sans_reinj(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
{	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  Emin[19];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  n_y_0 = 3.154*pow(10.,-30) ;
	//~ cout << " tau = " << tau << endl;
	double a=z;
	double B=(Z_x*n_y_0)/(E_0*H_r*tau);
	double Y;
	double Y_He_0=0.25,s_4He =0;
	double count_s4He_2d=0, tmp_z, s_4He_table = 0, tmp_s4He_2d;
	double count_s4He_d=0, tmp_s4He_d;
	double s=0;
	// for(int k=17;k<19;k++)
	// {
	// 	if(E_c>Emin[k])
	// 	{
	// 		if(k==17)s+=2*qsimp_E(func_sans_reinj,Emin[k],E_c,a,k,E_0,Z_x, z_0);
	// 		if(k==18)s+=qsimp_E(func_sans_reinj,Emin[k],E_c,a,k,E_0,Z_x, z_0);
	// 			// cout << " s = " << s << endl;
	// 	}
	// }

	//
	if(E_c>Emin[17]){if(z<=pow(10,5.04))s_4He+=pow(z,-2.51795)*pow(10,35.5388);
		else if(z>pow(10,5.04) && z<=pow(10,5.4))s_4He+=pow(z,-2.87902)*pow(10,37.3408);
		else if(z>pow(10,5.4) && z<=pow(10,5.86))s_4He+=pow(z,-4.3513)*pow(10,45.312);
		else if(z>pow(10,5.86) && z<=pow(10,6.11))s_4He+=pow(z,-7.06984)*pow(10,61.2727);
		else if(z>pow(10,6.11) && z<=pow(10,6.2))s_4He+=pow(z,-12.1155)*pow(10,92.1182);
		else if(z>pow(10,6.2) && z<=pow(10,6.3))s_4He+=pow(z,-30.3648)*pow(10,205.645);
		//~ else if(z>pow(10,6.19) && z<=pow(10,6.25))s_4He=pow(z,-7.50809)*pow(10,66.2643);
		//~ else {for(int k = 15; k<19;k++)if(E_c>Emin[k])s_4He=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);}
		else s_4He+=0;
		s_4He = E_0*s_4He/5000;}


		// cout << " s table = " << s_4He << " s_int = " << s << endl;
		Y =exp(-1/(2*H_r*pow(1+z,2)*tau))*Y_He_K_Perte_sans_reinj(a,z_0,i,z_0,E_0,Z_x, z_0)*s_4He*Y_He_0;
		//~ Y =-B*exp(-1/(2*H_r*pow(1+z,2)*tau))*s*Y_He_0;
		//~ Y = 1;
		//~ cout << " Y = " << Y << "z = " << z << endl;
		//~ cout << " Y_He_K_Perte(z,z_0,i,z_0,E_0,Z_x, z_0) = " << Y_He_K_Perte(z,z_0,i,z_0,E_0,Z_x, z_0) << endl;
		return Y;

	}
/*******EXTRA FONCTION : Limite inferieure via production de Deuterium uniquement**********/
	double  S_Gain_Check_De(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
{	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  Emin[19];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  n_y_0 = 3.154*pow(10.,-30) ;
	//~ cout << " tau = " << tau << endl;
	double a=z;
	double B=(Z_x*n_y_0)/(E_0*H_r*tau);
	double Y=0;
	double Y_He_0 = 0.25,s_4He=0;
	double s=0;
	// for(int k=17;k<19;k++)
	// 	{
	// 		if(E_c>Emin[k])
	// 			{
	// 				if(k==17)s+=2*qsimp_E(func_sans_reinj,Emin[k],E_c,a,k,E_0,Z_x, z_0);
	// 				else if(k==18)s+=qsimp_E(func_sans_reinj,Emin[k],E_c,a,k,E_0,Z_x, z_0);
	// 			}
	// 		}

	if(E_c>Emin[17]){	if(z<=pow(10,5.15))s_4He+=pow(z,-2.51665)*pow(10,37.5311);
		else if(z>pow(10,5.15) && z<=pow(10,5.56))s_4He+=pow(z,-3.36158)*pow(10,41.8865);
		else if(z>pow(10,5.56) && z<=pow(10,5.83))s_4He+=pow(z,-4.84291)*pow(10,50.1265 );
		else if(z>pow(10,5.83) && z<=pow(10,6.06))s_4He+=pow(z,-5.95335)*pow(10,56.6027);
		else if(z>pow(10,6.06) && z<=pow(10,6.19))s_4He+=pow(z,-10.4508)*pow(10,83.8837);
		else if(z>pow(10,6.19) && z<=pow(10,6.3))s_4He+=pow(z,-23.1571)*pow(10,162.642);
		//~ else if(z>pow(10,6.19) && z<=pow(10,6.25))s_4He=pow(z,-7.50809)*pow(10,66.2643);
		//~ else {for(int k = 15; k<19;k++)if(E_c>Emin[k])s_4He=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);}
		else s_4He+=0;}
		// if(E_c>Emin[17]){if(z<=pow(10,5.3))s=pow(z,-2.60454)*pow(10,36.0398);
		// 	else if(z>pow(10,5.3) && z<=pow(10,5.6))s=pow(z,-3.27568)*pow(10,39.5678);
		// 	else if(z>pow(10,5.6) && z<=pow(10,6.15))s=pow(z,-5.83282)*pow(10,53.9883);
		// 	else if(z>pow(10,6.2) && z<=pow(10,6.315))s=pow(z,-33.4332)*pow(10,224.703);
		// 	else s = 0;}
		Y =exp(-1/(2*H_r*pow(1+z,2)*tau))*s_4He*Y_He_0;

		return Y;
	}
	/*********************************Fin extra fonction**************************************/
	/*****************************************************************************************/
	/***********************Fin Fonction Deuterium Standard AVEC reinj************************/
	/*****************************************************************************************/













/************************************************************************************/
/***********************Debut Fonction Deuterium Non Standard************************/
/************************************************************************************/
double  func_non_standard(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  f;
	f=func_sigma(E_0,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(E_0,z,i,g,E_0,Z_x, z_0)+gamma_compton(E_0,z,i,g,E_0,Z_x, z_0)+gamma_phph(E_0,z,i,g,E_0,Z_x, z_0));
	// f=func_sigma(E_0,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0)*correction_phph(x,z,i,g,E_0,Z_x, z_0));
	return f;
}
double func_deux_interactions(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
return Spectre_deux_interactions(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

}
double func_trois_interactions(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
return Spectre_trois_interactions(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

}
double func_quatre_interactions(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
return Spectre_quatre_interactions(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

}
double func_quatre_interactions_corrected(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
double correction, dy;
linearint(vector_E_correction,vector_correction,vector_E_correction.size(),x,correction,dy);
return	(Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)+Spectre_deux_interactions(x,z,i,g,E_0,Z_x, z_0)+Spectre_trois_interactions(x,z,i,g,E_0,Z_x, z_0)+Spectre_quatre_interactions(x,z,i,g,E_0,Z_x, z_0))*func_sigma(x,z,i,g,E_0,Z_x, z_0)*correction/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

}
double func_quatre_interactions_produc2H(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
double correction, dy;
linearint(vector_E_correction,vector_correction,vector_E_correction.size(),x,correction,dy);
double f = (Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)+Spectre_deux_interactions(x,z,i,g,E_0,Z_x, z_0)+Spectre_trois_interactions(x,z,i,g,E_0,Z_x, z_0)+Spectre_quatre_interactions(x,z,i,g,E_0,Z_x, z_0))*correction*(2*func_sigma(x,z,17,g,E_0,Z_x, z_0)+func_sigma(x,z,18,g,E_0,Z_x, z_0))/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));
return f;
}
double func_cinq_interactions(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
return Spectre_cinq_interactions(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

}
double func_six_interactions(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
return Spectre_six_interactions(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

}
double func_sept_interactions(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
return Spectre_sept_interactions(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

}
double  func_une_interaction_monochromatique(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  f;
	f=Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
	// f=Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0)*correction_phph(x,z,i,g,E_0,Z_x, z_0));
	return f;
}

double func_spectre_compton(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0){
	double f;
	f = Spectre_Gamma_compton(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
	return f;
}
double  K_Perte_non_standard(double z, double x, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double E_min = 2.23;
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  n_y_0 = 3.154*pow(10.,-30) ;
	double s_2H_compton=0;
	double B= exp(-1/(2*H_r*tau*(z+1)*(z+1)));
	double K, s=0, s_2H = 0,s_2H_diffuse = 0,y,dy;
	//~ cout << " tau = " << tau << " z = " << z <<endl;
	if(v2==0){

				if(E_0<E_c)if(E_0>E_min)s_2H=func_non_standard(E_0,z,14,g,E_0,Z_x, z_0);
	}

	else{
					if(E_c>E_0){
								s_2H += func_non_standard(E_0,z,14,g,E_0,Z_x, z_0);
										if(avec_correction >=1){
											// cout << " E_0 " << E_0 << endl;
											if(E_0==20){if(z<pow(10,5.52))s_2H_diffuse = pow(z,-2.99869)*pow(10,37.0089);
											else if(z>pow(10,5.52) && z<pow(10,6.11)) s_2H_diffuse = pow(z,-2.68106)*pow(10,35.2528);
											else if(z>pow(10,6.11) && z<pow(10,6.36)) s_2H_diffuse = pow(z,-3.37479)*pow(10,39.4939);
											else if(z>pow(10,6.36) && z<pow(10,6.75)) s_2H_diffuse = pow(z,-4.17679)*pow(10,44.5983);
											else if(z>pow(10,6.75) && z<pow(10,7.20)) s_2H_diffuse = pow(z,-5.28615)*pow(10,52.0874);
											else if(z>pow(10,7.05) && z<pow(10,7.21)) s_2H_diffuse = pow(z,-7.89099)*pow(10,70.4448);
											else if(z>pow(10,7.21) && z<pow(10,7.32)) s_2H_diffuse = pow(z,-12.8303)*pow(10,106.033);
											else s_2H_diffuse=0;
											// cout << " s_2H_fit = " << s_2H_diffuse << endl;;
										}
										else if(E_0==4){if(z<pow(10,6.75))s_2H_diffuse = pow(z,-2.99563)*pow(10,36.6403);
										else if(z>pow(10,6.75) && z<pow(10,6.93)) s_2H_diffuse = pow(z,-3.50341)*pow(10,40.1183);
										else if(z>pow(10,6.93) && z<pow(10,7.1)) s_2H_diffuse = pow(z,-4.47742)*pow(10,46.8587);
										else if(z>pow(10,7.1) && z<pow(10,7.2)) s_2H_diffuse = pow(z,-8.69056)*pow(10,76.787);
										else if(z>pow(10,7.2) && z<pow(10,7.29)) s_2H_diffuse = pow(z,-13.2551)*pow(10,109.686);
										else if(z>pow(10,7.29) && z<pow(10,7.35)) s_2H_diffuse = pow(z,-30.332)*pow(10,234.321);
										else s_2H_diffuse=0;
									}
									//
										// else s_2H_diffuse +=qsimp_E(func_une_interaction_monochromatique,E_min,E_0,z,14,E_0,Z_x, z_0);
										s_2H_diffuse =qsimp_E(func_une_interaction_monochromatique,E_min,E_0,z,14,E_0,Z_x, z_0);
										// cout << " s_2H_vrai = " << s_2H_diffuse << endl;;

									// linearint(vector_z_s2H_destruc,vector_s2H_destruc,vector_z_s2H_destruc.size(),log10(z),y,dy);
									// s_2H_diffuse=pow(10,y);
							// cout << " s_2H_interpolation = " << s_2H_diffuse<< endl;;
										}
										if(avec_correction>=2){


																		if(E_0==70){if(z<pow(10,5.1))s_2H_diffuse += pow(z,-2.6718)*pow(10,35.105);
																		else if(z>pow(10,5.1) && z<pow(10,5.54)) s_2H_diffuse += pow(z,-1.88507)*pow(10,31.096);
																		else if(z>pow(10,5.54) && z<pow(10,5.84)) s_2H_diffuse += pow(z,-2.65906)*pow(10,35.3873);
																		else if(z>pow(10,5.84) && z<pow(10,6.57)) s_2H_diffuse += pow(z,-4.03745)*pow(10,43.4568);
																		else if(z>pow(10,6.57) && z<pow(10,6.92)) s_2H_diffuse += pow(z,-5.39545)*pow(10,52.3808);
																		else if(z>pow(10,6.92) && z<pow(10,7.11)) s_2H_diffuse += pow(z,-7.78715)*pow(10,68.9499);
																		else if(z>pow(10,7.11) && z<pow(10,7.33)) s_2H_diffuse += pow(z,-11.2867)*pow(10,93.8556);
																		else s_2H_diffuse+=0;
																		}
																		else if(E_0==30){if(z<pow(10,5.46))s_2H_diffuse += pow(z,-2.90147)*pow(10,36.2585);
																		else if(z>pow(10,5.46) && z<pow(10,6.13)) s_2H_diffuse += pow(z,-2.01308)*pow(10,31.4114);
																		else if(z>pow(10,6.13) && z<pow(10,6.74)) s_2H_diffuse += pow(z,-4.48436)*pow(10,46.5841);
																		else if(z>pow(10,6.74) && z<pow(10,7.05)) s_2H_diffuse += pow(z,-6.10206)*pow(10,57.5085);
																		else if(z>pow(10,7.05) && z<pow(10,7.18)) s_2H_diffuse += pow(z,-9.90638)*pow(10,84.3175);
																		else if(z>pow(10,7.18) && z<pow(10,7.25)) s_2H_diffuse += pow(z,-14.081)*pow(10,114.343);
																		else if(z>pow(10,7.25) && z<pow(10,7.33)) s_2H_diffuse += pow(z,-20.6108)*pow(10,161.7);
																		else s_2H_diffuse+=0;
																		// cout << " E_0 = " << E_0 << endl;

																		}
																		else {s_2H_diffuse += qsimp(func_deux_interactions,E_min,E_0,z,14,E_0,Z_x, z_0);
																			// cout << "s_2H_diffuse =  " << s_2H_diffuse << endl;
																			// cout << "s_2H_diffuse deux iterations=  " << s_2H_diffuse << endl;

																		}
																		// linearint(vector_z_s2H_destruc_deux_iterations,vector_s2H_destruc_deux_iterations,vector_s2H_destruc_deux_iterations.size(),log10(z),y,dy);
																		// s_2H_diffuse += pow(10,y);
																		// cout << " s_2H_diffuse interpolation = " << s_2H_diffuse << endl;
									}
									if(avec_correction>=3){

										linearint(vector_z_s2H_destruc_trois_iterations,vector_s2H_destruc_trois_iterations,vector_z_s2H_destruc_trois_iterations.size(),log10(z),y,dy);
										s_2H_diffuse+=pow(10,y);

										// else s_2H_diffuse += qsimp(func_trois_interactions,E_min,E_0,z,14,E_0,Z_x, z_0);

										// s_2H_compton = qsimp(func_spectre_compton,E_min,E_0,z,14,E_0,Z_x, z_0);
										// if(s_2H_diffuse <s_2H_compton) cout << "s_2H_compton =  " << s_2H_compton << endl;
										// cout << "s_2H_diffuse trois iterations=  " << s_2H_diffuse << endl;
									}
									if(avec_correction==0){

										linearint(vector_z_s2H_destruc_quatre_iterations,vector_s2H_destruc_quatre_iterations,vector_z_s2H_destruc_quatre_iterations.size(),log10(z),y,dy);
										s_2H_diffuse+=pow(10,y);

										// else s_2H_diffuse += qsimp(func_trois_interactions,E_min,E_0,z,14,E_0,Z_x, z_0);

										// s_2H_compton = qsimp(func_spectre_compton,E_min,E_0,z,14,E_0,Z_x, z_0);
										// if(s_2H_diffuse <s_2H_compton) cout << "s_2H_compton =  " << s_2H_compton << endl;
										// cout << "s_2H_diffuse quatre iterations=  " << s_2H_diffuse << endl;
									}




									// cout << " s_2H_diffuse deux interactions " << s_2H_diffuse << endl;

																// s_2H_bis=func_non_standard(E_0,z,14,g,E_0,Z_x, z_0);
																// cout << " s_2H_bis " << s_2H_bis << endl;
					}

					else {


						if(E_c>E_min){if(z<pow(10,5.91))s_2H = pow(z,-2.50882)*pow(10,36.8106);
						else if(z>pow(10,5.91) && z<pow(10,6.6)) s_2H = pow(z,-2.95593)*pow(10,39.4551);
						else if(z>pow(10,6.53) && z<pow(10,6.8)) s_2H = pow(z,-3.89315)*pow(10,45.5827);
						else if(z>pow(10,6.8) && z<pow(10,7.05)) s_2H = pow(z,-4.68502)*pow(10,50.966);
						else if(z>pow(10,7.05) && z<pow(10,7.20)) s_2H = pow(z,-7.04065)*pow(10,67.5667);
						else if(z>pow(10,7.2) && z<pow(10,7.31)) s_2H = pow(z,-11.7845)*pow(10,101.75);
						else s_2H=0;


						s_2H = E_0*s_2H/5000;
						// cout << "s_2H " << s_2H << endl;
					}
				}

		}
		// cout << "s_2H tot" << s_2H << endl;
		K = B*(s_2H+s_2H_diffuse+s_2H_compton);
		//~ cout << " s = " << s << " s_int = " << s_int << endl;
		//~ cout << " s = " << s << " B = " << B <<  " Y = " << Y << endl;
		return K;
	}


double  func_z_v2_non_standard(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
{


	double  Emin[19];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  A;
	double s_4He=0,s_4He_diffuse=0,s_2H=0,s_2H_diffuse=0;
	double s_4He_table=0,s_2H_table=0;
	double y, dy;
	//  cout << " z = " << z << endl;

	if(v2==0){

					if(E_0<E_c){
						for(int k = 15; k<19;k++)
						{
							if(E_0>Emin[k])s_4He+=func_non_standard(E_0,z,k,g,E_0,Z_x, z_0);
						}
						if(E_0>Emin[14])s_2H=func_non_standard(E_0,z,14,g,E_0,Z_x, z_0);
					}
	}




	else{
					if(E_c > E_0){
						for(int k = 15; k<19;k++){if(E_0>Emin[k])s_4He+=func_non_standard(E_0,z,k,g,E_0,Z_x, z_0);}
						if(avec_correction >=1){
							if(E_0==70){
								if(z<pow(10,4.92))s_4He_diffuse = pow(z,-2.94973)*pow(10,36.8374);
								else if(z>pow(10,4.92) && z<pow(10,5.47)) s_4He_diffuse = pow(z,-2.17361)*pow(10,33.0111);
								else if(z>pow(10,5.47) && z<pow(10,5.81)) s_4He_diffuse = pow(z,-3.76239)*pow(10,41.7045);
								else if(z>pow(10,5.81) && z<pow(10,6.12)) s_4He_diffuse = pow(z,-5.42405)*pow(10,51.3708);
								else if(z>pow(10,6.12) && z<pow(10,6.3)) s_4He_diffuse = pow(z,-7.23229)*pow(10,62.4409);
								else if(z>pow(10,6.3) && z<pow(10,6.37)) s_4He_diffuse = pow(z,-18.4285)*pow(10,133.054);
								else s_4He_diffuse=0;
							}
							else if(E_0==30){
								if(z<pow(10,4.92))s_4He_diffuse = pow(z,-2.98582)*pow(10,36.9087);
								else if(z>pow(10,4.92) && z<pow(10,5.47)) s_4He_diffuse = pow(z,-2.46896)*pow(10,34.1862);
								else if(z>pow(10,5.47) && z<pow(10,5.81)) s_4He_diffuse = pow(z,-3.88914)*pow(10,42.4291);
								else if(z>pow(10,5.81) && z<pow(10,6.12)) s_4He_diffuse = pow(z,-5.3442)*pow(10,51.1777);
								else if(z>pow(10,6.12) && z<pow(10,6.3)) s_4He_diffuse = pow(z,-8.59309)*pow(10,71.3725);
								else if(z>pow(10,6.3) && z<pow(10,6.37)) s_4He_diffuse = pow(z,-14.6491)*pow(10,109.499);
								else s_4He_diffuse=0;
							}
							//     	else {for(int k = 15; k<19;k++)if(E_0>Emin[k])s_4He_diffuse+=qsimp(func_une_interaction_monochromatique,Emin[k],E_0,z,k,E_0,Z_x, z_0);}
						}
						if(avec_correction >=2){
							// for(int k = 15; k<19;k++)s_4He_diffuse+=qsimp(func_deux_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
							if(E_0==70){
								if(z<pow(10,4.97))s_4He_diffuse += pow(z,-2.60856)*pow(10,34.5784);
								else if(z>pow(10,4.97) && z<pow(10,5.14)) s_4He_diffuse += pow(z,-1.57864)*pow(10,29.4573);
								else if(z>pow(10,5.14) && z<pow(10,5.42)) s_4He_diffuse += pow(z,-1.03726)*pow(10,26.6741);
								else if(z>pow(10,5.42) && z<pow(10,5.61)) s_4He_diffuse += pow(z,-2.0937)*pow(10,32.4017);
								else if(z>pow(10,5.61) && z<pow(10,5.86)) s_4He_diffuse += pow(z,-3.37311)*pow(10,39.5797);
								else if(z>pow(10,5.86) && z<pow(10,6.15)) s_4He_diffuse += pow(z,-7.14069)*pow(10,61.6933);
								else if(z>pow(10,6.15) && z<pow(10,6.35)) s_4He_diffuse += pow(z,-10.447)*pow(10,82.0597);
								else s_4He_diffuse+=0;
							}
							else if(E_0==30){
								if(z<pow(10,5.39))s_4He_diffuse += pow(z,-2.38224)*pow(10,32.6604);
								else if(z>pow(10,5.39) && z<pow(10,5.52)) s_4He_diffuse += pow(z,-1.65083)*pow(10,28.7143);
								else if(z>pow(10,5.52) && z<pow(10,5.72)) s_4He_diffuse += pow(z,-1.16574)*pow(10,26.0349);
								else if(z>pow(10,5.72) && z<pow(10,5.83)) s_4He_diffuse += pow(z,-1.72495)*pow(10,29.2353);
								else if(z>pow(10,5.83) && z<pow(10,5.99)) s_4He_diffuse += pow(z,-2.94351)*pow(10,36.3446);
								else if(z>pow(10,5.99) && z<pow(10,6.22)) s_4He_diffuse += pow(z,-4.82053)*pow(10,47.5996);
								else if(z>pow(10,6.22) && z<pow(10,6.35)) s_4He_diffuse += pow(z,-13.0298)*pow(10,98.687);
								else s_4He_diffuse+=0;
							}
							// else {for(int k = 15; k<19;k++)if(E_0>Emin[k])s_4He_diffuse+=qsimp(func_deux_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);}

							}
						// cout << " s_4He = " << s_4He << " s_4He_diffuse = " << s_4He_diffuse << endl;
						if(E_0>Emin[14])s_2H=func_non_standard(E_0,z,14,g,E_0,Z_x, z_0);
						if(avec_correction >=1){
							if(E_0>Emin[14])s_2H_diffuse += qsimp(func_une_interaction_monochromatique,Emin[14],E_0,z,14,E_0,Z_x, z_0);
						}
						if(avec_correction >=2){
							// if(E_0>Emin[14])s_2H_diffuse+=qsimp(func_deux_interactions,Emin[14],E_0,z,14,E_0,Z_x, z_0);
							// cout << " E_0 = " << E_0 << endl;

							if(E_0==70){if(z<pow(10,5.1))s_2H_diffuse += pow(z,-2.6718)*pow(10,35.105);
							else if(z>pow(10,5.1) && z<pow(10,5.54)) s_2H_diffuse += pow(z,-1.88507)*pow(10,31.096);
							else if(z>pow(10,5.54) && z<pow(10,5.84)) s_2H_diffuse += pow(z,-2.65906)*pow(10,35.3873);
							else if(z>pow(10,5.84) && z<pow(10,6.57)) s_2H_diffuse += pow(z,-4.03745)*pow(10,43.4568);
							else if(z>pow(10,6.57) && z<pow(10,6.92)) s_2H_diffuse += pow(z,-5.39545)*pow(10,52.3808);
							else if(z>pow(10,6.92) && z<pow(10,7.11)) s_2H_diffuse += pow(z,-7.78715)*pow(10,68.9499);
							else if(z>pow(10,7.11) && z<pow(10,7.33)) s_2H_diffuse += pow(z,-11.2867)*pow(10,93.8556);
							else s_2H_diffuse+=0;
							}
							else if(E_0==30){if(z<pow(10,5.46))s_2H_diffuse += pow(z,-2.90147)*pow(10,36.2585);
							else if(z>pow(10,5.46) && z<pow(10,6.13)) s_2H_diffuse += pow(z,-2.01308)*pow(10,31.4114);
							else if(z>pow(10,6.13) && z<pow(10,6.74)) s_2H_diffuse += pow(z,-4.48436)*pow(10,46.5841);
							else if(z>pow(10,6.74) && z<pow(10,7.05)) s_2H_diffuse += pow(z,-6.10206)*pow(10,57.5085);
							else if(z>pow(10,7.05) && z<pow(10,7.18)) s_2H_diffuse += pow(z,-9.90638)*pow(10,84.3175);
							else if(z>pow(10,7.18) && z<pow(10,7.25)) s_2H_diffuse += pow(z,-14.081)*pow(10,114.343);
							else if(z>pow(10,7.25) && z<pow(10,7.33)) s_2H_diffuse += pow(z,-20.6108)*pow(10,161.7);
							else s_2H_diffuse+=0;
							// cout << " E_0 = " << E_0 << endl;

							}
							//
							// else {if(E_0>Emin[14])s_2H_diffuse += qsimp(func_deux_interactions,Emin[14],E_0,z,14,E_0,Z_x, z_0);
							// 	cout << " E_0 = " << E_0 << endl;
							// }

						}
						if(avec_correction >=3){
							linearint(vector_z_s4He_destruc_trois_iterations,vector_s4He_destruc_trois_iterations,vector_z_s4He_destruc_trois_iterations.size(),log10(z),y,dy);
							s_4He_diffuse+=pow(10,y);
							linearint(vector_z_s2H_destruc_trois_iterations,vector_s2H_destruc_trois_iterations,vector_z_s2H_destruc_trois_iterations.size(),log10(z),y,dy);
							s_2H_diffuse+=pow(10,y);
							// cout << " s_2H = " << s_2H << " s_2H_diffuse 3 iterations = " << s_2H_diffuse << endl;
							// cout << " s_4He = " << s_4He << " s_4He_diffuse 3 iterations= " << s_4He_diffuse << endl;
						}

						if(avec_correction==0){
							linearint(vector_z_s4He_destruc_quatre_iterations,vector_s4He_destruc_quatre_iterations,vector_z_s4He_destruc_quatre_iterations.size(),log10(z),y,dy);
							s_4He_diffuse+=pow(10,y);
							linearint(vector_z_s2H_destruc_quatre_iterations,vector_s2H_destruc_quatre_iterations,vector_z_s2H_destruc_quatre_iterations.size(),log10(z),y,dy);
							s_2H_diffuse+=pow(10,y);
							// cout << " s_2H_diffuse 4 iterations = " << s_2H_diffuse << endl;
							// cout << " s_4He_diffuse 4 iterations= " << s_4He_diffuse << endl;
						}

					}

					else {
						// if(E_c>Emin[15]){
						// 	if(z<=pow(10,5.33))s_4He+=pow(z,-2.51233)*pow(10,36.7075);
						// 	else if(z>pow(10,5.33) && z<=pow(10,5.78))s_4He+=pow(z,-3.2278)*pow(10,40.5342);
						// 	else if(z>pow(10,5.78) && z<=pow(10,6.1))s_4He+=pow(z,-5.14742)*pow(10,51.6394);
						// 	else if(z>pow(10,6.1) && z<=pow(10,6.2))s_4He+=pow(z,-6.04727)*pow(10,57.132);
						// 	else if(z>pow(10,6.2) && z<=pow(10,6.3))s_4He+=pow(z,-7.78083)*pow(10,67.8923);
						// 	else if(z>pow(10,6.3) && z<=pow(10,6.37))s_4He+=pow(z,-13.083)*pow(10,101.287);
						// 	//~ else if(z>pow(10,6.19) && z<=pow(10,6.25))s_4He=pow(z,-7.50809)*pow(10,66.2643);
						// 	//~ else {for(int k = 15; k<19;k++)if(E_c>Emin[k])s_4He=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);}
						// 	else s_4He+=0;
						// 	s_4He = E_0*s_4He/5000;}
						if(E_c>Emin[15]){
							// for(int k=15;k<19;k++)s_4He_diffuse+=qsimp(func_trois_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
						linearint(vector_z_s4He_destruc_standard,vector_s4He_destruc_standard,vector_z_s4He_destruc_standard.size(),log10(z),y,dy);
							// output_Check << log10(z) << " "<< y << endl;
							s_4He=E_0*pow(10,y)/5000;
							// s_2H_compton = qsimp(func_spectre_compton,E_min,E_0,z,14,E_0,Z_x, z_0);
							// if(s_2H_diffuse <s_2H_compton) cout << "s_2H_compton =  " << s_2H_compton << endl;
							// cout << "s_4He_diffuse quatre iterations=  " << s_4He_diffuse << endl;
						}
							/****Fit realise avec E_0 = 5000 MeV soit M_X = 10000 MeV****/
							if(E_c>Emin[14]){if(z<pow(10,5.91))s_2H = pow(z,-2.50882)*pow(10,36.8106);
								else if(z>pow(10,5.91) && z<pow(10,6.6)) s_2H = pow(z,-2.95593)*pow(10,39.4551);
								else if(z>pow(10,6.53) && z<pow(10,6.8)) s_2H = pow(z,-3.89315)*pow(10,45.5827);
								else if(z>pow(10,6.8) && z<pow(10,7.05)) s_2H = pow(z,-4.68502)*pow(10,50.966);
								else if(z>pow(10,7.05) && z<pow(10,7.20)) s_2H = pow(z,-7.04065)*pow(10,67.5667);
								else if(z>pow(10,7.2) && z<pow(10,7.31)) s_2H = pow(z,-11.7845)*pow(10,101.75);
								else s_2H=0;
								s_2H = E_0*s_2H/5000;}
							}
			}

			A = exp(-1/(2*H_r*tau*(z+1)*(z+1)))*(s_4He+s_4He_diffuse-s_2H-s_2H_diffuse);
			//~ cout << " A = " << A << " z = " << z <<  endl;
			return A;
}

		double  Y_He_K_Perte_non_standard(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
	{
		double  T_0 = 2.7255*0.862*pow(10.,-10);
		double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
		double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
		double E_min = 2.23;
		double  tau =  5*pow(1+z_0,-2)/(2*H_r);
		double  n_y_0 = 3.154*pow(10.,-30) ;
		double B=(Z_x*n_y_0)/(E_0*H_r*tau);
		double Y_He,s;
		double y, dy;
		s=qsimp_E(func_z_v2_non_standard,z,z_0,z_0,14,E_0,Z_x, z_0);

		//~ cout << "z_0 " << z_0 << " z = " << z << endl;
		Y_He = exp(-B*s);
		// cout << " Y = " << Y_He << " B = " << B << " s = " << s << endl;
		return Y_He;
	}


	double  S_Gain_non_standard(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
{	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  Emin[19];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  n_y_0 = 3.154*pow(10.,-30) ;
	//~ cout << " tau = " << tau << endl;
	double a=z;
	double B=(Z_x*n_y_0)/(E_0*H_r*tau);
	double Y;
	double Y_He_0=0.25,s_4He =0, s_4He_diffuse=0;
	double s=0;
	double y, dy;
	if(v2==0){
					if(E_0<E_c){
						for(int k=17;k<19;k++){
							if(E_c>Emin[k]){
								if(k==17)s_4He+=2*func_non_standard(E_0,z,k,g,E_0,Z_x, z_0);
								if(k==18)s_4He+=func_non_standard(E_0,z,k,g,E_0,Z_x, z_0);
								if(avec_correction ==1){
									if(k==17)s_4He_diffuse+=2*qsimp_3(func_une_interaction_monochromatique,Emin[k],E_0,z,k,E_0,Z_x, z_0);
									if(k==18)s_4He_diffuse+=qsimp_3(func_une_interaction_monochromatique,Emin[k],E_0,z,k,E_0,Z_x, z_0);
								}
								// cout << " s_4He = " << s_4He << " s_4He_diffuse = " << s_4He_diffuse << endl;
							}
						}
					}
	}
	else{

					if(E_0<E_c){
						for(int k=17;k<19;k++){
							if(E_c>Emin[k]){
								if(k==17)s_4He+=2*func_non_standard(E_0,z,k,g,E_0,Z_x, z_0);
								if(k==18)s_4He+=func_non_standard(E_0,z,k,g,E_0,Z_x, z_0);
								if(avec_correction >=1){
									if(k==17)s_4He_diffuse+=2*qsimp_3(func_une_interaction_monochromatique,Emin[k],E_0,z,k,E_0,Z_x, z_0);
									if(k==18)s_4He_diffuse+=qsimp_3(func_une_interaction_monochromatique,Emin[k],E_0,z,k,E_0,Z_x, z_0);
								}

								// cout << " s_4He = " << s_4He << " s_4He_diffuse = " << s_4He_diffuse << endl;
							}
						}


						if(avec_correction >=2){
							// if(k==17)s_4He_diffuse+=2*qsimp(func_deux_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
							// if(k==18)s_4He_diffuse+=qsimp(func_deux_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
							//

							if(E_0==70){
								if(z<pow(10,4.97))s_4He_diffuse += pow(z,-2.79472)*pow(10,34.0765);
								else if(z>pow(10,4.97) && z<pow(10,5.41)) s_4He_diffuse += pow(z,-0.993372)*pow(10,25.1586);
								else if(z>pow(10,5.41) && z<pow(10,5.54)) s_4He_diffuse += pow(z,-2.22046)*pow(10,31.7923);
								else if(z>pow(10,5.54) && z<pow(10,5.69)) s_4He_diffuse += pow(z,-3.52412)*pow(10,39.0137);
								else if(z>pow(10,5.69) && z<pow(10,5.85)) s_4He_diffuse += pow(z,-4.60985)*pow(10,45.1994);
								else if(z>pow(10,5.85) && z<pow(10,5.99)) s_4He_diffuse += pow(z,-7.75257)*pow(10,63.5621);
								else if(z>pow(10,5.99) && z<pow(10,6.11)) s_4He_diffuse += pow(z,-10.898)*pow(10,82.3944);
								else if(z>pow(10,6.11) && z<pow(10,6.28)) s_4He_diffuse += pow(z,-15.7167)*pow(10,111.842);
								else s_4He_diffuse+=0;
							}
							if(E_0==30){
								if(z<pow(10,5.15))s_4He_diffuse += pow(z,-2.98899)*pow(10,32.5648);
								else if(z>pow(10,5.15) && z<pow(10,5.41)) s_4He_diffuse += pow(z,-2.20567)*pow(10,28.5335);
								else if(z>pow(10,5.41) && z<pow(10,5.75)) s_4He_diffuse += pow(z,-1.22)*pow(10,23.1918);
								else if(z>pow(10,5.75) && z<pow(10,5.99)) s_4He_diffuse += pow(z,-2.85291)*pow(10,32.5915);
								else if(z>pow(10,5.99) && z<pow(10,6.22)) s_4He_diffuse += pow(z,-5.28026)*pow(10,47.1239);
								else if(z>pow(10,6.22) && z<pow(10,6.31)) s_4He_diffuse += pow(z,-38.35)*pow(10,253.152);
								else s_4He_diffuse+=0;
							}
								// else{
								// 		for(int k = 17; k<19;k++){
								// 		if(k==17)s_4He_diffuse+=2*qsimp(func_deux_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
								// 		if(k==18)s_4He_diffuse+=qsimp(func_deux_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
								// 		}
								// }
						}



						if(avec_correction >=3){
							linearint(vector_z_s2H_produc_trois_iterations,vector_s2H_produc_trois_iterations,vector_z_s2H_produc_trois_iterations.size(),log10(z),y,dy);
							s_4He_diffuse+=pow(10,y);
							// cout << " s_4He = " << s_4He << " s_4He_diffuse produc 3 iterations = " << s_4He_diffuse << endl;
							// cout << " E_c = " << E_c << " E_0 " << E_0 <<endl;
						}
						if(avec_correction==0){
							linearint(vector_z_s2H_produc_quatre_iterations,vector_s2H_produc_quatre_iterations,vector_z_s2H_produc_quatre_iterations.size(),log10(z),y,dy);
							s_4He_diffuse+=pow(10,y);
							// cout << " s_4He = " << s_4He << " s_4He_diffuse produc 4 iterations= " << s_4He_diffuse << endl;
							// cout << " E_c = " << E_c << " E_0 " << E_0 <<endl;
						}


					}
					else{
						if(E_c>Emin[17]){if(z<=pow(10,5.04))s_4He+=pow(z,-2.51795)*pow(10,35.5388);
							else if(z>pow(10,5.04) && z<=pow(10,5.4))s_4He+=pow(z,-2.87902)*pow(10,37.3408);
							else if(z>pow(10,5.4) && z<=pow(10,5.86))s_4He+=pow(z,-4.3513)*pow(10,45.312);
							else if(z>pow(10,5.86) && z<=pow(10,6.11))s_4He+=pow(z,-7.06984)*pow(10,61.2727);
							else if(z>pow(10,6.11) && z<=pow(10,6.2))s_4He+=pow(z,-12.1155)*pow(10,92.1182);
							else if(z>pow(10,6.2) && z<=pow(10,6.3))s_4He+=pow(z,-30.3648)*pow(10,205.645);
							//~ else if(z>pow(10,6.19) && z<=pow(10,6.25))s_4He=pow(z,-7.50809)*pow(10,66.2643);
							//~ else {for(int k = 15; k<19;k++)if(E_c>Emin[k])s_4He=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);}
							else s_4He+=0;
							s_4He = E_0*s_4He/5000;}
					}
	}


	//

	// cout << " s table = " << s_4He << " s_int = " << s << endl;
	Y =exp(-1/(2*H_r*pow(1+z,2)*tau))*Y_He_K_Perte_non_standard(a,z_0,i,z_0,E_0,Z_x, z_0)*(s_4He+s_4He_diffuse)*Y_He_0;
	//~ Y =-B*exp(-1/(2*H_r*pow(1+z,2)*tau))*s*Y_He_0;
	//~ Y = 1;
	//~ cout << " Y = " << Y << "z = " << z << endl;
	//~ cout << " Y_He_K_Perte(z,z_0,i,z_0,E_0,Z_x, z_0) = " << Y_He_K_Perte(z,z_0,i,z_0,E_0,Z_x, z_0) << endl;
	return Y;

}
/************************************************************************************/
/************************Fin Fonctions Deuterium Non Standard************************/
/************************************************************************************/





/************************************************************************************/
/*******************Debut Fonctions Helium-4 (all case)*********************/
/************************************************************************************/
double  K_Perte_4He(double z, double x, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double  Emin[21];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	Emin[19] = 5.5 ; Emin[20] = 7.72;
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  n_y_0 = 3.154*pow(10.,-30) ;
	double y, dy;
	double B= exp(-1/(2*H_r*tau*(z+1)*(z+1)));
	double K, s=0, s_4He = 0,s_4He_diffuse=0;
	//~ cout << " tau = " << tau << " z = " << z <<endl;
	//~ if(E_c>E_min)s_int=qsimp(func,E_min,E_c,z,14,E_0,Z_x, z_0);
// 	if(avec_correction==1){
// 	for(int k=15;k<19;k++)
// 	{
// 		// if(E_c>Emin[k])s_4He+=qsimp_E(func_sans_reinj,Emin[k],E_c,z,k,E_0,Z_x, z_0);
// 		if(E_c>Emin[k])s_4He_diffuse+=qsimp_E(func_standard_helium,Emin[k],E_c,z,k,E_0,Z_x, z_0);
// 		if(s_4He_diffuse<0)cout << "s_4He  = " <<s_4He << "s_4He_diffuse = "<<s_4He_diffuse<<endl;
// 	}
// }
	// if(z<=pow(10,5.26))s_4He=pow(z,-2.5297)*pow(10,36.9046);
	// else if(z>pow(10,5.26) && z<=pow(10,5.62))s_4He=pow(z,-2.84119)*pow(10,38.5537);
	// else if(z>pow(10,5.62) && z<=pow(10,5.81))s_4He=pow(z,-3.8207)*pow(10,44.0669);
	// else if(z>pow(10,5.81) && z<=pow(10,6.))s_4He=pow(z,-4.78069 )*pow(10,49.6438);
	// else if(z>pow(10,6) && z<=pow(10,6.19))s_4He=pow(z,-6.05892)*pow(10,57.309);
	// else if(z>pow(10,6.19) && z<=pow(10,6.39))s_4He=pow(z,-7.50809)*pow(10,66.2643);
	// //~ else if(z>pow(10,6.19) && z<=pow(10,6.25))s_4He=pow(z,-7.50809)*pow(10,66.2643);
	// //~ else {for(int k = 15; k<19;k++)if(E_c>Emin[k])s_4He=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);}
	// else s_4He=0;
	// 	s_4He = E_0*s_4He/5000;
	if(E_c>Emin[15]){
		// for(int k=15;k<19;k++)s_4He_diffuse+=qsimp(func_trois_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
	linearint(vector_z_s4He_destruc_standard,vector_s4He_destruc_standard,vector_z_s4He_destruc_standard.size(),log10(z),y,dy);
		// output_Check << log10(z) << " "<< y << endl;
		s_4He=E_0*pow(10,y)/5000;
		// s_2H_compton = qsimp(func_spectre_compton,E_min,E_0,z,14,E_0,Z_x, z_0);
		// if(s_2H_diffuse <s_2H_compton) cout << "s_2H_compton =  " << s_2H_compton << endl;
		// cout << "s_4He_diffuse quatre iterations=  " << s_4He_diffuse << endl;
	}
	// if(E_c>Emin[15]){
	// 	if(z<=pow(10,5.33))s_4He+=pow(z,-2.51233)*pow(10,36.7075);
	// 	else if(z>pow(10,5.33) && z<=pow(10,5.78))s_4He+=pow(z,-3.2278)*pow(10,40.5342);
	// 	else if(z>pow(10,5.78) && z<=pow(10,6.1))s_4He+=pow(z,-5.14742)*pow(10,51.6394);
	// 	else if(z>pow(10,6.1) && z<=pow(10,6.2))s_4He+=pow(z,-6.04727)*pow(10,57.132);
	// 	else if(z>pow(10,6.2) && z<=pow(10,6.3))s_4He+=pow(z,-7.78083)*pow(10,67.8923);
	// 	else if(z>pow(10,6.3) && z<=pow(10,6.37))s_4He+=pow(z,-13.083)*pow(10,101.287);
	// 	//~ else if(z>pow(10,6.19) && z<=pow(10,6.25))s_4He=pow(z,-7.50809)*pow(10,66.2643);
	// 	//~ else {for(int k = 15; k<19;k++)if(E_c>Emin[k])s_4He=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);}
	// 	else s_4He+=0;
	// 	s_4He = E_0*s_4He/5000;}
	K = B*(s_4He);
	//~ cout << " s = " << s << " s_int = " << s_int << endl;
	//~ cout << " s = " << s << " B = " << B <<  " Y = " << Y << endl;
	return K;
}

double K_perte_non_standard_4He(double z, double x,int i, double g, double E_0,double Z_x,double z_0)
{
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double Emin[19];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double t_inj =  pow((1+z_0),-2)/(2*H_r);
	double tau = 5*t_inj;
	// cout << " tau = " << tau << endl;
	double y,dy;
	double f=0, s_4He =0, s_4He_diffuse =0, s_4He_compton=0;
	// cout << "E_0 " << E_0 << endl;
		if(E_c<=E_0){

			// if(E_c>Emin[15]){
			// 				if(z<=pow(10,5.26))s_4He=pow(z,-2.5297)*pow(10,36.9046);
			// 				else if(z>pow(10,5.26) && z<=pow(10,5.62))s_4He=pow(z,-2.84119)*pow(10,38.5537);
			// 				else if(z>pow(10,5.62) && z<=pow(10,5.81))s_4He=pow(z,-3.8207)*pow(10,44.0669);
			// 				else if(z>pow(10,5.81) && z<=pow(10,6.))s_4He=pow(z,-4.78069 )*pow(10,49.6438);
			// 				else if(z>pow(10,6) && z<=pow(10,6.19))s_4He=pow(z,-6.05892)*pow(10,57.309);
			// 				else if(z>pow(10,6.19) && z<=pow(10,6.39))s_4He=pow(z,-7.50809)*pow(10,66.2643);
			// 				//~ else if(z>pow(10,6.19) && z<=pow(10,6.25))s_4He=pow(z,-7.50809)*pow(10,66.2643);
			// 				//~ else {for(int k = 15; k<19;k++)if(E_c>Emin[k])s_4He=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);}
			// 				else s_4He=0;
			// 				s_4He = E_0*s_4He/5000;
			// 				// if(avec_correction==1){
			// 				// 	for(int k=15;k<19;k++){
			// 				// 		// if(E_c>Emin[k])s_4He+=qsimp_E(func_sans_reinj,Emin[k],E_c,z,k,E_0,Z_x, z_0);
			// 				// 		if(E_c>Emin[k])s_4He_diffuse+=qsimp_E(func_standard_helium,Emin[k],E_c,z,k,E_0,Z_x, z_0);
			// 				// 		// if(s_4He_diffuse>s_4He)	cout << "s_4He  = " <<s_4He << "s_4He_diffuse = "<<s_4He_diffuse<<" standard" <<endl;
			// 				// 	}
			// 				// }
			// }
			if(E_c>Emin[15]){
				// for(int k=15;k<19;k++)s_4He_diffuse+=qsimp(func_trois_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
			linearint(vector_z_s4He_destruc_standard,vector_s4He_destruc_standard,vector_z_s4He_destruc_standard.size(),log10(z),y,dy);
				// output_Check << log10(z) << " "<< y << endl;
				s_4He=E_0*pow(10,y)/5000;
				// s_2H_compton = qsimp(func_spectre_compton,E_min,E_0,z,14,E_0,Z_x, z_0);
				// if(s_2H_diffuse <s_2H_compton) cout << "s_2H_compton =  " << s_2H_compton << endl;
				// cout << "s_4He_diffuse quatre iterations=  " << s_4He_diffuse << endl;
			}
			// cout <<"spectre standard E_c = " << E_c << " z = " << z << endl;
		}
		else{
			// cout << "E_c " <<E_c<<" E_0 " << E_0 << "z " << z << " non universal " << endl;

					for(int k=15;k<19;k++){
						if(E_0>Emin[k]){
							s_4He+=func_sigma(E_0,z,k,g,E_0,Z_x,z_0)/(gamma_NPC(E_0,z,k,g,E_0,Z_x,z_0)+gamma_compton(E_0,z,k,g,E_0,Z_x,z_0)+gamma_phph(E_0,z,k,g,E_0,Z_x,z_0));
							// cout <<"spectre non standard E_c = " << E_c << " z = " << z << endl;
						}
					}
						if(avec_correction>=1){
						/******* FIT used before interpolation was implemented ******/
							if(E_0==70){
							if(z<pow(10,4.92))s_4He_diffuse = pow(z,-2.94973)*pow(10,36.8374);
							else if(z>pow(10,4.92) && z<pow(10,5.47)) s_4He_diffuse = pow(z,-2.17361)*pow(10,33.0111);
							else if(z>pow(10,5.47) && z<pow(10,5.81)) s_4He_diffuse = pow(z,-3.76239)*pow(10,41.7045);
							else if(z>pow(10,5.81) && z<pow(10,6.12)) s_4He_diffuse = pow(z,-5.42405)*pow(10,51.3708);
							else if(z>pow(10,6.12) && z<pow(10,6.3)) s_4He_diffuse = pow(z,-7.23229)*pow(10,62.4409);
							else if(z>pow(10,6.3) && z<pow(10,6.37)) s_4He_diffuse = pow(z,-18.4285)*pow(10,133.054);
							else s_4He_diffuse=0;
						}
						else if(E_0==30){
							if(z<pow(10,4.92))s_4He_diffuse = pow(z,-2.98582)*pow(10,36.9087);
							else if(z>pow(10,4.92) && z<pow(10,5.47)) s_4He_diffuse = pow(z,-2.46896)*pow(10,34.1862);
							else if(z>pow(10,5.47) && z<pow(10,5.81)) s_4He_diffuse = pow(z,-3.88914)*pow(10,42.4291);
							else if(z>pow(10,5.81) && z<pow(10,6.12)) s_4He_diffuse = pow(z,-5.3442)*pow(10,51.1777);
							else if(z>pow(10,6.12) && z<pow(10,6.3)) s_4He_diffuse = pow(z,-8.59309)*pow(10,71.3725);
							else if(z>pow(10,6.3) && z<pow(10,6.37)) s_4He_diffuse = pow(z,-14.6491)*pow(10,109.499);
							else s_4He_diffuse=0;
						}
						// // else {for(int k=15;k<19;k++)if(E_c>Emin[k])s_4He_diffuse+=qsimp_3(func_une_interaction_monochromatique,Emin[k],E_0,z,k,E_0,Z_x, z_0);}
						// 	// if(E_c>Emin[k])s_4He+=qsimp_E(func_sans_reinj,Emin[k],E_c,z,k,E_0,Z_x, z_0);
						// 	// for(int k=15;k<19;k++)s_4He_diffuse+=qsimp_E(func_une_interaction_monochromatique,Emin[k],E_0,z,k,E_0,Z_x, z_0);
						// 	// cout << "s_4He  = " <<s_4He << "s_4He_diffuse = "<<s_4He_diffuse<<" monochrom"<<endl;
						//
						//

						/*********** Interpolate with a polynomial of degree N-1 where N is the number of points**********/
							// polint(vector_z,vector_s4He_destruc,vector_z.size(),log10(z),y,dy);
							// s_4He_diffuse=pow(10,y);
						}
						if(avec_correction>=2){
							// if(E_c>Emin[k])s_4He+=qsimp_E(func_sans_reinj,Emin[k],E_c,z,k,E_0,Z_x, z_0);
							// for(int k=15;k<19;k++)s_4He_diffuse+=qsimp_E(func_deux_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
							linearint(vector_z_s4He_destruc_deux_iterations,vector_s4He_destruc_deux_iterations,vector_z_s4He_destruc_deux_iterations.size(),log10(z),y,dy);
							s_4He_diffuse+=pow(10,y);
							// if(E_0==30){
							// 	if(z<pow(10,5.39))s_4He_diffuse += pow(z,-2.38224)*pow(10,32.6604);
							// 	else if(z>pow(10,5.39) && z<pow(10,5.52)) s_4He_diffuse += pow(z,-1.65083)*pow(10,28.7143);
							// 	else if(z>pow(10,5.52) && z<pow(10,5.72)) s_4He_diffuse += pow(z,-1.16574)*pow(10,26.0349);
							// 	else if(z>pow(10,5.72) && z<pow(10,5.83)) s_4He_diffuse += pow(z,-1.72495)*pow(10,29.2353);
							// 	else if(z>pow(10,5.83) && z<pow(10,5.99)) s_4He_diffuse += pow(z,-2.94351)*pow(10,36.3446);
							// 	else if(z>pow(10,5.99) && z<pow(10,6.22)) s_4He_diffuse += pow(z,-4.82053)*pow(10,47.5996);
							// 	else if(z>pow(10,6.22) && z<pow(10,6.35)) s_4He_diffuse += pow(z,-13.0298)*pow(10,98.687);
							// 	else s_4He_diffuse+=0;
							// }
							// cout << "s_4He  = " <<s_4He << "s_4He_diffuse = "<<s_4He_diffuse<<" deux iterations"<<endl;
						}
						if(avec_correction>=3){
							// for(int k=15;k<19;k++)s_4He_diffuse+=qsimp(func_trois_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
						 linearint(vector_z_s4He_destruc_trois_iterations,vector_s4He_destruc_trois_iterations,vector_z_s4He_destruc_trois_iterations.size(),log10(z),y,dy);
							// output_Check << log10(z) << " "<< y << endl;
							s_4He_diffuse+=pow(10,y);
							// s_2H_compton = qsimp(func_spectre_compton,E_min,E_0,z,14,E_0,Z_x, z_0);
							// if(s_2H_diffuse <s_2H_compton) cout << "s_2H_compton =  " << s_2H_compton << endl;
							// cout << "s_4He_diffuse trois iterations=  " << s_4He_diffuse << endl;
						}
						if(avec_correction==0){
							// for(int k=15;k<19;k++)s_4He_diffuse+=qsimp(func_trois_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
						 linearint(vector_z_s4He_destruc_quatre_iterations,vector_s4He_destruc_quatre_iterations,vector_z_s4He_destruc_quatre_iterations.size(),log10(z),y,dy);
							// output_Check << log10(z) << " "<< y << endl;
							s_4He_diffuse+=pow(10,y);
							// s_2H_compton = qsimp(func_spectre_compton,E_min,E_0,z,14,E_0,Z_x, z_0);
							// if(s_2H_diffuse <s_2H_compton) cout << "s_2H_compton =  " << s_2H_compton << endl;
							// cout << "s_4He_diffuse quatre iterations=  " << s_4He_diffuse << endl;
						}
						if(avec_correction>=4){
							// for(int k=15;k<19;k++)s_4He_diffuse+=qsimp(func_trois_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
						linearint(vector_z_s4He_destruc_quatre_iterations,vector_s4He_destruc_quatre_iterations,vector_z_s4He_destruc_quatre_iterations.size(),log10(z),y,dy);
							// output_Check << log10(z) << " "<< y << endl;
							s_4He_diffuse+=pow(10,y);
							// s_2H_compton = qsimp(func_spectre_compton,E_min,E_0,z,14,E_0,Z_x, z_0);
							// if(s_2H_diffuse <s_2H_compton) cout << "s_2H_compton =  " << s_2H_compton << endl;
							// cout << "s_4He_diffuse quatre iterations=  " << s_4He_diffuse << endl;
						}


		}
				cout << exp(-1/(2*H_r*tau*(z+1)*(z+1))) << " " << s_4He<< endl;

				f = exp(-1/(2*H_r*tau*(z+1)*(z+1)))*(s_4He+s_4He_diffuse+s_4He_compton);
				if(f!=0)cout << "f = " <<f << endl;
				return f;
}
double  func_quatre_interactions_Destruc_Helium4(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double  f;
	double correction, dy;
	linearint(vector_E_correction,vector_correction,vector_E_correction.size(),x,correction,dy);
	f=(Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)+Spectre_deux_interactions(x,z,i,g,E_0,Z_x, z_0)+Spectre_trois_interactions(x,z,i,g,E_0,Z_x, z_0)+Spectre_quatre_interactions(x,z,i,g,E_0,Z_x, z_0))*correction*(func_sigma(x,z,15,g,E_0,Z_x, z_0)+func_sigma(x,z,16,g,E_0,Z_x, z_0)+func_sigma(x,z,17,g,E_0,Z_x, z_0)+func_sigma(x,z,18,g,E_0,Z_x, z_0))/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
	// f=Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0)*correction_phph(x,z,i,g,E_0,Z_x, z_0));
	return f;
}
/************************************************************************************/
/********************Fin Fonctions Helium-4 (all case)**********************/
/************************************************************************************/




/************************************************************************************/
/*************************Debut Fonctions Helium-4 SPECIAL***************************/
/************************************************************************************/

double func_spectre_special(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double T = T_0*(1+z);
	// cout << " T = " << T << endl;
	double f = 0 ;
	double m_e = 0.511;
	double E_gamma_b = 2.7*T_0*(1+z);
	double gamma_e = 4*E_gamma_b*E_0/(m_e*m_e);
	double S = 0, K = 0;
	double q = x/(gamma_e*(E_0-x));
	// if(q>1){
	// 	cout << " gamma_e = " << gamma_e << endl;
	// 	cout << " q  = " << q << " z = " << z << endl;}

	if(q<=1 && q>0)f = (2*q*log(q)+(1+2*q)*(1-q)+(gamma_e*q*gamma_e*q)/(2*(1+gamma_e*q))*(1-q));

		if(T<=pow(10,-4.))S=pow(T,1.84891)*pow(10,8.35556);
		else if(T>pow(10,-4) && T<=pow(10,-3.5))S=pow(T,1.58254)*pow(10,7.29116);
		else if(T>pow(10,-3.5) && T<=pow(10,-3))S=pow(T,1.29604)*pow(10,6.29365);
		else if(T>pow(10,-3) && T<=pow(10,-2.5))S=pow(T,0.949679)*pow(10,5.24496);
		else if(T>pow(10,-2.5) && T<=pow(10,-2))S=pow(T,0.719835)*pow(10,4.67385);
		else if(T>pow(10,-2) && T<=pow(10,-1.5))S=pow(T,0.504243)*pow(10,4.23143);
		else S=0;
		if(S!=0)K = E_0/S;
		// cout << " K = " << K << " S = " << S << endl;
		return x*f*K;
	}
	double func_special(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
{
	double f;
	/**Spectre SPECIAL dégradé SANS réinjection***/
	f=func_spectre_special(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));

	/**Spectre SPECIAL dégradé AVEC réinjection***/
	// f=func_spectre_special(x,z,i,g,E_0)*func_sigma(x,z,i,g,E_0)/(gamma_NPC(x,z,i,g,E_0)+spec_protheroe(x,z,i,g,E_0)*gamma_compton(x,z,i,g,E_0)+gamma_phph(x,z,i,g,E_0)*correction_phph(x,z,i,g,E_0));
	//	cout << f << endl;
	return f;
}

double  func_z_special(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
{

	double  Emin[19];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	double  T_0 = 2.7255*0.862*pow(10.,-10);
	double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double  tau =  5*pow(1+z_0,-2)/(2*H_r);
	double  A,s_4He,s_2H;
	s_4He=0,s_2H=0;
	//~ cout << "E_c = " << E_c <<endl;
	//~ if(E_0 < 3 && E_c>3)E_c == 3;
	double m_e = 0.511;
	double E_gamma_b = 2.7*T_0*(1+z);
	double E_max = 4*E_gamma_b*E_0*E_0/(m_e*m_e*(1+4*E_gamma_b*E_0/(m_e*m_e)));
	//  if(E_0 < 3 && E_c>3)E_c == 3;
	cout << " E_max = " << E_max << " tau " << tau << endl;
	if(E_max>E_c)E_max=E_c;

	for(i=15;i<19;i++)
	{
		if(i>=15 && E_max>Emin[i])
		{
			s_4He+=qsimp_E(func_special,Emin[i],E_max,z,i,E_0,Z_x, z_0);
		}

	}

			A = exp(-1/(2*H_r*tau*(z+1)*(z+1)))*(s_4He);
			//~ cout << " A = " << A << " z = " << z <<  endl;
			return A;
		}
		/************************************************************************************/
		/***************************Fin Fonctions Helium-4 SPECIAL***************************/
		/************************************************************************************/


		/************************************************************************************/
		/*******************Debut Fonctions Helium-3 Standard SANS reinj*********************/
		/************************************************************************************/
		double  func_Y_4He_sans_reinj(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
		{	double  T_0 = 2.7255*0.862*pow(10.,-10);
			double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
			double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
			double  Emin[21];
			Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;Emin[19] = 5.49 ; Emin[20] = 7.72 ;
			double  tau =  5*pow(1+z_0,-2)/(2*H_r);
			double  n_y_0 = 3.154*pow(10.,-30) ;
			double y, dy;

			//~ cout << " tau = " << tau << endl;
			double a=z;
			double B=(Z_x*n_y_0)/(E_0*H_r*tau);
			double Y;
			double Y_He_0 = 0.25,s_4He=0,s_3He=0;
			//~
			//~ for(int k=19;k<21;k++)
			//~ {
			//~ if(E_c>Emin[k])s_3He=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);
			//~ }
			//~ cout << " z = " << z <<endl;
			// if(E_c>Emin[15]){
			// 				if(z<=pow(10,5.26))s_4He=pow(z,-2.5297)*pow(10,36.9046);
			// 				else if(z>pow(10,5.26) && z<=pow(10,5.62))s_4He=pow(z,-2.84119)*pow(10,38.5537);
			// 				else if(z>pow(10,5.62) && z<=pow(10,5.81))s_4He=pow(z,-3.8207)*pow(10,44.0669);
			// 				else if(z>pow(10,5.81) && z<=pow(10,6.))s_4He=pow(z,-4.78069 )*pow(10,49.6438);
			// 				else if(z>pow(10,6) && z<=pow(10,6.19))s_4He=pow(z,-6.05892)*pow(10,57.309);
			// 				else if(z>pow(10,6.19) && z<=pow(10,6.39))s_4He=pow(z,-7.50809)*pow(10,66.2643);
			// 				//~ else if(z>pow(10,6.19) && z<=pow(10,6.25))s_4He=pow(z,-7.50809)*pow(10,66.2643);
			// 				//~ else {for(int k = 15; k<19;k++)if(E_c>Emin[k])s_4He=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);}
			// 				else s_4He=0;
			// 				s_4He = E_0*s_4He/5000;
			// 				// if(avec_correction==1){
			// 				// 	for(int k=15;k<19;k++){
			// 				// 		// if(E_c>Emin[k])s_4He+=qsimp_E(func_sans_reinj,Emin[k],E_c,z,k,E_0,Z_x, z_0);
			// 				// 		if(E_c>Emin[k])s_4He_diffuse+=qsimp_E(func_standard_helium,Emin[k],E_c,z,k,E_0,Z_x, z_0);
			// 				// 		// if(s_4He_diffuse>s_4He)	cout << "s_4He  = " <<s_4He << "s_4He_diffuse = "<<s_4He_diffuse<<" standard" <<endl;
			// 				// 	}
			// 				// }
			// }
			if(E_c>Emin[15]){
				// for(int k=15;k<19;k++)s_4He_diffuse+=qsimp(func_trois_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
			linearint(vector_z_s4He_destruc_standard,vector_s4He_destruc_standard,vector_z_s4He_destruc_standard.size(),log10(z),y,dy);
				// output_Check << log10(z) << " "<< y << endl;
				s_4He=E_0*pow(10,y)/5000;
				// s_2H_compton = qsimp(func_spectre_compton,E_min,E_0,z,14,E_0,Z_x, z_0);
				// if(s_2H_diffuse <s_2H_compton) cout << "s_2H_compton =  " << s_2H_compton << endl;
				// cout << "s_4He_diffuse quatre iterations=  " << s_4He_diffuse << endl;
			}

			if(E_c>Emin[19]){if(z<=pow(10,5.27))s_3He=pow(z,-2.51106)*pow(10,36.673);
				else if(z>pow(10,5.27) && z<=pow(10,5.91))s_3He=pow(z,-2.76477)*pow(10,38.0154);
				else if(z>pow(10,5.91) && z<=pow(10,6.26))s_3He=pow(z,-3.78462)*pow(10,44.0468);
				else if(z>pow(10,6.26) && z<=pow(10,6.5))s_3He=pow(z,-4.72676)*pow(10,49.9511);
				else if(z>pow(10,6.5) && z<=pow(10,6.68))s_3He=pow(z,-5.87075)*pow(10,57.3902);
				else if(z>pow(10,6.68) && z<=pow(10,6.9))s_3He=pow(z,-7.67094)*pow(10,69.4279);
				else s_3He=0;
				s_3He=E_0*s_3He/5000;}
			//~ cout << "s_4He = " << s_4He << " s_3He = " << s_3He << endl;


			Y =exp(-1/(2*H_r*pow(1+z,2)*tau))*(s_4He-s_3He);
			//~ Y =-B*exp(-1/(2*H_r*pow(1+z,2)*tau))*s*Y_He_0;
			//~ Y = 1;
			//~ cout << " Y = " << Y << "z = " << z << " tau =" << tau<< endl;
			//~ cout << " Y_He_K_Perte(z,z_0,i,z_0,E_0,Z_x, z_0) = " << Y_He_K_Perte(z,z_0,i,z_0,E_0,Z_x, z_0) << endl;
			return Y;
		}
		double  S_Gain_3He_sans_reinj(double z, double x,int i, double  g, double  E_0, double  Z_x, double z_0)
		{
			double  Emin[19];
			Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
			double  T_0 = 2.7255*0.862*pow(10.,-10);
			double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
			double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
			double  tau =  5*pow(1+z_0,-2)/(2*H_r);
			double  A,s,s2;
			double  n_y_0 = 3.154*pow(10.,-30) ;
			double B=(Z_x*n_y_0)/(E_0*H_r*tau);
				double s_4He=0;
			double Y,Y_He_0 = 0.25;
			//~ for(int k=16;k<17;k++)
			//~ {
			//~ if(E_c>Emin[k])
			//~ {
			//~ s+=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x, z_0);
			//~
			//~ }
			//~ }
			if(E_c>Emin[15]){if(z<=pow(10,5.34))s_4He=pow(z,-2.51051)*pow(10,36.6711);
				else if(z>pow(10,5.34) && z<=pow(10,5.72))s_4He=pow(z,-3.19543)*pow(10,40.3387);
				else if(z>pow(10,5.72) && z<=pow(10,5.99))s_4He=pow(z,-4.3813)*pow(10,47.1286);
				else if(z>pow(10,5.99) && z<=pow(10,6.2))s_4He=pow(z,-5.91205)*pow(10,56.297);
				else if(z>pow(10,6.2) && z<=pow(10,6.32))s_4He=pow(z,-7.77276)*pow(10,67.8416);
				else if(z>pow(10,6.32) && z<=pow(10,6.38))s_4He=pow(z,-15.9253)*pow(10,119.291);
				else s_4He=0;
				s_4He = E_0*s_4He/5000;}


			Y = -B*qsimp_E(func_Y_4He_sans_reinj,z,z_0,z_0,15,E_0,Z_x, z_0);
			A =exp(-1/(2*H_r*tau*(z+1)*(z+1)))*s_4He*exp(Y)*Y_He_0;
			//~ A =s*exp(Y);
			//~ cout << " A = " << A << " z = " << z << " Y = " << exp(Y) <<" s = " << s << endl;
			return A;
		}

		double  K_Perte_3He_sans_reinj(double z, double x, int i, double  g, double  E_0, double  Z_x, double z_0)
		{
			double  T_0 = 2.7255*0.862*pow(10.,-10);
			double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
			double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
			double  E_min[21];
			E_min[19] = 5.5 ; E_min[20] = 7.72;
			double  tau =  5*pow(1+z_0,-2)/(2*H_r);
			double  n_y_0 = 3.154*pow(10.,-30) ;

			double B= exp(-1/(2*H_r*tau*(z+1)*(z+1)));
			double K, s=0, s_3He = 0;
			//~ cout << " tau = " << tau << " z = " << z <<endl;
			//~ if(E_c>E_min)s_int=qsimp(func,E_min,E_c,z,14,E_0,Z_x, z_0);
			//~ for(int k=19;k<21;k++)
			//~ {
			//~ s+=qsimp(func,E_min[k],E_c,z,k,E_0,Z_x, z_0);
			//~ }
			if(E_c>E_min[19]){if(z<=pow(10,5.27))s_3He=pow(z,-2.51106)*pow(10,36.673);
				else if(z>pow(10,5.27) && z<=pow(10,5.91))s_3He=pow(z,-2.76477)*pow(10,38.0154);
				else if(z>pow(10,5.91) && z<=pow(10,6.26))s_3He=pow(z,-3.78462)*pow(10,44.0468);
				else if(z>pow(10,6.26) && z<=pow(10,6.5))s_3He=pow(z,-4.72676)*pow(10,49.9511);
				else if(z>pow(10,6.5) && z<=pow(10,6.68))s_3He=pow(z,-5.87075)*pow(10,57.3902);
				else if(z>pow(10,6.68) && z<=pow(10,6.9))s_3He=pow(z,-7.67094)*pow(10,69.4279);
				else s_3He=0;
				s_3He=E_0*s_3He/5000;}

			//~ else if(z>pow(10,6.76) && z<=pow(10,6.85))s_3He=pow(z,-7.50809)*pow(10,66.2643);
			//~ else {for(int k = 19; k<21;k++)if(E_c>E_min[k])s_3He=qsimp(func,E_min[k],E_c,z,k,E_0,Z_x, z_0);}

			K = B*s_3He;
			//~ cout << " s = " << s << " s_int = " << s_int << endl;
			//~ cout << " s = " << s << " B = " << B <<  " Y = " << Y << endl;
			return K;
		}
		/************************************************************************************/
		/********************Fin Fonctions Helium-3 Standard SANS reinj**********************/
		/************************************************************************************/



	/************************************************************************************/
	/*******************Debut Fonctions Helium-3 non standard*********************/
	/************************************************************************************/

	double  func_une_interaction_monochromatique_Helium3(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
	{
		double  f;
		f=Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)*(func_sigma(x,z,19,g,E_0,Z_x, z_0)+func_sigma(x,z,20,g,E_0,Z_x, z_0))/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		// f=Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0)*correction_phph(x,z,i,g,E_0,Z_x, z_0));
		return f;
	}
	double  func_quatre_interactions_Destruc_Helium3(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
	{
		double  f;
		double correction, dy;
		linearint(vector_E_correction,vector_correction,vector_E_correction.size(),x,correction,dy);
		f=(Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)+Spectre_deux_interactions(x,z,i,g,E_0,Z_x, z_0)+Spectre_trois_interactions(x,z,i,g,E_0,Z_x, z_0)+Spectre_quatre_interactions(x,z,i,g,E_0,Z_x, z_0))*correction*(func_sigma(x,z,19,g,E_0,Z_x, z_0)+func_sigma(x,z,20,g,E_0,Z_x, z_0))/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		// f=Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0)*correction_phph(x,z,i,g,E_0,Z_x, z_0));
		return f;
	}
	double  func_quatre_interactions_Produc_Helium3(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
	{
		double  f;
		double correction, dy;
		linearint(vector_E_correction,vector_correction,vector_E_correction.size(),x,correction,dy);
		f=(Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)+Spectre_deux_interactions(x,z,i,g,E_0,Z_x, z_0)+Spectre_trois_interactions(x,z,i,g,E_0,Z_x, z_0)+Spectre_quatre_interactions(x,z,i,g,E_0,Z_x, z_0))*correction*(func_sigma(x,z,15,g,E_0,Z_x, z_0)+func_sigma(x,z,16,g,E_0,Z_x, z_0))/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		// f=Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0)*correction_phph(x,z,i,g,E_0,Z_x, z_0));
		return f;
	}

	double  func_Y_4He_non_standard(double  z, double  x,int i, double  g, double  E_0, double  Z_x, double z_0)
	{	double  T_0 = 2.7255*0.862*pow(10.,-10);
		double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
		double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
		double  Emin[21];
		Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;Emin[19] = 5.49 ; Emin[20] = 7.72 ;
		double  tau =  5*pow(1+z_0,-2)/(2*H_r);
		double  n_y_0 = 3.154*pow(10.,-30) ;
		//~ cout << " tau = " << tau << endl;
		double a=z;
		double B=(Z_x*n_y_0)/(E_0*H_r*tau);
		double Y;
		double Y_He_0 = 0.25,s_4He=0,s_4He_diffuse=0,s_3He=0,s_3He_diffuse=0;
		double y, dy;
		//~
		// for(int k=19;k<21;k++)
		// {
		// 	if(E_c>Emin[k])s_3He+=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x,z_0);
		// }
		//~ cout << " z = " << z <<endl;



				if(v2==0){

					if(E_c > E_0){
						for(int k = 15; k<21;k++)
						{
							if(E_0>Emin[k]&& k<19)s_4He+=func_non_standard(E_0,z,k,g,E_0,Z_x,z_0);
							if(k==19 || k ==20)if(E_0>Emin[k])s_3He+=func_non_standard(E_0,z,k,g,E_0,Z_x,z_0);
						}
					}
				}




				else{
					if(E_c > E_0){
							for(int k = 15; k<21;k++)
							{
								if(E_0>Emin[k] && k<19)s_4He+=func_non_standard(E_0,z,k,g,E_0,Z_x,z_0);
								// if(avec_correction >=1 && k<19 && E_0>Emin[k])s_4He_diffuse+=qsimp_3(func_une_interaction_monochromatique,Emin[k],E_0,z,k,E_0,Z_x, z_0);
								if(k>19 && E_0>Emin[k])s_3He+=func_non_standard(E_0,z,k,g,E_0,Z_x,z_0);
							}
							if(avec_correction >=1){
									// s_3He_diffuse+=qsimp_3(func_une_interaction_monochromatique_Helium3,Emin[19],E_0,z,19,E_0,Z_x, z_0);
									if(E_0 == 70 || E_0==30){
									linearint(vector_z_s3He_destruc,vector_s3He_destruc,vector_z_s3He_destruc.size(),log10(z),y,dy);
									s_3He_diffuse=pow(10,y);
									}
								}

							if(avec_correction >=1){
								if(E_0==70){
									if(z<pow(10,4.92))s_4He_diffuse = pow(z,-2.94973)*pow(10,36.8374);
									else if(z>pow(10,4.92) && z<pow(10,5.47)) s_4He_diffuse = pow(z,-2.17361)*pow(10,33.0111);
									else if(z>pow(10,5.47) && z<pow(10,5.81)) s_4He_diffuse = pow(z,-3.76239)*pow(10,41.7045);
									else if(z>pow(10,5.81) && z<pow(10,6.12)) s_4He_diffuse = pow(z,-5.42405)*pow(10,51.3708);
									else if(z>pow(10,6.12) && z<pow(10,6.3)) s_4He_diffuse = pow(z,-7.23229)*pow(10,62.4409);
									else if(z>pow(10,6.3) && z<pow(10,6.37)) s_4He_diffuse = pow(z,-18.4285)*pow(10,133.054);
									else s_4He_diffuse=0;
								}
								else if(E_0==30){
									if(z<pow(10,4.92))s_4He_diffuse = pow(z,-2.98582)*pow(10,36.9087);
									else if(z>pow(10,4.92) && z<pow(10,5.47)) s_4He_diffuse = pow(z,-2.46896)*pow(10,34.1862);
									else if(z>pow(10,5.47) && z<pow(10,5.81)) s_4He_diffuse = pow(z,-3.88914)*pow(10,42.4291);
									else if(z>pow(10,5.81) && z<pow(10,6.12)) s_4He_diffuse = pow(z,-5.3442)*pow(10,51.1777);
									else if(z>pow(10,6.12) && z<pow(10,6.3)) s_4He_diffuse = pow(z,-8.59309)*pow(10,71.3725);
									else if(z>pow(10,6.3) && z<pow(10,6.37)) s_4He_diffuse = pow(z,-14.6491)*pow(10,109.499);
									else s_4He_diffuse=0;
								}
								// else {for(int k = 15; k<19;k++)if(E_0>Emin[k])s_4He_diffuse+=qsimp(func_une_interaction_monochromatique,Emin[k],E_0,z,k,E_0,Z_x, z_0);}
							}
							if(avec_correction >=2){

								if(E_0==70){
									if(z<pow(10,4.97))s_4He_diffuse += pow(z,-2.60856)*pow(10,34.5784);
									else if(z>pow(10,4.97) && z<pow(10,5.14)) s_4He_diffuse += pow(z,-1.57864)*pow(10,29.4573);
									else if(z>pow(10,5.14) && z<pow(10,5.42)) s_4He_diffuse += pow(z,-1.03726)*pow(10,26.6741);
									else if(z>pow(10,5.42) && z<pow(10,5.61)) s_4He_diffuse += pow(z,-2.0937)*pow(10,32.4017);
									else if(z>pow(10,5.61) && z<pow(10,5.86)) s_4He_diffuse += pow(z,-3.37311)*pow(10,39.5797);
									else if(z>pow(10,5.86) && z<pow(10,6.15)) s_4He_diffuse += pow(z,-7.14069)*pow(10,61.6933);
									else if(z>pow(10,6.15) && z<pow(10,6.35)) s_4He_diffuse += pow(z,-10.447)*pow(10,82.0597);
									else s_4He_diffuse+=0;
								}
								else if(E_0==30){
									if(z<pow(10,5.39))s_4He_diffuse += pow(z,-2.38224)*pow(10,32.6604);
									else if(z>pow(10,5.39) && z<pow(10,5.52)) s_4He_diffuse += pow(z,-1.65083)*pow(10,28.7143);
									else if(z>pow(10,5.52) && z<pow(10,5.72)) s_4He_diffuse += pow(z,-1.16574)*pow(10,26.0349);
									else if(z>pow(10,5.72) && z<pow(10,5.83)) s_4He_diffuse += pow(z,-1.72495)*pow(10,29.2353);
									else if(z>pow(10,5.83) && z<pow(10,5.99)) s_4He_diffuse += pow(z,-2.94351)*pow(10,36.3446);
									else if(z>pow(10,5.99) && z<pow(10,6.22)) s_4He_diffuse += pow(z,-4.82053)*pow(10,47.5996);
									else if(z>pow(10,6.22) && z<pow(10,6.35)) s_4He_diffuse += pow(z,-13.0298)*pow(10,98.687);
									else s_4He_diffuse+=0;
								}
								// else {for(int k = 15; k<19;k++)if(E_0>Emin[k])s_4He_diffuse+=qsimp(func_deux_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);}


								if(E_0==70){
									if(z<pow(10,5.06))s_3He_diffuse += pow(z,-2.60815)*pow(10,34.6418);
									else if(z>pow(10,5.06) && z<pow(10,5.54)) s_3He_diffuse += pow(z,-1.46172)*pow(10,28.8192);
									else if(z>pow(10,5.54) && z<pow(10,5.84)) s_3He_diffuse += pow(z,-2.81386)*pow(10,36.3061);
									else if(z>pow(10,5.84) && z<pow(10,6.23)) s_3He_diffuse += pow(z,-4.86971)*pow(10,48.3306);
									else if(z>pow(10,6.23) && z<pow(10,6.64)) s_3He_diffuse += pow(z,-5.96156)*pow(10,55.1417);
									else if(z>pow(10,6.64) && z<pow(10,6.83)) s_3He_diffuse += pow(z,-11.4043)*pow(10,91.2839);
									else if(z>pow(10,6.83) && z<pow(10,6.94)) s_3He_diffuse += pow(z,-20.6525)*pow(10,154.541);
									else s_3He_diffuse+=0;
								}
								else if(E_0==30){
									if(z<pow(10,5.24))s_3He_diffuse += pow(z,-3.02057)*pow(10,36.5253);
									else if(z>pow(10,5.24) && z<pow(10,5.52)) s_3He_diffuse += pow(z,-2.44919)*pow(10,33.5319);
									else if(z>pow(10,5.52) && z<pow(10,5.92)) s_3He_diffuse += pow(z,-1.48115)*pow(10,28.1824);
									else if(z>pow(10,5.92) && z<pow(10,6.21)) s_3He_diffuse += pow(z,-3.34608)*pow(10,39.2339);
									else if(z>pow(10,6.21) && z<pow(10,6.64)) s_3He_diffuse += pow(z,-6.26767)*pow(10,57.3954);
									else if(z>pow(10,6.64) && z<pow(10,6.83)) s_3He_diffuse += pow(z,-11.3226)*pow(10,91.0049);
									else if(z>pow(10,6.83) && z<pow(10,6.94)) s_3He_diffuse += pow(z,-20.5562)*pow(10,154.163);
									else s_3He_diffuse+=0;
								}
								// else {for(int k = 19; k<21;k++)if(E_0>Emin[k])s_3He_diffuse+=qsimp(func_deux_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);}
							}
							if(avec_correction >=3){
								linearint(vector_z_s4He_destruc_trois_iterations,vector_s4He_destruc_trois_iterations,vector_z_s4He_destruc_trois_iterations.size(),log10(z),y,dy);
								s_4He_diffuse+=pow(10,y);
								linearint(vector_z_s3He_destruc_trois_iterations,vector_s3He_destruc_trois_iterations,vector_z_s3He_destruc_trois_iterations.size(),log10(z),y,dy);
								s_3He_diffuse+=pow(10,y);
								// cout << "s_4He  = " <<s_4He << "s_4He_diffuse 3 iterations = "<<s_4He_diffuse<<endl;
								// cout << "s_3He  = " <<s_3He << "s_3He_diffuse 3 iterations= "<<s_3He_diffuse<<endl;
							}


					if(avec_correction==0){
						linearint(vector_z_s4He_destruc_quatre_iterations,vector_s4He_destruc_quatre_iterations,vector_z_s4He_destruc_quatre_iterations.size(),log10(z),y,dy);
						s_4He_diffuse+=pow(10,y);
						linearint(vector_z_s3He_destruc_quatre_iterations,vector_s3He_destruc_quatre_iterations,vector_z_s3He_destruc_quatre_iterations.size(),log10(z),y,dy);
						s_3He_diffuse+=pow(10,y);
						// cout << "s_4He  = " <<s_4He << "s_4He_diffuse 3 iterations= "<<s_4He_diffuse<<endl;
						// cout << "s_3He  = " <<s_3He << "s_3He_diffuse 3 iterations= "<<s_3He_diffuse<<endl;
					}
						}

					else {
						// if(E_c>Emin[15]){
						// 	if(z<=pow(10,5.33))s_4He+=pow(z,-2.51233)*pow(10,36.7075);
						// 	else if(z>pow(10,5.33) && z<=pow(10,5.78))s_4He+=pow(z,-3.2278)*pow(10,40.5342);
						// 	else if(z>pow(10,5.78) && z<=pow(10,6.1))s_4He+=pow(z,-5.14742)*pow(10,51.6394);
						// 	else if(z>pow(10,6.1) && z<=pow(10,6.2))s_4He+=pow(z,-6.04727)*pow(10,57.132);
						// 	else if(z>pow(10,6.2) && z<=pow(10,6.3))s_4He+=pow(z,-7.78083)*pow(10,67.8923);
						// 	else if(z>pow(10,6.3) && z<=pow(10,6.37))s_4He+=pow(z,-13.083)*pow(10,101.287);
						// 	//~ else if(z>pow(10,6.19) && z<=pow(10,6.25))s_4He=pow(z,-7.50809)*pow(10,66.2643);
						// 	//~ else {for(int k = 15; k<19;k++)if(E_c>Emin[k])s_4He=qsimp(func,Emin[k],E_c,z,k,E_0,Z_x,z_0);}
						// 	else s_4He+=0;
						// 	s_4He = E_0*s_4He/5000;}
						if(E_c>Emin[15]){
							// for(int k=15;k<19;k++)s_4He_diffuse+=qsimp(func_trois_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);
						linearint(vector_z_s4He_destruc_standard,vector_s4He_destruc_standard,vector_z_s4He_destruc_standard.size(),log10(z),y,dy);
							// output_Check << log10(z) << " "<< y << endl;
							s_4He=E_0*pow(10,y)/5000;
							// s_2H_compton = qsimp(func_spectre_compton,E_min,E_0,z,14,E_0,Z_x, z_0);
							// if(s_2H_diffuse <s_2H_compton) cout << "s_2H_compton =  " << s_2H_compton << endl;
							// cout << "s_4He_diffuse quatre iterations=  " << s_4He_diffuse << endl;
						}

							if(E_c>Emin[19]){if(z<=pow(10,5.27))s_3He=pow(z,-2.51106)*pow(10,36.673);
								else if(z>pow(10,5.27) && z<=pow(10,5.91))s_3He=pow(z,-2.76477)*pow(10,38.0154);
								else if(z>pow(10,5.91) && z<=pow(10,6.26))s_3He=pow(z,-3.78462)*pow(10,44.0468);
								else if(z>pow(10,6.26) && z<=pow(10,6.5))s_3He=pow(z,-4.72676)*pow(10,49.9511);
								else if(z>pow(10,6.5) && z<=pow(10,6.68))s_3He=pow(z,-5.87075)*pow(10,57.3902);
								else if(z>pow(10,6.68) && z<=pow(10,6.9))s_3He=pow(z,-7.67094)*pow(10,69.4279);
								else s_3He=0;
								s_3He=E_0*s_3He/5000;}
							}
						}
				Y =exp(-1/(2*H_r*pow(1+z,2)*tau))*(s_4He+s_4He_diffuse-s_3He-s_3He_diffuse);
				//~ Y =-B*exp(-1/(2*H_r*pow(1+z,2)*tau))*s*Y_He_0;
				//~ Y = 1;
				//~ cout << " Y = " << Y << "z = " << z << " tau =" << tau<< endl;
				//~ cout << " Y_He_K_Perte(z,z_0,i,z_0,E_0,Z_x) = " << Y_He_K_Perte(z,z_0,i,z_0,E_0,Z_x) << endl;
				return Y;
			}
			double  S_Gain_3He_non_standard(double  z, double  x,int i, double  g, double  E_0, double  Z_x, double z_0)
			{
			double  Emin[19];
			Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
			double  T_0 = 2.7255*0.862*pow(10.,-10);
			double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
			double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
			double  tau =  5*pow(1+z_0,-2)/(2*H_r);
			double  A,s,s2;
			double  n_y_0 = 3.154*pow(10.,-30) ;
			double B=(Z_x*n_y_0)/(E_0*H_r*tau);
			double s_4He=0,s_4He_diffuse=0;
			double Y,Y_He_0 = 0.25;
			double y, dy;


			if(v2==0){
				if(E_0<E_c){
					for(int k=15;k<17;k++){
						if(E_0>Emin[k]){
							s_4He+=func_non_standard(E_0,z,k,g,E_0,Z_x,z_0);
							// cout << " s = " << s << endl;
						}
					}
				}
			}
			else{
			if(E_0<E_c){
			for(int k=15;k<17;k++){
			// 	if(E_0>Emin[k]){
					s_4He+=func_non_standard(E_0,z,k,g,E_0,Z_x,z_0);
			// 		if(avec_correction >=1){s_4He_diffuse+=qsimp_3(func_une_interaction_monochromatique,Emin[k],E_0,z,k,E_0,Z_x, z_0);}
				}
			// }
			// 		cout << "s_4He_diffuse vraie production = " << s_4He_diffuse << endl;
					if(avec_correction >=1){
						if(E_0 == 70 || E_0==30){
						linearint(vector_z_s3He_produc,vector_s3He_produc,vector_z_s3He_produc.size(),log10(z),y,dy);
						s_4He_diffuse=pow(10,y);
						}
					// cout << "s_4He_diffuse interpolation production= " << s_4He_diffuse << endl;
					}

			// cout << "s_4He  = " <<s_4He << "s_4He_diffuse 1 iteration = "<<s_4He_diffuse<<endl;
			if(avec_correction >=2){
				// for(int k=15;k<17;k++)s_4He_diffuse+=qsimp_3(func_deux_interactions,Emin[k],E_0,z,k,E_0,Z_x, z_0);

											if(E_0==70){
												if(z<pow(10,4.97))s_4He_diffuse += pow(z,-2.85052)*pow(10,35.7447);
												else if(z>pow(10,4.97) && z<pow(10,5.16)) s_4He_diffuse += pow(z,-1.65059)*pow(10,29.8121);
												else if(z>pow(10,5.16) && z<pow(10,5.44)) s_4He_diffuse += pow(z,-1.03743)*pow(10,26.6457);
												else if(z>pow(10,5.44) && z<pow(10,5.67)) s_4He_diffuse += pow(z,-2.37186)*pow(10,33.9069);
												else if(z>pow(10,5.67) && z<pow(10,5.85)) s_4He_diffuse += pow(z,-3.72351)*pow(10,41.5778);
												else if(z>pow(10,5.85) && z<pow(10,6.1)) s_4He_diffuse += pow(z,-6.64222)*pow(10,58.6734);
												else if(z>pow(10,6.1) && z<pow(10,6.2)) s_4He_diffuse += pow(z,-8.49589)*pow(10,69.9916);
												else if(z>pow(10,6.2) && z<pow(10,6.38)) s_4He_diffuse += pow(z,-12.9154)*pow(10,97.501);
												else s_4He_diffuse+=0;
											}
											if(E_0==30){
												if(z<pow(10,5.11))s_4He_diffuse += pow(z,-2.99831)*pow(10,35.8277);
												else if(z>pow(10,5.11) && z<pow(10,5.41)) s_4He_diffuse += pow(z,-2.45589)*pow(10,33.0502);
												else if(z>pow(10,5.41) && z<pow(10,5.79)) s_4He_diffuse += pow(z,-1.37872)*pow(10,27.2266);
												else if(z>pow(10,5.79) && z<pow(10,5.99)) s_4He_diffuse += pow(z,-2.50366)*pow(10,33.738);
												else if(z>pow(10,5.99) && z<pow(10,6.22)) s_4He_diffuse += pow(z,-5.00479)*pow(10,48.7279);
												else if(z>pow(10,6.22) && z<pow(10,6.35)) s_4He_diffuse += pow(z,-13.1222)*pow(10,99.2742);
												else s_4He_diffuse+=0;
											}

			}
			if(avec_correction>=3){
											if(E_0 == 70 || E_0==30){
											linearint(vector_z_s3He_produc_trois_iterations,vector_s3He_produc_trois_iterations,vector_z_s3He_produc_trois_iterations.size(),log10(z),y,dy);
											s_4He_diffuse+=pow(10,y);

											}
											// cout << "s_4He  = " <<s_4He << "s_4He_diffuse produc 3 iterations = "<<s_4He_diffuse<<endl;

			}
			if(avec_correction==0){
											if(E_0 == 70 || E_0==30){
											linearint(vector_z_s3He_produc_quatre_iterations,vector_s3He_produc_quatre_iterations,vector_z_s3He_produc_quatre_iterations.size(),log10(z),y,dy);
											s_4He_diffuse+=pow(10,y);

											}
											// cout << "s_4He  = " <<s_4He << "s_4He_diffuse produc 4 iterations = "<<s_4He_diffuse<<endl;

			}
			}
								else{

										/***Taking only 3He production into account***/
										// if(E_c>Emin[16]){if(z<=pow(10,5.34))s_4He=pow(z,-2.51107)*pow(10,36.3434);
										// 	else if(z>pow(10,5.34) && z<=pow(10,5.72))s_4He=pow(z,-3.22559)*pow(10,40.1692);
										// 	else if(z>pow(10,5.72) && z<=pow(10,5.99))s_4He=pow(z,-4.43413)*pow(10,47.0873);
										// 	else if(z>pow(10,5.99) && z<=pow(10,6.2))s_4He=pow(z,-5.99072)*pow(10,56.4142);
										// 	else if(z>pow(10,6.2) && z<=pow(10,6.32))s_4He=pow(z,-8.1021)*pow(10,69.5151);
										// 	else if(z>pow(10,6.32) && z<=pow(10,6.38))s_4He=pow(z,-18.9133)*pow(10,137.759);
										// 	else s_4He=0;
										// 	s_4He = E_0*s_4He/5000;}
										/***Taking also t production into account***/
										if(E_c>Emin[15]){if(z<=pow(10,5.34))s_4He=pow(z,-2.51051)*pow(10,36.6711);
											else if(z>pow(10,5.34) && z<=pow(10,5.72))s_4He=pow(z,-3.19543)*pow(10,40.3387);
											else if(z>pow(10,5.72) && z<=pow(10,5.99))s_4He=pow(z,-4.3813)*pow(10,47.1286);
											else if(z>pow(10,5.99) && z<=pow(10,6.2))s_4He=pow(z,-5.91205)*pow(10,56.297);
											else if(z>pow(10,6.2) && z<=pow(10,6.32))s_4He=pow(z,-7.77276)*pow(10,67.8416);
											else if(z>pow(10,6.32) && z<=pow(10,6.38))s_4He=pow(z,-15.9253)*pow(10,119.291);
											else s_4He=0;
											s_4He = E_0*s_4He/5000;}


								}
			}
			Y = -B*qsimp_E(func_Y_4He_non_standard,z,z_0,z_0,15,E_0,Z_x,z_0);
			// A =exp(-1/(2*H_r*tau*(z+1)*(z+1)))*s_4He*Y_He_0;
			A =exp(-1/(2*H_r*tau*(z+1)*(z+1)))*(s_4He+s_4He_diffuse)*exp(Y)*Y_He_0;
			//~ A =s*exp(Y);
			cout << " A = " << A << " z = " << z << " Y = " << exp(Y) <<" s_4He = " << s_4He << endl;
			return A;
			}

		double  K_Perte_3He_non_standard(double  z, double  x, int i, double  g, double  E_0, double  Z_x, double z_0)
	{
		double  T_0 = 2.7255*0.862*pow(10.,-10);
		double  H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
		double  E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
		double  E_min[21];
		E_min[19] = 5.5 ; E_min[20] = 7.72;
		double  tau =  5*pow(1+z_0,-2)/(2*H_r);
		double  n_y_0 = 3.154*pow(10.,-30) ;
		double y, dy;

		double B= exp(-1/(2*H_r*tau*(z+1)*(z+1)));
		double K, s=0, s_3He = 0, s_3He_diffuse = 0;

			if(v2==0){
				if(E_0<E_c){
					for(int k=19;k<21;k++){
						if(E_0>E_min[k]){
							s_3He+=func_non_standard(E_0,z,k,g,E_0,Z_x,z_0);
							// cout << " s = " << s << endl;
						}
					}
				}
			}


			else{
				if(E_0<E_c){
					for(int k=19;k<21;k++){
						if(E_0>E_min[k]){
							s_3He+=func_non_standard(E_0,z,k,g,E_0,Z_x,z_0);
							// if(avec_correction >=1){
							// 	s_3He_diffuse+=qsimp_E(func_une_interaction_monochromatique,E_min[k],E_0,z,k,E_0,Z_x, z_0);
							// 	cout << "s_3He_diffuse vraie = " << s_3He_diffuse << endl;
							// 	}
							// cout << " s = " << s << endl;
						}
					}

					if(avec_correction >=1){
						if(E_0 == 70 || E_0==30){
						linearint(vector_z_s3He_destruc,vector_s3He_destruc,vector_z_s3He_destruc.size(),log10(z),y,dy);
						s_3He_diffuse=pow(10,y);
						}
						// cout << "s_3He_diffuse interpolation= " << s_3He_diffuse << endl;
					}
					if(avec_correction >=2){
						// for(int k=19;k<21;k++) s_3He_diffuse+=qsimp_E(func_deux_interactions,E_min[k],E_0,z,k,E_0,Z_x, z_0);

					if(E_0==70){
						if(z<pow(10,5.06))s_3He_diffuse += pow(z,-2.60815)*pow(10,34.6418);
						else if(z>pow(10,5.06) && z<pow(10,5.54)) s_3He_diffuse += pow(z,-1.46172)*pow(10,28.8192);
						else if(z>pow(10,5.54) && z<pow(10,5.84)) s_3He_diffuse += pow(z,-2.81386)*pow(10,36.3061);
						else if(z>pow(10,5.84) && z<pow(10,6.23)) s_3He_diffuse += pow(z,-4.86971)*pow(10,48.3306);
						else if(z>pow(10,6.23) && z<pow(10,6.64)) s_3He_diffuse += pow(z,-5.96156)*pow(10,55.1417);
						else if(z>pow(10,6.64) && z<pow(10,6.83)) s_3He_diffuse += pow(z,-11.4043)*pow(10,91.2839);
						else if(z>pow(10,6.83) && z<pow(10,6.94)) s_3He_diffuse += pow(z,-20.6525)*pow(10,154.541);
						else s_3He_diffuse+=0;
					}
					else if(E_0==30){
						if(z<pow(10,5.24))s_3He_diffuse += pow(z,-3.02057)*pow(10,36.5253);
						else if(z>pow(10,5.24) && z<pow(10,5.52)) s_3He_diffuse += pow(z,-2.44919)*pow(10,33.5319);
						else if(z>pow(10,5.52) && z<pow(10,5.92)) s_3He_diffuse += pow(z,-1.48115)*pow(10,28.1824);
						else if(z>pow(10,5.92) && z<pow(10,6.21)) s_3He_diffuse += pow(z,-3.34608)*pow(10,39.2339);
						else if(z>pow(10,6.21) && z<pow(10,6.64)) s_3He_diffuse += pow(z,-6.26767)*pow(10,57.3954);
						else if(z>pow(10,6.64) && z<pow(10,6.83)) s_3He_diffuse += pow(z,-11.3226)*pow(10,91.0049);
						else if(z>pow(10,6.83) && z<pow(10,6.94)) s_3He_diffuse += pow(z,-20.5562)*pow(10,154.163);
						else s_3He_diffuse+=0;
					}
				}

				if(avec_correction>=3){
					if(E_0 == 70 || E_0==30){
					linearint(vector_z_s3He_destruc_trois_iterations,vector_s3He_destruc_trois_iterations,vector_z_s3He_destruc_trois_iterations.size(),log10(z),y,dy);
					s_3He_diffuse+=pow(10,y);
					}
				}
					// cout << "s_3He  = " <<s_3He << "s_3He_diffuse 3 iterations = "<<s_3He_diffuse<<endl;
					if(avec_correction==0){
						if(E_0 == 70 || E_0==30){
						linearint(vector_z_s3He_destruc_quatre_iterations,vector_s3He_destruc_quatre_iterations,vector_z_s3He_destruc_quatre_iterations.size(),log10(z),y,dy);
						s_3He_diffuse+=pow(10,y);
						}
					}
						// cout << "s_3He  = " <<s_3He << "s_3He_diffuse 4 iterations = "<<s_3He_diffuse<<endl;
				}
				else{
					if(E_c>E_min[19]){if(z<=pow(10,5.27))s_3He=pow(z,-2.51106)*pow(10,36.673);
						else if(z>pow(10,5.27) && z<=pow(10,5.91))s_3He=pow(z,-2.76477)*pow(10,38.0154);
						else if(z>pow(10,5.91) && z<=pow(10,6.26))s_3He=pow(z,-3.78462)*pow(10,44.0468);
						else if(z>pow(10,6.26) && z<=pow(10,6.5))s_3He=pow(z,-4.72676)*pow(10,49.9511);
						else if(z>pow(10,6.5) && z<=pow(10,6.68))s_3He=pow(z,-5.87075)*pow(10,57.3902);
						else if(z>pow(10,6.68) && z<=pow(10,6.9))s_3He=pow(z,-7.67094)*pow(10,69.4279);
						else s_3He=0;
						s_3He=E_0*s_3He/5000;}
					}
				}
			K = B*(s_3He+s_3He_diffuse);
			//~ cout << " s = " << s << " s_int = " << s_int << endl;
			//~ cout << " s = " << s << " B = " << B <<  " Y = " << Y << endl;
			return K;
		}
		/************************************************************************************/
		/********************Fin Fonctions Helium-3 non standard**********************/
		/************************************************************************************/

		/****************************************************************************************/
		/******************Debut Fonctions pour etudes des differents spectres*******************/
		/****************************************************************************************/
		double integrale_deux_interactions(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
		{
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		return Spectre_deux_interactions(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		// return x*Spectre_deux_interactions(x,z,i,g,E_0,Z_x, z_0);
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

		}
		double integrale_trois_interactions(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
		{
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		return Spectre_trois_interactions(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		// return x*Spectre_trois_interactions(x,z,i,g,E_0,Z_x, z_0);
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

		}
		double integrale_quatre_interactions(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
		{
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
			double y,dy;
		linearint(vector_E,vector_Spectre,vector_E.size(),x,y,dy);
		double Spectre_quatre_interactions = y;
		return Spectre_quatre_interactions;
		// return Spectre_quatre_interactions/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		// return Spectre_quatre_interactions(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

		}
		double integrale_cinq_interactions(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
		{
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		double y,dy;
		linearint(vector_E_cinq_iterations,vector_Spectre_cinq_iterations,vector_E_cinq_iterations.size(),x,y,dy);
	double Spectre_cinq_interactions = y;
	return Spectre_cinq_interactions;
	// 	return Spectre_cinq_interactions/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		// return Spectre_cinq_interactions(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

		}
		double integrale_six_interactions(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
		{
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		return Spectre_six_interactions(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

		}
		double integrale_sept_interactions(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
		{
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		double y,dy;
		// linearint(vector_E_sept_iterations,vector_Spectre_sept_iterations,vector_E_sept_iterations.size(),x,y,dy);
		// double Spectre_sept_iterations = y;
		// return Spectre_sept_iterations;
		return Spectre_sept_interactions(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

		}
		double integrale_sept_interactions_table(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
		{
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		double y,dy;
		linearint(vector_E_sept_iterations,vector_Spectre_sept_iterations,vector_E_sept_iterations.size(),x,y,dy);
		double Spectre_sept_iterations = y;
		return Spectre_sept_iterations;
		// return Spectre_sept_interactions(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
		// return func_spectre(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_tot_phph(x,z,i,g,E_0,Z_x, z_0));

		}
		double  integrale_une_interaction_monochromatique(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
		{
			double  f;
			f= Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
			// f= x*Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0);
			// f=Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0)*correction_phph(x,z,i,g,E_0,Z_x, z_0));
			return f;
		}
		double  integrale_spectre_standard(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
		{
			double  f;
			f= func_spectre(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
			// f=Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0)*correction_phph(x,z,i,g,E_0,Z_x, z_0));
			return f;
		}
		double  integrale_spectre_electron_compton(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
		{
				double  f=0;
				double T_0 = 2.7255*0.862*pow(10.,-10);
				double T = T_0*(1+z);
				// cout << " T " << T << endl;
				double m_e=0.511;
				double alpha = 1./137;
				double r_e = alpha/m_e;
				double pi = 3.14159;
				double spectre_compton = 0;
				double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
				double Gamma_electron;
				double E_gamma_b = 2.701*T_0*(1+z);
				Gamma_electron = 2*pi*r_e*r_e*m_e*m_e/(x*x)*qsimp_E(Integrand_gamma_e,0.0001,1,z,i,E_0,x, z_0)*int_bb/E_gamma_b;
				cout << " Gamma_electron " << Gamma_electron << endl;
				if(Gamma_electron>0)f= Spectre_electron_compton(x,z,i,g,E_0,Z_x, z_0)/(Gamma_electron);
			// f=Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0)*correction_phph(x,z,i,g,E_0,Z_x, z_0));
			return f;
		}
		double  integrale_spectre_gamma_compton(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
		{
			double  f;
			f= Spectre_Gamma_compton(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
			// f=Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0)*correction_phph(x,z,i,g,E_0,Z_x, z_0));
			return f;
		}
		double  spectre_dirac(double  x, double  z, int i, double  g, double  E_0, double  Z_x, double z_0)
		{
			double  f=0;
			int A=E_0;
			int X = x ;
			// cout << " x = " << x << " E_0 = " << E_0 << endl;
			if(E_0==x){f= 1./(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0));
				cout << " f = " << f << endl;
			}
			// f=Spectre_une_interaction_monochromatique(x,z,i,g,E_0,Z_x, z_0)*func_sigma(x,z,i,g,E_0,Z_x, z_0)/(gamma_NPC(x,z,i,g,E_0,Z_x, z_0)+gamma_compton(x,z,i,g,E_0,Z_x, z_0)*spec_protheroe(x,z,i,g,E_0,Z_x, z_0)+gamma_phph(x,z,i,g,E_0,Z_x, z_0)*correction_phph(x,z,i,g,E_0,Z_x, z_0));
			return f;
		}
		/****************************************************************************************/
		/******************Fin Fonctions pour etudes des differents spectres*******************/
		/****************************************************************************************/
