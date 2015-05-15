#include "../Include/EM_cascade.h"
#include "../Include/Injected_spectrum.h"
#include "../Include/tools.h"
double  dsigma_compton(double  x, double  z, double g){

	double  dsigma = pi*pow(r_e,2)*m_e/pow(x,2)*(x/g+g/x+pow(m_e/g-m_e/x,2)-2*m_e*(1/g-1/x))*n_e*pow(1+z,3);
	return dsigma;

}
double  dsigma_phph(double  x, double  z,  double g){

	double  dsigma = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(a*r_e,2)*pow(m_e,-6)*pow(x,2)*pow(1-g/x+pow(g/x,2),2)/(63*10125*pi);
	return dsigma;

}
double  gamma_compton(double  x, double  z){

	double X = 2*x/m_e;
	double sigma_cs = 2*pi*pow(r_e,2)/X*((1-4/X-8/pow(X,2))*log(1+X)+1/2+8/X+1/(2*pow(1+X,2)));
	double Gamma = sigma_cs*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3);
	return Gamma;

}
double  gamma_NPC(double  x, double  z){

	double  k = x/m_e;
	double  rho = (2*k-4)/(k+2+2*pow(2*k,0.5));
	double  sigma_PCN;
	sigma_PCN = a*pow(r_e,2)*(28/9*log(2*k)-218/27) ;
	double  Gamma_2 = sigma_PCN*pow(1+z,3)*eta*n_y_0;
	return Gamma_2;

}

double  gamma_phph(double  x, double  z){

	double  Gamma_3;
	Gamma_3 = 0.1513*pow(a,4)*m_e*pow(x/m_e,3)*pow(T_0*(1+z)/m_e,6);
	return Gamma_3;

}

double Dirac_Spectrum_After_One_Iteration(double  x, double  z, double E_0){


	double Gamma_tot = gamma_NPC(E_0,z)+gamma_compton(E_0,z)+gamma_phph(E_0,z);

	double T = T_0*(1+z);
	double int_BB = 8./63.*pow(pi,4)*pow(T_0*(1+z),6);

	double spectre_gamma_gamma = 1112./(10125*pi)*pow(a*r_e,2)*pow(m_e,-6)*pow(E_0,2)*pow(1-x/E_0+pow(x/E_0,2),2)*int_BB;
	double spectre_compton = pi*r_e*r_e*m_e*pow(E_0,-2)*(x/E_0+E_0/x+pow(m_e/x-m_e/E_0,2)-2*m_e*(1/x-1/E_0))*n_e*pow(1+z,3);
	double f = (spectre_gamma_gamma+spectre_compton)/(Gamma_tot);
	if(x>E_0)f=0;
	return f;

}
/*
void  Cascade_Spectrum_Calculation_From_File(ifstream &file, struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum, long n_step)
{
	double dE, h;
	double E1, E2, E3;
	double f1, f2, f3;

    dE = (pt_Injected_Spectrum->E_0 - pt_Injected_Spectrum->E_min)/ (double) n_step;
    h = dE/2.;
		//Initialization
    pt_Cascade_Spectrum->Gamma_Spectrum.resize(n_step);
    for(int i=n_step;i>=0;i--){

      if(i==n_step){
        E1=pt_Injected_Spectrum->E_0;
        }
      else{
        E1=E3;
      }

      E2=E1 - h;
      E3=E2 - h;
      linearint(pt_Injected_Spectrum->Gamma_Energy, pt_Injected_Spectrum->Gamma_Spectrum, pt_Injected_Spectrum->Gamma_Energy.size(), E1, f1)
      *(dsigma_phph(E1,z,k,g,pt_Injected_Spectrum)+dsigma_compton(E1,z,k,g,pt_Injected_Spectrum))/(gamma_NPC(E1,z,k,g,pt_Injected_Spectrum)+gamma_compton(E1,z,k,g,pt_Injected_Spectrum)+gamma_phph(E1,z,k,g,pt_Injected_Spectrum));
      linearint(pt_Injected_Spectrum->Gamma_Energy, pt_Injected_Spectrum->Gamma_Spectrum, pt_Injected_Spectrum->Gamma_Energy.size(), E2, f2)
      *(dsigma_phph(E2,z,k,g,pt_Injected_Spectrum)+dsigma_compton(E2,z,k,g,pt_Injected_Spectrum))/(gamma_NPC(E2,z,k,g,pt_Injected_Spectrum)+gamma_compton(E2,z,k,g,pt_Injected_Spectrum)+gamma_phph(E2,z,k,g,pt_Injected_Spectrum));
      linearint(pt_Injected_Spectrum->Gamma_Energy, pt_Injected_Spectrum->Gamma_Spectrum, pt_Injected_Spectrum->Gamma_Energy.size(), E3, f3)
      *(dsigma_phph(E3,z,k,g,pt_Injected_Spectrum)+dsigma_compton(E3,z,k,g,pt_Injected_Spectrum))/(gamma_NPC(E3,z,k,g,pt_Injected_Spectrum)+gamma_compton(E3,z,k,g,pt_Injected_Spectrum)+gamma_phph(E3,z,k,g,pt_Injected_Spectrum))
      resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
      pt_Cascade_Spectrum->Gamma_Spectrum[i]+=f1;
      for(int j =0; j<=i)pt_Cascade_Spectrum->Gamma_Spectrum[j]+=resultat;
    }
}
*/
void  Cascade_Spectrum_Calculation_From_Function(double (*func)(double,double,double),double z, struct Structure_Particle_Physics_Model * pt_Particle_Model, struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum, long n_step, int iterations){
	double dE, h;
	double E1, E2, E3, E_gamma;
	double f1, f2, f3;
	double resultat;
	pt_Cascade_Spectrum->redshift = z;
    dE = (pt_Particle_Model->E_0 - E_min)/ (double) n_step;
    h = dE/2.;

    /*****Initialization****/
		pt_Cascade_Spectrum->Gamma_Energy.resize(Gamma_Table_Size);
		pt_Cascade_Spectrum->Gamma_Spectrum.resize(Gamma_Table_Size);
		cout<<"iteration : 0 " << endl;
		for(int i=0;i<Gamma_Table_Size;i++){
			E1=E_min+i*dE;
			pt_Cascade_Spectrum->Gamma_Energy[i]=E1;
			pt_Cascade_Spectrum->Gamma_Spectrum[i]=(*func)(E1,z,pt_Particle_Model->E_0);
			// cout << " E1 = " << E1 << "spectre = " << pt_Cascade_Spectrum->Gamma_Spectrum[i] << endl;
			print_spectrum_automatic_names(0, pt_Cascade_Spectrum, pt_Particle_Model);

		}

    // pt_Cascade_Spectrum->Gamma_Spectrum.resize(n_step);
		for(int k = 0; k<iterations;k++){
						cout<<"iteration : " << k+1 << endl;
						for(int j =0; j<Gamma_Table_Size;j++){
							resultat=0;

					    for(int i=j;i<=n_step;i++){
					      if(i==j){
					        E1=E_min+i*dE;
									E_gamma = E1;
					        }
					      else{
					        E1=E3;
					      }

					      E2=E1 + h;
					      E3=E2 + h;

								linearint(pt_Cascade_Spectrum->Gamma_Energy, pt_Cascade_Spectrum->Gamma_Spectrum, pt_Cascade_Spectrum->Gamma_Energy.size(), E1, f1);
								linearint(pt_Cascade_Spectrum->Gamma_Energy, pt_Cascade_Spectrum->Gamma_Spectrum, pt_Cascade_Spectrum->Gamma_Energy.size(), E2, f2);
								linearint(pt_Cascade_Spectrum->Gamma_Energy, pt_Cascade_Spectrum->Gamma_Spectrum, pt_Cascade_Spectrum->Gamma_Energy.size(), E3, f3);
								f1 *=(dsigma_phph(E1,z,E_gamma)+dsigma_compton(E1,z,E_gamma))/(gamma_NPC(E1,z)+gamma_compton(E1,z)+gamma_phph(E1,z));
								f2 *=(dsigma_phph(E2,z,E_gamma)+dsigma_compton(E2,z,E_gamma))/(gamma_NPC(E2,z)+gamma_compton(E2,z)+gamma_phph(E2,z));
								f3 *=(dsigma_phph(E3,z,E_gamma)+dsigma_compton(E3,z,E_gamma))/(gamma_NPC(E3,z)+gamma_compton(E3,z)+gamma_phph(E3,z));

					      resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
					    }
							pt_Cascade_Spectrum->Gamma_Spectrum[j]+=resultat;
						// cout << " E1 = " << E_min+j*dE << "spectre = " << pt_Cascade_Spectrum->Gamma_Spectrum[j] << endl;
						}
						print_spectrum_automatic_names(k+1, pt_Cascade_Spectrum, pt_Particle_Model);

	}
}




double  cross_section(double  x, int i)
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
