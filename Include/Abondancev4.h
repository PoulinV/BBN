double gamma_tot_compton(double x, double z, int i, double g, double E_0);
double gamma_tot_compton_2(double x, double z, int i, double g, double E_0);
double gamma_tot_compton_3(double x, double z, int i, double g, double E_0);
double gamma_tot_phph(double x, double z, int i, double g, double E_0);
double gamma_tot_phph_2(double x, double z, int i, double g, double E_0);
double gamma_tot_phph_3(double x, double z, int i, double g, double E_0);

//g correspond au stockage de x_0 pour la prise en compte des gains dans le taux d'interaction relatif au photon
/*******Routine pour integrer, issu de Numerical Recipes*******/
double v2=1;
double spec_protheroe(double x, double z, int i, double g, double E_0)
{
double m_e=0.511;
double k;
k=1-4/(3*(log(2*x/(m_e))+1/2));
return k;
}
double trapzd(double (*func)(double,double,int,double,double), double a, double b, double d, int n, int i, double E_0)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;
	if (n == 1) {
	return (s = 0.5*(b-a)*(FUNC(a,d,i,a,E_0)+FUNC(b,d,i,a,E_0)));
	}

	else
	{
		for (it = 1,j = 1;j<n-1;j++) it <<=  1;
		tnm = it;
		del = (b-a)/tnm;
		x = a+0.5*del;
		for (sum = 0.0,j = 1;j<= it;j++,x+= del) sum +=  FUNC(x,d,i,a,E_0);
		s = 0.5*(s+(b-a)*sum/tnm);
		return s;
	}

}

double qsimp(double (*func)(double,double,int,double,double), double a, double b, double d, int i, double E_0)
{
	double trapzd(double (*func)(double,double,int,double,double), double a, double b, double d, int n, int i, double E_0);
	int j;
	double s,st,ost = 0.0,os = 0.0;
	for (j = 1;j<= JMAX;j++)
	{
		//~ if(j%10==0) cout << j*100/JMAX << "%" << endl;
		st = trapzd(func,a,b,d,j,i,E_0);
		s = (4.0*st-ost)/3.0;
		if (j > 5)
		if (fabs(s-os) < EPS*fabs(os) ||
		(s == 0.0 && os == 0.0)) return s;
		os = s;
		ost = st;
	}
	return 0.0;
}
double trapzd2(double (*func)(double,double,int,double,double), double a, double b, double d, int n, int i, double E_0)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;
	if (n == 1) {
	return (s = 0.5*(b-a)*(FUNC(a,d,i,a,E_0)+FUNC(b,d,i,a,E_0)));
	}

	else
	{
		for (it = 1,j = 1;j<n-1;j++) it <<=  1;
		tnm = it;
		del = (b-a)/tnm;
		x = a+0.5*del;
		for (sum = 0.0,j = 1;j<= it;j++,x+= del) sum +=  FUNC(x,d,i,a,E_0);
		s = 0.5*(s+(b-a)*sum/tnm);
		return s;
	}

}

double qsimp2(double (*func)(double,double,int,double,double), double a, double b, double d, int i, double E_0)
{
	double trapzd2(double (*func)(double,double,int,double,double), double a, double b, double d, int n, int i, double E_0);
	int j;
	double s,st,ost = 0.0,os = 0.0;
	for (j = 1;j<= JMAX;j++)
	{
		//~ if(j%10==0) cout << j*100/JMAX << "%" << endl;
		st = trapzd2(func,a,b,d,j,i,E_0);
		s = (4.0*st-ost)/3.0;
		if (j > 5)
		if (fabs(s-os) < EPS*fabs(os) ||
		(s == 0.0 && os == 0.0)) return s;
		os = s;
		ost = st;
	}
	return 0.0;
}
double dsigma_compton(double x, double z, int i, double g, double E_0)
{
	double pi = 3.14159;
	double m_e = 0.511;
	double r_e = 1.42*pow(10.,-2);
	double dsigma = pi*pow(r_e,2)*m_e/pow(x,2)*(x/g+g/x+pow(m_e/g-m_e/x,2)-2*m_e*(1/g-1/x));
	//~ cout << "dsigma = " << dsigma << endl;
	return dsigma;
}
double dsigma_phph(double x, double z, int i, double g, double E_0)
{
	double dsigma = pow(x,2)*pow(1-g/x+pow(g/x,2),2);
	return dsigma;

}

double func_gain_cs(double x, double z, int i, double g, double E_0)
{
	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double Y = 0.24;
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double m_e = 0.511;
	double r_e = 1.42*pow(10.,-2);
	double a = 0.0073;
	double K_0 = E_0/(pow(E_x,2)*(2+log(E_c/E_x)));
	double pi = 3.15159;
		double f_gain_cs;
	//~ double f_gain_cs =K_0*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3)*(pow(E_x,2)/(60*x*x)*(((15*x*E_x*(x*x+2*x*m_e-2*m_e*m_e)+12*x*x*m_e*m_e+20*m_e*E_x*E_x*(m_e-2*x)+30*x*pow(E_x,3))/pow(E_x,5) - 15*x*E_c*(x*x+2*x*m_e-2*m_e*m_e)+12*x*x*m_e*m_e+20*m_e*E_c*E_c*(m_e-2*x)+30*x*pow(E_c,3))/pow(E_c,5)) + 2*pow(E_x,3/2)/(315*pow(x,9/2))*((-pow(x,5/2)*(45*x*E_x*(x*x+2*x*m_e-2*m_e*m_e)+35*x*x*m_e*m_e+63*m_e*E_x*E_x*(m_e-2*x)-105*x*pow(E_x,3)))/pow(E_x,9/2)+150*pow(x,2)-36*x*m_e+8*pow(m_e,2)));
	//~ double f_gain_cs = K_0*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3)*pi*pow(r_e,2)*m_e/(1260*pow(x,5)*pow(E_x,3)*pow(E_c,7))*(-30*pow(x,6)*(7*pow(E_x,7)*E_c+5*E_x*pow(E_c,7))-20*pow(x,5)*m_e*(3*pow(E_x,7)*(7*E_c+3*m_e)+15*E_x*pow(E_c,7)+5*pow(E_c,7)*m_e)+3*pow(x,4)*E_x*E_c*(20*pow(m_e,2)*(7*pow(E_x,6)+5*pow(E_c,6))+168*E_x*E_c*m_e*(pow(E_x,5)+pow(E_c,5))-35*pow(E_x,2)*pow(E_c,2)*(3*pow(E_x,4)+5*pow(E_c,4)))-252*pow(x,3)*pow(E_x,2)*pow(E_c,2)*pow(m_e,2)*(pow(E_x,5)+pow(E_c,5))+1200*pow(E_x,7)*pow(E_c,7)*pow(a/E_x,5/2)-288*pow(E_x,6)*pow(E_c,7)*m_e*pow(a/E_x,3/2)+64*pow(E_x,5)*pow(E_c,7)*pow(m_e,2)*pow(x/E_x,0.5));

	//~ f_gain_cs = K_0*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3)*pi*pow(r_e,2)*m_e/(1260*pow(x,5)*pow(E_x,3)*pow(E_c,7))*(-30*pow(x,6)*(7*pow(E_x,7)*E_c+5*E_x*pow(E_c,7))-20*pow(x,5)*m_e*(3*pow(E_x,7)*(7*E_c+3*m_e)+15*E_x*pow(E_c,7)+5*pow(E_c,7)*m_e)+3*pow(x,4)*E_x*E_c*(20*pow(m_e,2)*(7*pow(E_x,6)+5*pow(E_c,6))+168*E_x*E_c*m_e*(pow(E_x,5)+pow(E_c,5))-35*pow(E_x,2)*pow(E_c,2)*(3*pow(E_x,4)+5*pow(E_c,4)))-252*pow(x,3)*pow(E_x,2)*pow(E_c,2)*pow(m_e,2)*(pow(E_x,5)+pow(E_c,5))+1200*pow(E_x,7)*pow(E_c,7)*pow(x/E_x,5/2)-288*pow(E_x,6)*pow(E_c,7)*m_e*pow(x/E_x,3/2)+64*pow(E_x,5)*pow(E_c,7)*pow(m_e,2)*pow(x/E_x,0.5));

	//~ f_gain_cs = K_0*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3)*pi*pow(r_e,2)*m_e/(1260*pow(x,3))*(-45*pow(x,4)*(7*pow(E_x,4)+pow(E_c,4))/(pow(E_x,2)*pow(E_c,4))-2*pow(x,3)*m_e*(63*pow(E_x,5)*(5*E_c+2*m_e)+45*E_x*pow(E_c,5)+14*pow(E_c,5)*m_e)/(pow(E_x,3)*pow(E_c,5))+6*pow(x,2)*(-35*pow(E_x,2)*(3*pow(E_c,2)-4*E_c*m_e-3*pow(m_e,2))/pow(E_c,4)+15*m_e*m_e/pow(E_c,2)+28*m_e/E_x-35)+64*m_e*m_e*pow(E_x/x,3/2)-84*x*m_e*m_e*(5*pow(E_x,3)+pow(E_c,3))/(E_x*pow(E_c,3))-288*E_x*m_e*pow(E_x/x,0.5)+1200*x*E_x*pow(E_x/x,0.5));

	return f_gain_cs;
}

double func_gain_phph(double x, double z, int i, double g, double E_0)
{
	double pi = 3.14159;
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double m_e = 0.511;
	double r_e = 1.42*pow(10.,-2);
	double a = 0.0073;
	double K_0 = E_0/(pow(E_x,2)*(2+log(E_c/E_x)));
	double f_gain_phph ;
	//~ f_gain_phph = K_0*35584/(10125*pi)*pow(a*r_e,2)*pow(m_e,-6)*8*pow(pi,4)*pow(T_0*(1+z),6)/63*(pow(x,2)/3*(-pow(x,4)/pow(E_c,3)-3*pow(x,3)/pow(E_c,2)+3*pow(x,2)/E_c+(pow(x,4)+3*pow(x,3)*E_x-3*pow(x,2)*pow(E_x,2)-3*pow(E_x,4))/pow(E_x,3)+6*x*log(E_x/E_c)+3*E_c)+2/15*(23*pow(x,2)+(-3*pow(x,4)-10*pow(x,3)*E_x+15*pow(x,2)*pow(E_x,2)-30*x*pow(E_x,3)+5*pow(E_x,4))*(x/pow(E_x,3/2))/E_x));
	//~ f_gain_phph = 1/(30*E_x*pow(E_c,5))*(-6*pow(x,4)*(pow(E_x,5)+pow(E_c,5))+5*pow(x,3)*(3*pow(E_x,5)*E_c+5*E_x*pow(E_c,5))-30*pow(x,2)*pow(E_x,2)*pow(E_c,2)*(pow(E_x,3)+5*pow(E_c,3))+6*x*pow(E_x,3)*pow(E_c,3)*(pow(E_c,2)*(42*pow(x/E_x,0.5)-25)+5*pow(E_x,2))+10*pow(E_x,4)*pow(E_c,4)*(5*E_c-3*E_x));
	//~ f_gain_phph=K_0*35584/(10125*pi)*pow(a*r_e,2)*pow(m_e,-6)*8*pow(pi,4)*pow(T_0*(1+z),6)/63*f_gain_phph;}
	//~ f_gain_phph=-pow(x,4)*pow(E_x,2)/(3*pow(E_c,3))-pow(x,4)/(15*E_x)+pow(x,3)*pow(E_x,2)/pow(E_c,2)+42*pow(x,3)*pow(E_x/x,3/2)/5+pow(x,3)/3-3*pow(x,2)*pow(E_x,2)/E_c-3*pow(x,2)*E_x-2*x*pow(E_x,2)*log(E_c)-4*x*pow(E_x,2)+2*x*pow(E_x,2)*log(E_x)-pow(E_x,3)+pow(E_x,2)*E_c;
	//~ f_gain_phph=K_0*35584/(10125*pi)*pow(a*r_e,2)*pow(m_e,-6)*8*pow(pi,4)*pow(T_0*(1+z),6)/63*f_gain_phph;
	f_gain_phph=0;
	if(f_gain_phph<0) {f_gain_phph=0;}
	return f_gain_phph;
}

double gamma_compton(double x, double z, int i, double g, double E_0)
{
	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double Y = 0.24;
	double m_e = 0.511;
	double X = 2*x/m_e;
	double pi = 3.14159;
	double a = 0.0073;
	double r_e = 1.42*pow(10.,-2);
	double sigma_cs = 2*pi*pow(r_e,2)/X*((1-4/X-8/pow(X,2))*log(1+X)+1/2+8/X+1/(2*pow(1+X,2)));
	//~ double sigma_cs = 2*pi*pow(r_e,2)/X*log(X);
	double Gamma = sigma_cs*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3);
	return Gamma;
}
double gamma_NPC(double x, double z, int i, double g, double E_0)
{
	double n_y_0 = 3.154*pow(10.,-30);
	double eta = 6.05*pow(10.,-10);
	double Y = 0.24;
	double m_e = 0.511;
	double k = x/m_e;
	 //~ cout << "x = " << x << endl;
	double r_e = 1.42*pow(10.,-2);
	double rho = (2*k-4)/(k+2+2*pow(2*k,0.5));
	double sigma_PCN;
	double a = 0.0073;
	double pi = 3.14159;
	//~ if(k<4) sigma_PCN = a*pow(r_e,2)*2*pi/3*pow((k-2)/k,3)*(1+rho/2+23*pow(rho,2)/40+11*pow(rho,3)/60+20*pow(rho,4)/960);

	//~ else sigma_PCN = a*pow(r_e,2)*(28/9*log(2*k)-218/27+pow(2/k,2)*(2/3*pow(log(2*k),3)-pow(log(2*k),2)+(6-pi*pi/3)*log(2*k)+2*1.2021+pi*pi/6-7/2)-pow(2/k,4)*(3/16*log(2*k)+1/8)-pow(2/k,6)*(29/2304*log(2*k)-77/13824)) ;
	sigma_PCN = a*pow(r_e,2)*(28/9*log(2*k)-218/27) ;
	double Gamma_2 = sigma_PCN*pow(1+z,3)*eta*n_y_0;
	//~ cout << "sigma_PCN = " << sigma_PCN << endl;
 //~ if(Gamma_2!=0) cout << "Gamma NPC = "<< Gamma_2<<endl;
	return Gamma_2;
}

double gamma_phph(double x, double z, int i, double g, double E_0)
{
	double m_e = 0.511;
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double a = 0.0073;
double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	//~ double Gamma = 0.1513*pow(a,4)*m_e*pow(x/m_e,3)*pow(T_0/m_e,6)*pow(1+z,6);
	double Gamma_3;
	Gamma_3 = 0.1513*pow(a,4)*m_e*pow(x/m_e,3)*pow(T_0*(1+z)/m_e,6);

	return Gamma_3;
}
double correction_phph(double x, double z, int k, double g, double E_0)
{
	double T_0 = 2.7255*0.862*pow(10.,-10),a1,b1,a2,b2,A,B;
	double T = T_0*(1+z);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double y=0;

		if(T*pow(10,6)<=100) a1 = 0.00963738 ,  b1 = 1.22003;

	else if(T*pow(10,6)>100 && T*pow(10,6) <=1000)a1 =  0.0257008 , b1 = 1.1504;

	else if(T*pow(10,6)>1000 && T*pow(10,6) <=2000) a1 = 0.0114048, b1 = 1.23979;
	else if(T*pow(10,6)>2000 && T*pow(10,6) <3000) a1 = 0.0129126, b1 = 1.23948;
	else if(T*pow(10,6)>=3000) a1 = 0.0152741, b1 = 1.22058;

	A = pow(T*pow(10,6),a1)*exp(b1);
	B = log10(gamma_phph(E_c,z,k,g,E_0)*pow(10,-3))-A*log10(E_c*pow(10,-3));

	//~ cout << " T = " <<T*pow(10,6) << endl;
	//~ if(T*pow(10,6)>=1000 && T*pow(10,6) <=1100)A = 3.78, B = -17.53;
	//~ cout << "A = " << A << "  " << " B = " << B << endl;

	y=pow(x*pow(10,-3),A)*pow(10,B)/(gamma_phph(x,z,k,g,E_0)*pow(10,-3));

	return y;
}
double func_spectre(double x, double z, int i, double g, double E_0)
{
//  cout << "z spectre" << z << endl;
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double K_0 = E_0/(pow(E_x,2)*(2+log(E_c/E_x)));
	double f = 0;

	if(x < E_x) f = K_0*pow(E_x/x,1.5);
	else if(x > E_x && x < E_c) f =  K_0*pow(E_x/x,2);
	else {f = 0;}
	return f;
}
double func_spectre_special(double x, double z, int i, double g, double E_0)
{
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double f = 0 ;
	double m_e = 0.511;
	double E_gamma_b = 2.7*T_0*(1+z);
	double gamma_e = 4*E_gamma_b*E_0/(m_e*m_e);

	double q = x/(gamma_e*(E_0-x));
	if(q>1){
	cout << " gamma_e = " << gamma_e << endl;
	cout << " q  = " << q << " z = " << z << endl;}
	if(q<=1 && q>0)f = x*(2*q*log(q)+(1+2*q)*(1-q)+(gamma_e*q*gamma_e*q)/(2*(1+gamma_e*q))*(1-q));
	return f;
}
double func_spectre2(double x, double z, int i, double g, double E_0)
{

	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double K_0 = E_0/(pow(E_x,2)*(2+log(E_c/E_x)));
	double K_1 = pow(E_0,1/2)/(2*pow(E_x,3/2));
	double K_2 = E_0/(pow(E_x,2)*(2+log(E_0/E_x)));
	double f;
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
double func_spectre_gauss(double x, double z, int i, double g, double E_0)
{
double T_0 = 2.7255*0.862*pow(10.,-10);
double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
double pi = 3.14159;
double sigma = E_0/100;
double E_mean = E_0;
double f;
//~ cout << E_mean << endl;
//~ double K = E_0/(exp(-E_mean*E_mean/(2*pow(sigma,2)))/pow(2*pi,0.5) + E_mean/2*(1+erf(E_mean/(pow(2,0.5)*sigma))));
//~ double K = ((-exp(-(pow(E_mean - E_c,2)/(2*pow(sigma,2)))) + exp(-(pow(E_mean,2)/(2*pow(sigma,2)))))*sigma)/pow(2*pi,0.5) + 1/2*E_mean*(erf(E_mean/(pow(2,0.5)*sigma)) - erf((-E_c + E_mean)/(pow(2,0.5)*sigma)));

//~ double K = 2*E_0/(exp(-pow(E_0,2)/(2*pow(sigma,2)))*pow(2/pi,0.5)*sigma+E_0*(1+erf(E_0/(pow(2,0.5)*sigma))));
double K = E_0/(exp(-pow(E_0,2)/(2*pow(sigma,2)))*sigma/pow(2*pi,0.5)+E_0*(1/2+erf(E_0/(pow(2,0.5)*sigma))));
//~ cout << " K " << K << endl;
f = K*exp(-pow(x-E_mean,2)/(2*pow(sigma,2)))/(sigma*sqrt(2*pi));
//~ cout << " K = " << K << endl;
return f;

}
double func_spectre_new(double x, double z, int i, double g, double E_0)
{
double T_0 = 2.7255*0.862*pow(10.,-10);
double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
double pi = 3.14159;
double sigma = E_0/1000;
double E_mean = E_0;
double K = 2*E_0/(exp(-pow(E_0,2)/(2*pow(sigma,2)))*pow(2/pi,0.5)*sigma+E_0*(1+erf(E_0/(pow(2,0.5)*sigma))));
//~ cout << " K " << K << endl;
return K*exp(-pow(x-E_mean,2)/(2*pow(sigma,2)))/(sigma*sqrt(2*pi));

}


double func_int_phph(double x, double z, int k, double g, double E_0)
{

double integrand ;
integrand = func_spectre_gauss(x,z,k,g,E_0)*dsigma_phph(x,z,k,g,E_0)/(gamma_NPC(x,z,k,g,E_0)+spec_protheroe(x,z,k,g,E_0)*gamma_compton(x,z,k,g,E_0)+gamma_phph(x,z,k,g,E_0));

//~ else integrand = func_spectre(x,z,k,g,E_0)*dsigma_phph(x,z,k,g,E_0)/(gamma_tot_compton(x,z,k-1,g,E_0)+gamma_NPC(x,z,k-1,g,E_0)+gamma_tot_phph(x,z,k-1,g,E_0));
return integrand;

}double func_int_phph2(double x, double z, int k, double g, double E_0)
{

double integrand ;
integrand = func_spectre_gauss(x,z,k,g,E_0)*dsigma_phph(x,z,k,g,E_0)/(gamma_NPC(x,z,k,g,E_0)+spec_protheroe(x,z,k,g,E_0)*gamma_compton(x,z,k,g,E_0)+gamma_tot_phph(x,z,k,g,E_0));

//~ else integrand = func_spectre(x,z,k,g,E_0)*dsigma_phph(x,z,k,g,E_0)/(gamma_tot_compton(x,z,k-1,g,E_0)+gamma_NPC(x,z,k-1,g,E_0)+gamma_tot_phph(x,z,k-1,g,E_0));
return integrand;

}double func_int_phph3(double x, double z, int k, double g, double E_0)
{

double integrand ;
integrand = func_spectre(x,z,k,g,E_0)*dsigma_phph(x,z,k,g,E_0)/(gamma_NPC(x,z,k,g,E_0)+spec_protheroe(x,z,k,g,E_0)*gamma_compton(x,z,k,g,E_0)+gamma_tot_phph_2(x,z,k,g,E_0));

//~ else integrand = func_spectre(x,z,k,g,E_0)*dsigma_phph(x,z,k,g,E_0)/(gamma_tot_compton(x,z,k-1,g,E_0)+gamma_NPC(x,z,k-1,g,E_0)+gamma_tot_phph(x,z,k-1,g,E_0));
return integrand;

}
double func_int_compton(double x, double z, int k, double g, double E_0)
{

double integrand ;
double Gamma = (gamma_NPC(x,z,k,g,E_0)+gamma_compton(x,z,k,g,E_0)+gamma_phph(x,z,k,g,E_0)*correction_phph(x,z,k,g,E_0));

integrand = func_spectre(x,z,k,g,E_0)*dsigma_compton(x,z,k,g,E_0)/Gamma;

//~ else integrand = func_spectre(x,z,k,g,E_0)*dsigma_compton(x,z,k,g,E_0)/(gamma_tot_compton(x,z,k-1,g,E_0)+gamma_NPC(x,z,k-1,g,E_0)+gamma_tot_phph(x,z,k-1,g,E_0));
return integrand;

}double func_int_compton_2(double x, double z, int k, double g, double E_0)
{

double integrand ;
integrand = func_spectre_gauss(x,z,k,g,E_0)*dsigma_compton(x,z,k,g,E_0)/(gamma_NPC(x,z,k,g,E_0)+gamma_tot_compton(x,z,k,g,E_0)+gamma_phph(x,z,k,g,E_0)*correction_phph(x,z,k,g,E_0));

//~ else integrand = func_spectre(x,z,k,g,E_0)*dsigma_compton(x,z,k,g,E_0)/(gamma_tot_compton(x,z,k-1,g,E_0)+gamma_NPC(x,z,k-1,g,E_0)+gamma_tot_phph(x,z,k-1,g,E_0));
return integrand;

}double func_int_compton_3(double x, double z, int k, double g, double E_0)
{

double integrand ;
integrand = func_spectre(x,z,k,g,E_0)*dsigma_compton(x,z,k,g,E_0)/(gamma_NPC(x,z,k,g,E_0)+gamma_tot_compton_2(x,z,k,g,E_0)+gamma_phph(x,z,k,g,E_0)*correction_phph(x,z,k,g,E_0));

//~ else integrand = func_spectre(x,z,k,g,E_0)*dsigma_compton(x,z,k,g,E_0)/(gamma_tot_compton(x,z,k-1,g,E_0)+gamma_NPC(x,z,k-1,g,E_0)+gamma_tot_phph(x,z,k-1,g,E_0));
return integrand;

}
double gamma_tot_phph(double x, double z, int k, double g, double E_0)
{	double a = 0.0073;
	double pi = 3.14159;
	double r_e = 1.42*pow(10.,-2);
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));

	double m_e = 0.511;
	double Gamma;

	//~ double GammaTot = gamma_phph(x,z,k,g,E_0)/(1+pow(T_0*(1+z),6)*8*pow(pi,4)*35584*pow(a*r_e,2)*pow(m_e,-6)*qsimp(func_int_phph,x,E_c,z,k)/(63*10125*pi*func_spectre(x,z,k,g,E_0)));
	Gamma = qsimp(func_int_phph,x,E_c,z,k,E_0);

	double gain = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(a*r_e,2)*pow(m_e,-6)*Gamma*(gamma_NPC(x,z,k,g,E_0)+gamma_compton(x,z,k,g,E_0)*spec_protheroe(x,z,k,g,E_0)+gamma_phph(x,z,k,g,E_0))/(63*10125*pi*func_spectre_gauss(x,z,k,g,E_0));
	//~ double gain = 0;
		//~ cout << " gain num = " << gain << endl;
	double GammaTot = gamma_phph(x,z,k,g,E_0)-gain;
	//~ cout << "GammaTot = " << GammaTot << " Gamma phph = " <<gamma_phph(x,z,k,g,E_0) << endl;
	return GammaTot;
}
double gamma_tot_phph_2(double x, double z, int k, double g, double E_0)
{	double a = 0.0073;
	double pi = 3.14159;
	double r_e = 1.42*pow(10.,-2);
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));

	double m_e = 0.511;
	double Gamma;

	//~ double GammaTot = gamma_phph(x,z,k,g,E_0)/(1+pow(T_0*(1+z),6)*8*pow(pi,4)*35584*pow(a*r_e,2)*pow(m_e,-6)*qsimp(func_int_phph,x,E_c,z,k)/(63*10125*pi*func_spectre(x,z,k,g,E_0)));
	Gamma = qsimp(func_int_phph2,x,E_c,z,k,E_0);

	double gain = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(a*r_e,2)*pow(m_e,-6)*Gamma*(gamma_NPC(x,z,k,g,E_0)+gamma_compton(x,z,k,g,E_0)*spec_protheroe(x,z,k,g,E_0)+gamma_tot_phph(x,z,k,g,E_0))/(63*10125*pi*func_spectre_gauss(x,z,k,g,E_0));
	//~ double gain = 0;

	double GammaTot = gamma_phph(x,z,k,g,E_0)-gain;
	//~ cout << "GammaTot = " << GammaTot << " Gamma phph = " <<gamma_phph(x,z,k,g,E_0) << endl;
	return GammaTot;
}
double gamma_tot_phph_3(double x, double z, int k, double g, double E_0)
{	double a = 0.0073;
	double pi = 3.14159;
	double r_e = 1.42*pow(10.,-2);
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));

	double m_e = 0.511;
	double Gamma;

	//~ double GammaTot = gamma_phph(x,z,k,g,E_0)/(1+pow(T_0*(1+z),6)*8*pow(pi,4)*35584*pow(a*r_e,2)*pow(m_e,-6)*qsimp(func_int_phph,x,E_c,z,k)/(63*10125*pi*func_spectre(x,z,k,g,E_0)));
	Gamma = qsimp(func_int_phph3,x,E_c,z,k,E_0);

	double gain = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(a*r_e,2)*pow(m_e,-6)*Gamma*(gamma_NPC(x,z,k,g,E_0)+gamma_compton(x,z,k,g,E_0)*spec_protheroe(x,z,k,g,E_0)+gamma_tot_phph_2(x,z,k,g,E_0))/(63*10125*pi*func_spectre(x,z,k,g,E_0));
	//~ double gain = 0;
		//~ cout << " gain num = " << gain << endl;
	double GammaTot = gamma_phph(x,z,k,g,E_0)-gain;
	//~ cout << "GammaTot = " << GammaTot << " Gamma phph = " <<gamma_phph(x,z,k,g,E_0) << endl;
	return GammaTot;
}

double gamma_tot_compton(double x, double z, int k, double g, double E_0)
{
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));

	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double Y = 0.24;
	double Gamma;

	//~ double GammaTot = gamma_compton(x,z,k,g,E_0)/(1+(eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3)*qsimp(func_int_compton,x,E_c,z,k))/func_spectre(x,z,k,g,E_0));
	Gamma=(gamma_NPC(x,z,k,g,E_0)+gamma_compton(x,z,k,g,E_0)+gamma_phph(x,z,k,g,E_0)*correction_phph(x,z,k,g,E_0));
	//~ cout << " Gamma = " << Gamma << endl;
	double gain = eta*n_y_0*(1+Y/2)*pow(1+z,3)*qsimp(func_int_compton,x,E_c,z,k,E_0)*Gamma/(func_spectre(x,z,k,g,E_0)*(1+Y));
	//~ cout << "gain = " << gain << endl;
	double GammaTot = gamma_compton(x,z,k,g,E_0)-gain;
	//~ cout << "Gamma tot = " << GammaTot<< " Gamma compton = " << gamma_compton(x,z,k,g,E_0)<<endl;
	return GammaTot;

}
double gamma_tot_compton_2(double x, double z, int k, double g, double E_0)
{
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));

	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double Y = 0.24;
	double Gamma;

	//~ double GammaTot = gamma_compton(x,z,k,g,E_0)/(1+(eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3)*qsimp(func_int_compton,x,E_c,z,k))/func_spectre(x,z,k,g,E_0));

	double gain = eta*n_y_0*(1+Y/2)*pow(1+z,3)*qsimp(func_int_compton_2,x,E_c,z,k,E_0)*(gamma_NPC(x,z,k,g,E_0)+gamma_tot_compton(x,z,k,g,E_0)+gamma_phph(x,z,k,g,E_0)*correction_phph(x,z,k,g,E_0))/(func_spectre(x,z,k,g,E_0)*(1+Y));
	double GammaTot = gamma_compton(x,z,k,g,E_0)-gain;
	//~ cout << "Gamma tot = " << GammaTot<< " Gamma compton = " << gamma_compton(x,z,k,g,E_0)<<endl;
	return GammaTot;

}double gamma_tot_compton_3(double x, double z, int k, double g, double E_0)
{
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));

	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double Y = 0.24;
	double Gamma;

	//~ double GammaTot = gamma_compton(x,z,k,g,E_0)/(1+(eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3)*qsimp(func_int_compton,x,E_c,z,k))/func_spectre(x,z,k,g,E_0));

	double gain = eta*n_y_0*(1+Y/2)*pow(1+z,3)*qsimp(func_int_compton_3,x,E_c,z,k,E_0)*(gamma_NPC(x,z,k,g,E_0)+gamma_tot_compton_2(x,z,k,g,E_0)+gamma_phph(x,z,k,g,E_0)*correction_phph(x,z,k,g,E_0))/(func_spectre(x,z,k,g,E_0)*(1+Y));
	double GammaTot = gamma_compton(x,z,k,g,E_0)-gain;
	//~ cout << "Gamma tot = " << GammaTot<< " Gamma compton = " << gamma_compton(x,z,k,g,E_0)<<endl;
	return GammaTot;

}
double func_kawmor(double x, double z, int h, double E_0)
{	//Attention résultat en GEV
	double a_pp[3],a_low[3],N_pp[3],N_low[3];

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

	double f;
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));

	double T = T_0*(1+z);


	if(x < E_x) f = E_0*N_low[h]*pow(T*pow(10.,-3),-3)*pow(x*pow(10.,-3),a_low[h])*pow(10.,-3);
	else if(x > E_x && x < E_c) f = E_0*N_pp[h]*pow(T*pow(10.,-3),-6)*pow(x*pow(10.,-3),a_pp[h])*pow(10.,-3);
	else f = 0;

	return f;
}



double func_sigma(double x, double z, int i, double g, double E_0)
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

	double y;

	if(i == 0)
	{
	double Q = 2.467032;
	//y = (0.105*2371*pow(x,-2)*exp(-2.5954*pow(x-Q,-0.5))*exp(-2.056*(x-Q))*(1+2.2875*pow(x-Q,2)-1.1798*pow(x-Q,3)+2.5279*pow(x-Q,4)))*pow(x,d);
	y = 0.057*931.434*pow(x,-2)*exp(-2.59*pow(x-Q,-0.5));
	}
	else if(i == 1)
	{
	double Q = 7.249962;
	y = (0.176*pow(Q,1.51)*pow(x-Q,0.49)*pow(x,-2)+1205*pow(Q,5.5)*pow(x-Q,5)*pow(x,-10.5)+0.06/(1+pow((x-7.46)/0.188,2)));
	}
	else if(i == 2)
	{
	double Q = 10.948850;
	y = 122*pow(Q,4)*pow(x-Q,3)*pow(x,-7);
	}
	else if(i == 3)
	{
	double Q0 = 8.725;
	double Q1 = 23;
	y = 3.8 * pow(Q0,2.3)*(x-Q0)/pow(x,3.3);
	if(x>=Q1) y += 2.1*pow(Q1,1.5)*(x-Q1)/pow(x,2.5);

	}
	else if(i == 4)
	{
	double Q = 9.98;
	y = 10.8 *pow(Q,2)*pow(x-Q,1.2)*(x-Q)/pow(x,3.2);
	}
	else if(i == 5)
	{
	double Q = 22.28;
	y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
	}
	else if(i == 6)
	{
	double Q = 23.05;
	y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
	}
	else if(i == 7)
	{
		double Q = 1.586627;
		// if(x>1.59)y =  (0.504*2371*pow(x,-2)*exp(-5.1909*pow(x-Q,-0.5))*exp(-0.548*(x-Q))*(1-0.428*pow(x-Q,2)+0.543*pow(x-Q,3)-0.115*pow(x-Q,4)));
		//if(x>1.59)y =  (409*pow(x,-2)*exp(-5.1909*pow(x-Q,-0.5))*exp(-0.548*(x-Q))*(1-0.428*pow(x-Q,2)+0.543*pow(x-Q,3)-0.115*pow(x-Q,4)));
		if(x>1.59)y = 0.26*931.434*pow(x,-2)*exp(-5.19*pow(x-Q,-0.5));
		else y=0;
	}
	else if(i == 8)
	{
		double Q = 5.605794;
		y = (32.6*pow(Q,10)*pow((x-Q),2)*pow(x,-12)+2.27*pow(10.,6)*pow(Q,8.8335)*pow((x-Q),13)*pow(x,-21.8335));
	}
	else if(i == 9)
	{
		double Q = 9.30468;
		y = 133*pow(Q,4)*pow(x-Q,3)*pow(x,-7);
	}
	else if(i == 10)
	{
	double Q0 = 7.08;
	double Q1 = 23;
	y = 3.8 * pow(Q0,2.3)*(x-Q0)/pow(x,3.3);
	if(x>=Q1) y += 2.1*pow(Q1,1.5)*(x-Q1)/pow(x,2.5);
	}
	else if(i == 11)
	{
	double Q = 10.68;
	y = 10.8 *pow(Q,2)*pow(x-Q,1.2)*(x-Q)/pow(x,3.2);
	}
	else if(i == 12)
	{
	double Q = 22.17;
	y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
	}
	else if(i == 13)
	{
	double Q = 21.4;
	y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
	}
	else if(i == 14)
	{
	double Q = 2.224573;
	y = 18.75*(pow(pow(Q*(x-Q),0.5)/x,3)+0.007947*pow(pow(Q*(x-Q),0.5)/x,2)*pow(pow(Q,0.5)-pow(0.037,0.5),2)/(x-Q+0.037));
	}
	else if(i == 15)
	{
	double Q = 19.813852;
	//~ y = 128.9*pow(Q,4.524)*pow(x-Q,2.512)/pow(x,4.524+2.512);
	y = 19.5*pow(Q,3.5)*pow(x-Q,1)/pow(x,4.5);
	}
	else if(i == 16)
	{
	double Q = 20.577615;
	//~ y = 31.68*pow(Q,3.663)*pow(x-Q,1.580)/pow(x,3.663+1.580);
	y = 17.1*pow(Q,3.5)*pow(x-Q,1.)/pow(x,4.5);
	}
	else if(i == 17)
	{
	double Q = 23.846527;
	y = 10.7*pow(Q,10.2)*pow(x-Q,3.4)/pow(x,13.6);
	}
	else if(i == 18)
	{
	double Q = 26.0711;
	y = 21.7*pow(Q,4.0)*pow(x-Q,3.0)/pow(x,7.0);
	}
	else if(i == 19)
	{
	double Q = 5.483485;
	y = 8.88*pow(Q,1.75)*pow(x-Q,1.65)/pow(x,3.4);
	}
	else if(i == 20)
	{
	double Q = 7.718058;
	y = 16.7*pow(Q,1.95)*pow(x-Q,2.3)/pow(x,4.25);
	}
	y=y*2.569*pow(10.,-6);//Conversion mb en Mev^-2
	return y;

}

double func_z(double z,double T, int i, double g, double E_0)
{
	double H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double tau =  pow(T,-2);
	double A;
	A = exp(-1/(2*H_r*tau*(z+1)*(z+1)));
	return A;
}

double func_gauss(double x, double z, int i, double g, double E_0)
{
	double f;
	//~ /**Spectre GAUSSIEN dégradé AVEC réinjection**/
	//~ f=func_spectre_gauss(x,z,i,g,E_0)*func_sigma(x,z,i,g,E_0)/(gamma_NPC(x,z,i,g,E_0)+spec_protheroe(x,z,i,g,E_0)*gamma_compton(x,z,i,g,E_0)+gamma_phph(x,z,i,g,E_0)*correction_phph(x,z,i,g,E_0));
	/**Spectre GAUSSIEN dégradé SANS réinjection**/
	f=func_spectre_gauss(x,z,i,g,E_0)*func_sigma(x,z,i,g,E_0)/(gamma_NPC(x,z,i,g,E_0)+gamma_compton(x,z,i,g,E_0)+gamma_phph(x,z,i,g,E_0));
	//~ /**Spectre STANDARD dégradé AVEC réinjection***/
	//~ f=func_spectre(x,z,i,g,E_0)*func_sigma(x,z,i,g,E_0)/(gamma_NPC(x,z,i,g,E_0)+spec_protheroe(x,z,i,g,E_0)*gamma_compton(x,z,i,g,E_0)+gamma_phph(x,z,i,g,E_0)*correction_phph(x,z,i,g,E_0));

	return f;
}
double func_standard(double x, double z, int i, double g, double E_0)
{
	double f;
	/**Spectre STANDARD dégradé SANS réinjection***/
	f=func_spectre(x,z,i,g,E_0)*func_sigma(x,z,i,g,E_0)/(gamma_NPC(x,z,i,g,E_0)+gamma_compton(x,z,i,g,E_0)+gamma_phph(x,z,i,g,E_0));

	 /**Spectre STANDARD dégradé AVEC réinjection***/
	// f=func_spectre(x,z,i,g,E_0)*func_sigma(x,z,i,g,E_0)/(gamma_NPC(x,z,i,g,E_0)+spec_protheroe(x,z,i,g,E_0)*gamma_compton(x,z,i,g,E_0)+gamma_phph(x,z,i,g,E_0)*correction_phph(x,z,i,g,E_0));
//	cout << f << endl;
	return f;
}
double func_spectre_degraded(double x, double z, int i, double g, double E_0)
{
	double f;
	/**Spectre STANDARD dégradé SANS réinjection***/
	f=func_spectre(x,z,i,g,E_0)/(gamma_NPC(x,z,i,g,E_0)+gamma_compton(x,z,i,g,E_0)+gamma_phph(x,z,i,g,E_0));

	/**Spectre STANDARD dégradé AVEC réinjection***/
	// f=func_spectre(x,z,i,g,E_0)*func_sigma(x,z,i,g,E_0)/(gamma_NPC(x,z,i,g,E_0)+spec_protheroe(x,z,i,g,E_0)*gamma_compton(x,z,i,g,E_0)+gamma_phph(x,z,i,g,E_0)*correction_phph(x,z,i,g,E_0));
//	cout << f << endl;
	return f;
}
double Spectre_une_interaction_monochromatique(double  x, double  z, int i, double  g, double  E_0)
{

	double Gamma_tot = gamma_NPC(x,z,i,g,E_0)+gamma_compton(x,z,i,g,E_0)+gamma_phph(x,z,i,g,E_0);
	// cout << "Gamma_tot = " << Gamma_tot<< endl;
	double proba_ph = gamma_phph(x,z,i,g,E_0)/Gamma_tot;
	double proba_compton = gamma_compton(x,z,i,g,E_0)/Gamma_tot;
	// cout << " proba ph = " << proba_ph << " proba compton = " << proba_compton << endl;
	double E = x/0.511;
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double m_e=0.511;
	double alpha = 1./137;
	double r_e = alpha/m_e;
	double E_gamma = E_0/m_e;
	double n_y_0 = 3.154*pow(10.,-30) ;
	double eta = 6.05*pow(10.,-10) ;
	double Y = 0.25;
	double pi = 3.14157;
	double n_e = eta*n_y_0*(1+Y/2);
	double int_BB = 8./63*pow(pi,4)*pow(T_0*(1+z),6);
	// double R = 1.83*pow(10,-27)*50*pow(1+z,6)*pow(E_gamma,3);
	// double spectre_gamma_gamma = R*(20/7)/E_gamma*pow(1-E/E_gamma+pow(E/E_gamma,2),2);
	double spectre_gamma_gamma = 35584./(10125*pi)*pow(alpha*r_e,2)*pow(m_e,-6)*pow(E_0,2)*pow(1-x/E_0+pow(x/E_0,2),2)*int_BB;
	// double spectre_gamma_gamma = 0;
	double spectre_compton = pi*r_e*r_e*m_e*pow(E_0,-2)*(x/E_0+E_0/x+pow(m_e/x-m_e/E_0,2)-2*m_e*(1/x-1/E_0))*n_e*pow(1+z,3);
	double f = (proba_ph*spectre_gamma_gamma/gamma_phph(x,z,i,g,E_0) + proba_compton*spectre_compton/gamma_compton(x,z,i,g,E_0));
	// double S=pow(T,1.83987)*pow(10,9.17361);
	// cout << " f " << f << endl;
	return x*f;
}
double  func_une_interaction_monochromatique(double  x, double  z, int i, double  g, double  E_0)
{
	double  f;
	f=Spectre_une_interaction_monochromatique(x,z,i,g,E_0)*func_sigma(x,z,i,g,E_0)/(gamma_NPC(x,z,i,g,E_0)+gamma_compton(x,z,i,g,E_0)+gamma_phph(x,z,i,g,E_0));

	return f;
}
double func_special(double x, double z, int i, double g, double E_0)
{
	double f;
	/**Spectre SPECIAL dégradé SANS réinjection***/
	f=func_spectre_special(x,z,i,g,E_0)*func_sigma(x,z,i,g,E_0)/(gamma_NPC(x,z,i,g,E_0)+gamma_compton(x,z,i,g,E_0)+gamma_phph(x,z,i,g,E_0));

	/**Spectre SPECIAL dégradé AVEC réinjection***/
	// f=func_spectre_special(x,z,i,g,E_0)*func_sigma(x,z,i,g,E_0)/(gamma_NPC(x,z,i,g,E_0)+spec_protheroe(x,z,i,g,E_0)*gamma_compton(x,z,i,g,E_0)+gamma_phph(x,z,i,g,E_0)*correction_phph(x,z,i,g,E_0));
	//	cout << f << endl;
	return f;
}
double func_theo(double z, double T, int i, double g, double E_0)
{
double a = (1+z)/(1+g);
//~ cout << "E = " << E_0*a;
double y = func_z(z,T,i,g,E_0)*func_sigma(E_0*a, z, i, E_0,E_0)/(gamma_NPC(E_0*a,z,i,E_0,E_0)+gamma_compton(E_0*a,z,i,E_0,E_0)+gamma_phph(E_0*a,z,i,E_0,E_0));
//~ cout << " y = " << y << endl;
return y;
}
double func_z_v2(double z, double z_0,int i, double g, double E_0)
{

	double Emin[19];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double tau =  5*pow((1+z_0),-2)/(2*H_r);
	double A,s=0,s_int;
	double m_e = 0.511;
	double E_gamma_b = 2.7*T_0*(1+z);
	double E_max = 4*E_gamma_b*E_0*E_0/(m_e*m_e*(1+4*E_gamma_b*E_0/(m_e*m_e)));
	//  if(E_0 < 3 && E_c>3)E_c == 3;
	//  if(E_max>E_c)E_max=E_c;
	// cout << " E_max = " << E_max << endl;
	//  if(i>=14 && E_max>Emin[i])
	//  {
	// 	 s+=qsimp2(func_special,Emin[i],E_max,z,i,E_0);
	//  }
	 if(E_c>Emin[i])
	 {
		 s+=qsimp2(func_standard,Emin[i],E_c,z,i,E_0);
	 }

	A = exp(-1/(2*H_r*tau*(z+1)*(z+1)))*s;
	//~ cout << " A = " << A << endl;
	return A;
}

double func_Be_non_standard(double z, double z_0,int i, double g, double E_0)
{
		double T_0 = 2.7255*0.862*pow(10.,-10);
		double Emin[19];
		Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
		double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
		double H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
		double t_inj =  pow((1+z_0),-2)/(2*H_r);
		double tau = 5*t_inj;
		double f=0, s =0;
		if(v2==1){
		if(E_c<=E_0 && E_c>Emin[i]){
			s = qsimp2(func_standard,Emin[i],E_c,z,i,E_0);
			// cout <<"spectre standard E_c = " << E_c << " z = " << z << endl;
		}
		if(E_c>E_0 && E_0>Emin[i]){
			s=func_sigma(E_0,z,i,g,E_0)/(gamma_NPC(E_0,z,i,g,E_0)+gamma_compton(E_0,z,i,g,E_0)+gamma_phph(E_0,z,i,g,E_0));
			// cout <<"spectre non standard E_c = " << E_c << " z = " << z << endl;
		}}
		else{
		if(E_0<E_c)s=func_sigma(E_0,z,i,g,E_0)/(gamma_NPC(E_0,z,i,g,E_0)+gamma_compton(E_0,z,i,g,E_0)+gamma_phph(E_0,z,i,g,E_0));}


			f = exp(-1/(2*H_r*tau*(z+1)*(z+1)))*s;

		return f;
}
double func_Anal_avec_reinj(double z, double z_0,int i, double g, double E_0)
{
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double Emin[19];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double t_inj =  pow((1+z_0),-2)/(2*H_r);
	double tau = 5*t_inj;
	double f=0, s =0;
	if(v2==1){
		if(E_c<=E_0 && E_c>Emin[i]){
			s = qsimp2(func_standard,Emin[i],E_c,z,i,E_0);
			// cout <<"spectre standard E_c = " << E_c << " z = " << z << endl;
		}
		if(E_c>E_0 && E_0>Emin[i]){
			s=func_sigma(E_0,z,i,g,E_0)/(gamma_NPC(E_0,z,i,g,E_0)+gamma_compton(E_0,z,i,g,E_0)+correction_phph(E_0,z,i,g,E_0));
			// cout <<"spectre non standard E_c = " << E_c << " z = " << z << endl;
		}}
		else{
			if(E_0<E_c)s=func_sigma(E_0,z,i,g,E_0)/(gamma_NPC(E_0,z,i,g,E_0)+spec_protheroe(E_0,z,i,g,E_0)*gamma_compton(E_0,z,i,g,E_0)+correction_phph(E_0,z,i,g,E_0)*gamma_phph(E_0,z,i,g,E_0));}


				f = exp(-1/(2*H_r*tau*(z+1)*(z+1)))*s;

				return f;
}
/*double func_z_v3(double z, double z_0,int i, double g, double E_0)
{

	double Emin[19];
	Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
	double T_0 = 2.7255*0.862*pow(10.,-10);
	double E_x = 0.261121/(80*T_0*(1+z)), E_c = 0.261121/(22*T_0*(1+z));
	double H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5);
	double tau =  pow(T_0*(1+z_0),-2);
	double A,s,s_He_perte=0,s_De_gain=0,s_De_perte=0;

	s_De_perte=exp(-E_0*H_r*tau/(Zeta*n_y_0)*qsimp(func_z_v2,z_max,z,z,14,E_0));
	for(int j=15;j<19;j++
	{
		if(E_c>Emin[j])s_He_perte+=qsimp(func_z_v2,z_0,z,z,j,E_0);
		if(j==17 || j==18)
		{
			s_De_gain+=qsimp(func,Emin[j],E_c,z,j,E_0);;
		}
	}
	s_He_perte=exp(-E_0*H_r*tau/(Zeta*n_y_0)*s_He_perte);
	A=-E_0*H_r*tau/(Zeta*n_y_0)*exp(-1/(2*H_r*tau*(z+1)*(z+1)))*(s_De_gain*s_He_perte/S_De_perte);
}*/
