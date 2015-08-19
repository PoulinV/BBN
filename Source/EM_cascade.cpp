#include "../Include/EM_cascade.h"
#include "../Include/Injected_spectrum.h"
#include "../Include/tools.h"
double  dsigma_compton(double  x, double  z, double g){

	double  dsigma =0;
	if(g<x && g>(x/(1+2*x)))dsigma = pi*pow(r_e,2)*m_e/pow(x,2)*(x/g+g/x+pow(m_e/g-m_e/x,2)-2*m_e*(1/g-1/x))*n_e*pow(1+z,3);
 	return dsigma;

}
double  dsigma_phph(double  x, double  z,  double g){

	double  dsigma = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(a*r_e,2)*pow(m_e,-6)*pow(x,2)*pow(1-g/x+pow(g/x,2),2)/(63*10125*pi);
	return dsigma;

}

// double dsigma_NPC(double x, double z, double g){
//
// 	double p = sqrt(x*x-m_e*m_e);
// 	double L = log((x*x+p*p+m_e*m_e)/(x*x-p*p+m_e*m_e));
// 	double l = log((x+p)/(x-p));
//
// 	return a*pow(r_e,2)*pow(1+z,3)*eta*n_y_0*
//
// }

double  rate_compton(double  x, double  z){

	double X = 2*x/m_e;
	double sigma_cs = 2*pi*pow(r_e,2)/X*((1-4/X-8/pow(X,2))*log(1+X)+1/2+8/X+1/(2*pow(1+X,2)));
	double Gamma = sigma_cs*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3);
	return Gamma;

}
double  rate_NPC(double  x, double  z){

	double  k = x/m_e;
	double  rho = (2*k-4)/(k+2+2*pow(2*k,0.5));
	double  sigma_PCN;
	sigma_PCN = a*pow(r_e,2)*(28/9*log(2*k)-218/27) ;
	double  Gamma_2 = sigma_PCN*pow(1+z,3)*eta*n_y_0;
	return Gamma_2;

}

double  rate_gg_scattering(double  x, double  z){

	double  Gamma_3;
	Gamma_3 = 0.1513*pow(a,4)*m_e*pow(x/m_e,3)*pow(T_0*(1+z)/m_e,6);
	return Gamma_3;

}
double Function_Integrand_Spectre_Compton(double E_gamma, double E_e, double E_gamma_bar){
	double f;
	double Gamma_e = 4*E_gamma_bar*E_e/(m_e*m_e);
	double q = E_gamma/(Gamma_e*(E_e-E_gamma));

	if(q>=0 && q<=1) {f = 2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q);
}
	else f=0;
	return f;
}

double Function_Integrand_Spectre_Compton_version_q(double q, double E_e, double E_gamma_bar){
	double f;
	double Gamma_e = 4*E_gamma_bar*E_e/(m_e*m_e);
	double E_gamma=Gamma_e*E_e*q/(1+Gamma_e*q);

	if(q>=0 && q<=1) {f = (2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q))*pow(E_e-E_gamma,2)*Gamma_e/E_e;
}
	else f=0;
	if(f<=0)f=0;
	// cout <<"q = "<< q << " Gamma_e = "<<Gamma_e<<" E_gamma = "<<E_gamma<<" E_e ="<<E_e<<" f = " << f << endl;

	return f;
}
void  Spectre_electron_compton(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
													 struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
													 struct Structure_Spectrum * pt_Gamma_Spectrum,
													 struct Structure_Spectrum * pt_Electron_Spectrum){

 	double z = pt_Electron_Spectrum->redshift;
	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
	double E_0 = pt_Particle_Physics_Model->E_0;
	double E_gamma_bb = 2.701*T_0*(1+z);
	double E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	double resultat = 0, Gamma_electron = 0;
	double dE, dE_2, dq, h;
	double E1, E2, E3, f1, f2, f3, F1, F2, F3;
	double E_gamma, E_e;
	double E_gamma_1, E_gamma_2, E_gamma_3;
	double q_min=0.0001, q1, q2, q3;
	double resultat_monochromatique_1, resultat_monochromatique_2;
	dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Electron_Table_Size-1);
	dq = (1-q_min)/ (double) (n_step-1);
	for(int j =0; j<Electron_Table_Size-1;j++){
		resultat=0;
		Gamma_electron=0;
		dE_2 = (pt_Particle_Physics_Model->E_0 - (E_min+j*dE))/ (double) (n_step-1);
		h = dE_2/2;
		for(int i=0;i<n_step-1;i++){
			if(i==0){
				E1=E_min+i*dE_2+j*dE;
				q1=q_min+i*dq;
				E_e = E1;
				}
			else{
				E1=E3;
				E_gamma_1=E_gamma_3;
				q1=q3;
			}

			E2=E1 + h;
			E3=E2 + h;
			q2=q1+dq/2.;
			q3=q2+dq/2.;
			F1 = Function_Integrand_Spectre_Compton_version_q(q1,E_e, E_gamma_bb);
			F2 = Function_Integrand_Spectre_Compton_version_q(q2,E_e, E_gamma_bb);
			F3 = Function_Integrand_Spectre_Compton_version_q(q3,E_e, E_gamma_bb);
			linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E1, f1);
			linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E2, f2);
			linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E3, f3);
			f1*=(dsigma_compton(E1,z,(E1+m_e-E_e)));
			// cout << "f1 = " << f1<<endl;
			f1/=(rate_NPC(E1,z)+rate_compton(E1,z)+rate_gg_scattering(E1,z));
			f2*=(dsigma_compton(E2,z,(E2+m_e-E_e)));
			f2/=(rate_NPC(E2,z)+rate_compton(E2,z)+rate_gg_scattering(E2,z));
			f3*=(dsigma_compton(E3,z,(E3+m_e-E_e)));
			f3/=(rate_NPC(E3,z)+rate_compton(E3,z)+rate_gg_scattering(E3,z));

			// cout << "E = " << E_min+j*dE<<" E1 = " << E1 << endl;
			resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
			Gamma_electron += dq/2*(F1/3+4.*F2/3+F3/3.);
		}
		E_e = (E_min+j*dE);
		if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice=="Dirac")resultat += (dsigma_compton(E_0,z,(E_0+m_e-E_e)))/(rate_NPC(E_0,z)+rate_compton(E_0,z)+rate_gg_scattering(E_0,z)) ;
		pt_Electron_Spectrum->Energy[j]=E_e;
		Gamma_electron*=2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb*E_e*E_e);
		pt_Electron_Spectrum->Spectrum[j]=resultat;
		pt_Electron_Spectrum->Spectrum[j]/=(Gamma_electron);
		// cout << " E = " << pt_Electron_Spectrum->Energy[j] << " resultat = " << pt_Electron_Spectrum->Spectrum[j] << " Gamma_electron = " << Gamma_electron << endl;

	}

}
void Spectrum_electron_scattered(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																 struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
																 struct Structure_Spectrum * pt_Input_Electron_Spectrum,
															 	 struct Structure_Spectrum * pt_Output_Electron_Spectrum){

	  	double z = pt_Input_Electron_Spectrum->redshift;
		 	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
		 	double E_0 = pt_Particle_Physics_Model->E_0;
		 	double E_gamma_bb = 2.701*T_0*(1+z);
		 	double E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
		 	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
		 	double resultat = 0, Gamma_electron = 0, Gamma_electron_Dirac=0;
		 	double dE, dE_2, dq, h;
		 	double E1, E2, E3, f1, f2, f3, F1, F2, F3;
		 	double E_gamma, E_e;
		 	double E_gamma_1, E_gamma_2, E_gamma_3;
		 	double q_min=0.0001, q1, q2, q3;
		 	double resultat_monochromatique_1, resultat_monochromatique_2;
		 	dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Electron_Table_Size-1);
		 	dq = (1-q_min)/ (double) (n_step-1);

			for(int i=0;i<n_step-1;i++){
				if(i==0){
					q1=q_min+i*dq;
					}
				else{
					q1=q3;
				}
				q2=q1+dq/2.;
				q3=q2+dq/2.;
				F1 = Function_Integrand_Spectre_Compton_version_q(q1,E_0, E_gamma_bb);
				F2 = Function_Integrand_Spectre_Compton_version_q(q2,E_0, E_gamma_bb);
				F3 = Function_Integrand_Spectre_Compton_version_q(q3,E_0, E_gamma_bb);
				Gamma_electron_Dirac += dq/2*(F1/3+4.*F2/3+F3/3.);
			}
			Gamma_electron_Dirac*=2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb*E_0*E_0);
			// cout << " Gamma_electron = " << Gamma_electron<<endl;

		 	for(int j =0; j<Electron_Table_Size-1;j++){
		 		resultat=0;
		 		Gamma_electron=0;
		 		dE_2 = (pt_Particle_Physics_Model->E_0 - (E_min+j*dE))/ (double) (n_step-1);
		 		h = dE_2/2;
		 		for(int i=0;i<n_step-1;i++){
		 			if(i==0){
		 				E1=E_min+i*dE_2+j*dE;
		 				q1=q_min+i*dq;
		 				E_e = E1;
		 				}
		 			else{
		 				E1=E3;
		 				E_gamma_1=E_gamma_3;
		 				q1=q3;
		 			}

		 			E2=E1 + h;
		 			E3=E2 + h;
		 			q2=q1+dq/2.;
		 			q3=q2+dq/2.;
		 			F1 = Function_Integrand_Spectre_Compton(E1+E_gamma_bb-E_e,E1, E_gamma_bb);
		 			F2 = Function_Integrand_Spectre_Compton(E2+E_gamma_bb-E_e,E2, E_gamma_bb);
		 			F3 = Function_Integrand_Spectre_Compton(E3+E_gamma_bb-E_e,E3, E_gamma_bb);
					if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice!="Dirac"){
					if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E1, f1);
		 		 	else f1=0;
		 		 	if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E2, f2);
		 		 	else f2=0;
		 		 	if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E3, f3);
		 		 	else f3=0;
					// cout << " E1 = " << E1 << " f1 = " << f1 << endl;
					f1*=F1/(E1*E1);
					f2*=F2/(E2*E2);
					f3*=F3/(E3*E3);
		 			resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
					}
					F1 = Function_Integrand_Spectre_Compton_version_q(q1,E_e, E_gamma_bb);
					F2 = Function_Integrand_Spectre_Compton_version_q(q2,E_e, E_gamma_bb);
					F3 = Function_Integrand_Spectre_Compton_version_q(q3,E_e, E_gamma_bb);
		 			Gamma_electron += dq/2*(F1/3+4.*F2/3+F3/3.);
		 		}
		 		E_e = (E_min+j*dE);
		 		pt_Output_Electron_Spectrum->Energy[j]=E_e;
		 		Gamma_electron*=2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb*E_e*E_e);

		 		if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice=="Dirac"){
					pt_Output_Electron_Spectrum->Spectrum[j]=2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb*E_0*E_0*Gamma_electron_Dirac)*Function_Integrand_Spectre_Compton(E_0+E_gamma_bb-E_e,E_0,E_gamma_bb);
				}
				else pt_Output_Electron_Spectrum->Spectrum[j]=resultat*2*pi*r_e*r_e*m_e*m_e;
				// cout << " E = " << pt_Output_Electron_Spectrum->Energy[j] << " resultat spectre diffuse avant = " << pt_Output_Electron_Spectrum->Spectrum[j] << " Gamma_electron = " << Gamma_electron << endl;

					pt_Output_Electron_Spectrum->Spectrum[j]/=(Gamma_electron);
				// cout << " E = " << pt_Output_Electron_Spectrum->Energy[j] << " resultat spectre diffuse = " << pt_Output_Electron_Spectrum->Spectrum[j] << endl;
		}

}
void Spectre_gamma_compton(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
													 struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
													 struct Structure_Spectrum * pt_Electron_Spectrum,
													 struct Structure_Spectrum * pt_Gamma_Spectrum){

	double z = pt_Gamma_Spectrum->redshift;

	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
	double E_gamma_bb = 2.701*T_0*(1+z);
	double E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	double resultat = 0;
	double dE, dE_2, h;
	double E1, E2, E3, f1, f2, f3, F1, F2, F3;
	double E_gamma, E_e;
	double q_min=0.0001, q1, q2, q3, dq;
	double Gamma_electron=0;
	double E_0=pt_Particle_Physics_Model->E_0;
	dq = (1-q_min)/ (double) (n_step-1);
	dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Gamma_Table_Size-1);


		for(int j =0; j<Gamma_Table_Size;j++){
			resultat=0;
			dE_2 = (pt_Particle_Physics_Model->E_0 - (E_min+j*dE+m_e))/ (double) (n_step-1);
			h = dE_2/2;
			for(int i=0;i<n_step-1;i++){
				if(i==0){
					E1=E_min+i*dE_2+j*dE+m_e;
					E_gamma = E1-m_e;
					}
				else{
					E1=E3;
				}

				E2=E1 + h;
				E3=E2 + h;
				if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E1, f1);
				else f1=0;
				if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E2, f2);
				else f2=0;
				if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E3, f3);
				else f3=0;
				F1 = Function_Integrand_Spectre_Compton(E_gamma,E1, E_gamma_bb);
				F2 = Function_Integrand_Spectre_Compton(E_gamma,E2, E_gamma_bb);
				F3 = Function_Integrand_Spectre_Compton(E_gamma,E3, E_gamma_bb);
				// if(f1!=0)cout << "E1 = " << E1<< " f1 = " << f1 << endl;

				f1*=2*F1/(E1*E1);
				f2*=2*F2/(E2*E2);
				f3*=2*F3/(E3*E3);
				// cout << "E = " << E_min+j*dE<<"F1 = "<< F1 << " f1 apres = " << f1 << endl;

				resultat += h * (f1/3. + 4.*f2/3. + f3/3.);

			}
			pt_Gamma_Spectrum->Energy[j]=E_min+j*dE;
			pt_Gamma_Spectrum->Spectrum[j]=resultat*2*pi*r_e*r_e*m_e*m_e*int_bb/E_gamma_bb;
			if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice=="Dirac"){
				for(int i=0;i<n_step-1;i++){
					if(i==0){
						q1=q_min+i*dq;
						}
					else{
						q1=q3;
					}
					q2=q1+dq/2.;
					q3=q2+dq/2.;
					F1 = Function_Integrand_Spectre_Compton_version_q(q1,E_0, E_gamma_bb);
					F2 = Function_Integrand_Spectre_Compton_version_q(q2,E_0, E_gamma_bb);
					F3 = Function_Integrand_Spectre_Compton_version_q(q3,E_0, E_gamma_bb);
					Gamma_electron += dq/2*(F1/3+4.*F2/3+F3/3.);
				}
				Gamma_electron*=2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb*E_0*E_0);
				// cout << " Gamma_electron = " << Gamma_electron<<endl;
				pt_Gamma_Spectrum->Spectrum[j]+=2*pi*r_e*r_e*m_e*m_e*2*int_bb*Function_Integrand_Spectre_Compton(E_min+j*dE,E_0,E_gamma_bb)/(E_gamma_bb*Gamma_electron*E_0*E_0);
			}
			// cout << "E = "  << E_min+j*dE<< "resultat = " << 	pt_Gamma_Spectrum->Spectrum[j] << endl;

		}


}

void Cascade_Spectrum_Reading_From_File(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																				struct Structure_Spectrum * pt_Cascade_Spectrum,
																				double z,
																				int iterations){
  ostringstream os;
  string name;
	double tmp_Energy, tmp_Spectrum;
  pt_Cascade_Spectrum->redshift = z;

    os << "Cascade_Spectrum_Folder/Spectrum_m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << iterations <<"iterations.dat";
    name = os.str();
    ifstream file(name);
    if(file)cout << "Importing file " << name << " in structure Cascade_Spectrum." << endl;
		else{
			cout << "I couldn't recognize cascade spectrum file. Please check that it is in present in the folder Cascade_Specrum_Folder with proper name : 'Spectrum_mXXX_zXXX_XXXiterations.dat' corresponding to the value of m, z and iterations you are using."<<endl;
			return;
		}
    while(file){
			string line;
	    getline(file, line);
	    // stringstream is(ligne);
	    if(line[0] == '#' or line[0] == '\0') continue;
      file >> tmp_Energy >> tmp_Spectrum ;
			pt_Cascade_Spectrum->Energy.push_back(tmp_Energy);
			pt_Cascade_Spectrum->Spectrum.push_back(tmp_Spectrum*(rate_NPC(tmp_Energy,z)+rate_compton(tmp_Energy,z)+rate_gg_scattering(tmp_Energy,z)));
			  // cout << z << "  " << tmp_Energy << "  " << tmp_Spectrum<<endl;

		}

		file.close();


}

void  Cascade_Spectrum_Calculation(double z,
																	 struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																	 struct Structure_Spectrum * pt_Cascade_Spectrum,
																	 struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){
	double dE, dE_2, h, h2;
	// double E1, E2, E3, E_gamma, Old_spectrum;
	// double f1, f2, f3;
	double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7,E_gamma, Old_spectrum;
	double resultat, integrale;
	double E_c = E_c_0/(1+z);
	struct Structure_Spectrum Electron_Spectrum;
	struct Structure_Spectrum Inverse_Compton_Spectrum;
	struct Structure_Spectrum Diffused_Electron_Spectrum;

	double E_0 = pt_Particle_Physics_Model->E_0;
	dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Gamma_Table_Size-1);
	/*****Initialization****/
	pt_Cascade_Spectrum->Energy.resize(Gamma_Table_Size);
	pt_Cascade_Spectrum->Spectrum.resize(Gamma_Table_Size);
	pt_Cascade_Spectrum->species = "photon";
	Inverse_Compton_Spectrum.Energy.resize(Gamma_Table_Size);
	Inverse_Compton_Spectrum.Spectrum.resize(Gamma_Table_Size);
	Inverse_Compton_Spectrum.species = "photon";
	Electron_Spectrum.Energy.resize(Electron_Table_Size);
	Electron_Spectrum.Spectrum.resize(Electron_Table_Size);
	Electron_Spectrum.species="electron";
	Diffused_Electron_Spectrum.Energy.resize(Electron_Table_Size);
	Diffused_Electron_Spectrum.Spectrum.resize(Electron_Table_Size);
	Diffused_Electron_Spectrum.species = "electron";
	Diffused_Electron_Spectrum.spectrum_name = "diffused_electron_";

	Inverse_Compton_Spectrum.redshift=z;
	Electron_Spectrum.redshift=z;
	pt_Cascade_Spectrum->redshift = z;
	Diffused_Electron_Spectrum.redshift=z;

	if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "universal"){
					for(int i=0;i<Gamma_Table_Size;i++){
						E1=E_min+i*dE;
						pt_Cascade_Spectrum->Energy[i]=E1;
						pt_Cascade_Spectrum->Spectrum[i]=Universal_Spectrum(E1,z,pt_Particle_Physics_Model->E_0);
					}
					check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,integrale);
					for(int i=0;i<Gamma_Table_Size;i++){
						pt_Cascade_Spectrum->Spectrum[i]=pt_Particle_Physics_Model->E_0/integrale*pt_Cascade_Spectrum->Spectrum[i];
					}
					check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,integrale);

	}
	else{
					if(E_c <= pt_Particle_Physics_Model->E_0){
								for(int i=0;i<Gamma_Table_Size;i++){
									E1=E_min+i*dE;
									pt_Cascade_Spectrum->Energy[i]=E1;
									pt_Cascade_Spectrum->Spectrum[i]=Universal_Spectrum(E1,z,pt_Particle_Physics_Model->E_0);
								}
								check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,integrale);
								for(int i=0;i<Gamma_Table_Size;i++){
									pt_Cascade_Spectrum->Spectrum[i]=pt_Particle_Physics_Model->E_0/integrale*pt_Cascade_Spectrum->Spectrum[i];
								}
								check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,integrale);
								if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
										if(verbose>1)cout <<" I will now print the spectrum in files." << endl;
										print_spectrum_automatic_names(0, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
								}


					}


					else{
									/********First step : compute initial ICS spectrum from the electon spectrum injected**********/
									if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "none"){
											for(int i=0;i<Gamma_Table_Size;i++)Inverse_Compton_Spectrum.Spectrum[i]=0;
									}
								  else{
											for(int i=0;i<Electron_Table_Size;i++){
													E1=E_min+i*dE;
													Electron_Spectrum.Energy[i]=E1;
													if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice=="Dirac")Electron_Spectrum.Spectrum[i]=0;
													else Electron_Spectrum.Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum(E1,z,pt_Particle_Physics_Model->E_0);
											}
									Spectrum_electron_scattered(pt_Particle_Physics_Model,
 																					 	pt_Spectrum_and_Precision_Parameters,
 																						&Electron_Spectrum,
 																						&Diffused_Electron_Spectrum);
									if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
										Diffused_Electron_Spectrum.spectrum_name = "electron_diffused_from_dirac_";
										print_spectrum_automatic_names(101, &Diffused_Electron_Spectrum, pt_Particle_Physics_Model);
									}
									Spectre_gamma_compton(pt_Particle_Physics_Model,
																				pt_Spectrum_and_Precision_Parameters,
																				&Electron_Spectrum,
																				&Inverse_Compton_Spectrum);
									if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
											Inverse_Compton_Spectrum.spectrum_name = "ICS_from_e_injection_";
											if(verbose>1)cout <<" I will now print the spectrum in files." << endl;
											print_spectrum_automatic_names(0, &Inverse_Compton_Spectrum, pt_Particle_Physics_Model);
									}
								}
								/*********Second step : compute the gamma spectrum from gg->gg, ge->ge, gN->Nee for a certain number of iterations.*********/

									for(int i=0;i<Gamma_Table_Size;i++){
												E1=E_min+i*dE;
												pt_Cascade_Spectrum->Energy[i]=E1;
												if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice=="none")pt_Cascade_Spectrum->Spectrum[i]=0;
												else pt_Cascade_Spectrum->Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E1,z,pt_Particle_Physics_Model->E_0)/(rate_NPC(E1,z)+rate_compton(E1,z)+rate_gg_scattering(E1,z));
												// pt_Cascade_Spectrum->Spectrum[i]+=Inverse_Compton_Spectrum.Spectrum[i];
										}
										// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,integrale);
										// for(int i=0;i<Gamma_Table_Size;i++){
										// 	pt_Cascade_Spectrum->Spectrum[i]=pt_Particle_Physics_Model->E_0/integrale*pt_Cascade_Spectrum->Spectrum[i];
										// }
										// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,integrale);

										if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
												pt_Cascade_Spectrum->spectrum_name = "Cascade_";
												if(verbose>1)cout <<" I will now print the spectrum in files." << endl;
												print_spectrum_automatic_names(0, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
										}
										/*double E_gamma_minus_1, E_gamma_plus_1, E_j, E_j_minus_1, E_j_plus_1, dE_j;

												for(double i = (Gamma_Table_Size-1); i>=0 ;i--){
												resultat=0;
												E_gamma = E_min*pow(E_0/E_min,i/(Gamma_Table_Size-1));
												E_gamma_plus_1 = E_min*pow(E_0/E_min,(i+1)/(Gamma_Table_Size-1));
												E_gamma_minus_1 = E_min*pow(E_0/E_min,(i-1)/(Gamma_Table_Size-1));
												dE = (E_gamma_plus_1 - E_gamma_minus_1)/2.;
												pt_Cascade_Spectrum->Energy[i] = E_gamma;
												pt_Cascade_Spectrum->Spectrum[i] = 0;
												if(i==(Gamma_Table_Size-1)){
													pt_Cascade_Spectrum->Spectrum[i]=1/(dE*((rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z))));
												}
												else{
													for (int j = i+1 ; j < Gamma_Table_Size ; j ++){
														E_j = E_min*pow(E_0/E_min,j/(Gamma_Table_Size-1));
														E_j_plus_1 = E_min*pow(E_0/E_min,(j+1)/(Gamma_Table_Size-1));
														E_j_minus_1 = E_min*pow(E_0/E_min,(j-1)/(Gamma_Table_Size-1));
														dE_j = (E_j_plus_1 - E_j_minus_1)/2.;
														if(j==i+1)pt_Cascade_Spectrum->Energy[i] = E_gamma;
														pt_Cascade_Spectrum->Spectrum[i] += dE_j * pt_Cascade_Spectrum->Spectrum[j] * (dsigma_phph(E_j,z,E_gamma)+dsigma_compton(E_j,z,E_gamma))/(rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z));
													}
												}
												cout << " E = "<< pt_Cascade_Spectrum->Energy[i] << " pt_Cascade_Spectrum->Spectrum[i] = " << pt_Cascade_Spectrum->Spectrum[i] << endl;
												}

											check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,integrale);
											// for(int i=0;i<Gamma_Table_Size;i++){
											// 	pt_Cascade_Spectrum->Spectrum[i]=pt_Particle_Physics_Model->E_0/integrale*pt_Cascade_Spectrum->Spectrum[i];
											// }
											// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,integrale);

											if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
												pt_Cascade_Spectrum->spectrum_name = "Cascade_";
												print_spectrum_automatic_names(1, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
											}*/


									for(int k = 0; k<pt_Spectrum_and_Precision_Parameters->iterations;k++){
												if(verbose>1)cout<<"iteration : " << k+1 << endl;
												/*for(int j =0; j<Gamma_Table_Size;j++){
												resultat=0;
												dE_2 = (pt_Particle_Physics_Model->E_0 - (E_min+j*dE))/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
												h = dE_2/2;
												for(int i=0;i<pt_Spectrum_and_Precision_Parameters->n_step-1;i++){
													if(i==0){
														E1=E_min+j*dE+i*dE_2;
														E_gamma = E1;
														}
													else{
														E1=E3;
													}

													E2=E1 + h;
													E3=E2 + h;

													linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E1, f1);
													linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E2, f2);
													linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E3, f3);
													//
													// f1 *=(dsigma_phph(E1,z,E_gamma)+dsigma_compton(E1,z,E_gamma))/(rate_NPC(E1,z)+rate_compton(E1,z)+rate_gg_scattering(E1,z));
													// f2 *=(dsigma_phph(E2,z,E_gamma)+dsigma_compton(E2,z,E_gamma))/(rate_NPC(E2,z)+rate_compton(E2,z)+rate_gg_scattering(E2,z));
													// f3 *=(dsigma_phph(E3,z,E_gamma)+dsigma_compton(E3,z,E_gamma))/(rate_NPC(E3,z)+rate_compton(E3,z)+rate_gg_scattering(E3,z));
													// //
													f1 *= (dsigma_phph(E1,z,E_gamma)+dsigma_compton(E1,z,E_gamma));
													f2 *= (dsigma_phph(E2,z,E_gamma)+dsigma_compton(E2,z,E_gamma));
													f3 *= (dsigma_phph(E3,z,E_gamma)+dsigma_compton(E3,z,E_gamma));
													//
													// f1 *= (dsigma_compton(E1,z,E_gamma));
													// f2 *= (dsigma_compton(E2,z,E_gamma));
													// f3 *= (dsigma_compton(E3,z,E_gamma));

													resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
												}
												pt_Cascade_Spectrum->Energy[j]=E_min+j*dE;
												pt_Cascade_Spectrum->Spectrum[j]=resultat/(rate_NPC(E_min+j*dE,z)+rate_compton(E_min+j*dE,z)+rate_gg_scattering(E_min+j*dE,z));
												// pt_Cascade_Spectrum->Spectrum[j]+=Inverse_Compton_Spectrum.Spectrum[j];
												if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice!="none")pt_Cascade_Spectrum->Spectrum[j]+=pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E_min+j*dE,z,pt_Particle_Physics_Model->E_0)/(rate_NPC(E_min+j*dE,z)+rate_compton(E_min+j*dE,z)+rate_gg_scattering(E_min+j*dE,z));

											}*/



											for(int i=0;i<Gamma_Table_Size;i++){
											  E_gamma = E_min+i*dE;
											  resultat = 0;
											        for(int j=0; j<pt_Spectrum_and_Precision_Parameters->n_step-1;j++){
											          dE_2 = (pt_Particle_Physics_Model->E_0 - (E_gamma))/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
											          h2 = dE_2/6.;

											          if(j==0){
											            E1=E_gamma;
											            }
											          else{
											            E1=E7;
											          }

											          E2=E1 + h2;
											          E3=E1 + 2*h2;
											          E4=E1 + 3*h2;
											          E5=E1 + 4*h2;
											          E6=E1 + 5*h2;
											          E7=E1 + 6*h2;
											          if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E1, f1);
											          else f1=0;
											          if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E2, f2);
											          else f2=0;
											          if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E3, f3);
											          else f3=0;
											          if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E4, f4);
											          else f4=0;
											          if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E5, f5);
											          else f5=0;
											          if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E6, f6);
											          else f6=0;
											          if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E7, f7);
											          else f7=0;
											          // f1*=(dsigma_compton(E1,z,E_gamma));
											          // f2*=(dsigma_compton(E2,z,E_gamma));
											          // f3*=(dsigma_compton(E3,z,E_gamma));
											          // f4*=(dsigma_compton(E4,z,E_gamma));
											          // f5*=(dsigma_compton(E5,z,E_gamma));
											          // f6*=(dsigma_compton(E6,z,E_gamma));
											          // f7*=(dsigma_compton(E7,z,E_gamma));
											          f1*=(dsigma_phph(E1,z,E_gamma)+dsigma_compton(E1,z,E_gamma));
											          f2*=(dsigma_phph(E2,z,E_gamma)+dsigma_compton(E2,z,E_gamma));
											          f3*=(dsigma_phph(E3,z,E_gamma)+dsigma_compton(E3,z,E_gamma));
											          f4*=(dsigma_phph(E4,z,E_gamma)+dsigma_compton(E4,z,E_gamma));
											          f5*=(dsigma_phph(E5,z,E_gamma)+dsigma_compton(E5,z,E_gamma));
											          f6*=(dsigma_phph(E6,z,E_gamma)+dsigma_compton(E6,z,E_gamma));
											          f7*=(dsigma_phph(E7,z,E_gamma)+dsigma_compton(E7,z,E_gamma));


											          resultat += dE_2/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);


											          // cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << " i = " << i << endl;
											        }
															pt_Cascade_Spectrum->Energy[i]=E_gamma;
															pt_Cascade_Spectrum->Spectrum[i]=resultat/(rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z));
															// pt_Cascade_Spectrum->Spectrum[j]+=Inverse_Compton_Spectrum.Spectrum[j];
															if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice!="none")pt_Cascade_Spectrum->Spectrum[i]+=pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E_gamma,z,pt_Particle_Physics_Model->E_0)/(rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z));

											  }
											// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,integrale);
											// for(int i=0;i<Gamma_Table_Size;i++){
											// 	pt_Cascade_Spectrum->Spectrum[i]=pt_Particle_Physics_Model->E_0/integrale*pt_Cascade_Spectrum->Spectrum[i];
											// }
											// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,integrale);

											if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
												pt_Cascade_Spectrum->spectrum_name = "Cascade_";
												print_spectrum_automatic_names(k+1, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
											}

								}
								/**********Compute the associated electron spectrum and the gamma spectrum from ICS.**********/
								if(pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering == "yes"){
									pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice = "from_cascade";
									Spectre_electron_compton(pt_Particle_Physics_Model,
																					 pt_Spectrum_and_Precision_Parameters,
																					 pt_Cascade_Spectrum,
																					 &Electron_Spectrum);
								//  check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,&Electron_Spectrum,integrale);

								//  Spectrum_electron_scattered(pt_Particle_Physics_Model,
								// 													 	pt_Spectrum_and_Precision_Parameters,
								// 														&Electron_Spectrum,
								// 														&Diffused_Electron_Spectrum);
								 //
								//  check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,&Diffused_Electron_,&integraleSpectrum);
								 //
								// 	for(int j=0; j<Electron_Table_Size; j++){
								// 		Electron_Spectrum.Spectrum[j]+=Diffused_Electron_Spectrum.Spectrum[j];
								// 	}
									Spectre_gamma_compton(pt_Particle_Physics_Model,
																				pt_Spectrum_and_Precision_Parameters,
																				&Electron_Spectrum,
																				&Inverse_Compton_Spectrum);
									for(int i = 0; i<Gamma_Table_Size ; i++){
										pt_Cascade_Spectrum->Spectrum[i]+=Inverse_Compton_Spectrum.Spectrum[i];
										// cout << " Inverse_Compton_Spectrum.Spectrum[i] = " << Inverse_Compton_Spectrum.Spectrum[i] << endl;
									}
									if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
										Electron_Spectrum.spectrum_name = "electron_";
										print_spectrum_automatic_names(0, &Electron_Spectrum, pt_Particle_Physics_Model);
									}
									if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
										pt_Cascade_Spectrum->spectrum_name = "total_";
										print_spectrum_automatic_names(100, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
									}
									check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,integrale);

									if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
										Inverse_Compton_Spectrum.spectrum_name = "ICS_from_cascade_";
										print_spectrum_automatic_names(101, &Inverse_Compton_Spectrum, pt_Particle_Physics_Model);
									}
								}
					}
	}



}
