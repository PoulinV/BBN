#include "bbn/EM_cascade.h"
#include "bbn/injected_spectrum.h"
#include "bbn/structures.h"
#include "bbn/BBN_constraints.h"
#include "bbn/tools.h"
#include "bbn/test_functions.h"

using namespace std;


int main(int argc, char** argv){
  if(argc!=2){
    cout <<"Error : you have not given an input file! \nPlease, restart cBBNFast with one '.ini' file. If you want to use default parameters, just give 'default_param.ini'." << endl;
    return 0;
  }


    time_t t1,t2;
    double duree;
    t1 = time(NULL);

    map_parameters map_parameters;
    ifstream file_default("default_param.ini");
    ifstream file_input(argv[1]);
    fill_default_parameters(file_default,map_parameters);
    string task = "task_test";
    if(argc==2)get_parameter_from_file(file_input,task);
    if(task=="default" || argc==1) task=map_parameters["task_test"];



struct Structure_Spectrum_and_Precision_Parameters Spectrum_and_Precision_Parameters;
fill_structure_spectrum_and_precision_parameters(file_input, map_parameters, &Spectrum_and_Precision_Parameters);
struct Structure_Particle_Physics_Model Particle_Physics_Model;
fill_structure_particle_physics_model(file_input, map_parameters, &Particle_Physics_Model);
struct Structure_Output_Options Output_Options;
fill_structure_output_options(file_input, map_parameters, &Output_Options);

if(task == "print_polylog"){
  string name_file = "test_files";
  if(argc==2)get_parameter_from_file(file_input,name_file);
  if(name_file=="default" || argc==1) task=map_parameters["name_file"];
  ostringstream os;
  if(name_file=="automatic"){
    os << "test/polylog_2.dat";
  }
  else os << name_file << ".dat";
  name_file = os.str();
  ofstream file_output(name_file);

  for(double z = -50;z<1;z+=0.1){
    file_output << " z = " << z << " " << polylog_2(z, &Spectrum_and_Precision_Parameters) << endl;
  }

}

if(task == "print_interaction_rate"){

  string redshift = "redshift";
  double redshift_d;
  string temperature = "temperature";
  if(argc==2)get_parameter_from_file(file_input,redshift);
  if(redshift=="default" && argc!=1){
      get_parameter_from_file(file_input,temperature);
      if(temperature!="default")redshift_d=atof(temperature.c_str())/T_0-1;
      else redshift_d=atof(map_parameters["redshift"].c_str());
    }
    else redshift_d=atof(map_parameters["redshift"].c_str());


  print_interaction_rate(redshift_d,
                         atof(map_parameters["E_min_table"].c_str()),
                         atof(map_parameters["m_x"].c_str())/2.,
                         &Output_Options,
                         &Spectrum_and_Precision_Parameters);
}

if(task == "print_func_kawmor"){
  string redshift = "redshift";
  double redshift_d;
  string temperature = "temperature";
  if(argc==2)get_parameter_from_file(file_input,redshift);
  if(redshift=="default" && argc!=1){
      get_parameter_from_file(file_input,temperature);
      if(temperature!="default")redshift_d=atof(temperature.c_str())/T_0-1;
      else redshift_d=atof(map_parameters["redshift"].c_str());
    }
    else redshift_d=atof(map_parameters["redshift"].c_str());

    print_func_kawmor( redshift_d, atof(map_parameters["m_x"].c_str())/2., &Spectrum_and_Precision_Parameters);
}
if(task == "integrate_dsigma_compton" || task =="integrate_dsigma_phph" || task == "integrate_dsigma_pair_creation" || task == "integrate_dsigma_NPC"){

  string redshift = "redshift";
  double redshift_d, E_g;
  string temperature = "temperature";
  if(argc==2)get_parameter_from_file(file_input,redshift);
  if(redshift=="default" && argc!=1){
      get_parameter_from_file(file_input,temperature);
      if(temperature!="default")redshift_d=atof(temperature.c_str())/T_0-1;
      else redshift_d=atof(map_parameters["redshift"].c_str());
    }
    else redshift_d=atof(map_parameters["redshift"].c_str());
    for(int i = (Spectrum_and_Precision_Parameters.Energy_Table_Size-1); i>=0 ; i--){

    E_g = Spectrum_and_Precision_Parameters.E_min_table*pow((Particle_Physics_Model.E_0)/Spectrum_and_Precision_Parameters.E_min_table,(double) i/(Spectrum_and_Precision_Parameters.Energy_Table_Size-1));

    if(task == "integrate_dsigma_compton"){
        integrate_dsigma_compton(E_g/(m_e+2*E_g),E_g,redshift_d,&Spectrum_and_Precision_Parameters,&Output_Options);
      }
    if(task == "integrate_dsigma_pair_creation"){
      integrate_dsigma_pair_creation(m_e,E_g,redshift_d,&Spectrum_and_Precision_Parameters,&Output_Options);
    }
    if(task == "integrate_dsigma_phph")integrate_dsigma_phph(E_g*pow(10,-3),E_g,redshift_d,&Spectrum_and_Precision_Parameters,&Output_Options);
    if(task == "integrate_dsigma_NPC")integrate_dsigma_phph(E_g*pow(10,-3),E_g,redshift_d,&Spectrum_and_Precision_Parameters,&Output_Options);

    }

}

file_input.close();
file_default.close();
t2 = time(NULL);
duree = difftime(t2,t1);
cout << " It has last : "<< duree << " s." << endl;
return 0;
}
