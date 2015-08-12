#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include "../Include/common.h"
void Spectrum_and_cross_sections_convolution(struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum,
                                            struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                            int i_min,
                                            int i_max,
                                            double &resultat,
                                            double z,
                                            long n_step,
                                            string spectrum_choice);

// void Check_model_from_destruction_only(string nuclei,
//                                        struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum,
//                                        struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
//                                        double &Abundance,
//                                        double z_initial,
//                                        double z_final,
//                                        int z_step,
//                                        long n_step,
//                                        int iterations,
//                                        string spectrum_choice,
//                                        string spectrum_mode);

void Compute_Constraints_from_destruction_only(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                               struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                               struct Structure_Scan_Parameters * pt_Scan_Parameters,
                                               struct Structure_Output_Options * pt_Output_Options);

void Compute_constraints_from_destruction_and_production(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                                         struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                                         struct Structure_Scan_Parameters * pt_Scan_Parameters,
                                                         struct Structure_Output_Options * pt_Output_Options);
// 
// void Check_model_from_destruction_and_production(string nuclei,
//                                                  struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum,
//                                                  struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
//                                                  double &Abundance,
//                                                  double z_initial,
//                                                  double z_final,
//                                                  int z_step,
//                                                  long n_step,
//                                                  int iterations,
//                                                  string spectrum_choice,
//                                                  string spectrum_mode);

double cross_section(double  x, int i);
void Check_nuclei(string nuclei,
                  int &i_min,
                  int &i_max,
                  int &k_min,
                  int &k_max,
                  double &Y_min,
                  double &Y_max,
                  double &Y_0);
