set logscale
set xrange [2:71]
plot "Cascade_Spectrum_Folder/Spectrum_total_electron_triangular_m140_z425643_0iterations_500.dat" u 1:2 w l, \
"Cascade_Spectrum_Folder/Spectrum_total_photon_triangular_m140_z425643_0iterations_500.dat" u 1:2 w l, \
"Cascade_Spectrum_Folder/Spectrum_total_electron_triangular_m140_z425643_0iterations.dat" u 1:2 w l, \
"Cascade_Spectrum_Folder/Spectrum_total_photon_triangular_m140_z425643_0iterations.dat" u 1:2 w l, \
"Cascade_Spectrum_Folder/Spectrum_total_electron_triangular_m140_z425643_0iterations_2000.dat" u 1:2 w l, \
"Cascade_Spectrum_Folder/Spectrum_total_photon_triangular_m140_z425643_0iterations_2000.dat" u 1:2 w l , \
"Cascade_Spectrum_Folder/Spectrum_total_electron_triangular_m140_z648069_0iterations.dat" u 1:2 w l , \
"Cascade_Spectrum_Folder/Spectrum_total_photon_triangular_m140_z648069_0iterations.dat" u 1:2 w l
