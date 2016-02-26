set logscale
set xrange [2:80]
plot "Output/Cascade_Spectrum_Folder/Spectrum_total_photon_triangular_m140_z648069_0iterations.dat" u 1:2 w l , \
"Output/Cascade_Spectrum_Folder/Spectrum_total_electron_triangular_m140_z648069_0iterations.dat" u 1:2 w l , \
"Output/Cascade_Spectrum_Folder/Spectrum_total_photon_triangular_m140_z648069_triangular.dat" u 1:2 , \
"Output/Cascade_Spectrum_Folder/Spectrum_total_electron_triangular_m140_z648069_triangular.dat" u 1:2, \
"Output/Cascade_Spectrum_Folder/Spectrum_universal_photon_m140_z648069.dat" u 1:2
