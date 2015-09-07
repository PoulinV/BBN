set logscale
set xrange [2:80]
plot "Output/Cascade_Spectrum_Folder/Spectrum_total_photon_triangular_m140_z648069_0iterations.dat" u 1:2 w l , \
"Output/Cascade_Spectrum_Folder/Spectrum_total_electron_triangular_m140_z648069_0iterations.dat" u 1:2 w l , \
"test/spectrum_photon_dirac_70MeV_photon.dat" u 1:2 , \
"test/spectrum_photon_dirac_70MeV_electron.dat" u 1:2
