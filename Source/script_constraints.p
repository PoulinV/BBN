set logscale
set yrange [1e-12:1e-3]
set xrange [1e4:1e10]
plot 'Output/Results_destruc_only_4He_m5000MeV_Universal_Spectrum.dat' u 1:2 w l ,\
'Output/Results_destruc_only_4He_m60MeV_Dirac_Spectrum_7iterations.dat' u 1:2 w l ,\
'Output/Results_destruc_only_4He_m140MeV_Dirac_Spectrum_7iterations.dat' u 1:2 w l
