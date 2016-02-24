set logscale x
set xrange [1:1e6]
#f1(x)=log10(1.4e2*(1e-12)*(0.001*x)**(-1.56))
#f2(x)=log10(6.2e-12*(1e-12)*(x*10**-3)**(-5.07)*1e30)

plot "Output/Cascade_Spectrum_Folder/Spectrum_total_photons_m2000_z425643_triangular.dat" u 1:(log10($2)) w l t 'photons', \
"Output/Cascade_Spectrum_Folder/Spectrum_universal_photon_m2000_z425643.dat" u 1:(log10($2)) w l t 'universal', \
"Output/Cascade_Spectrum_Folder/Spectrum_total_electrons_m2000_z425643_triangular.dat" u 1:(log10($2)) w l t 'electrons'
#f1(x) w l t 'Kawasaki & Moroi', \
#f2(x) w l
