set logscale x
set xrange [1:1e5]
set yrange [0:40]
#f1(x)=log10(1.4e2*(1e-12)*(0.001*x)**(-1.56))
#f2(x)=log10(6.2e-12*(1e-12)*(x*10**-3)**(-5.07)*1e30)

plot "Output/Cascade_Spectrum_Folder/Spectrum_total_photons_m200000_z2000_triangular.dat" u 1:(log10($2)) w l t 'total photons', \
"Output/Cascade_Spectrum_Folder/Spectrum_total_electrons_m200000_z2000_triangular.dat" u 1:(log10($2)) w l t 'total electrons', \
"Output/Cascade_Spectrum_Folder/Spectrum_universal_photon_m200000_z2000.dat" u 1:(log10($2)) w l t 'universal' lc 'black'
#f1(x) w l t 'Kawasaki & Moroi', \
#f2(x) w l
