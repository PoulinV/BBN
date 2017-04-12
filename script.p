#set logscale
plot 'test/5000steps_T10eV_100GeV_photon.datphoton.dat' u (log10($1)):(log10($2)) w l ,\
'test/5000steps_T10eV_100GeV_photon.datelectron.dat' u (log10($1)):(log10($2)) w l ,\
'test/5000steps_T10eV_100GeV_checkElectronRate_photon.dat' u (log10($1)):(log10($2)) w l ,\
'test/5000steps_T10eV_100GeV_checkElectronRate_electron.dat' u (log10($1)):(log10($2)) w l ,\
'output/Cascade_Spectrum_Folder/Spectrum_universal_photon_m200000_z42563.4.dat' u (log10($1)):(log10($2))
