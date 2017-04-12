plot 'test/3000steps_T10eV_100GeV_noICS_photon.dat' u (log10($1)):(log10($2)) w l ,\
'test/3000steps_T10eV_100GeV_noICS_electron.dat' u (log10($1)):(log10($2)) w l ,\
'test/5000steps_T10eV_100GeV_noICS_photon.dat' u (log10($1)):(log10($2)) w l ,\
'test/5000steps_T10eV_100GeV_noICS_electron.dat' u (log10($1)):(log10($2)) w l, \
'test/10000steps_T10eV_100GeV_noICS_photon.dat' u (log10($1)):(log10($2)) w l ,\
'test/10000steps_T10eV_100GeV_noICS_electron.dat' u (log10($1)):(log10($2)) w l
