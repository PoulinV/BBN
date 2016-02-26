set logscale
set yrange [1e-12:1e-3]
set xrange [1e4:1e12]
plot 'Output/Result_Scan_Folder/Results_destruc_only_4He_m200000MeV_universal.dat' u 1:2  t '4He analytic', \
'Output/Result_Scan_Folder/Results_destruc_only_4He_m200000MeV_Dirac.dat' u 1:2  t '4He full calculation', \
'Output/Result_Scan_Folder/Results_destruc_only_2H_m200000MeV_universal.dat' u 1:2  t '2H analytic', \
'2H_test.dat' u 1:2  t '2H test production and destruction', \
'Output/Result_Scan_Folder/Results_destruc_only_2H_m200000MeV_Dirac.dat' u 1:2  t '2H full calculation'
