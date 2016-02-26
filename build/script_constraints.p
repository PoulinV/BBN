set logscale
set yrange [1e-12:1e-3]
set xrange [1e4:1e12]
plot 'Output/Result_Scan_Folder/Results_destruc_only_4He_m200000MeV_universal.dat' u 1:2  t '4He analytic'
