set logscale
plot 'Output/Result_Scan_Folder/Results_destruc_only_4He_m200000MeV_Dirac.dat' u 1:2 w l t 'Full calculation' , \
'Output/Result_Scan_Folder/Results_destruc_only_4He_m200000MeV_universal.dat' u 1:2 w l t 'analytical' , \
'Output/Result_Scan_Folder/Results_destruc_only_2H_m200000MeV_Dirac.dat' u 1:2 w l t 'Full calculation' , \
'Output/Result_Scan_Folder/Results_destruc_only_2H_m200000MeV_universal.dat' u 1:2 w l t 'analytical'
