set logscale
set xrange [2:80]
plot "Cascade_Spectrum_Folder/Spectrum_Cascade_m140_z425643_0iterations.dat" u 1:2 w l,\
"Cascade_Spectrum_Folder/Spectrum_Cascade_m140_z425643_1iterations.dat" u 1:2 w l,\
"Cascade_Spectrum_Folder/Spectrum_Cascade_m140_z425643_2iterations.dat" u 1:2 w l,\
"Cascade_Spectrum_Folder/Spectrum_total_m140_z425643_100iterations.dat" u 1:2 w l , \
"Cascade_Spectrum_Folder/Spectrum_ICS_from_cascade_m140_z425643_101iterations.dat" u 1:2 
