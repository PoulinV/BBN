set logscale
plot "Output/Interaction_Rate_Folder/Interaction_Rate_at_z425643.dat" u 1:2 w l t 'yp->ee', \
"Output/Interaction_Rate_Folder/Interaction_Rate_at_z425643.dat" u 1:3 w l t 'ye->ye', \
"Output/Interaction_Rate_Folder/Interaction_Rate_at_z425643.dat" u 1:4 w l t 'yy->yy', \
"Output/Interaction_Rate_Folder/Interaction_Rate_at_z425643.dat" u 1:5 w l t 'yy->ee', \
"Output/Interaction_Rate_Folder/Interaction_Rate_at_z425643.dat" u 1:6 w l t 'tot', \
"Output/Interaction_Rate_Folder/Interaction_Rate_at_z425643.dat" u 1:7 w l t 'ey->ye', \
"Output/Interaction_Rate_Folder/Interaction_Rate_at_z425643.dat" u 1:8 w l t 'ey->ye v2', \
"Output/Interaction_Rate_Folder/Interaction_Rate_at_z425643.dat" u 1:9 w l t 'ey->ye v3'
