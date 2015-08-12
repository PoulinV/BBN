set logscale
set xrange [2:80]
plot "Cascade_Spectrum_Folder/Spectrum_m140_z425643_0iterations.dat" u 1:2 w l notitle, \
"Cascade_Spectrum_Folder/Spectrum_m140_z425643_1iterations.dat" u 1:2 w l notitle, \
"Cascade_Spectrum_Folder/Spectrum_m140_z425643_2iterations.dat" u 1:2 w l notitle, \
"Cascade_Spectrum_Folder/Spectrum_m140_z425643_3iterations.dat" u 1:2 w l notitle, \
"Cascade_Spectrum_Folder/Spectrum_m140_z425643_4iterations.dat" u 1:2 w l notitle, \
"Cascade_Spectrum_Folder/Spectrum_m140_z425643_5iterations.dat" u 1:2 w l notitle, \
"Cascade_Spectrum_Folder/Spectrum_m140_z425643_6iterations.dat" u 1:2 w l notitle, \
"Output/Universal_spectrum.dat" u 1:2 w l notitle, \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:9 notitle, \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:3 notitle, \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:4 notitle, \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:5 notitle, \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:6 notitle, \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:7 notitle, \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:8 notitle, \
"Cascade_Spectrum_Folder/Spectrum_m140_z425643_100iterations.dat" u 1:2 w l , \
"Cascade_Spectrum_Folder/Spectrum_m140_z425643_101iterations.dat" u 1:2 w l 
