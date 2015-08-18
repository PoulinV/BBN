set logscale
set xrange [2:80]
plot "Cascade_Spectrum_Folder/Spectrum_Cascade_m140_z425643_0iterations.dat" u 1:2 w l , \
"Cascade_Spectrum_Folder/Spectrum_Cascade_m140_z425643_1iterations.dat" u 1:2 w l , \
"Cascade_Spectrum_Folder/Spectrum_Cascade_m140_z425643_2iterations.dat" u 1:2 w l , \
"Cascade_Spectrum_Folder/Spectrum_Cascade_m140_z425643_3iterations.dat" u 1:2 w l , \
"Cascade_Spectrum_Folder/Spectrum_Cascade_m140_z425643_4iterations.dat" u 1:2 w l , \
"Cascade_Spectrum_Folder/Spectrum_Cascade_m140_z425643_5iterations.dat" u 1:2 w l , \
"Cascade_Spectrum_Folder/Spectrum_Cascade_m140_z425643_6iterations.dat" u 1:2 w l , \
"Output/Universal_spectrum.dat" u 1:2 w l , \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:9 , \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:3 , \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:4 , \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:5 , \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:6 , \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:7 , \
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:8 ,\
"/Users/poulin/Documents/Doc_Labo/ProgrammeBBN/ExtraProgram/OutputExtra/SpectreTotaliterations_high_precision.dat" u 1:11
