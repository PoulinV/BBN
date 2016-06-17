import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec

x1,x2 = np.loadtxt("output/Cascade_Spectrum_Folder/Spectrum_universal_photon_m200000_z425643.dat",unpack=True)
y1,y2 = np.loadtxt("spectre_kawasaki_moroi_z425643_m200000.dat",unpack=True)

# # z1,z2,z3,z4,z5 = np.loadtxt("class_public-2.4.2/output_cl/Standard_Cl_79_cl.dat",unpack=True);
# #

plt.figure()
plt.axis([1,1e4,1e-4,50])
plt.tick_params(axis='x',labelsize=16)
plt.tick_params(axis='y',labelsize=16)
# plt.axis([2,2000,-0.1,0.1])

# line,=ax1.plot(x1,x2,'r',label=r'Semi-analytical model',lw=2)
dashpoint = [20,3,2,3]
line1,=plt.plot(y1,y2,'g',lw=2,label='Kawasaki and Moroi')
line1.set_dashes(dashpoint)
line2,=plt.plot(x1,x2,'b',lw=2,label='Universal')


# f.subplots_adjust(hspace=0)
# line5,=plt.plot([],[],'k',lw=2);
# line6,=plt.plot([],[],'k-',lw=2);
# dashpoint = [20,3,2,3]
# line6.set_dashes(dashpoint)
# leg_zreio=ax1.legend([line5,line6],['Standard','Semi-analytical model'],fontsize=15,handlelength=2.5,loc='upper right',numpoints=1, frameon=0,title='Stars')
# leg_zreio.get_title().set_fontsize('15');
# plt.gca().add_artist(leg_zreio)
#
# line4,=plt.plot([],[],'g',lw=2)
# line5,=plt.plot([],[],'b',lw=2);
# # line6,=plt.plot([],[],'r',lw=2);
# leg_DM=ax1.legend([line4,line5,line6],['$z_\mathrm{h} = 30$, $\mathcal{B}(z=0)=10^6$','$z_\mathrm{h} = 30$, $\mathcal{B}(z=0)=10^7$','$z_\mathrm{h} = 30$, $\mathcal{B}(z=0)=10^8$'],fontsize=15,handlelength=2,loc='upper left',numpoints=1, frameon=0,title='Halos')
# leg_DM.get_title().set_fontsize('15');
# plt.gca().add_artist(leg_DM)

plt.semilogx()
plt.semilogy()
# ax2.semilogx()
# ax1.semilogy()
plt.legend(prop={'size':18},loc='upper left',numpoints=1,frameon=False,handlelength=3)
plt.title(r'$T = 10$ keV, $E_0 = 100$ GeV',fontsize=18,y=1.01)
plt.savefig('universal_vs_moroi_T10keV_E0100GeV.pdf', bbox_inches='tight')


plt.show()
