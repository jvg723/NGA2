from cmath import pi
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
# from scipy.optimize import curve_fit
import string
import sys
import os
# from scipy import stats

mpl.rcParams['figure.constrained_layout.use'] = True

plot_name = "volume_weighted_pdf"

dimless_x_label = r'$d/D$'

# Parameters
lz=20 # ligament length in z direction
D=1.0 # initial ligament diameter
initial_volume=np.pi*(D/2.0)**2.0*lz
logbins=np.logspace(np.log10(0.01),np.log10(1.0), num=30, base=10.0) # binning

# # Start plot
# fig1, ax1 = plt.subplots(dpi=300,figsize=[6,5]) 
# fig2, ax2 = plt.subplots(dpi=300,figsize=[6,5])
# fig3, ax3 = plt.subplots(dpi=300,figsize=[6,5]) 

# os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/ligament/pdf_data')

# #01_lig (viscr=50,Case 1)
# diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
# diam,u,v,w,vel_total = np.genfromtxt('01_lig.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
# weights = np.pi/6.0*diam**3.0 # Weight by volume
# vals1,edges1 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
# edgemed1 = (edges1[:-1]+edges1[1:])/2
# ax1.plot(edgemed1,vals1,'-x',linewidth=2.0,label='Case 1')
# ax3.plot(edgemed1,vals1,'k-x',linewidth=2.0,label=r'$Oh_l=0.05$')
# # Figure 1: visc_r=50 comparison
# ax1.set_xlabel(dimless_x_label)
# ax1.set_ylabel(r'volume weighted PDF')
# ax1.legend()
# fig1.suptitle(r'$\mu_l/\mu_g=50$',fontsize=20)
# fig1.savefig(plot_name+'_viscr_50.png') 
# area = sum(np.diff(edges1)*vals1)
# print(area)

# #02_lig (viscr=500,Case 1)
# diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
# diam,u,v,w,vel_total = np.genfromtxt('02_lig.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
# weights = np.pi/6.0*diam**3.0 # Weight by volume
# vals2,edges2 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
# edgemed2 = (edges2[:-1]+edges2[1:])/2
# ax2.plot(edgemed2,vals2,'-x',linewidth=2.0,label='Case 1')
# ax3.plot(edgemed2,vals2,'r-o',linewidth=2.0,label=r'$Oh_l=0.5$')

# # Figure 2: visc_r=50 comparison
# ax2.set_xlabel(dimless_x_label)
# ax2.set_ylabel(r'volume weighted PDF')
# ax2.legend()
# fig2.suptitle(r'$\mu_l/\mu_g=500$',fontsize=20)
# fig2.savefig(plot_name+'_viscr_500.png') 

# # Figure 3: Ratio comparison
# ax3.set_xlabel(dimless_x_label)
# ax3.set_ylabel(r'volume weighted PDF')
# ax3.legend()
# fig3.suptitle(r'$\mu_l/\mu_g-comparison$',fontsize=20)
# fig3.savefig(plot_name+'_viscr_comparison.png') 


# # ax3.plot(edgemed2,vals2*np.sum(weights)/initial_volume,'r-o',linewidth=2.0,label=r'$\mu_l/\mu_g=500$')



# Start plot (for DONTDELETE_pdf_data)
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/ligament/DONTDELETE_pdf_data')
fig1, ax1 = plt.subplots(dpi=300,figsize=[6,5]) 
fig2, ax2 = plt.subplots(dpi=300,figsize=[6,5]) 
fig3, ax3 = plt.subplots(dpi=300,figsize=[6,5]) 
fig4, ax4 = plt.subplots(dpi=300,figsize=[6,5])
fig5, ax5 = plt.subplots(dpi=300,figsize=[6,5])

# arrays for summing total
# diam_005=np.empty([1,2],dtype=float) #Column 0 = diameters Column 1 = weights
# # print(diam_005)
# # diam_005[0,0]=23
# # print(diam_005)
# # diam_005[0,1]=46
# # print(diam_005)
# # diam=np.array([3,4,5,6])
# # print(diam)
# # weigths=np.array([0.9,0.5,16.2,87.4])
# # print(weigths)
# # diam_005=np.append(diam_005,[[3,4,5,6],[0.9,0.5,16.2,87.4]],axis = 1)
# # print(diam_005)

# # arr = np.array([0]) 
# # print(arr)
# # result = np.append(arr,[[5],[7],[9]],axis = 1)
# # print(result)


# filenames = ['lig_02.txt', 'lig_05.txt', 'lig_08.txt']
# with open('/Users/josephgiliberto/Builds/NGA2/examples/ligament/DONTDELETE_pdf_data/out.txt', 'w') as outfile:
#     for fname in filenames:
#         with open(fname) as infile:
#             for line in infile:
#                 outfile.write(line)

# filenames = ['lig_03.txt', 'lig_06.txt', 'lig_09.txt']
# with open('/Users/josephgiliberto/Builds/NGA2/examples/ligament/DONTDELETE_pdf_data/Oh05.txt', 'w') as outfile:
#     for fname in filenames:
#         with open(fname) as infile:
#             for line in infile:
#                 outfile.write(line)

#Oh=0.05
diam=0.0; weights=0.0; vals1=0.0; edges1=0.0; edgemed1=0.0
diam,u,v,w,vel_total = np.genfromtxt('Oh005.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
weights = np.pi/6.0*diam**3.0 # Weight by volume
vals1,edges1 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
edgemed1 = (edges1[:-1]+edges1[1:])/2
ax5.plot(edgemed1,vals1,'k-x',linewidth=3.0,label=r'$Oh_l=0.05$')
#Oh=0.5
diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
diam,u,v,w,vel_total = np.genfromtxt('Oh05.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
weights = np.pi/6.0*diam**3.0 # Weight by volume
vals2,edges2 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
edgemed2 = (edges2[:-1]+edges2[1:])/2
ax5.plot(edgemed2,vals2,'r-o',linewidth=3.0,label=r'$Oh_l=0.5$')

# Figure 3: Ratio comparison
ax5.set_xlabel(dimless_x_label,fontsize='20')
ax5.set_ylabel(r'volume weighted PDF',fontsize='20')
ax5.legend(fontsize='20')
ax5.set_xticks([0,0.2,0.4,0.6,0.8]) 
ax5.set_xticklabels([0,0.2,0.4,0.6,0.8], fontsize=20)
ax5.set_yticks([0,1,2,3,4,5,6]) 
ax5.set_yticklabels([0,1,2,3,4,5,6], fontsize=20)
ax5.set_xlim(0,0.8)
ax5.set_ylim(0,6)
fig5.savefig(plot_name+'_viscr_average.png')



# #02_lig (viscr=50,Case 1)
# diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
# diam,u,v,w,vel_total = np.genfromtxt('lig_02.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
# weights = np.pi/6.0*diam**3.0 # Weight by volume
# vals1,edges1 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
# edgemed1 = (edges1[:-1]+edges1[1:])/2
# ax1.plot(edgemed1,vals1,'-x',linewidth=2.0,label='Case 1')
# ax3.plot(edgemed1,vals1,'k-x',linewidth=2.0,label=r'$Oh_l=0.05$')
# #05_lig (viscr=50, Case 2)
# diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
# diam,u,v,w,vel_total = np.genfromtxt('lig_05.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
# # appending array

# weights = np.pi/6.0*diam**3.0 # Weight by volume
# vals1,edges1 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
# edgemed1 = (edges1[:-1]+edges1[1:])/2
# ax1.plot(edgemed1,vals1,'-x',linewidth=2.0,label='Case 2')
# ax3.plot(edgemed1,vals1,'k-x',linewidth=2.0)
# #08_lig (viscr=50, Case 3)
# diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
# diam,u,v,w,vel_total = np.genfromtxt('lig_08.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
# weights = np.pi/6.0*diam**3.0 # Weight by volume
# vals1,edges1 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
# edgemed1 = (edges1[:-1]+edges1[1:])/2
# ax1.plot(edgemed1,vals1,'-x',linewidth=2.0,label='Case 3')
# ax3.plot(edgemed1,vals1,'k-x',linewidth=2.0)

# # Figure 1: visc_r=50 comparison
# ax1.set_xlabel(dimless_x_label)
# ax1.set_ylabel(r'volume weighted PDF')
# ax1.legend()
# fig1.suptitle(r'$\mu_l/\mu_g=50$',fontsize=20)
# fig1.savefig(plot_name+'_viscr_50.png') 

# #####

# #03_lig (viscr=500,Case 1)
# diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
# diam,u,v,w,vel_total = np.genfromtxt('lig_03.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
# weights = np.pi/6.0*diam**3.0 # Weight by volume
# vals2,edges2 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
# edgemed2 = (edges2[:-1]+edges2[1:])/2
# ax2.plot(edgemed2,vals2,'-x',linewidth=2.0,label='Case 1')
# ax3.plot(edgemed2,vals2,'r-o',linewidth=2.0,label=r'$Oh_l=0.5$')
# #06_lig (viscr=500, Case 2)
# diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
# diam,u,v,w,vel_total = np.genfromtxt('lig_06.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
# weights = np.pi/6.0*diam**3.0 # Weight by volume
# vals2,edges2 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
# edgemed2 = (edges2[:-1]+edges2[1:])/2
# ax2.plot(edgemed2,vals2,'-x',linewidth=2.0,label='Case 2')
# ax3.plot(edgemed2,vals2,'r-o',linewidth=2.0)
# #09_lig (viscr=500, Case 3)
# diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
# diam,u,v,w,vel_total = np.genfromtxt('lig_09.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
# weights = np.pi/6.0*diam**3.0 # Weight by volume
# vals2,edges2 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
# edgemed2 = (edges2[:-1]+edges2[1:])/2
# ax2.plot(edgemed2,vals2,'-x',linewidth=2.0,label='Case 3')
# ax3.plot(edgemed2,vals2,'r-o',linewidth=2.0)
# ax4.plot(edgemed2,vals2,'k-o',linewidth=3.0,label=r'$n=1.0$')

# # Figure 2: visc_r=50 comparison
# ax2.set_xlabel(dimless_x_label)
# ax2.set_ylabel(r'volume weighted PDF')
# ax2.legend()
# fig2.suptitle(r'$\mu_l/\mu_g=500$',fontsize=20)
# fig2.savefig(plot_name+'_viscr_500.png') 

# # Figure 3: Ratio comparison
# ax3.set_xlabel(dimless_x_label,fontsize='20')
# ax3.set_ylabel(r'volume weighted PDF',fontsize='20')
# ax3.legend(fontsize='20')
# ax3.set_xticks([0,0.2,0.4,0.6,0.8]) 
# ax3.set_xticklabels([0,0.2,0.4,0.6,0.8], fontsize=20)
# ax3.set_yticks([0,1,2,3,4,5,6]) 
# ax3.set_yticklabels([0,1,2,3,4,5,6], fontsize=20)
# ax3.set_xlim(0,0.8)
# ax3.set_ylim(0,6)
# fig3.savefig(plot_name+'_viscr_comparison.png')

# ### Shear thinning data

# #11_lig (viscr=500,n=0.8)
# diam=0.0; weights=0.0; vals4=0.0; edges4=0.0; edgemed4=0.0
# diam,u,v,w,vel_total = np.genfromtxt('lig_11_st.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
# weights = np.pi/6.0*diam**3.0 # Weight by volume
# vals4,edges4 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
# edgemed4 = (edges4[:-1]+edges4[1:])/2
# ax4.plot(edgemed4,vals4,'r-o',linewidth=3.0,label=r'$n=0.8$')

# #12_lig (viscr=500,n=0.5)
# diam=0.0; weights=0.0; vals4=0.0; edges4=0.0; edgemed4=0.0
# diam,u,v,w,vel_total = np.genfromtxt('lig_12_st.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
# weights = np.pi/6.0*diam**3.0 # Weight by volume
# vals4,edges4 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
# edgemed4 = (edges4[:-1]+edges4[1:])/2
# ax4.plot(edgemed4,vals4,'b-o',linewidth=3.0,label=r'$n=0.5$')

# #13_lig (viscr=500,n=0.1)
# diam=0.0; weights=0.0; vals4=0.0; edges4=0.0; edgemed4=0.0
# diam,u,v,w,vel_total = np.genfromtxt('lig_13_st.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
# weights = np.pi/6.0*diam**3.0 # Weight by volume
# vals4,edges4 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
# edgemed4 = (edges4[:-1]+edges4[1:])/2
# ax4.plot(edgemed4,vals4,'g-o',linewidth=3.0,label=r'$n=0.1$')

# #14_lig (viscr=500,n=0.01)
# diam=0.0; weights=0.0; vals4=0.0; edges4=0.0; edgemed4=0.0
# diam,u,v,w,vel_total = np.genfromtxt('lig_14_st.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
# weights = np.pi/6.0*diam**3.0 # Weight by volume
# vals4,edges4 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
# edgemed4 = (edges4[:-1]+edges4[1:])/2
# ax4.plot(edgemed4,vals4,'c-o',linewidth=3.0,label=r'$n=0.01$')

# # Figure 4: Shear-thinning comparison
# ax4.set_xlabel(dimless_x_label,fontsize='20')
# ax4.set_ylabel(r'volume weighted PDF',fontsize='20')
# ax4.set_xticks([0,0.2,0.4,0.6,0.8]) 
# ax4.set_xticklabels([0,0.2,0.4,0.6,0.8], fontsize=20)
# ax4.set_yticks([0,1,2,3,4,5,6]) 
# ax4.set_yticklabels([0,1,2,3,4,5,6], fontsize=20)    
# ax4.legend(fontsize='20')
# ax4.set_xlim(0,0.8)
# ax4.set_ylim(0,6)
# fig4.savefig(plot_name+'_n_comparison.png')
