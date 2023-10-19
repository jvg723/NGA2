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

# Start plot
fig1, ax1 = plt.subplots(dpi=300,figsize=[6,5]) 
fig2, ax2 = plt.subplots(dpi=300,figsize=[6,5]) 
fig3, ax3 = plt.subplots(dpi=300,figsize=[6,5]) 

os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/ligament/pdf_data')
#02_lig (viscr=50,Case 1)
diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
diam,u,v,w,vel_total = np.genfromtxt('lig_02.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
weights = np.pi/6.0*diam**3.0 # Weight by volume
vals1,edges1 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
edgemed1 = (edges1[:-1]+edges1[1:])/2
ax1.plot(edgemed1,vals1*np.sum(weights)/initial_volume,'-x',linewidth=2.0,label='Case 1')
ax3.plot(edgemed1,vals1*np.sum(weights)/initial_volume,'k-x',linewidth=2.0,label=r'$\mu_l/\mu_g=50$')
#05_lig (viscr=50, Case 2)
diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
diam,u,v,w,vel_total = np.genfromtxt('lig_05.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
weights = np.pi/6.0*diam**3.0 # Weight by volume
vals1,edges1 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
edgemed1 = (edges1[:-1]+edges1[1:])/2
ax1.plot(edgemed1,vals1*np.sum(weights)/initial_volume,'-x',linewidth=2.0,label='Case 2')
ax3.plot(edgemed1,vals1*np.sum(weights)/initial_volume,'k-x',linewidth=2.0)
#05_lig (viscr=50, Case 3)
diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
diam,u,v,w,vel_total = np.genfromtxt('lig_08.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
weights = np.pi/6.0*diam**3.0 # Weight by volume
vals1,edges1 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
edgemed1 = (edges1[:-1]+edges1[1:])/2
ax1.plot(edgemed1,vals1*np.sum(weights)/initial_volume,'-x',linewidth=2.0,label='Case 3')
ax3.plot(edgemed1,vals1*np.sum(weights)/initial_volume,'k-x',linewidth=2.0)

# Figure 1: visc_r=50 comparison
ax1.set_xlabel(dimless_x_label)
ax1.set_ylabel(r'volume weighted probability')
ax1.legend()
fig1.suptitle(r'$\mu_l/\mu_g=50$',fontsize=20)
fig1.savefig(plot_name+'_viscr_50.png') 

#####

#03_lig (viscr=500,Case 1)
diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
diam,u,v,w,vel_total = np.genfromtxt('lig_03.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
weights = np.pi/6.0*diam**3.0 # Weight by volume
vals2,edges2 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
edgemed2 = (edges2[:-1]+edges2[1:])/2
ax2.plot(edgemed2,vals2*np.sum(weights)/initial_volume,'-x',linewidth=2.0,label='Case 1')
ax3.plot(edgemed2,vals2*np.sum(weights)/initial_volume,'r-o',linewidth=2.0,label=r'$\mu_l/\mu_g=500$')
#06_lig (viscr=500, Case 2)
diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
diam,u,v,w,vel_total = np.genfromtxt('lig_06.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
weights = np.pi/6.0*diam**3.0 # Weight by volume
vals2,edges2 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
edgemed2 = (edges2[:-1]+edges2[1:])/2
ax2.plot(edgemed2,vals2*np.sum(weights)/initial_volume,'-x',linewidth=2.0,label='Case 2')
ax3.plot(edgemed2,vals2*np.sum(weights)/initial_volume,'r-o',linewidth=2.0)
#09_lig (viscr=500, Case 3)
diam=0.0; weights=0.0; vals2=0.0; edges2=0.0; edgemed2=0.0
diam,u,v,w,vel_total = np.genfromtxt('lig_09.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
weights = np.pi/6.0*diam**3.0 # Weight by volume
vals2,edges2 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
edgemed2 = (edges2[:-1]+edges2[1:])/2
ax2.plot(edgemed2,vals2*np.sum(weights)/initial_volume,'-x',linewidth=2.0,label='Case 3')
ax3.plot(edgemed2,vals2*np.sum(weights)/initial_volume,'r-o',linewidth=2.0)

# Figure 2: visc_r=50 comparison
ax2.set_xlabel(dimless_x_label)
ax2.set_ylabel(r'volume weighted probability')
ax2.legend()
fig2.suptitle(r'$\mu_l/\mu_g=500$',fontsize=20)
fig2.savefig(plot_name+'_viscr_500.png') 

# Figure 3: Ratio comparison
ax3.set_xlabel(dimless_x_label)
ax3.set_ylabel(r'volume weighted probability')
ax3.legend()
fig3.suptitle(r'$\mu_l/\mu_g-comparison$',fontsize=20)
fig3.savefig(plot_name+'_viscr_comparison.png') 
