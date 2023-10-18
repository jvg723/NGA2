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

plot_name = "npdf"

dimless_x_label = r'$d/D$'

# Parameters
lz=20 # ligament length in z direction
D=1.0 # initial ligament diameter
initial_volume=np.pi*(D/2.0)**2.0*lz
logbins=np.logspace(np.log10(0.01),np.log10(1.0), num=30, base=10.0) # binning

# Start plot
fig2, ax2 = plt.subplots(dpi=300,figsize=[6,5]) 

os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/ligament/pdf_data')
#02_lig ( Case 1)
diam,u,v,w,vel_total = np.genfromtxt('lig_02.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)
weights = np.pi/6.0*diam**3.0 # Weight by volume
vals2,edges2 = np.histogram(diam/D,bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
edgemed2 = (edges2[:-1]+edges2[1:])/2
ax2.plot(edgemed2,vals2*np.sum(weights)/initial_volume,'-x',linewidth=2.0,label='Case 1')

# Figure 1: vpdf from film, rim, and nodes

ax2.set_xlabel(dimless_x_label)
ax2.set_ylabel(r'volume weighted probability')
ax2.legend()
fig2.suptitle(r'$\mu_l/\mu_g=50$',fontsize=20)
fig2.savefig(plot_name+'_weighted_f10.png') 

