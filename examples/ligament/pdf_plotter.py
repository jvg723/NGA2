# Importing packages
import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import re
import pandas as pd


# base_dir=os.getcwd() 

# ## Define functions

# # Log-normal curve
# def log_norm(x,mu,sigma):
#     return 1.0/(x*sigma*np.sqrt(2.0*np.pi))*np.exp(-(np.log(x)-mu)**2/(2*sigma**2))

# # Ligament diameter
# Diam=1.0




#02_lig files
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/ligament/pdf_data')
diam,u,v,w,vel_total = np.loadtxt('lig_02.txt',skiprows=1,usecols=(0,1,2,3,4),unpack=True)
# Weight by volume
weights = np.pi/6.0*diam**3.0
logbins=np.logspace(1.0,11.0,num=30,base=2.0)
vals3,edges3 = np.histogram(diam, bins=logbins,range=(0,270),weights=weights,density=True)  # arguments are passed to np.histogram
edgemed3 = (edges3[:-1]+edges3[1:])/2
fig3, ax3 = plt.subplots(dpi=300,figsize=[6,3])
ax3.plot(edgemed3,vals3,'-x',linewidth=1.0)

# Show the plot
plt.show()

# df_02_lig=pd.read_csv('droplets.002262', delim_whitespace=True, header=None, skiprows=1, usecols=[0], names=['Diameter'])
# mean=np.mean(np.log(df_02_lig['Diameter']))                                                                       # mean of ln(df[Diameter])   
# std=np.std(np.log(df_02_lig['Diameter']))                                                                         # standard deviation of ln(df[Diameter])                                                                     
# df_02_lig['PDF']=log_norm(df_02_lig['Diameter'],mean,std) 

# #02_lig files
# os.chdir('/home/fs01/jvg36/Builds/NGA2/examples/ligament/cases/02_lig/spray')
# df_02_lig=pd.read_csv('droplets.002262', delim_whitespace=True, header=None, skiprows=1, usecols=[0], names=['Diameter'])
# mean=np.mean(np.log(df_02_lig['Diameter']))                                                                       # mean of ln(df[Diameter])   
# std=np.std(np.log(df_02_lig['Diameter']))                                                                         # standard deviation of ln(df[Diameter])                                                                     
# df_02_lig['PDF']=log_norm(df_02_lig['Diameter'],mean,std) 

# # resest
# mean=0.0; std=0.0

# #03_lig files
# os.chdir('/home/fs01/jvg36/Builds/NGA2/examples/ligament/cases/03_lig/spray')
# df_03_lig=pd.read_csv('droplets.002700', delim_whitespace=True, header=None, skiprows=1, usecols=[0], names=['Diameter'])
# mean=np.mean(np.log(df_03_lig['Diameter']))                                                                       # mean of ln(df[Diameter])   
# std=np.std(np.log(df_03_lig['Diameter']))                                                                         # standard deviation of ln(df[Diameter])                                                                     
# df_03_lig['PDF']=log_norm(df_03_lig['Diameter'],mean,std) 

# # resest
# mean=0.0; std=0.0

# #04_lig files
# os.chdir('/home/fs01/jvg36/Builds/NGA2/examples/ligament/cases/04_lig/spray')
# df_04_lig=pd.read_csv('droplets.016682', delim_whitespace=True, header=None, skiprows=1, usecols=[0], names=['Diameter'])
# mean=np.mean(np.log(df_04_lig['Diameter']))                                                                       # mean of ln(df[Diameter])   
# std=np.std(np.log(df_04_lig['Diameter']))                                                                         # standard deviation of ln(df[Diameter])                                                                     
# df_04_lig['PDF']=log_norm(df_04_lig['Diameter'],mean,std) 

# # resest
# mean=0.0; std=0.0

# #08_lig files
# os.chdir('/home/fs01/jvg36/Builds/NGA2/examples/ligament/cases/08_lig/spray')
# df_08_lig=pd.read_csv('droplets.002290', delim_whitespace=True, header=None, skiprows=1, usecols=[0], names=['Diameter'])
# mean=np.mean(np.log(df_08_lig['Diameter']))                                                                       # mean of ln(df[Diameter])   
# std=np.std(np.log(df_08_lig['Diameter']))                                                                         # standard deviation of ln(df[Diameter])                                                                     
# df_08_lig['PDF']=log_norm(df_08_lig['Diameter'],mean,std) 

# # resest
# mean=0.0; std=0.0

# #11_lig files
# os.chdir('/home/fs01/jvg36/Builds/NGA2/examples/ligament/cases/11_lig/spray')
# df_11_lig=pd.read_csv('droplets.015380', delim_whitespace=True, header=None, skiprows=1, usecols=[0], names=['Diameter'])
# mean=np.mean(np.log(df_11_lig['Diameter']))                                                                       # mean of ln(df[Diameter])   
# std=np.std(np.log(df_11_lig['Diameter']))                                                                         # standard deviation of ln(df[Diameter])                                                                     
# df_11_lig['PDF']=log_norm(df_11_lig['Diameter'],mean,std) 

# # resest
# mean=0.0; std=0.0

# #12_lig files
# os.chdir('/home/fs01/jvg36/Builds/NGA2/examples/ligament/cases/12_lig/spray')
# df_12_lig=pd.read_csv('droplets.017861', delim_whitespace=True, header=None, skiprows=1, usecols=[0], names=['Diameter'])
# mean=np.mean(np.log(df_12_lig['Diameter']))                                                                       # mean of ln(df[Diameter])   
# std=np.std(np.log(df_12_lig['Diameter']))                                                                         # standard deviation of ln(df[Diameter])                                                                     
# df_12_lig['PDF']=log_norm(df_12_lig['Diameter'],mean,std) 

# # resest
# mean=0.0; std=0.0

# #13_lig files
# os.chdir('/home/fs01/jvg36/Builds/NGA2/examples/ligament/cases/13_lig/spray')
# df_13_lig=pd.read_csv('droplets.015364', delim_whitespace=True, header=None, skiprows=1, usecols=[0], names=['Diameter'])
# mean=np.mean(np.log(df_13_lig['Diameter']))                                                                       # mean of ln(df[Diameter])   
# std=np.std(np.log(df_13_lig['Diameter']))                                                                         # standard deviation of ln(df[Diameter])                                                                     
# df_13_lig['PDF']=log_norm(df_13_lig['Diameter'],mean,std) 

# #15_lig files
# os.chdir('/home/fs01/jvg36/Builds/NGA2/examples/ligament/cases/15_lig/spray')
# df_13_lig=pd.read_csv('droplets.017289', delim_whitespace=True, header=None, skiprows=1, usecols=[0], names=['Diameter'])

# os.chdir(base_dir)          

# # Plot fonts
# mpl.rcParams['text.usetex'] = False
# mpl.rcParams['font.family'] = 'DejaVu Sans'
# mpl.rcParams["mathtext.fontset"]
# mpl.rcParams["axes.edgecolor"] = "black"
# mpl.rcParams["axes.linewidth"] = 2.50

# # Plot PDF comparing impact of viscosity on drop size 
# fig1,ax1=plt.subplots(figsize=(30, 25))
# ax1.scatter(df_02_lig['Diameter']/Diam, df_02_lig['PDF']*Diam,marker="o",s=600,color='black',label=r'$\mu_l/\mu_g=50, We=10$')
# ax1.scatter(df_08_lig['Diameter']/Diam, df_08_lig['PDF']*Diam,marker="o",s=600,color='green',label=r'$\mu_l/\mu_g=50, We=20$')
# ax1.scatter(df_04_lig['Diameter']/Diam, df_04_lig['PDF']*Diam,marker="o",s=600,color='red',  label=r'$\mu_l/\mu_g=500, We=10$')

# ax1.set_xscale('log')
# ax1.set_yscale('log')
# ax1.set_xlabel(r'$d/D$',fontsize='50')
# ax1.set_xticks(1e-4,1)
# ax1.set_ylabel(r'$PDF$',fontsize='50')
# ax1.tick_params(axis='x',labelsize=60)
# ax1.tick_params(axis='y',labelsize=60)
# ax1.legend(loc='lower left',fontsize='50')
# # ax.set_ylabel('\\textit{Velocity (\N{DEGREE SIGN}/sec)}', fontsize=16)
# fig1.savefig('/home/fs01/jvg36/Builds/NGA2/examples/ligament/visc_ratio.png')

# Plot PDF comparing impact of shear thinining on drop size 
# os.chdir(base_dir)          # move back to base director
# # Set plot properies
# # plt.rcParams['text.usetex'] = True

# # fig, ax = plt.subplots(figsize=(6, 4))

# # plt.scatter(df_02_lig['Diameter']/Diam, df_02_lig['PDF']*Diam,marker="o",s=80,color='black')
# plt.scatter(df_04_lig['Diameter']/Diam, df_04_lig['PDF']*Diam,marker="o",s=80,color='red')
# # plt.scatter(df_08_lig['Diameter']/Diam, df_08_lig['PDF']*Diam,marker="o",s=80,color='green')
# plt.scatter(df_11_lig['Diameter']/Diam, df_11_lig['PDF']*Diam,marker="o",s=80,color='blue')
# plt.scatter(df_13_lig['Diameter']/Diam, df_13_lig['PDF']*Diam,marker="o",s=80,color='orange')
# plt.xscale('log')
# plt.yscale('log')
# # ax.set_xlabel(r'\textbf{d/D}')
# # ax.set_ylabel('\\textit{Velocity (\N{DEGREE SIGN}/sec)}', fontsize=16)
# plt.savefig('/home/fs01/jvg36/Builds/NGA2/examples/ligament/ST_comparsion.png')







