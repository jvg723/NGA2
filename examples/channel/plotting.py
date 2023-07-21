# Importing packages
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import re
import pandas as pd


# Set plot properies
plt.rc('text', usetex=True)

#Beta=1.0 files
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_1/velocity')
df_b1_vel=pd.read_csv('Uavg_1.00000E+02', delim_whitespace=True, header=None, skiprows=1, usecols=[0,1], names=['Height','Uavg'])
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_1/plots')
df_b1_thy=pd.read_csv('theory', delim_whitespace=True, header=None, skiprows=1, usecols=[0,1,2,3], names=['Height','u','Txx','Txy'])
# Beta=0.8 files
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_08/velocity')
df_b08_vel=pd.read_csv('Uavg_1.00000E+02', delim_whitespace=True, header=None, skiprows=1, usecols=[0,1], names=['Height','Uavg'])
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_08/stress')
df_b08_str=pd.read_csv('Cavg_1.00000E+02', delim_whitespace=True, header=None, skiprows=1, usecols=[0,4,5], names=['Height','Txx_avg','Txy_avg'])
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_08/plots')
df_b08_thy=pd.read_csv('theory', delim_whitespace=True, header=None, skiprows=1, usecols=[0,1,2,3], names=['Height','u','Txx','Txy'])
# Beta=0.6 files
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_06/velocity')
df_b06_vel=pd.read_csv('Uavg_1.00016E+02', delim_whitespace=True, header=None, skiprows=1, usecols=[0,1], names=['Height','Uavg'])
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_06/stress')
df_b06_str=pd.read_csv('Cavg_1.00016E+02', delim_whitespace=True, header=None, skiprows=1, usecols=[0,4,5], names=['Height','Txx_avg','Txy_avg'])
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_06/plots')
df_b06_thy=pd.read_csv('theory', delim_whitespace=True, header=None, skiprows=1, usecols=[0,1,2,3], names=['Height','u','Txx','Txy'])
# Beta=0.4 files
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_04/velocity')
df_b04_vel=pd.read_csv('Uavg_1.00000E+02', delim_whitespace=True, header=None, skiprows=1, usecols=[0,1], names=['Height','Uavg'])
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_04/stress')
df_b04_str=pd.read_csv('Cavg_1.00000E+02', delim_whitespace=True, header=None, skiprows=1, usecols=[0,4,5], names=['Height','Txx_avg','Txy_avg'])
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_04/plots')
df_b04_thy=pd.read_csv('theory', delim_whitespace=True, header=None, skiprows=1, usecols=[0,1,2,3], names=['Height','u','Txx','Txy'])
# Beta=0.2 files
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_02/velocity')
df_b02_vel=pd.read_csv('Uavg_1.00000E+02', delim_whitespace=True, header=None, skiprows=1, usecols=[0,1], names=['Height','Uavg'])
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_02/stress')
df_b02_str=pd.read_csv('Cavg_1.00000E+02', delim_whitespace=True, header=None, skiprows=1, usecols=[0,4,5], names=['Height','Txx_avg','Txy_avg'])
os.chdir('/Users/josephgiliberto/Builds/NGA2/examples/channel/cases/beta_02/plots')
df_b02_thy=pd.read_csv('theory', delim_whitespace=True, header=None, skiprows=1, usecols=[0,1,2,3], names=['Height','u','Txx','Txy'])

#Normalization parameters
H=1.0
Uc=1.5
mu_s=0.01
non_dim_stress=(mu_s*Uc/H)

# Subplot configuration and size
fig,axs=plt.subplots(1,3,figsize=(20,8))

# Plotting velocity
ax1=plt.subplot(131)
ax1.plot(df_b1_thy['u']/Uc,df_b1_thy['Height']/H,linewidth='5',color='black')
ax1.scatter(df_b1_vel['Uavg']/Uc,df_b1_vel['Height']/H,marker="o",s=80,color='black')
ax1.plot(df_b08_thy['u']/Uc,df_b08_thy['Height']/H,linewidth='5',color='blue')
ax1.scatter(df_b08_vel['Uavg']/Uc,df_b08_vel['Height']/H,marker="o",s=80,color='blue')
ax1.plot(df_b06_thy['u']/Uc,df_b06_thy['Height']/H,linewidth='5',color='red')
ax1.scatter(df_b06_vel['Uavg']/Uc,df_b06_vel['Height']/H,marker="o",s=80,color='red')
ax1.plot(df_b04_thy['u']/Uc,df_b04_thy['Height']/H,linewidth='5',color='green')
ax1.scatter(df_b04_vel['Uavg']/Uc,df_b04_vel['Height']/H,marker="o",s=80,color='green')
ax1.plot(df_b02_thy['u']/Uc,df_b02_thy['Height']/H,linewidth='5',color='orange')
ax1.scatter(df_b02_vel['Uavg']/Uc,df_b02_vel['Height']/H,marker="o",s=80,color='orange')
# Setting velocity labels
ax1.set_xlabel(r'$u/U_c$',fontsize='30')
ax1.set_ylabel(r'$y/H$',fontsize='30')
# Setting velocity limits
ax1.set_xlim(0.0,3.0)
ax1.set_ylim(-0.5,0.5)
# # Setting velocity tick locations
# ax1.set_xticks(np.arange(0.0,np.ceil(np.max(u))))
# ax1.set_yticks(np.arange(np.min(y),np.max(y)+H/2.0,H/2.0))
# Setting velocity tick label size
ax1.tick_params(axis='x',labelsize=20)
ax1.tick_params(axis='y',labelsize=20)

# Plotting normal stress
ax2=plt.subplot(132)
ax2.plot(df_b08_thy['Txx']/non_dim_stress,df_b08_thy['Height']/H,linewidth='5',color='blue')
ax2.scatter(df_b08_str['Txx_avg']/non_dim_stress,df_b08_str['Height']/H,marker="o",s=80,color='blue')
ax2.plot(df_b06_thy['Txx']/non_dim_stress,df_b06_thy['Height']/H,linewidth='5',color='red')
ax2.scatter(df_b06_str['Txx_avg']/non_dim_stress,df_b06_str['Height']/H,marker="o",s=80,color='red')
ax2.plot(df_b04_thy['Txx']/non_dim_stress,df_b04_thy['Height']/H,linewidth='5',color='green')
ax2.scatter(df_b04_str['Txx_avg']/non_dim_stress,df_b04_str['Height']/H,marker="o",s=80,color='green')
ax2.plot(df_b02_thy['Txx']/non_dim_stress,df_b02_thy['Height']/H,linewidth='5',color='orange')
ax2.scatter(df_b02_str['Txx_avg']/non_dim_stress,df_b02_str['Height']/H,marker="o",s=80,color='orange')
# Setting normal stress label
ax2.set_xlabel(r'$\tau_{xx}/(\mu_sU_c/H)$',fontsize='30')
# # Setting normal stress limits
ax2.set_xlim(0.0,300.0)
ax2.set_ylim(-0.5,0.5)
# Setting normal stress tick locations
ax2.set_yticks([])
# Setting velocity tick label size
ax2.tick_params(axis='x',labelsize=20)
ax2.tick_params(axis='y',labelsize=20)

# Plotting shear stress
ax3=plt.subplot(133)
ax3.plot(df_b08_thy['Txy']/non_dim_stress,df_b08_thy['Height']/H,linewidth='5',color='blue')
ax3.scatter(df_b08_str['Txy_avg']/non_dim_stress,df_b08_str['Height']/H,marker="o",s=80,color='blue')
ax3.plot(df_b06_thy['Txy']/non_dim_stress,df_b06_thy['Height']/H,linewidth='5',color='red')
ax3.scatter(df_b06_str['Txy_avg']/non_dim_stress,df_b06_str['Height']/H,marker="o",s=80,color='red')
ax3.plot(df_b04_thy['Txy']/non_dim_stress,df_b04_thy['Height']/H,linewidth='5',color='green')
ax3.scatter(df_b04_str['Txy_avg']/non_dim_stress,df_b04_str['Height']/H,marker="o",s=80,color='green')
ax3.plot(df_b02_thy['Txy']/non_dim_stress,df_b02_thy['Height']/H,linewidth='5',color='orange')
ax3.scatter(df_b02_str['Txy_avg']/non_dim_stress,df_b02_str['Height']/H,marker="o",s=80,color='orange')
# Setting shear stress label
ax3.set_xlabel(r'$\tau_{xy}/(\mu_sU_c/H)$',fontsize='30')
# # Setting shear stress limits
ax3.set_xlim(-10.0,10.0)
ax3.set_ylim(-0.5,0.5)
# Setting normal stress tick locations
ax3.set_yticks([])
# Setting velocity tick label size
ax3.tick_params(axis='x',labelsize=20)
ax3.tick_params(axis='y',labelsize=20)

# Show the plot
plt.show()




