# Importing packages
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import re

## Functions ##
# F values at a given x position
def fvalues(x,A,C):
    # F^+(x)
    Fplus=(C*x+math.sqrt(A**3.0+(C*x)**2.0))**(1.0/3.0)
    # F^-(x)
    Fminus=-abs(C*x-math.sqrt(A**3.0+(C*x)**2.0))**(1.0/3.0)
    # Pass out values
    return Fplus, Fminus

# G values at a given x position
def gvalues(x,A,C):
    # G^+(x)
    Gplus=3.0*C*x+math.sqrt(A**3.0+(C*x)**2.0)
    # G^-(x)
    Gminus=3.0*C*x-math.sqrt(A**3.0+(C*x)**2.0)
    # Pass out values
    return Gplus, Gminus 

# Read parameters from input
def search_str(file_path,word):
    with open(file_path, 'r') as fp:
        # read all lines in a list
        lines=fp.readlines()
        # pattern for numbers
        pattern="[0-9]+.?[0-9]*"
        for line in lines:
            # check if string present on a current line
            if line.find(word) != -1:
                param=(re.findall(pattern,line))
                variable=float(param[0])
    # Pass out values
    return variable 
                
# Path to input file
path=os.getcwd()+'/input'

# Read input
H     =search_str(path,'Ly')                                 # Channel height
ny    =search_str(path,'ny')                                 # Number of grid points
visc_s=search_str(path,'Solvent dynamic viscosity')          # Solvent viscosity
rho   =search_str(path,'Density')                            # Fluid density
Uc    =search_str(path,'Ubulk')                              # Average channel flow
L     =search_str(path,'Maximum polymer extensibility')      # Maximum extension of polymer
lam   =search_str(path,'Polymer relaxation time')            # Polymer relaxation time
visc_p=search_str(path,'Polymer viscosity')                  # Polymer viscosity 

# Generate grid
y=np.linspace(-H/2.0,H/2.0,int(ny)+1)

# Charateristic length scale
Lc=H

# Flow rate
Q=Uc*H

# Polymer terms
b=L**2.0-3.0
epsilon=1.0/(3.0*(b+5.0))
lam=((b+2.0)/(b+5.0))*lam

# Fluid viscosity
beta=visc_s/(visc_s+visc_p)
visc  =visc_s+visc_p
print('beta=',beta)

# Pressure gradient
px=-12.0*(visc/H**3.0)*Q
print('px=',px)

# Constant coefficents 
A=(visc_p**2.0/(6.0*epsilon*lam**2.0))*(1.0+visc_p/visc_s)
C=(visc_p**2.0/(4.0*epsilon*lam**2.0))*(visc_p/visc_s)*px

# Initalize arrays for stress and velocity
Txx=np.zeros(len(y))
Txy=np.zeros(len(y))
u  =np.zeros(len(y))

# Loop over channel length to calculate velocity and stress
for i in range(0,len(y)):
    # Posistion varrying coefficent
    B=C*y[i]
    # Shear stress
    Txy[i]=(B+math.sqrt((A**3.0)+(B**2.0)))**(1.0/3.0)-abs(B-math.sqrt((A**3.0)+(B**2.0)))**(1.0/3.0)
    # Normal stress
    Txx[i]=2.0*(lam/visc_p)*Txy[i]**2.0
    # FENE-P curve 
    if y[i]<0.0:
        # @ -H/2
        FH=fvalues(-H/2.0,A,C);Fp_H=FH[0];Fm_H=FH[1]
        GH=gvalues(-H/2.0,A,C);Gp_H=GH[0];Gm_H=GH[1]
    elif y[i]>=0.0:
    # @ H/2
        FH=fvalues(H/2.0,A,C);Fp_H=FH[0];Fm_H=FH[1]
        GH=gvalues(H/2.0,A,C);Gp_H=GH[0];Gm_H=GH[1]
    # @ y[i]
    Fy=fvalues(y[i],A,C);Fp_y=Fy[0];Fm_y=Fy[1]
    Gy=gvalues(y[i],A,C);Gp_y=Gy[0];Gm_y=Gy[1]
    # Velocity
    if beta==1.0:
        u[i]=(-(px*(H/2.0)**2.0)/(2.0*visc_s))*(1.0-(2.0*y[i]/H)**2.0)
    else:
        u[i]=(-(px*(H/2.0)**2.0)/(2.0*visc_s))*(1.0-(2.0*y[i]/H)**2.0)+(3.0/(8.0*C*visc_s))*(Fp_H*Gm_H-Fp_y*Gm_y+Fm_H*Gp_H-Fm_y*Gp_y)

# Set plot properies
plt.rc('text', usetex=True)

# Subplot configuration and size
fig,axs=plt.subplots(1,3,figsize=(20,8))

# Plotting velocity
ax1=plt.subplot(131)
ax1.plot(u,y,linewidth='5',color='black')
# Setting velocity labels
ax1.set_xlabel(r'$u$',fontsize='30')
ax1.set_ylabel(r'$y$',fontsize='30')
# Setting velocity limits
ax1.set_xlim(left=0.0)
ax1.set_ylim([np.min(y),np.max(y)])
# Setting velocity tick locations
ax1.set_xticks(np.arange(0.0,np.ceil(np.max(u))))
ax1.set_yticks(np.arange(np.min(y),np.max(y)+H/2.0,H/2.0))
# Setting velocity tick label size
ax1.tick_params(axis='x',labelsize=20)
ax1.tick_params(axis='y',labelsize=20)

# Plotting normal stress
ax2=plt.subplot(132)
ax2.plot(Txx,y,linewidth='5',color='black')
# Setting normal stress label
ax2.set_xlabel(r'$\tau_{xx}$',fontsize='30')
# Setting normal stress limits
ax2.set_xlim(left=0.0)
ax2.set_ylim([np.min(y),np.max(y)])
# Setting normal stress tick locations
ax2.set_yticks([])
# Setting velocity tick label size
ax2.tick_params(axis='x',labelsize=20)
ax2.tick_params(axis='y',labelsize=20)

# Plotting shear stress
ax3=plt.subplot(133)
ax3.plot(Txy,y,linewidth='5',color='black')
# Setting shear stress label
ax3.set_xlabel(r'$\tau_{xy}$',fontsize='30')
# Setting normal stress limits
ax3.set_ylim([np.min(y),np.max(y)])
# Setting normal stress tick locations
ax3.set_yticks([])
# Setting velocity tick label size
ax3.tick_params(axis='x',labelsize=20)
ax3.tick_params(axis='y',labelsize=20)

# # Plotting velocity
# ax1=plt.subplot(131)
# ax1.plot(u/Uc,y/H,linewidth='5',color='black')
# # Setting velocity labels
# ax1.set_xlabel(r'$u/U_c$',fontsize='30')
# ax1.set_ylabel(r'$y/H$',fontsize='30')
# # Setting velocity limits
# ax1.set_xlim(left=0.0)
# ax1.set_ylim([np.min(y/H),np.max(y/H)])
# # Setting velocity tick locations
# ax1.set_xticks(np.arange(0.0,np.ceil(np.max(u/Uc))+1))
# ax1.set_yticks(np.arange(np.min(y/H),np.max(y/H)+H/2.0,H/2.0))
# # Setting velocity tick label size
# ax1.tick_params(axis='x',labelsize=20)
# ax1.tick_params(axis='y',labelsize=20)

# # Plotting normal stress
# ax2=plt.subplot(132)
# ax2.plot(Txx/(visc*Uc/H),y/H,linewidth='5',color='black')
# # Setting normal stress label
# ax2.set_xlabel(r'$\tau_{xx}/(\mu_0U_c/H)$',fontsize='30')
# # Setting normal stress limits
# ax2.set_xlim(left=0.0)
# ax2.set_ylim([np.min(y/H),np.max(y/H)])
# # Setting normal stress tick locations
# ax2.set_yticks([])
# # Setting velocity tick label size
# ax2.tick_params(axis='x',labelsize=20)
# ax2.tick_params(axis='y',labelsize=20)

# # Plotting shear stress
# ax3=plt.subplot(133)
# ax3.plot(Txy/(visc*Uc/H),y/H,linewidth='5',color='black')
# # Setting shear stress label
# ax3.set_xlabel(r'$\tau_{xy}/(\mu_0U_c/H)$',fontsize='30')
# # Setting normal stress limits
# ax3.set_ylim([np.min(y/H),np.max(y/H)])
# # Setting normal stress tick locations
# ax3.set_yticks([])
# # Setting velocity tick label size
# ax3.tick_params(axis='x',labelsize=20)
# ax3.tick_params(axis='y',labelsize=20)

# Add figure title
fig.suptitle(r'$\beta=%0.2f$' % beta, fontsize=40)

# Show the plot
plt.show()

# Saving the figure.
# plt.savefig(os.getcwd()+'/theory.png')


