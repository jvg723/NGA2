# Importing packages
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import re

## Functions ##

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

# Case Input
visc_l=search_str(path,'Liquid dynamic viscosity')      # (N*s/m^2) Liquid viscosity
rho_l  =search_str(path,'Liquid density')               # (kg/m^3)  Liquid density
visc_g=search_str(path,'Gas dynamic viscosity')         # (N*s/m^2) Gas viscosity
rho_g  =search_str(path,'Gas density')                  # (kg/m^3)  Gas density
sigma  =search_str(path,'Surface tension coefficient')  # (N/m)     Surface tension
g      =9.81                                            # (m/s^2)   gravity
D      =search_str(path,'Bubble diameter')              # (m)       bubble diameter


# Dimensionless numbers
Mo=(g*visc_l**4.0)/(rho_l*sigma**3.0)   # Morton number
print('Mo=',Mo)
#Terminal velocity correlation (Angelino 1965)
m=0.167*(1.0+0.34*Mo**0.24)
K=25.0/(1.0+0.33*Mo**0.29)
V=(4.0/3.0)*math.pi*(D/2.0)**3.0*1e6
print('Vol=',V)
print
vt=K*V**m
vt=vt/100.0
print('vt=',vt)
Eo=(g*D**2.0*rho_l)/sigma               # Eotvos number
print('Eo=',Eo)
Re=(rho_l*D*vt)/visc_l                  # Reynolds number
print('Re=',Re)


