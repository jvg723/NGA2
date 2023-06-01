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
Reg   =search_str(path,'Reynolds number')  # Gas Reynolds number
Weg   =search_str(path,'Weber number')     # Gas Weber number
visc_r=search_str(path,'Viscosity ratio')  # Viscosity ratio (visc_l/visc_g)
rho_r =search_str(path,'Density ratio')    # Density ratio (rho_l/rho_g)
D     =1.0                                 # Ligament diameter
Cd    =0.42

# Flow properties
visc_g=1.0/Reg;      rho_g=1.0
visc_l=visc_r*visc_g; rho_l=rho_r
sigma=1.0/Weg

#Ohnesorge number for case
Oh=visc_l/math.sqrt(rho_l*sigma*D)

lambda_RT=D*math.sqrt((6.0*math.pi**3.0*rho_l)/(Cd*(rho_l-rho_g)))*Weg**-0.5
print(lambda_RT)
