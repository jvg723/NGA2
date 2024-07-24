# Importing packages
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import re

# Experimental data from Pilz & Brenn (2007) for aqueous 0.8% weight P2500 solution
exp_data=np.zeros((30,2))
####### strain rate (s^-1) ######## ####### viscosity (Pa*s) ########
exp_data[0][0] =0.03320264167012359; exp_data[0][1] =1.555542634556734
exp_data[1][0] =0.04759095611380148; exp_data[1][1] =1.555542634556734
exp_data[2][0] =0.06821442481378769; exp_data[2][1] =1.529332243485468
exp_data[3][0] =0.09777504241665053; exp_data[3][1] =1.529332243485468
exp_data[4][0] =0.1422638188714286 ; exp_data[4][1] =1.529332243485468
exp_data[5][0] =0.19494059178649895; exp_data[5][1] =1.503563489039807
exp_data[6][0] =0.28364081830232807; exp_data[6][1] =1.478228929785222
exp_data[7][0] =0.41270067495909846; exp_data[7][1] =1.381088160370716
exp_data[8][0] =0.5740604865125165 ; exp_data[8][1] =1.334938488311054
exp_data[9][0] =0.8352646550276838 ; exp_data[9][1] =1.1652541657741733
exp_data[10][0]=1.1618406833946238 ; exp_data[10][1]=0.9831503229233347
exp_data[11][0]=1.640528249962994  ; exp_data[11][1]=0.8437217834532673
exp_data[12][0]=2.35144868059379   ; exp_data[12][1]=0.6880791154915206
exp_data[13][0]=3.4730945968949705 ; exp_data[13][1]=0.551692928104595 
exp_data[14][0]=4.831022808952545  ; exp_data[14][1]=0.4499212767026571
exp_data[15][0]=7.029194126327743  ; exp_data[15][1]=0.3486867337294896
exp_data[16][0]=10.075284747070821 ; exp_data[16][1]=0.28436395033707657
exp_data[17][0]=14.441394121461085 ; exp_data[17][1]=0.23190689071079063
exp_data[18][0]=20.087763937159924 ; exp_data[18][1]=0.1923680051224494
exp_data[19][0]=28.364081830232834 ; exp_data[19][1]=0.15423817643968538
exp_data[20][0]=40.65561370092818  ; exp_data[20][1]=0.12158243827894
exp_data[21][0]=59.15438869368025  ; exp_data[21][1]=0.09748322521355639
exp_data[22][0]=84.78885337587096  ; exp_data[22][1]=0.08224872063395836
exp_data[23][0]=119.72253273028251 ; exp_data[23][1]=0.06939503726003293
exp_data[24][0]=174.1976714117908  ; exp_data[24][1]=0.05564004260551747
exp_data[25][0]=245.9684921635237  ; exp_data[25][1]=0.05024663830396288
exp_data[26][0]=357.88700421234427 ; exp_data[26][1]=0.04097756288702974
exp_data[27][0]=512.976795051085   ; exp_data[27][1]=0.0363819195224116
exp_data[28][0]=702.9194126327736  ; exp_data[28][1]=0.032301678647516555
exp_data[29][0]=1007.5284747070853 ; exp_data[29][1]= 0.02967049031542662

# # Bisection root finder function
# def bisection(f,a,b)
#     # method parameters
#     k_max=200 # max number of iterations
#     tol=1e-10 # tolerance

#     # root estimate at first iteration
#     c=(a+b)/2.0

#     # return root estimate at first iteration if it is a root of f(x)
#     if f(c)


# Plotting
plt.figure(figsize=(10, 6))
plt.scatter(exp_data[:, 0], exp_data[:, 1], marker='o', color='k', label='Experimental Data (Pilz & Brenn, 2007)')
plt.xscale('log')  # Set logarithmic scale for the x-axis
plt.yscale('log')  # Set logarithmic scale for the y-axis
plt.xlabel('Strain Rate (s^-1)')
plt.ylabel('Viscosity (Pa-s)')
plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()
