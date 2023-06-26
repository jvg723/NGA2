# -*- coding: utf-8 -*-

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_daq as daq
from dash import dcc
from dash import html
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import matplotlib.pyplot as plt
import math
import re
import os
import numpy as np
from scipy import stats
pd.options.plotting.backend = "plotly"

## Define functions

# Log-normal curve
def log_norm(x,mu,sigma):
    return 1.0/(x*sigma*np.sqrt(2.0*np.pi))*np.exp(-(np.log(x)-mu)**2/(2*sigma**2))

# Process the input file to obtain simulation parameters
with open('input') as fobj:
    for line in fobj:
        line_data = re.split(':',line)
        if line_data[0].rstrip()=='Reynolds number':
            Reynolds=float(line_data[1])
        if line_data[0].rstrip()=='Weber number':
            Weber=float(line_data[1])
        if line_data[0].rstrip()=='Viscosity ratio':
            visc_r=float(line_data[1])
        if line_data[0].rstrip()=='Density ratio':
            rho_r=float(line_data[1])
        if line_data[0].rstrip()=='Target Re_lambda':
            Re_lambda=float(line_data[1])
        if line_data[0].rstrip()=='Turbulence intensity':
            TI=float(line_data[1])

# Gas density
rho_g=1
# Compute gas viscosity
visc_g=1/Reynolds
# Liquid density
rho_l=rho_r
# Liquid viscosity
visc_l=visc_g*visc_r
# Ligament diameter
Diam=1.0

# Directories
case_dir='04_lig'                                                                         # case directory
base_dir=os.getcwd()                                                                      # base directroy
spray_dir='/Users/josephgiliberto/Builds/NGA2/examples/ligament/cases/'+case_dir+'/spray' # spray data directory

# Prepare PDF of generated drops
os.chdir(spray_dir)                                                                                        # move into spray directory
sprayfile='droplets.005051'                                                                                # current timestep
df=pd.read_csv(sprayfile, delim_whitespace=True, header=None, skiprows=1, usecols=[0], names=['Diameter']) # Load data
# df=pd.DataFrame(columns=['time','diameters','PDF','mean','std'])    # Empty struct
#Loop through output files 
# for dropfile in os.listdir(spray_dir):
#     diam=pd.read_csv(dropfile, delim_whitespace=True, header=None, skiprows=1, usecols=[0], names=['Diameter'])   # read in data file
#     mean=np.mean(np.log(diam['Diameter']))                                                                        # mean of ln(diameter)
#     mean_print=np.mean(diam['Diameter'])                                                                          # mean of diameter
#     std=np.std(np.log(diam['Diameter']))                                                                          # standard deviation of ln(diameter)
#     std_print=np.std(diam['Diameter'])                                                                            # standard deivation of diameter
#     diam['PDF']=log_norm(diam['Diameter'],mean,std)                                                               # pdf of drop_data
#     time_stamp=re.split('s',dropfile)                                                                             # time stamp
#     df=df._append({'time':float(time_stamp[1]),'diameters':[diam['Diameter']],'PDF':[diam['PDF']],'mean':mean_print,'std':std_print}, ignore_index = True)
    
# df=df.sort_values('time')   # sort data frame by time
# os.chdir(base_dir)          # move back to base director

mean=np.mean(np.log(df['Diameter']))                                                                       # mean of ln(df[Diameter])
mean_print=np.mean(df['Diameter'])    
std=np.std(np.log(df['Diameter']))                                                                         # standard deviation of ln(df[Diameter])
std_print=np.std(df['Diameter'])                                                                        
df['PDF']=log_norm(df['Diameter'],mean,std)                                                                # pdf of drop_data


# Plot PDF
os.chdir(base_dir)          # move back to base director
fig1=go.Figure()
fig1.add_trace(go.Scatter(mode='markers',x=df['Diameter']/Diam, y=df['PDF']*Diam, showlegend=False))
fig1.update_layout(width=800,height=600)
fig1.update_xaxes(type='log',title_text='$\Large{d/D}$',exponentformat="power",title_font_size=48,tickfont_size=24)
fig1.update_yaxes(type='log',title_text='$\Large{\mathrm{{PDF}}(d)D}$', exponentformat="power",title_font_size=48,tickfont_size=24)

# This is where we define the dashboard layout
app = dash.Dash(__name__)
app.layout = html.Div(style={"margin-left": "15px"},children=[

    # Title of doc
    dcc.Markdown('''# Break-up of a liquid ligament'''),
    dcc.Markdown('''*NGA2 Dashboard written by J. Giliberto, last updated 06/22/2023*'''),
    
    # Intro
    dcc.Markdown('''
    ## Overview

    In this dashboard, we post-process the raw data generated by NGA2's ligament
    case. This simulation is based on the experimental work of Nasser Ashgriz to study the break-up behavior of an isolated liquid 
    ligament in a turbulent cross flow.
    '''),
    
    # Simulation parameters
    dcc.Markdown(f'''
    ---
    ## Simulation parameters

    By analyzing the input file, we have detected the following parameters:
    - Reynolds number: $\mathrm{{Re}}$={Reynolds}
    - Weber number: $\mathrm{{We}}$={Weber}
    - Viscosity ratio: $\mu_l / \mu_g$={visc_r}
    - Density ratio: $\\rho_l / \\rho_g$={rho_r}
    - Taylor-scale Reynolds number: $\mathrm{{R}}_{{\lambda}}$={Re_lambda}
    - Turbulence intensity: $\mathrm{{TI}}$={TI}

    ''',mathjax=True),
    
    # Droplet size distribution 
    # first row
    html.Div([
            dcc.Markdown(f'''
            ---
            ## Droplet size distribution
            
            The graph below shows the probability density function of the drops generated during the break-up of a liquid liqament. The droplet 
                        diameter, $d$, and probability density function, $\mathrm{{PDF}}(d)$, have been normalized by the ligament 
                        diameter, $D$. The slider allows for the $\mathrm{{PDF}}(d)$ to be viewed at different output times. 
            ''',mathjax=True),

    ],className='row'),

    
    # second row
    html.Div(
    [
        # first column
        html.Div(
            [
                dcc.Graph(id="drop_size", figure=fig1,mathjax=True)],
            style={"width": "65%", "float": "left", "display": "inline-block"},
        ),
        # second column
        html.Div(
            [
                dcc.Markdown(f'''
                ## Statistical analysis
                - Log-normal probability density function:
                                 
                $\qquad \qquad \large{{\mathrm{{PDF}} = \\frac{{1}}{{x\sigma\sqrt{{2\pi}}}} \exp(-\\frac{{\ln x -\mu^2}}{{2\sigma^2}})}}$   
                             
                - Mean diameter: $\mu/D$ = {round(mean_print,4)}
                - Standard deviation of diameter: $\sigma/D$ = {round(std_print,4)}
                
                ## Export raw $\mathrm{{PDF}}$ data
                ''',mathjax=True),
            ],
            style={
                "width": "35%","display": "inline-block", "margin-top": "60px",
            },
        ),

    ],
    style={"margin": "40px"},className='row'
),
    
    
]) # End app

# @app.callback(
#     dash.dependencies.Output('Radius-slider','figure'),
#     [dash.dependencies.Input('time-slider','value')])
# def update_figure(value):
#     myind=value
#     myrf=rf0.copy()
#     myrf['y']=rf['y']/Rinit
#     myrf['R']=rf[ftime[myind]]/Rinit
#     myrf=myrf[myrf['y']>=0]
#     myrf=myrf[myrf['R']>0]

#     fig2=go.Figure()
#     fig2.add_trace(go.Scatter(x=myrf['R'], y=myrf['y'], showlegend=False, mode='lines+markers', line_color='red', line_shape='spline'))
#     fig2.update_layout(width=700,height=800)
#     fig2.update_xaxes(title_text='Normalized distance from centerline axis of drop',title_font_size=24,tickfont_size=24,range=[0,2.5])
#     fig2.update_yaxes(title_text='Normalized height above plate',title_font_size=24,tickfont_size=24,range=[0,3])
#     fig2.add_annotation(x=1.5, y=2.75,text='Normalized time={:.3f}'.format(rtime[myind]/Tref),showarrow=False, font_size=24, font_color='darkslategrey')

#     return fig2


if __name__ == '__main__':
    app.run_server(debug=True)
