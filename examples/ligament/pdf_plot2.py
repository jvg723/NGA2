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

# Arguments
# to run 'python $python_script $data_to_plot'
home_directory = os.getcwd()
os.chdir(sys.argv[1])
base_directory = os.getcwd()
plot_lig_only=False
plot_film_only=False
plot_det_only=False

if len(sys.argv)>=3:
    if sys.argv[2]=='ligament':
        plot_lig_only=True
        plot_directory  = base_directory+'/pdf-ligament'
        os.makedirs(plot_directory,exist_ok=True)
    elif sys.argv[2]=='film':
        plot_film_only=True
        plot_directory  = base_directory+'/pdf-film'
        os.makedirs(plot_directory,exist_ok=True)
    elif sys.argv[2]=='detached':
        plot_det_only=True
        plot_directory  = base_directory+'/pdf-detached'
        os.makedirs(plot_directory,exist_ok=True)
    elif sys.argv[2]=='all':
        plot_directory  = base_directory+'/pdf-all'
        os.makedirs(plot_directory,exist_ok=True) 
else:
    plot_directory  = base_directory+'/pdf-all'
    os.makedirs(plot_directory,exist_ok=True)

# spray_ts=sys.argv[2]
# sprayfile="droplets.001680" # 128
# sprayfile="droplets.002820" # 256
# plot_name="velvsdiam"
plot_name_loglog = "npdfloglog"
plot_name_semilogy = "npdfsemilogy"
plot_name_2 = "npdfHR"
plot_name = "npdf"
plot_name_film = "npdf-film"

label_0 = "$D_{10}$"
label_1 = "$D_{32}$"
label_2 = "MMD"
label_3 = "Volume fit $5^3$"
label_4 = "Width=10"

dimensional_x_label = r'diameter $(d)$ [$\mathrm{\mu}$m]'
dimless_x_label = r'$d/d_0$'
dimensional_y_label = r'$pd(d)$ [1/$\mathrm{\mu}$m]'
dimless_y_label = r'$pd(d/d_0)$'

## Parameters
lz=10 #0.02
D=0.00254 # initial diameter
Dm=D/1e-6 # initial diameter in microns
# film_fraction=1 #0.13 #0.0855 #0.1235
initial_volume=np.pi/6.0*Dm**3.0

## Data folders


# # turb-inflow
# ndata=2
# order=[0,1]
# nx_list=[200,256,128,128]
# # folders=['bag-only/160-smallts','bag-only/256']
# folders=['ligament-breakup/160','ligament-breakup/256']
# elvira=[False,False]
# doplot=[True,True]
# labels=[r'$D_0/\Delta x=$'+f"{1/(lz/nx_list[0]):.1f}",r'$D_0/\Delta x=$'+f"{1/(lz/nx_list[1]):.1f}",r'$D_0/\Delta x=$'+f"{D/(lz/nx_list[2]):.1f}"+', ELVIRA',r'$D_0/\Delta x=$'+f"{D/(lz/nx_list[3]):.1f}"+', R2P']
# colors=["#5790fc", "#f89c20", "#e42536", "#964a8b"]

# # 128-batch
# ndata=5
# order=np.arange(24)
# nx_list=128*np.ones(24)
# # folders=['bag-only/160-smallts','bag-only/256']
# folders=list(map(str,np.arange(24)))
# labels=list(map(str,np.arange(24)))
# # colors=["#5790fc", "#f89c20", "#e42536", "#964a8b"]
# maxrows=1e5

# # 160-batch
# ndata=12
# order=np.arange(24)
# nx_list=160*np.ones(24)
# # folders=['bag-only/160-smallts','bag-only/256']
# folders=list(map(str,np.arange(24)))
# labels=list(map(str,np.arange(24)))
# # colors=["#5790fc", "#f89c20", "#e42536", "#964a8b"]
# maxrows=1e5

# # 200-batch
# ndata=12
# order=np.arange(24)
# nx_list=200*np.ones(24)
# # folders=['bag-only/160-smallts','bag-only/256']
# folders=list(map(str,np.arange(24)))
# labels=list(map(str,np.arange(24)))
# # colors=["#5790fc", "#f89c20", "#e42536", "#964a8b"]
# maxrows=1e5

# restarts
ndata=3
order=[0,1,2,3,4,5]
nx_list=[200,256,300]
folders=['200/0','256/0','300/0']
# folders=['160-laminar','200','200-1','200-1-0','256','256-1']
# elvira=[False,False,False,False]
# doplot=[True,True,True,True]
# labels=list(map(str,nx_list))
labels=folders
# labels=[r'$D_0/\Delta x=$'+f"{1/(lz/nx_list[0]):.1f}",r'$D_0/\Delta x=$'+f"{1/(lz/nx_list[1]):.1f}",r'$D_0/\Delta x=$'+f"{1/(lz/nx_list[2]):.1f}"]
# colors=["#5790fc", "#f89c20", "#e42536", "#964a8b"]
maxrows=1e5

if len(sys.argv)>=4:
    ndata=1
    order=[0]
    nx_list=[int(sys.argv[3])]
    folders=['.']
    elvira=[False]
    doplot=[True]
    labels=list(map(str,nx_list))
    maxrows=1e5

## Log-normal curve
def func(x,mu,sigma):
    return 1.0/(x*sigma*np.sqrt(2.0*np.pi))*np.exp(-(np.log(x)-mu)**2/(2*sigma**2))

# Experimental data
# e_npdf = np.loadtxt('/home/fs01/ah2262/cases/drop_breakup/guildenbecher-bag/guildenbecher_7b.txt',skiprows=1,usecols=(1),unpack=True)
# e_npdf_highres_bins,e_npdf_highres = np.loadtxt('/home/fs01/ah2262/cases/drop_breakup/guildenbecher-bag/guildenbecher_9.txt',skiprows=1,usecols=(0,1),delimiter=',',unpack=True)
# e_vpdf_bins,e_vpdf_vals = np.loadtxt('/home/fs01/ah2262/cases/drop_breakup/guildenbecher-bag/guildenbecher_10.txt',skiprows=1,usecols=(0,1),delimiter=',',unpack=True)
# e_vpdf_bins,e_vpdf_vals = np.loadtxt('/home/fs01/ah2262/cases/drop_breakup/guildenbecher-bag/guildenbecher_10_27micron.txt',skiprows=1,usecols=(0,1),delimiter=',',unpack=True)
e_npdf = np.loadtxt(home_directory+'/guildenbecher-bag/guildenbecher_7b.txt',skiprows=1,usecols=(1),unpack=True)
e_npdf_highres_bins,e_npdf_highres = np.loadtxt(home_directory+'/guildenbecher-bag/guildenbecher_9.txt',skiprows=1,usecols=(0,1),delimiter=',',unpack=True)
e_vpdf_bins,e_vpdf_vals = np.loadtxt(home_directory+'/guildenbecher-bag/guildenbecher_10.txt',skiprows=1,usecols=(0,1),delimiter=',',unpack=True)

# Start plot
fig,ax = plt.subplots(dpi=600,figsize=[5,3]) # Fig. 9 high-res experiment
fig1, ax1 = plt.subplots(dpi=300,figsize=[6,3]) # Fig. 7 low-res experiment
fig2, ax2 = plt.subplots(dpi=300,figsize=[6,5]) # Fig. 10 volume-weighted experiment, later time
fig3, ax3 = plt.subplots(dpi=300,figsize=[6,3]) # Fig. 9 volume-weighted high-res experiment

y_data = {}
for j in range(0,ndata):
    i=order[j]
    # if doplot[i]==False:
    #     continue
    os.chdir(folders[i])
    # sprayfile="droplets.00"+spray_ts[i]
    # filmfile="cells.00"+bagcells_ts[i]
    nx=nx_list[i]
    print(os.getcwd())
    # skipping 5 rows because of mistake in formatting
    # diam_m,u,v,w,vel_total = np.loadtxt('spray/'+sprayfile,skiprows=5,usecols=(0,1,2,3,4),unpack=True)
    diam_m,u,v,w,vel_total = np.genfromtxt('spray-all/droplets',skip_header=1,usecols=(0,1,2,3,4),unpack=True,max_rows=maxrows)
    origin = np.genfromtxt('spray-all/droplets',skip_header=1,usecols=(8),dtype='str',max_rows=maxrows)
    # diam_m,u,v,w,vel_total,x,y,z = np.loadtxt('spray/'+sprayfile,skiprows=5,usecols=(0,1,2,3,4,5,6,7),unpack=True)
    # timestep,time_s,ts_size,n_drops,d10_m,d32_m,mmd_m = np.loadtxt('monitor/dropstats',skiprows=2,usecols=(0,1,2,3,4,5,6),unpack=True)
    # id,min_thickness,thickness_temp,volume_temp = np.loadtxt('film/'+filmfile,skiprows=1,usecols=(0,1,2,3),unpack=True)
    # id_mode = stats.mode(id)[0]
    # thickness = thickness_temp[id==id_mode]/1e-6
    # volume = volume_temp[id==id_mode]/1e-6

    diam = diam_m*D/1e-6 # droplet diameter in microns

    ## for dropstats
    # d10 = d10_m/1e-6
    # d32 = d32_m/1e-6
    # mmd = mmd_m/1e-6
    # time = time_s/1e-3
    diam_cutoff = 27 # Minimum diameter detectable in experiments of Guildenbecher et al.
    flotsam_cutoff = 1e6*lz/nx_list[i]/2 # microns flotsam
    elvira_cutoff = max(diam_cutoff,flotsam_cutoff)
    # large_diam = diam[diam>diam_cutoff]
    # if elvira[i]==True:
    #     large_diam=diam[diam>elvira_cutoff]
    #     diam=diam[diam>flotsam_cutoff]
    # else:
        # large_diam = diam[diam>diam_cutoff]
    if plot_lig_only:
        diam=diam[origin=='ligament']
    elif plot_film_only:
        diam=diam[origin=='film']
    elif plot_det_only:
        diam=diam[origin=='detached']
        
    large_diam = diam[diam>diam_cutoff]    
    # Weight by volume
    weights = np.pi/6.0*diam**3.0
    # large_diam_vel_total = vel_total[diam>diam_cutoff]
    os.chdir(base_directory)
    # if no droplets skip
    if (np.sum(diam_m) <= 0.0):
        continue
    # logbins=np.logspace(1.0,11.0,num=30,base=2.0)
    logbins=e_vpdf_bins/Dm
    logbins_hr=e_npdf_highres_bins
    vals1,edges1 = np.histogram(large_diam/Dm, bins=24,range=(0,600/Dm),density=True)  # arguments are passed to np.histogram
    vals,edges = np.histogram(diam/Dm, bins=27,range=(0,270/Dm),density=True)  # arguments are passed to np.histogram
    # vals2,edges2 = np.histogram(diam, bins=logbins,range=(0,2000),weights=weights,density=True)  # arguments are passed to np.histogram
    vals2,edges2 = np.histogram(diam/Dm, bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
    # vals,edges= np.histogram(diam, bins=27,range=(0,270),density=False)  # arguments are passed to np.histogram
    # vals3,edges3 = np.histogram(thickness, bins=25,range=(0,50),weights=volume,density=True)  # arguments are passed to np.histogram
    # vals3,edges3 = np.histogram(thickness, bins=25,range=(0,50),density=True)  # arguments are passed to np.histogram
    # vals3,edges3 = np.histogram(diam, bins=logbins_hr,weights=weights,density=True)  # arguments are passed to np.histogram
    # vals3,edges3 = np.histogram(diam, bins=27,range=(0,270),weights=weights,density=True)  # arguments are passed to np.histogram
    vals3,edges3 = np.histogram(diam/Dm, bins=27,range=(0,270/Dm),weights=weights,density=True)  # arguments are passed to np.histogram

    # Standard res
    buf1=np.copy(vals1)
    buf1[:]=0.0
    buf1[0:np.size(e_npdf)]=e_npdf*Dm

    # High res
    buf=e_npdf_highres*Dm

    edgemed = (edges[:-1]+edges[1:])/2
    edgemed1 = (edges1[:-1]+edges1[1:])/2
    edgemed2 = (edges2[:-1]+edges2[1:])/2
    edgemed3 = (edges3[:-1]+edges3[1:])/2
    # ax3.plot(edgemed3,vals3,'-x',linewidth=1.0,label=labels[i],color=colors[i])
    # ax2.plot(edgemed2,vals2*film_fraction,'-x',linewidth=2.0,label=labels[i],color=colors[i])
    # ax1.plot(edgemed1,vals1,'-x',linewidth=2.0,label=labels[i],color=colors[i])
    # ax.plot(edgemed,vals,'-x',linewidth=2.0,label=labels[i],color=colors[i])
    ax3.plot(edgemed3,vals3,'-x',linewidth=1.0,label=labels[i])
    ax2.plot(edgemed2,vals2*np.sum(weights)/initial_volume,'-x',linewidth=2.0,label=labels[i])
    y_data[j]=vals2*np.sum(weights)/initial_volume
    ax1.plot(edgemed1,vals1,'-x',linewidth=2.0,label=labels[i])
    ax.plot(edgemed,vals,'-x',linewidth=2.0,label=labels[i])
# vpdf_bins

# Standard res
ax1.bar(edges1[:-2],buf1[:-1],width=np.diff(edges1[:-1]),color='0.7',edgecolor="black",align="edge",label='Experiment')

# High res
ax.bar(edges[:-2],buf[1:],width=np.diff(edges[:-1]),color='0.7',edgecolor="black",align="edge",label='Experiment')

# Volume-weighted
ax2.bar(logbins[:-1],e_vpdf_vals[:-1]*Dm,width=np.diff(logbins),color='0.7',edgecolor="black",align="edge",label='Experiment')

# Volume-weighted HR - need to normalize
area=sum((edges3[1:-1]-edges3[:-2])*buf[1:]*np.pi/6.0*(edgemed3[:-1])**3.0)
ax3.bar(edges3[:-2],buf[1:]*np.pi/6.0*(edgemed3[:-1])**3.0/area,width=np.diff(edges3[:-1]),color='0.7',edgecolor="black",align="edge",label='Experiment')


ax.set_xlabel(dimless_x_label)
ax.set_ylabel(dimless_y_label)
# ax.set_title(my_title)
# ax.set_ylim(0,.03)
# ax.set_ylim(0,.06)
ax.legend(bbox_to_anchor=(1.05, 1),
                         loc='upper left', borderaxespad=0.)
ax1.legend()
ax2.legend()
# ax2.set_ylim(top=0.00125)
# ax2.legend(bbox_to_anchor=(0, 1.05),
#                          loc='lower left')
ax3.legend()
# fig.tight_layout()
os.chdir(plot_directory)
if plot_film_only:
    fig.savefig(plot_name_2+'_f9.png') # fig 9 in Guildenbecher et al. 2017
    ax1.set_xlabel(dimless_x_label)
    ax1.set_ylabel(dimless_y_label)
    fig1.savefig(plot_name+'_LR_f7.png') # fig 7 in Guildenbecher et al. 2017
    # ax3.set_xlabel(dimless_x_label)
    ax3.set_xlabel(dimless_x_label)
    ax3.set_ylabel(r'volume weighted density')
    fig3.savefig(plot_name_2+'_weighted_f9.png')
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.savefig(plot_name_loglog+'_HR_f9.png')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    fig1.savefig(plot_name_loglog+'_LR_f7.png')

# Figure 10: vpdf from film, rim, and nodes
ax2.set_xlabel(dimless_x_label)
ax2.set_ylabel(r'volume weighted probability [1/$\mathrm{\mu}$m]')
fig2.savefig(plot_name+'_weighted_f10.png') # fig 10 in Guildenbecher et al. 2017
# ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylim(bottom=5e-5)
fig2.savefig(plot_name_semilogy+'_weighted_f10.png')
# fig.savefig(plot_name_loglog+'count.png')

## output tables
np.savetxt("results/convtable.dat", np.column_stack([edgemed2,y_data[0]
                                            # ,list(y_data["zalesak"]["r2p"].values())
                                            # ,list(y_data["def2d"]["elvira"].values())
                                            # ,list(y_data["def2d"]["r2p"].values())
                                            ,y_data[1]
                                            ,y_data[2]])
                                            , fmt='%.8e', delimiter="  "
                                            , header="nx              zalesak_elvira" 
                                            + "  zalesak_r2p     def2d_elvira   "
                                            + " def2d_r2p       def3d_elvira   "
                                            + " def3d_r2p", comments="")

############################
# PLOTS
############################
# Number PDF

# fig2, ax2 = plt.subplots(1, 1)

# vals,edges,patches = ax2.hist(diam, bins=60,range=(0,300),edgecolor="black",density=True)  # arguments are passed to np.histogram
# mids=(edges[0:-1]+edges[1:])/2
# popt, pcov = curve_fit(func, mids, vals)
# plotx = np.linspace(0, 300, 1000)
# # print('mu sigma')
# # print(popt)
# ax2.plot(plotx,func(plotx,popt[0],popt[1]))
# ax2.set_xlabel(dimensional_x_label)
# ax2.set_ylabel(r'$pd(d)$ [1/$\mathrm{\mu}$s]')
# # ax.set_ylim(0,.016)
# # ax.set_ylim(0,.05)
# fig2.savefig(plot_name_2+'.png')

# # Number PDF
# # fig = plt.figure(dpi=300)
# # ax = plt.subplot(111)
# vals,edges,patches = ax.hist(large_diam, bins=24,range=(0,600),edgecolor="black",density=True)  # arguments are passed to np.histogram
# buf=np.copy(vals)
# buf[:]=0.0
# buf[0:np.size(e_npdf)]=e_npdf
# # ax.bar(edges[:-1],buf,width=np.diff(edges),edgecolor="black",align="edge")
# ax.set_xlabel(dimensional_x_label)
# ax.set_ylabel(r'$pd(d)$ [1/$\mathrm{\mu}$m]')
# ax.set_title('Simulation')
# # ax.set_ylim(0,.016)
# ax.set_ylim(0,.05)
# # fig.savefig(plot_name_2+'largediam.png')

# # Number PDF
# # fig = plt.figure(dpi=300)
# # ax = plt.subplot(111)
# edgemed = (edges[:-1]+edges[1:])/2
# ax.bar(edges[:-1],buf,width=np.diff(edges),edgecolor="black",align="edge",label='Experiment')
# ax.plot(edgemed,vals,'-xr',linewidth=2.0,label='Simulation')
# # ax.bar(edges[:-1],vals,width=np.diff(edges),edgecolor="black",align="edge",label='Simulation')
# # ax.plot(edgemed,buf,'-xr',linewidth=2.0,label='Experiment')
# ax.set_xlabel(dimensional_x_label)
# ax.set_ylabel(r'$pd(d)$ [1/$\mathrm{\mu}$m]')
# # ax.set_title('Experiment')
# # ax.set_ylim(0,.016)
# ax.set_ylim(0,.05)
# ax.legend()
# # fig.savefig(plot_name_2+'largediam_exp_swap.png')
