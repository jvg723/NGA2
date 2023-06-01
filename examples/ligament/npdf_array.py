#Importing packages
import os
import sys

# Set path to case folder(s)
base_directory = os.getcwd()
sprayfile="droplets.004368" # 128
path=base_directory+'/cases/01_lig/spray'+sprayfile
print(path)

# Arguments
to run 'python $python_script $data_to_plot'
os.chdir(sys.argv[1])
base_directory = os.getcwd()
if len(sys.argv)==1:
    my_data = 'guildenbecher-bag/256'
else:
    my_data = sys.argv[1]
if len(sys.argv)==2:
    spray_ts="1680"
else:
    spray_ts = sys.argv[2]
spray_ts=sys.argv[2]
sprayfile="droplets.001680" # 128
sprayfile="droplets.002820" # 256
plot_name="velvsdiam"
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

## Parameters
# lz=0.02
# D=0.00254

## Data folders
# ndata=6
# nx_list=[50,64,100,128,200,256]
# folders=['50','64','100','128','200','256']
# spray_ts=['0980','1120','1500','1310','2640','2700'] # eth-air-we-17-highRe

# ndata=6
# nx_list=[50,64,100,128,200,256]
# folders=['50','64','100','128','200','256']
# spray_ts=['0781','1000','1230','1220','2320','4520'] # eth-air-we-20-highRe
# folders=['50','64','100','128','200','256-0']
# spray_ts=['0781','1000','1230','1220','2320','2533'] # eth-air-we-20-highRe using 256-0

# # guildenbecher-bag
# ndata=5
# nx_list=[64,100,128,200,256]
# folders=['64','100','128','200','256']
# spray_ts=['1310','1800','1680','2940','2820']
# # bagcells_ts=['1300','1790','1670','2930','2810']
# # label=r'$D_0/\Delta x=$'
# labels=folders

# # thickness-threshold
# ndata=4
# nx_list=[64,64,64,128]
# # folders=['64/1e-6','64/1e-7','64/1e-8','128/1e-8']
# # spray_ts=['2320','2420','0920','1680']#'2490']
# # bagcells_ts=['2310','2410','0910','1670']#'2490']
# folders=['64/1e-6','64/5e-7-lessfilm','64/1e-8','128/1e-8']
# spray_ts=['2320','0350','0920','1680']#'2490']
# bagcells_ts=['2310','0340','0910','1670']#'2490']
# labels=folders

# # film-thickness-to-diameter
# ndata=1
# nx_list=[100,64,64,128,64]
# folders=['power-law/100','64/5e-7-lessfilm','64/1e-8','128/1e-8','64/1e-9']
# spray_ts=['0230','0350','0920','1680','2490']
# bagcells_ts=['0220','0340','0910','1670','2490']
# labels=folders

# # r2p-ST
# ndata=3
# nx_list=[100,100,100,128,64]
# folders=['100','100-1e-6','100-csf','128/1e-8','64/1e-9']
# spray_ts=['2360','0750','0510','1680','2490']
# bagcells_ts=['2350','0740','0500','1670','2490']
# labels=folders

# movframe
# ndata=1
# # order=[0,2,1,3]
# # order=[3]
# order=[0]
# nx_list=[100,100,256,256]
# folders=['retraction/100-1','100-9','256-elvira-stats','256-2']
# elvira=[False,False,True,False]
# doplot=[True,True,True,True]
# # doplot=[False,False,False,True]
# labels=[r'$D_0/\Delta x=$'+f"{D/(lz/nx_list[0]):.1f}"+', ELVIRA',r'$D_0/\Delta x=$'+f"{D/(lz/nx_list[1]):.1f}"+', R2P',r'$D_0/\Delta x=$'+f"{D/(lz/nx_list[2]):.1f}"+', ELVIRA',r'$D_0/\Delta x=$'+f"{D/(lz/nx_list[3]):.1f}"+', R2P']
# colors=["#5790fc", "#f89c20", "#e42536", "#964a8b"]

# ## Log-normal curve
# def func(x,mu,sigma):
#     return 1.0/(x*sigma*np.sqrt(2.0*np.pi))*np.exp(-(np.log(x)-mu)**2/(2*sigma**2))

# Experimental data
# e_npdf = np.loadtxt('/home/fs01/ah2262/cases/drop_breakup/guildenbecher-bag/guildenbecher_7b.txt',skiprows=1,usecols=(1),unpack=True)
# e_npdf_highres = np.loadtxt('/home/fs01/ah2262/cases/drop_breakup/guildenbecher-bag/guildenbecher_9.txt',skiprows=1,usecols=(1),unpack=True)
# e_vpdf_bins,e_vpdf_vals = np.loadtxt('/home/fs01/ah2262/cases/drop_breakup/guildenbecher-bag/guildenbecher_10.txt',skiprows=1,usecols=(0,1),delimiter=',',unpack=True)
# e_vpdf_bins,e_vpdf_vals = np.loadtxt('/home/fs01/ah2262/cases/drop_breakup/guildenbecher-bag/guildenbecher_10_27micron.txt',skiprows=1,usecols=(0,1),delimiter=',',unpack=True)
# e_npdf = np.loadtxt('guildenbecher-bag/guildenbecher_7b.txt',skiprows=1,usecols=(1),unpack=True)
# e_npdf_highres_bins,e_npdf_highres = np.loadtxt('guildenbecher-bag/guildenbecher_9.txt',skiprows=1,usecols=(0,1),delimiter=',',unpack=True)
# e_vpdf_bins,e_vpdf_vals = np.loadtxt('guildenbecher-bag/guildenbecher_10.txt',skiprows=1,usecols=(0,1),delimiter=',',unpack=True)

# Start plot
# fig,ax = plt.subplots(dpi=600,figsize=[5,3]) # high-res experiment
# fig1, ax1 = plt.subplots(dpi=300,figsize=[6,3]) # low-res experiment
# fig2, ax2 = plt.subplots(dpi=300,figsize=[6,5]) # volume-weighted experiment, later time
# fig3, ax3 = plt.subplots(dpi=300,figsize=[6,3]) # volume-weighted high-res experiment

# for j in range(0,ndata):
#     i=order[j]
#     if doplot[i]==False:
#         continue
#     os.chdir(folders[i])
#     # sprayfile="droplets.00"+spray_ts[i]
#     # filmfile="cells.00"+bagcells_ts[i]
#     nx=nx_list[i]

#     # skipping 5 rows because of mistake in formatting
#     # diam_m,u,v,w,vel_total = np.loadtxt('spray/'+sprayfile,skiprows=5,usecols=(0,1,2,3,4),unpack=True)
#     diam_m,u,v,w,vel_total = np.loadtxt('spray-all/droplets',skiprows=5,usecols=(0,1,2,3,4),unpack=True)
#     # diam_m,u,v,w,vel_total,x,y,z = np.loadtxt('spray/'+sprayfile,skiprows=5,usecols=(0,1,2,3,4,5,6,7),unpack=True)
#     # timestep,time_s,ts_size,n_drops,d10_m,d32_m,mmd_m = np.loadtxt('monitor/dropstats',skiprows=2,usecols=(0,1,2,3,4,5,6),unpack=True)
#     # id,min_thickness,thickness_temp,volume_temp = np.loadtxt('film/'+filmfile,skiprows=1,usecols=(0,1,2,3),unpack=True)
#     # id_mode = stats.mode(id)[0]
#     # thickness = thickness_temp[id==id_mode]/1e-6
#     # volume = volume_temp[id==id_mode]/1e-6

#     diam = diam_m/1e-6

#     ## for dropstats
#     # d10 = d10_m/1e-6
#     # d32 = d32_m/1e-6
#     # mmd = mmd_m/1e-6
#     # time = time_s/1e-3
#     diam_cutoff = 27 # Minimum diameter detectable in experiments of Guildenbecher et al.
#     flotsam_cutoff = 1e6*lz/nx_list[i]/2 # microns flotsam
#     elvira_cutoff = max(diam_cutoff,flotsam_cutoff)
#     # large_diam = diam[diam>diam_cutoff]
#     if elvira[i]==True:
#         large_diam=diam[diam>elvira_cutoff]
#         diam=diam[diam>flotsam_cutoff]
#     else:
#         large_diam = diam[diam>diam_cutoff]
#     # Weight by volume
#     weights = np.pi/6.0*diam**3.0
#     # large_diam_vel_total = vel_total[diam>diam_cutoff]
#     os.chdir(base_directory)

#     # logbins=np.logspace(1.0,11.0,num=30,base=2.0)
#     logbins=e_vpdf_bins
#     logbins_hr=e_npdf_highres_bins
#     vals1,edges1 = np.histogram(large_diam, bins=24,range=(0,600),density=True)  # arguments are passed to np.histogram
#     vals,edges = np.histogram(diam, bins=27,range=(0,270),density=True)  # arguments are passed to np.histogram
#     # vals2,edges2 = np.histogram(diam, bins=logbins,range=(0,2000),weights=weights,density=True)  # arguments are passed to np.histogram
#     vals2,edges2 = np.histogram(diam, bins=logbins,weights=weights,density=True)  # arguments are passed to np.histogram
#     # vals,edges= np.histogram(diam, bins=27,range=(0,270),density=False)  # arguments are passed to np.histogram
#     # vals3,edges3 = np.histogram(thickness, bins=25,range=(0,50),weights=volume,density=True)  # arguments are passed to np.histogram
#     # vals3,edges3 = np.histogram(thickness, bins=25,range=(0,50),density=True)  # arguments are passed to np.histogram
#     # vals3,edges3 = np.histogram(diam, bins=logbins_hr,weights=weights,density=True)  # arguments are passed to np.histogram
#     vals3,edges3 = np.histogram(diam, bins=27,range=(0,270),weights=weights,density=True)  # arguments are passed to np.histogram

#     # Standard res
#     buf1=np.copy(vals1)
#     buf1[:]=0.0
#     buf1[0:np.size(e_npdf)]=e_npdf

#     # High res
#     buf=e_npdf_highres

#     edgemed = (edges[:-1]+edges[1:])/2
#     edgemed1 = (edges1[:-1]+edges1[1:])/2
#     edgemed2 = (edges2[:-1]+edges2[1:])/2
#     edgemed3 = (edges3[:-1]+edges3[1:])/2
#     ax3.plot(edgemed3,vals3,'-x',linewidth=1.0,label=labels[i],color=colors[i])
#     ax2.plot(edgemed2,vals2,'-x',linewidth=2.0,label=labels[i],color=colors[i])
#     ax1.plot(edgemed1,vals1,'-x',linewidth=2.0,label=labels[i],color=colors[i])
#     ax.plot(edgemed,vals,'-x',linewidth=2.0,label=labels[i],color=colors[i])
# # vpdf_bins

# # Standard res
# ax1.bar(edges1[:-2],buf1[:-1],width=np.diff(edges1[:-1]),color='0.7',edgecolor="black",align="edge",label='Experiment')

# # High res
# ax.bar(edges[:-2],buf[1:],width=np.diff(edges[:-1]),color='0.7',edgecolor="black",align="edge",label='Experiment')

# # Volume-weighted
# ax2.bar(e_vpdf_bins[:-1],e_vpdf_vals[:-1],width=np.diff(e_vpdf_bins),color='0.7',edgecolor="black",align="edge",label='Experiment')

# # Volume-weighted HR
# ax3.bar(edges3[:-2],buf[1:]*np.pi/6.0*(edgemed3[:-1])**3.0*1e-6,width=np.diff(edges3[:-1]),color='0.7',edgecolor="black",align="edge",label='Experiment')


# ax.set_xlabel(r'diameter $(d)$ [$\mathrm{\mu}$m]')
# ax.set_ylabel(r'$pd(d)$ [1/$\mathrm{\mu}$m]')
# # ax.set_title(my_title)
# # ax.set_ylim(0,.03)
# # ax.set_ylim(0,.06)
# # ax.legend(bbox_to_anchor=(1.05, 1),
# #                          loc='upper left', borderaxespad=0.)
# # ax1.legend()
# # # ax2.legend()
# # ax2.legend(bbox_to_anchor=(0, 1.05),
# #                          loc='lower left')
# # ax3.legend()
# # fig.tight_layout()

# fig.savefig(plot_name_2+'.png') # fig 7 in Guildenbecher et al. 2017
# ax1.set_xlabel(r'diameter $(d)$ [$\mathrm{\mu}$m]')
# ax1.set_ylabel(r'$pd(d)$ [1/$\mathrm{\mu}$m]')
# fig1.savefig(plot_name_2+'_LR.png') # fig 9 in Guildenbecher et al. 2017
# ax2.set_xlabel(r'diameter $(d)$ [$\mathrm{\mu}$m]')
# ax2.set_ylabel(r'volume weighted probability [1/$\mathrm{\mu}$m]')
# fig2.savefig(plot_name+'_weighted.png') # fig 10 in Guildenbecher et al. 2017
# ax3.set_xlabel(r'diameter $(d)$ [$\mathrm{\mu}$m]')
# ax3.set_ylabel(r'volume weighted probability [1/$\mathrm{\mu}$m]')
# fig3.savefig(plot_name_2+'_weighted.png')
# ax.set_xscale('log')
# ax.set_yscale('log')
# fig.savefig(plot_name_loglog+'.png')
# ax1.set_xscale('log')
# ax1.set_yscale('log')
# fig1.savefig(plot_name_loglog+'_LR.png')
# # ax2.set_xscale('log')
# ax2.set_yscale('log')
# ax2.set_ylim(bottom=5e-5)
# fig2.savefig(plot_name_semilogy+'_weighted.png')
# fig.savefig(plot_name_loglog+'count.png')

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
# ax2.set_xlabel(r'diameter $(d)$ [$\mathrm{\mu}$m]')
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
# ax.set_xlabel(r'diameter $(d)$ [$\mathrm{\mu}$m]')
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
# ax.set_xlabel(r'diameter $(d)$ [$\mathrm{\mu}$m]')
# ax.set_ylabel(r'$pd(d)$ [1/$\mathrm{\mu}$m]')
# # ax.set_title('Experiment')
# # ax.set_ylim(0,.016)
# ax.set_ylim(0,.05)
# ax.legend()
# # fig.savefig(plot_name_2+'largediam_exp_swap.png')
