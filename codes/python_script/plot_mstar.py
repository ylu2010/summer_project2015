#!/Users/luyu/Library/Enthought/Canopy_64bit/User/bin/python

# plot_mstar.py

import numpy as np
import os, sys
import math
#import pylab as plt
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from os.path import getsize as getFileSize
import os.path
from scipy import stats
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
from matplotlib.ticker import MaxNLocator


mname = 'm12'
rundir = "/Users/luyu/project/disc_metal/data/"+mname+"/ihs_run2/"
rundir = "/Users/luyu/project/summer_project2015/results/model_ej_lhc/"+mname+'/'
nfiles = 300
ndim = 3
Hubble_h =0.7

##########################################
fname_range = rundir + 'param_range.txt'
data = np.loadtxt(fname_range, comments = '#')
range0 = data[1,]
range1 = data[2,]
range2 = data[3,]


prior_set_array = [range0, range1, range2]

##########################################
fname_ilhs = rundir+'ILHS.txt'
data = np.loadtxt(fname_ilhs, comments = '#')
p0 = data[:,0]
p1 = data[:,1]
p2 = data[:,2]

##########################################
fname_list = rundir+'list.dat'
p = np.zeros(shape=(ndim,nfiles))
data = np.loadtxt(fname_list, comments = '#')
p0 = data[:,1]
p1 = data[:,2]
p2 = data[:,3]

p = np.transpose(data[:,1:])
##########################################
#### meaning from main_run_pre_sample.c
#     Par.SNLoadingFactor = params[0];
#     Par.PreheatEntropy = params[1];
#     Par.ReincorporationTimeScale = params[2];
#     Par.SNLoadingFractionToHot = params[3];
#     Par.ZFractionYieldToEject = params[4];
#    Par.ZFractionYieldToHot = params[5];
pname0 = r'$\alpha_{LD}$'
pname1 = r'$\beta_{LD}$'
pname2 = r'$\mu_{\lambda}$'
pname3 = r'$\alpha_{LD, HOT}$'
pname4 = r'$\eta_{Z, EJ}$'
pname5 = r'$\eta_{Z, HOT}$'
pname6 = r'$\beta_{LD}$'

param_names = [pname0, pname1, pname2, pname3, pname4, pname5, pname6]

#######################################
filebase = 'sample_z0.0_m'


mhalo_vec = np.zeros(nfiles)
mstar_vec = np.zeros(nfiles)
mcold_vec = np.zeros(nfiles)
mhot_vec = np.zeros(nfiles)
for i in range(nfiles):
    fname = rundir+filebase+str(i)+'.dat'
    #print i, fname
    line = np.loadtxt(fname, comments = '#')
    #print line
    mh = line[0]
    mst = line[3]
    mcg = line[5]
    mca = line[6]
    sfr = line[7]
    mho = line[8]
    mej = line[9]
    mzho = line[10]
    mzcg = line[11]
    mzst = line[12]
    mzej = line[13]
    mhalo_vec[i] = mh
    mstar_vec[i] = mst
    mcold_vec[i] = mcg
    mhot_vec[i] = mho
################################################
# observation
# ;; Behroozi 
zlabel='z0.10'
obsfile='/Users/luyu/project_data/obsdata/behroozi/release-sfh_z0_z8_052913/smmr/c_smmr_'+zlabel+'_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat'
data=np.loadtxt(obsfile, comments='#')
lghm_beh=data[:,0]
lgsm_beh=data[:,1]+lghm_beh
lgsmup_beh=lgsm_beh+data[:,2]
lgsmlo_beh=lgsm_beh-data[:,3]


mh_beh=10.**lghm_beh
ms_beh=10.**lgsm_beh
msup_beh=10.**lgsmup_beh
mslo_beh=10.**lgsmlo_beh
########
#;; Lv 2013
obsfile='/Users/luyu/project_data/obsdata/MW_lv/xy_m3.dat'
data=np.loadtxt(obsfile, comments='#')
lghm_lv=data[:,0]-np.log10(Hubble_h)
lgsm_lv=data[:,1]+lghm_lv
lgsmlo_lv=data[:,2]+lghm_lv
lgsmup_lv=data[:,3]+lghm_lv

mh_lv=10.**lghm_lv
ms_lv=10.**lgsm_lv
mslo_lv=10.**lgsmlo_lv
msup_lv=10.**lgsmup_lv
#########
# combine the constraints
lghm_obs = lghm_beh[3:32]
lgsm1 = interp1d(lghm_beh, lgsmlo_beh)(lghm_obs)
lgsm2 = interp1d(lghm_lv, lgsmlo_lv)(lghm_obs)
lgsmlo_obs = np.min([lgsm1, lgsm2], axis=0)

lgsm1 = interp1d(lghm_beh, lgsmup_beh)(lghm_obs)
lgsm2 = interp1d(lghm_lv, lgsmup_lv)(lghm_obs)
lgsmup_obs = np.max([lgsm1, lgsm2], axis=0)

mh_obs = 10.**lghm_obs
msup_obs = 10.**lgsmup_obs
mslo_obs = 10.**lgsmlo_obs
################################################
# selection
if mname == 'm10':
    mstar_min, mstar_max = [1e5, 1e6]
    mcold_min, mcold_max = [1e6, 1e7]
    mhot_min, mhot_max = [0, 0.17*1e10]

if mname == 'm11':
    m1 = 10.**interp1d(lghm_obs, lgsmup_obs)(11.)
    m2 =  10.**interp1d(lghm_obs, lgsmlo_obs)(11.)
    mstar_min, mstar_max = m2, m1
    #mstar_min, mstar_max = [3e8, 1.1e9]
    mcold_min, mcold_max = [7e8, 1.5e9]
    mhot_min, mhot_max = [0, 0.17*1e11]

if mname == 'm12':
    m1 = 10.**interp1d(lghm_obs, lgsmup_obs)(12.)
    m2 =  10.**interp1d(lghm_obs, lgsmlo_obs)(12.)
    #mstar_min, mstar_max = [1.9e10, 6e10]
    mstar_min, mstar_max = m2, m1
    mcold_min, mcold_max = [3.2e9, 1.1e10]
    mhot_min, mhot_max = [0, 1e11]
    
x_vec = mhalo_vec
x_range = [1e10, 1e13]
y_vec = mstar_vec
y_min, y_max = mstar_min, mstar_max
y_range = [1e8, 1e11]

#isel = np.where(mstar_vec >= mstar_min) & (mstar_vec <= mstar_max))[0]   
#isel = np.where((mcold_vec >= mcold_min) & (mcold_vec <= mcold_max))[0] 
#isel = np.where((mhot_vec >= mhot_min) & (mhot_vec <= mhot_max))[0]   
isel = np.where((y_vec >= y_min) & (y_vec <= y_max))[0]   

if len(isel) > 1:
    mhalo_sel = mhalo_vec[isel]
    mstar_sel = mstar_vec[isel]
    mcold_sel = mcold_vec[isel]
    mhot_sel = mhot_vec[isel]
    p0sel = p0[isel]
    p1sel = p1[isel]
    p2sel = p2[isel]
    p3sel = p3[isel]
    p4sel = p4[isel]
    p5sel = p5[isel]
else:
    print 'no model falls into the selection!\n'
    exit
################################################
# make distributions
bins = 15
bin_edges = np.zeros(shape=(ndim, bins))
hists=np.zeros(shape = (ndim, bins))
binsizes = np.zeros(ndim)
for i in range(ndim):
    hist0, bin_edge0 = np.histogram(p[i,isel], bins=bins, range=prior_set_array[i], density=True)
    bin_edges[i,:] = bin_edge0[:-1]
    hists[i,:] = hist0 * len(isel)/nfiles
    binsizes[i] = bin_edge0[1]-bin_edge0[0]
    
################################################
#ndim = 6
fpdf = rundir+'/figs/test_mstar_'+mname+'.pdf'
################################################
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)

alpha=0.5
ms=4
plt.close()
plt.figure(figsize(12,12))
nrow=ndim
ncol=ndim
################################################
color='blue'
levels = [0.33, 0.67, 0.95]
################################################
for i in range(ndim):
    ip = i*nrow + i+1

    ax = plt.subplot(ncol, nrow, ip)
    if i < ndim-1:
        ax.get_xaxis().set_visible(False)
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()

    #a = den1d[:,i]
    #b = function_array_scale_to_cdf.function_array_scale_to_cdf(a)
    #cset = ax.plot(x1d[:,i], den1d[:,i],'g-')
    #cset = ax.plot(x1d[:,i], b,'g-')
    
    xmin, xmax = prior_set_array[i]
    
    
    ax.set_xlim(xmin, xmax)
    #ax.set_ylim(0.0, 1.1*ymax)
    plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))
    plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
    plt.locator_params(nbins=3)

for i in range(nrow):
    ip = i*nrow + i +1
    ax = plt.subplot(ncol,nrow,ip)
    ax.bar(bin_edges[i,:], hists[i,:], width = binsizes[i], color='red', alpha=alpha)
    xmin, xmax = prior_set_array[i]
    ax.set_ylim(0, 1./(xmax-xmin))
    #ax.plot(bin_edges[i,:], hists[i,:], '-', color='red')
    #ax.set_ylim(0,max(hists[i,:]))
    
    if i== nrow-1:
        ax.set_xlabel(param_names[i])
    for j in range(ncol):
        ip = i*nrow + j + 1

        if j<i:
            ax = plt.subplot(ncol,nrow,ip)
            if i < nrow-1:
                ax.get_xaxis().set_visible(False)
            if j > 0:
                ax.get_yaxis().set_visible(False)
            #a = den2d[:,:,i,j]
            #b = function_array_scale_to_cdf.function_array_scale_to_cdf(a)
            #cset = ax.contour(b, extent=extent2d[:,i,j], linewidths = 2, colors=color, levels = levels)
            #ax.plot(p_maxlikelihood[j], p_maxlikelihood[i], "*", color='red', ms=ms)
            ax.plot(p[j,:], p[i,:], 'o', ms = ms, alpha=alpha, color='blue')
            ax.plot(p[j,isel], p[i,isel], 'o', color='red', ms = ms, alpha = alpha)
            xmin, xmax = prior_set_array[j]
            ymin, ymax = prior_set_array[i]
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_xlabel(param_names[j])
            ax.set_ylabel(param_names[i])
            plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))
            plt.gca().yaxis.set_major_locator(MaxNLocator(prune='upper'))
            plt.locator_params(nbins=3)


################################################
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0, wspace=0)
################################################
#axz = plt.subplot(ncol, nrow, ncol-1)

left, width = 0.6, 0.3
bottom, height = 0.6, 0.3

square_box = [left, bottom, width, height]

box = plt.axes(square_box)

#zoom_effect02(box, ax_zoom, color='grey')
################################

box.plot(x_vec, y_vec, 'o', color=color, ms=ms, alpha = alpha)
box.plot(x_vec[isel], y_vec[isel], 'o', color = 'red', ms=ms, alpha=alpha)
box.plot(x_range, [y_min, y_min], '--', color='m')
box.plot(x_range, [y_max, y_max], '--', color='m')

box.fill_between(mh_obs, mslo_obs, msup_obs, color='m', alpha = 0.3)
#box.fill_between(mh_beh, mslo_beh, msup_beh, color='m', alpha = 0.3)
#box.fill_between(mh_lv, mslo_lv, msup_lv, color='g', alpha = 0.3)
#box.plot(mh_beh, msup_beh, ':', color='m')
#box.plot(mh_beh, ms_beh, ':', color='m')
#box.plot(mh_beh, mslo_beh, ':', color='m')


#box.get_xaxis().set_visible(True)
box.set_xlabel(r'$M_{halo}  [M_{\odot}]$')
box.set_ylabel(r'$M_{*}  [M_{\odot}]$')
box.set_xlim(x_range)
box.set_xscale('log')
box.set_ylim(y_range)
box.set_yscale('log')
########

#plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))
#plt.gca().yaxis.set_major_locator(MaxNLocator(prune='upper'))
#plt.locator_params(nbins=3)



################################################
plt.show()

plt.savefig(fpdf)
