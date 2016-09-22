#!/Users/luyu/Library/Enthought/Canopy_64bit/User/bin/python

# plot_mzcold.py

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


rname = 'pr'

mname = 'm12'
rundir = "/Users/luyu/project/disc_metal/data/"+mname+"/ihs_run2/"
#rundir = "/Users/luyu/project/summer_project2015/results/model_ej_lhc/"+mname+'/'
rundir = "/Users/luyu/project/summer_project2015/results/model_"+rname+"_lhc/"+mname+'/'
nfiles = 600
ndim = 4
Hubble_h = 0.7
oxyenyield=0.44
zsolar = 0.0134

##########################################
fname_range = rundir + 'param_range.txt'
data = np.loadtxt(fname_range, comments = '#')
range0 = data[1,]
range1 = data[2,]
range2 = data[3,]
range3 = data[4,]

prior_set_array = [range0, range1, range2, range3, range4, range5, range6]
##########################################
fname_ilhs = rundir+'ILHS.txt'
data = np.loadtxt(fname_ilhs, comments = '#')
p0 = data[:,0]
p1 = data[:,1]
p2 = data[:,2]
p3 = data[:,3]
##########################################
fname_list = rundir+'list.dat'
p = np.zeros(shape=(ndim,nfiles))
data = np.loadtxt(fname_list, comments = '#')
p0 = data[:,1]
p1 = data[:,2]
p2 = data[:,3]
p3 = data[:,4]
p = np.transpose(data[:,1:])
##########################################
#### meaning from main_run_pre_sample.c
#     Par.SNLoadingFactor = params[0];
#     Par.PreheatEntropy = params[1];
#     Par.ReincorporationTimeScale = params[2];
#     Par.SNLoadingFractionToHot = params[3];
#     Par.ZFractionYieldToEject = params[4];
#    Par.ZFractionYieldToHot = params[5];
if rname == 'ej':
    pname0 = r'$\alpha_{LD}$'
    pname1 = r'$\beta_{LD}$'
if rname == 'pr':
    pname0 = r'log $M_{PR}$'
    pname1 = r'$\gamma_{PR}$'
pname2 = r'$\mu_{\lambda}$'
pname3 = r'$\tau_{RI}$'
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
mzstar_vec = np.zeros(nfiles)
mzcold_vec = np.zeros(nfiles)
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
    mzstar_vec[i] = mzst
    mzcold_vec[i] = mzcg
zstar_vec = mzstar_vec/mstar_vec/zsolar
zcold_vec = mzcold_vec/mcold_vec/zsolar
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

lghms_obs = lghm_obs
mhs_obs = 10.**lghm_obs
msup_obs = 10.**lgsmup_obs
mslo_obs = 10.**lgsmlo_obs
################################################
# observation
#;; papastergis
obsfile='/Users/luyu/project_data/obsdata/papastergis_baryon_fraction.dat'
data=np.loadtxt(obsfile, comments='#')
lghm_obs=data[:,0]
lgsm_obs=data[:,1]
lgcbmlo_obs=data[:,2]
lgcbmup_obs=data[:,3]

lghmc_obs = lghm_obs
mhc_obs = 10.**lghm_obs
mcup_obs = 10.**lgcbmup_obs - 10.**lgsm_obs
mclo_obs = 10.**lgcbmlo_obs - 10.**lgsm_obs
lgmcup_obs = np.log10(mcup_obs)
lgmclo_obs = np.log10(mclo_obs)
################################################
#;;; Tremonti et al. 2004 (fitting model from Kewley & Ellison 2008 and fitting coefficients from Peeples & Shankar (2011) 
lgohsun = 8.69 #;; Allende Prieto et al. 2001; ;Asplund et al. 2009
lgmstar_obs=np.arange(5,11.5,0.1)
lgoh12_tremonti = -0.759210 + 1.30177  * lgmstar_obs + 0.003261 * lgmstar_obs**2 - 0.00364112 * lgmstar_obs**3

mstar_tremonti = 10.**lgmstar_obs
lgmstar_tremonti=lgmstar_obs
lgzcold_tremonti = lgoh12_tremonti-lgohsun
zcold_tremonti=10.**(lgzcold_tremonti)

#;;; Andrews & Martini (2013)
gam = 0.64
lgmt0 = 8.901
lgohsun_andrews=8.89
lgoh12_andrews = 8.798 - np.log10(1.+ 10.**(gam * (lgmt0 - lgmstar_obs)))

mstar_andrews = 10.**lgmstar_obs
lgmstar_andrews=lgmstar_obs
lgzcold_andrews = lgoh12_andrews-lgohsun #-lgohsun_andrews

lgmszcold_obs = lgmstar_obs
lgzcold_obs = lgzcold_andrews
lgzcoldup_obs = lgzcold_andrews + 0.25
lgzcoldlo_obs = lgzcold_andrews - 0.25

mszcold_obs = 10.**lgmszcold_obs
zcold_obs=10.**(lgzcold_obs)
zcoldup_obs = 10.**lgzcoldup_obs
zcoldlo_obs = 10.**lgzcoldlo_obs

################################################
# selection
if mname == 'm10':
    mstar_min, mstar_max = [1e5, 1e6]
    mcold_min, mcold_max = [1e6, 1e7]
    mhot_min, mhot_max = [0, 0.17*1e10]

if mname == 'm11':
    m1 = 10.**interp1d(lghms_obs, lgsmup_obs)(11.)
    m2 =  10.**interp1d(lghms_obs, lgsmlo_obs)(11.)
    mstar_min, mstar_max = m2, m1
    #mstar_min, mstar_max = [10**(8.75)/3, 10**(8.75)*3]
    m1 = 10.**interp1d(lghm_obs, lgmclo_obs)(11.)
    m2 =  10.**interp1d(lghm_obs, lgmcup_obs)(11.)
    mcold_min, mcold_max = [m1, m2]
    
    mzcold_min, mzcold_max = [1e6, 2e7]
    mhot_min, mhot_max = [0, 0.17*1e11]
    mzstar_min, mzstar_max = [8e5, 4e7]
    m1 = 10.**interp1d(lgmszcold_obs, lgzcoldlo_obs)(np.log10(mstar_min))
    m2 = 10.**interp1d(lgmszcold_obs, lgzcoldup_obs)(np.log10(mstar_max))
    zcold_min, zcold_max = [m1, m2]

if mname == 'm12':
    m1 = 10.**interp1d(lghms_obs, lgsmup_obs)(12.)
    m2 =  10.**interp1d(lghms_obs, lgsmlo_obs)(12.)
    mstar_min, mstar_max = m2, m1
    
    m1 = 10.**interp1d(lghmc_obs, lgmclo_obs)(12.)
    m2 =  10.**interp1d(lghmc_obs, lgmcup_obs)(12.)
    mcold_min, mcold_max = [m1, m2]
    
    mzcold_min, mzcold_max = [1e6, 2e7]
    mhot_min, mhot_max = [0, 1e11]
    mzstar_min, mzstar_max = [1e8, 4e8]
    mzhot_min, mzhot_max = [2e7, 1e9]
    m1 = 10.**interp1d(lgmszcold_obs, lgzcoldlo_obs)(np.log10(mstar_min))
    m2 = 10.**interp1d(lgmszcold_obs, lgzcoldup_obs)(np.log10(mstar_max))
    zcold_min, zcold_max = [m1, m2]
    

x_vec = mstar_vec
x_min, x_max = mstar_min, mstar_max
x_range = [1e7, 1e11]
y_vec = zcold_vec
y_min, y_max = zcold_min, zcold_max
y_range = [1e-2, 10]

bool_mstar = (mstar_vec >= mstar_min) & (mstar_vec <= mstar_max)
#bool_zstar = (interp1d(lgms_vec, lgzslo_vec)(np.log10(mstar_vec)) < np.log10(zstar_vec)) & (interp1d(lgms_vec, lgzsup_vec)(np.log10(mstar_vec)) > np.log10(zstar_vec))
bool_zcold = (interp1d(lgmszcold_obs, lgzcoldlo_obs, fill_value='extrapolate')(np.log10(mstar_vec)) < np.log10(zcold_vec)) & (interp1d(lgmszcold_obs, lgzcoldup_obs, fill_value='extrapolate')(np.log10(mstar_vec)) > np.log10(zcold_vec))
bool_mcold = (mcold_vec >= mcold_min) & (mcold_vec <= mcold_max)
bool_mzcold = (mzcold_vec >= mzcold_min) & (mzcold_vec <= mzcold_max)
bool_mhot = (mhot_vec >= mhot_min) & (mhot_vec <= mhot_max)
#bool_mzhot = (mzhot_vec >= mzhot_min) & (mzhot_vec <= mzhot_max)
#isel = np.where(bool_mstar & bool_mcold & bool_mhot)[0]  
isel = np.where(bool_mstar & bool_zcold)[0] 
isel_mzcold = isel

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
    #bin_centers[i,:] =  bin_edges[i,:]+
        
    
################################################
#ndim = 6
fpdf = rundir+'/figs/test_mzcold_'+mname+'.pdf'
################################################
font = {'family' : 'serif',
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
box.plot(x_range, [y_min, y_min], '--', color='green')
box.plot(x_range, [y_max, y_max], '--', color='green')
box.plot([x_min, x_min], y_range, '--', color='green')
box.plot([x_max, x_max], y_range, '--', color='green')

#box.plot(mstar_tremonti, zcold_tremonti, '-', color='green')
#box.plot(mstar_andrews, zcold_andrews, '--', color='green')
box.fill_between(mszcold_obs, zcoldlo_obs, zcoldup_obs, color='green', alpha = 0.3)
#box.fill_between(ms_vec, zclo_vec, zcup_vec, color='m', alpha = 0.3)

#box.get_xaxis().set_visible(True)
box.set_xlabel(r'$M_{*}  [M_{\odot}]$')
box.set_ylabel(r'$Z_{gas}/Z_{\odot}$')
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
