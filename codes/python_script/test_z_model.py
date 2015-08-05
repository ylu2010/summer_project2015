import numpy as np
import scipy as sp
from scipy import special as sf
import mpmath
import matplotlib.pyplot as plt
from matplotlib import cm
from pylab import *
from mpl_toolkits.mplot3d.axes3d import Axes3D

#########################################################################
def metallicity_open_box(eta, lam, mratio, yield0=0.03, freturn=0.42, zsolar=0.014):
    yield1 = yield0/zsolar
    pindex = -lam / (1.-freturn + eta - lam)
    z = yield1 / lam *(1.- (1.+ (1.+(eta - lam)/(1.-freturn)) * mratio)**pindex)
    return z
    
#########################################################################    
ngrid=40
mratio_vec=np.logspace(-3,3, ngrid, base=10)
lgmratio_vec = np.log10(mratio_vec)
eta_vec = np.arange(0., 5, (5.-0.)/ngrid)
lam_vec = np.arange(0., 0.6, (0.6-0)/ngrid)

z_arr = np.zeros((ngrid,ngrid))

#########################################################################
X1, X2 = np.meshgrid(lgmratio_vec, eta_vec)
lam = 0.3
for i in range(ngrid):
    mratio = mratio_vec[i]
    for j in range(ngrid):
        eta = eta_vec[j]
        z = metallicity_open_box(eta, lam, mratio)
        z_arr[j, i] = z
        
lgz_arr1 = np.log10(z_arr)
#########################################################################  
Y1, Y2 = np.meshgrid(lgmratio_vec, lam_vec)
eta = 1.0
for i in range(ngrid):
    mratio = mratio_vec[i]
    for j in range(ngrid):
        lam = lam_vec[j]
        z = metallicity_open_box(eta, lam, mratio)
        z_arr[j, i] = z
        
lgz_arr2 = np.log10(z_arr)
#########################################################################    
fpdf = 'z_model.pdf'
feps = 'z_model.eps'

plt.close()
fig = plt.figure(figsize=(16,6))

#########################################################################
# `ax` is a 3D-aware axis instance because of the projection='3d' keyword argument to add_subplot
ax = fig.add_subplot(1, 2, 1, projection='3d')
#ax = fig.gca(projection='3d')
ax.set_xlabel(r'$log M_*/M_{cold}$', fontsize=18)
ax.set_ylabel(r'$\eta$', fontsize=18)
ax.set_zlabel(r'$log Z$', fontsize=18)
#cm.coolwarm
cc=cm.gist_rainbow
#cc=cm.gray
#cc=cm.winter
p = ax.plot_surface(X1, X2, lgz_arr1,  rstride=1, cstride=1, cmap=cc, linewidth=0.0, alpha=0.7, antialiased=False)
#cb = fig.colorbar(p, shrink=0.5)
cc=cm.gist_rainbow
#p = ax.plot_wireframe(X1, X2, lgz_arr, rstride=5, cstride=5, color='white')
#p = ax.plot_surface(X1, X2, lgz_arr,  rstride=2, cstride=2, cmap=cm.coolwarm, linewidth=0.5, alpha=0.5, antialiased=False)
ax.view_init(20, 60)

#########################################################################
ax = fig.add_subplot(1, 2, 2, projection='3d')
ax.set_xlabel(r'$log M_*/M_{cold}$', fontsize=18)
ax.set_ylabel(r'$\Lambda$', fontsize=18)
ax.set_zlabel(r'$log Z$', fontsize=18)
#cm.coolwarm
#cc=cm.gist_rainbow
#cc=cm.gray
#cc=cm.winter
p = ax.plot_surface(Y1, Y2, lgz_arr2,  rstride=1, cstride=1, cmap=cc, linewidth=0.0, alpha=0.7, antialiased=False)
#cb = fig.colorbar(p, shrink=0.5)
#cc=cm.gist_rainbow

#p = ax.plot_wireframe(X1, X2, lgz_arr, rstride=5, cstride=5, color='white')
#p = ax.plot_surface(X1, X2, lgz_arr,  rstride=2, cstride=2, cmap=cm.coolwarm, linewidth=0.5, alpha=0.5, antialiased=False)
ax.view_init(20, 60)

#########################################################################

plt.savefig(fpdf)
plt.savefig(feps, format='eps', dpi=50)
