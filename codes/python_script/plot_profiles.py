import numpy as np
import os, sys
import matplotlib
import matplotlib.pyplot as plt
from os.path import getsize as getFileSize

data_dir = '/Users/luyu/project/test/run/'
file1 = data_dir+'disc.dat'
arr = np.loadtxt(file1,comments='#',skiprows=1,dtype=None)
iradius = arr[:,0]
r_inner = arr[:,1]
sd_cold = arr[:,2]
sd_star = arr[:,3]
sd_mole = arr[:,4]
sd_atom = arr[:,5]
r_outer = arr[:,6]
r_mid = 0.5*(r_inner+r_outer)
sd_sfr = arr[:,15]/1e12
sd_ofr = arr[:,16]/1e12

################################################
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)

fpdf = 'profile.pdf'
feps = 'profile.eps'

alpha=0.7
ms=10
plt.close()
#fig = plt.figure(figsize(12.5,6))
nrow=1
ncol=1
fig, ((ax11)) = plt.subplots(nrow, ncol, sharex=False, sharey=False)
fig.set_size_inches(6,5)
################################################
ax11.plot(r_mid, sd_star, '-', color='r', alpha=alpha, ms=ms)
ax11.plot(r_mid, sd_cold, '-', color='b', alpha=alpha, ms=ms)
ax11.plot(r_mid, sd_atom, '-', color='g', alpha=alpha, ms=ms)
ax11.plot(r_mid, sd_sfr, '--', color='m',  alpha=alpha, ms=ms)
#ax11.plot(r_mid, sd_ofr/sd_cold, '--', color='m',  alpha=alpha, ms=ms)
#ax11.plot(mstar, mag_nuv, 'bo', alpha=alpha, ms=ms)

ax11.set_xlim(0.1,200)
ax11.set_ylim(0.01, 1000)
ax11.set_xscale('log')
ax11.set_yscale('log')

ax11.set_xlabel(u'r [kpc]')
ax11.set_ylabel(u'Sigma')
################################################
fig.tight_layout()
plt.show()

plt.savefig(fpdf)
plt.savefig(feps, format='eps', dpi=50)