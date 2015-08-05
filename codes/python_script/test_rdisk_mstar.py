import numpy as np
import os, sys
import matplotlib
import matplotlib.pyplot as plt
from os.path import getsize as getFileSize

# test_rdisk_mstar.py

zsolar=0.0134
rundir='/Users/luyu/project/metalgrad/data/results/'

model_name = 'EJ'

dirname = rundir + '/edit10_'+model_name

filename = dirname + '/sample_z0.0.dat'

data = np.loadtxt(filename, comments='#')

mstar_ej = data[:,3]
rstar_ej = data[:,4]
################################################
model_name = 'PR'

dirname = rundir + '/edit10_'+model_name

filename = dirname + '/sample_z0.0.dat'

data = np.loadtxt(filename, comments='#')

mstar_pr = data[:,3]
rstar_pr = data[:,4]
################################################
model_name = 'RI'

dirname = rundir + '/edit10_'+model_name

filename = dirname + '/sample_z0.0.dat'

data = np.loadtxt(filename, comments='#')

mstar_ri = data[:,3]
rstar_ri = data[:,4]
################################################
ngrid=10
ms_vec=np.logspace(8,12, ngrid, base=10)
rd_z0_dutton11=10.**0.72*(ms_vec/10**10.44)**0.18*(0.5+0.5*(ms_vec/10**10.44)**1.8)**((0.52-0.18)/1.8)
################################################
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)

fpdf = 'rdisk.pdf'
feps = 'rdisk.eps'
 
alpha=0.7
ms=10
plt.close()
nrow=1
ncol=1
fig, ax = plt.subplots(nrow, ncol, sharex=False, sharey=False)
fig.set_size_inches(6,5)
################################################
ax.plot(mstar_ej, rstar_ej, 'ro')
ax.plot(mstar_pr, rstar_pr, 'bo')
ax.plot(mstar_ri, rstar_ri, 'go')
ax.plot(ms_vec, rd_z0_dutton11, '-', color='m', alpha=alpha, ms=ms)


ax.set_xlim(1e8,1e12)
ax.set_ylim(5e-1,30.)
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel(u'$M_* [M_\u2609]$')
ax.set_ylabel('$r_d$ [kpc]')
################################################
fig.tight_layout()
plt.show()

plt.savefig(fpdf)
plt.savefig(feps, format='eps', dpi=50)
