#!/usr/bin/env python

# extract_host_mah.py

import numpy as np
import os, sys, shutil
import math
import matplotlib
import matplotlib.pyplot as plt
from os.path import getsize as getFileSize

import read_subhalotrees
reload(read_subhalotrees)

readit = 1
tree_dir= '/Users/luyu/project/milkyways/trees/'
#tree_dir = '/Users/luyu/project/milkyways/trees/trees_chinL125_2048/'

# read snaphsot_alist file
snapshot_file=tree_dir+'Snapshot_alist'
alist =np.array(np.genfromtxt(snapshot_file, dtype=None))

NSnapshots = alist.__len__()
zlist = 1./alist -1.
snapnum0 = NSnapshots - 1

Hubble_h = 0.7
tree_dir1 = 'first_trees/'
mah_dir1 = 'first_trees_mah/'
#tree_dir1 = 'trees_chinL125_2048_m12.1_single_sequence/'
#mah_dir1 = 'trees_chinL125_2048_m12.1_single_sequence_mah/'

inputdir = tree_dir + tree_dir1
fin_base = tree_dir+tree_dir1+'/Tree_first'
#fin_base = tree_dir+tree_dir1+'/Tree_sel'

plotdir = tree_dir + mah_dir1
ndig = 3

if readit :
    nsnap = NSnapshots
    nfile = len([name for name in os.listdir(inputdir) if os.path.isfile(os.path.join(inputdir, name))])
    
    lgmah = np.zeros(shape=(nfile, nsnap))

    for iFile in range(nfile):
        mahfile = tree_dir + mah_dir1 + '/mah_'+str(iFile).zfill(ndig)+'.txt'
        halofile = tree_dir + mah_dir1 + '/halo_'+str(iFile).zfill(ndig)+'.txt'
        fmah = open(mahfile, 'w')
        fhalo = open(halofile, 'w')
        
        HalosPerTree = []
        H , HalosPerTree = read_subhalotrees.read_subhalotrees(fin_base, iFile, ndig=ndig)

        # going back along the main branch
        k = 0L
        ih = 0
        j = ih
        while H.FirstProgenitor[j] > 0:
            #print j, H.FirstHaloInGroup[j]
            #if H.FirstHaloInGroup[j] == j:
            lgmah[iFile, NSnapshots-k-1] = np.log10(H.Mvir[j]/Hubble_h)+10.
            fmah.writelines("%d %f %g\n" % (NSnapshots-k-1, zlist[NSnapshots-k-1], H.Mvir[j]/Hubble_h*1e10))

            j = H.FirstProgenitor[j]
            k+=1
        fhalo.writelines("%g %g\n" % (H.Mvir[ih]/Hubble_h*1e10, H.Concen[ih]))
        fhalo.close()
        fmah.close()

vec_snapshot = np.arange(NSnapshots)+1

################################################
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)

fpdf = plotdir+'mah.pdf'
feps = plotdir+'mah.eps'
xspread = 0.3+0.8*np.random.rand(nfile)*0.0
yspread = 2+3.*np.random.rand(nfile)
alpha=0.7
ms=15
plt.close()
#fig = plt.figure(figsize(12.5,6))
nrow=1
ncol=1
fig, ((ax11)) = plt.subplots(nrow, ncol, sharex=False, sharey=False)
fig.set_size_inches(5,5)
################################################

for i in range(nfile):
    ax11.plot(zlist+1, lgmah[i,:])
    
ax11.set_xscale('log')
ax11.set_xlim([1,14])
ax11.set_ylim([7,13])
ax11.set_ylabel(u'log $M_{host}$    [$M_\u2609$]')
ax11.set_xlabel('$1+z$')

################################################
#fig.subplots_adjust(hspace=0.15, wspace=0.3)
fig.tight_layout()
plt.show()

plt.savefig(fpdf)
plt.savefig(feps, format='eps', dpi=50)
