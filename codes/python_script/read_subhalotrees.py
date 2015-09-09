#!/usr/bin/env python

import numpy as np
import os, sys
from os.path import getsize as getFileSize


def read_subhalotrees(tree_name, fnr, ndig=3):

    # The input galaxy structure:
    Halodesc_full = [
        ('SnapNum'                      , np.int32),                    # 
        ('ID'                           , np.int32),                    # 
        ('Descendant'                   , np.int32),                    #
        ('CountProgenitor'              , np.int32),                    #
        ('FirstProgenitor'              , np.int32),                    #
        ('NextProgenitor'               , np.int32),                    #
        ('FirstHaloInGroup'             , np.int32),                    #
        ('NextHaloInGroup'              , np.int32),                    #
        ('Mvir'                         , np.float32),                  # 
        ('Rvir'                         , np.float32),                  # all radii in Mpc/h (physical)
        ('Concen'                       , np.float32),                  # concentration
        ('Pos'                          , (np.float32, 3)),             # all positions in Mpc/h
        ('Vel'                          , (np.float32, 3))              # all velocities in km/s
        ]

    names = [Halodesc_full[i][0] for i in xrange(len(Halodesc_full))]
    formats = [Halodesc_full[i][1] for i in xrange(len(Halodesc_full))]
    Halodesc = np.dtype({'names':names, 'formats':formats}, align=True)

    print "Determining array storage requirements"
                
    # Read the file and determine the total number of halos to be read in
    fname = tree_name+'_'+str(fnr).zfill(ndig)  # Complete filename
        
    if getFileSize(fname) == 0:
        print "File\t%s  \tis empty! " % (fname)
        
    fin = open(fname, 'rb')  # Open the file
    Ntrees = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file
    Nhalos = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of halos in file.

    print
    print "Input files contain:\t%d trees ;\t%d halos ." % (Ntrees, Nhalos)
    print

    buf = np.fromfile(fin,np.dtype(np.int32),1)   # skipping a buffer by reading it 
        
    HalosPerTree = np.fromfile(fin, np.dtype((np.int32, Ntrees)),1) # Read the number of halos in each tree
    #print "HalosPerTree=", HalosPerTree
    print "Reading N=", Nhalos, "   \thalos from file: ", fname
    HH = np.fromfile(fin, Halodesc, Nhalos)  # Read in the halo structures
    fin.close()  # Close the file


    # Convert the Halo array into a recarray
    HH = HH.view(np.recarray)

    return HH, HalosPerTree
    
#tree_name='/Users/luyu/project_data/mergertrees/MilkyWay/trees/Tree_sub'

#H, HalosPerTree = read_subhalotrees_single(tree_name, 188)
