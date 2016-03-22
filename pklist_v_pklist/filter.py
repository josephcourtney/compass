#!/usr/bin/env python
#-*- coding:utf-8 -*-

# Import peaks from SPARKY peaklist with columns:
# Assignment    w1  w2
# filter the peaks to exclude peaks that are:
#   not within the aliphatic region (10 ppm, 80 ppm)
#   not matched on either side of the diagonal withint cross_tol
#   within diag_tol of the diagonal
# and export them as a csv file
# Usage:
#   filter.py sparky_peaks.list > filtered.pks


import sys
import re
import numpy as np
from scipy.spatial import cKDTree

def aliph_filt(pks):
    filt_pks = []
    for pk in pks:
        pk_1, pk_2 = pk
        if pk_1 > 10 and pk_1 < 80 and pk_2 > 10 and pk_2 < 80:
            filt_pks.append((pk_1, pk_2))
    return np.array(filt_pks)

def diag_filt(pks, diag_tol = 0.5):
    filt_pks = []
    for pk in pks:
        pk_1, pk_2 = pk
        if np.abs(pk_1-pk_2) > diag_tol:
            filt_pks.append((pk_1, pk_2))
    return np.array(filt_pks)

def cross_filt(pks, cross_tol = 0.4):
    pktree = cKDTree(pks)
    filt_pks = []
    for pk in pks:
        d, i = pktree.query(pk[::-1])
        if d < cross_tol:
            if tuple(pk) not in filt_pks:
                filt_pks.append(tuple(pk))
            if tuple(pks[i]) not in filt_pks:
                filt_pks.append(tuple(pks[i]))
    return np.array(filt_pks)


pks = []
with open(sys.argv[1], 'r') as f:
    for line in f.readlines()[2:]:
        mtch = re.match("\s+\?\-\?\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
        _pk_1, _pk_2 = mtch.group(1,2)
        pk_1, pk_2 = float(_pk_1), float(_pk_2)
        pks.append((pk_1, pk_2))
pks = np.array(pks)
print len(pks)


filt_pks = pks
filt_pks = aliph_filt(pks)
filt_pks = diag_filt(filt_pks)
filt_pks = cross_filt(filt_pks)

with open(sys.argv[2],'w') as f:
    for pk in filt_pks:
        f.write('%.3f,%.3f\n'%(pk[0], pk[1]))
print len(filt_pks)
