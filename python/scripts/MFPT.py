#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
V1.0 05/04/2020 Compute the mean first passage time for
                a defined barrier

Authors: L. Sagresti
"""

import numpy as np
import matplotlib.pyplot as plt; import matplotlib as mpl
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
from scipy.linalg import logm
import datetime
import sys; from sys import argv
import MDAnalysis as mda


#error if argv less than 10 arguments
if len(sys.argv) != 5:
    sys.exit("\nYou did not provide the correct number\
of arguments in the command line.\nUSAGE IS:\n\
python3 script.py file.pdb file.xtc start_pos end_pos")

##########################################################

# File reading

# Path to the pdb file
PDB = argv[1]

# Path to the wrapped xtc
XTC = argv[2]

u = mda.Universe(PDB, XTC)

z_start = np.float(argv[3])

z_end   = np.float(argv[4])

#if (z_start<z_end):
#
#    print("Script written to assign start position the greater value: \
#swap between end and start positions. Run continue...")
#
#    z_start, z_end = z_end, z_start

delta_t = u.trajectory[1].time-u.trajectory[0].time

delta_z = 0.5

NAR_ion = u.select_atoms("resname NAR")

z_coordinate_NAR = np.zeros(len(u.trajectory))

m=0
k=0

# load trajectory
for frm in u.trajectory[:]:

    if ( ( m % int(len(u.trajectory)/10)  ) == 0 ) :

        print(" %i / 10 of file loaded- %s" %(k, datetime.datetime.now()))

        k+=1

    z_coordinate_NAR[m] = NAR_ion.positions[:,2]

    m+=1

def fitfunction(t, tau, A):

    # tau is in ps only if t is

    initial_size = A

    y = initial_size * np.exp(-t / tau)

    return y

def non_increasing(L):

    return all(x>=y for x, y in zip(L, L[1:]))

def non_decreasing(L):

    return all(x<=y for x, y in zip(L, L[1:]))

#Compute MFPT from z_start to z_end FROM MD

mfpt  = []

idz   = 0

for i in range(0,len(u.trajectory)):

    if i < idz:

        continue

    if z_start - delta_z < z_coordinate_NAR[i] < z_start + delta_z:

        for j in range(i+1,len(u.trajectory)):

            if z_coordinate_NAR[j] > (z_start + delta_z):

                break 

            if z_end - delta_z < z_coordinate_NAR[j] < z_end + delta_z:

                mfpt.append(j - i)

                idz = j

                break

mfpt = np.array([t * delta_t for t in mfpt])

t_nbins_fit = 20 #to check

max_t = np.amax(mfpt)

min_t = np.amin(mfpt)

#t_bins = np.linspace(0.1,max_t,t_nbins_fit)
#
#hist, _ = np.histogram(mfpt, bins=t_bins, density=True)
#
#init_vals = [1.0, hist[0]]
#
#best_vals, covar = curve_fit(fitfunction, (t_bins[1:]+t_bins[:-1])/2, hist, p0=init_vals)

print("max time obtained for passing the barrier: %8.3f ps \n" %(max_t))
print("min time obtained for passing the barrier: %8.3f ps \n" %(min_t))

#comment here below to mute the plot part
#Plot

#fig = plt.figure()
#
#axes = fig.add_subplot(121)
#
#axes.set_yscale('log')
#axes.bar(t_bins[:-1], hist, alpha = 0.5, width=t_bins[0]-t_bins[1] ) #, label = "$P_{MD}(x)$ " + str(interval) )
#axes.plot((t_bins[1:]+t_bins[:-1])/2, fitfunction((t_bins[1:]+t_bins[:-1])/2, *best_vals),'r')
#axes.set_title('MFPT- from (%.3f Ang) to (%.3f Ang)' %(z_start,z_end))
#axes.set_ylabel('probability')
#axes.set_xlabel('time [ps]')

#fig.savefig('MFPT-MD.jpg')
print("We had %d events from passing the barrier, with mean time: %.3f ps +/- %.3f ps" %(len(mfpt), np.mean(mfpt),np.std(mfpt)))
#print("We had %d events from passing the barrier, with mean time: %.3f ps \n tau fitted: %.2f +/- %.2f ps" %(len(mfpt), np.mean(mfpt),best_vals[0],np.sqrt(covar[0,0])))
