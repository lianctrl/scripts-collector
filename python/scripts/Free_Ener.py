#!/usr/bin/env python
"""
Author:Luca Sagresti

Date: V1.0 04/01/20

Details:Script computing Free energy surface for the NAR ion starting from its positional distribution along the z axis

"""
######################################
#               USAGE
# python script.py file.pdb file.xtc bin
#
#    SCRIPT CURRENT VERSION IS V1.0
######################################
import sys; from sys import argv

import numpy as np

from scipy.interpolate import CubicSpline

import matplotlib.pyplot as plt

import matplotlib as mpl

import MDAnalysis as mda

import datetime

#mpl.rcParams['figure.dpi']=100
#mpl.rcParams['figure.titlesize']=20
#mpl.rcParams['axes.facecolor']='white'
#mpl.rcParams['lines.linewidth']=2.0
#mpl.rcParams['axes.linewidth']=2.0
#mpl.rcParams['xtick.major.pad']=6
#mpl.rcParams['ytick.major.pad']=6
#mpl.rcParams['xtick.labelsize']=14
#mpl.rcParams['ytick.labelsize']=14
#mpl.rcParams['axes.titlesize']=18
#mpl.rcParams['axes.labelsize']=18
#mpl.rcParams['axes.grid']='True'
#mpl.rcParams['axes.axisbelow']='line'
#mpl.rcParams['legend.fontsize']=12


#error if argv less than 3 arguments
if len(sys.argv) != 4:
    sys.exit("\nYou did not provide the correct number\
of arguments in the command line.\n\
USAGE IS: python script.py file.pdb file.xtc bin")

# File reading

# Path to the pdb file
PDB = argv[1]

# Path to the xtc file
XTC = argv[2]

# Loading info
u = mda.Universe(PDB, XTC)

ref_atom = u.select_atoms("resname NAR")

z_coord = np.zeros((len(ref_atom),len(u.trajectory)))

box_lenghts = u.dimensions

z_dimension = box_lenghts[2]

m=0
k=0

for frm in u.trajectory[:]:

	if ( ( m % int(len(u.trajectory)/10)  ) == 0 ) :

        	print(" %i / 10 of file loaded- %s" %(k, datetime.datetime.now()))
        	k+=1

	z_coord[:,m] = ref_atom.positions[:,2]

	m+=1

md_z_inf = np.min(z_coord)

md_z_sup = np.max(z_coord)

bins = np.int(argv[3])

z_binning = np.linspace(md_z_inf,md_z_sup,bins+1)

bin_width = (md_z_sup - md_z_inf) / bins

print("Calculating FES for %d bins of width %3.2f" %(bins,bin_width))
print("Range of z coordinate obtained in MD:  [ %.3f  :  %.3f ]" % (md_z_inf, md_z_sup))

# Create an histogram where pile up the z positions and normalize the histogram 
hist, _ = np.histogram(z_coord, bins=z_binning, density=True)

#binning centers
z_free_ener = (z_binning[1:] + z_binning[:-1]) / 2

# simulation at T=298 K 
KbT = 0.0083144621 * 298   #KJoule/mol

delta_bin = z_binning[1] - z_binning[0] #equal to binwidth! 

F_z = -KbT * ( np.log( hist / delta_bin ))

# Put the minimum at 0 energy
#F_z2 = F_z - abs(np.amin(F_z))

# compute spline of free energy
cs = CubicSpline(z_free_ener, F_z)

# print to file the results

plot_data = np.array([z_free_ener,F_z])

np.savetxt("Free_Energy.dat", plot_data.T)

# Starting to plot
#fig = plt.figure()
#
#xs = np.arange(md_z_inf, md_z_sup, 0.1)
#ys = cs(xs) + abs(np.amin(cs(xs)))
#
#ax0 = fig.add_subplot(111) 
#ax0.plot(z_free_ener, F_z2, label="Discrete")  # Free energy from MD, discrete points
#ax0.plot(xs, ys, label="Spline")  # Free energy from MD, continue view
#
#ax0.set_title('Free Energy')
#ax0.legend(loc='best')
#ax0.set_ylim([-np.max(free_ener)*0.1,np.max(free_ener)*1.5])
#
#ax0.set_xlabel('Z coordinate')
#ax0.set_ylabel('Free Energy [KJoule/mol]')
#
#plt.show()
