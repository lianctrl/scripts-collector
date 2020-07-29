#!/usr/bin/env python
"""
Author:Luca Sagresti

Date: V1.0 11/24/19
      V2.0 03/19/20

Details:Script computing 1D local MSD along a speicfied axis
	for reference or all the sodium atoms. It is also possible to
	tune the window of time and the binning specifics of the
	chosen coordinate. Here adapted for computing other

"""
######################################
#               USAGE
# python script.py file.dat bin_num window lag 
#
#    SCRIPT CURRENT VERSION IS V2.0
######################################

import sys; from sys import argv

import numpy as np

import matplotlib.pyplot as plt

import matplotlib as mpl

import MDAnalysis as mda

mpl.rcParams['figure.dpi']=100
mpl.rcParams['figure.titlesize']=20
mpl.rcParams['axes.facecolor']='white'
mpl.rcParams['lines.linewidth']=2.0
mpl.rcParams['axes.linewidth']=2.0
mpl.rcParams['xtick.major.pad']=6
mpl.rcParams['ytick.major.pad']=6
mpl.rcParams['xtick.labelsize']=14
mpl.rcParams['ytick.labelsize']=14
mpl.rcParams['axes.titlesize']=18
mpl.rcParams['axes.labelsize']=18
mpl.rcParams['axes.grid']='True'
mpl.rcParams['axes.axisbelow']='line'
mpl.rcParams['legend.fontsize']=12


#error if argv less than 8 arguments
if len(sys.argv) != 5:
    sys.exit("\nYou did not provide the correct number\
of arguments in the command line.\nUSAGE IS:\n\
python script.py file.dat bin_num window lag")

##########################################################

# File reading

# N.B. IT IS IMPORTANT TO IMPORT THE TRAJECTORY (XTC) UNWRAPPED WITH
# NO PBC JUMPS OTHERWISE LOCAL_MSD COULD BE WRONG COMPUTED

# Path to the coord file
file_coord = argv[1]

####
#Space discretization

binning_numbers = np.int(argv[2])

md_x_inf = 5.723 
md_x_sup = 9.131

bin_coord = np.linspace(md_x_inf,md_x_sup,binning_numbers)

####
#Time discretization

# Time between frames (ps)
delta_t = 0.004 

# Window of time recording MSD
T   = np.float(argv[3]) #300 # ps

#how many picoseconds you want to skip (crucial param)
lag = np.float(argv[4]) # ps

#check if user provided a credible lag time
if (lag < delta_t):
    sys.exit("\n lag time less than time between frames!!")

shift = np.int(lag/delta_t)

N = np.int(T/(delta_t*shift))

# time array used to plot msd
t = np.linspace(0,T,np.int(T/delta_t))


# Reading data from trajectories

data = np.loadtxt(file_coord)

skipped_index = 1000 

dynamics = data[skipped_index:]

npoints = len(dynamics)

# Definition of the function MSD call in the main (now with shift option)

def compute_msd(traj, starting, steps):

    msds = np.zeros(steps)

    for n in range(0,steps):

        diffs = traj[starting] - traj[starting+n*shift]

        sqdist = np.square(diffs)

        if n==0:

            msds[n] = sqdist

        else:
            
            # here compute msd divide by two times the time recorded (Einstein Smoluchowsky equation)
            # get out the estimated diffusion
            msds[n] = sqdist / (2*delta_t*(n+shift))

    return msds


# main loop on the binning, then on the atoms and finally on the single atom positions

D_temp = np.zeros(N)

Dloc = np.zeros((N, len(bin_coord)-1))

for j in range (0,len(bin_coord)-1):

    count=0

    for i in range (0,npoints):

    # check if position is inside the current bin
        if (bin_coord[j]<=dynamics[i]<bin_coord[j+1] and i<(npoints - (N * shift))):

            D_temp[:]+=compute_msd(dynamics, i, N)

            count+=1

#    insert some statistic here to give an idea of the goodness of the algorithm

#    print (count)

    Dloc[:,j] = D_temp[:]/float(count)

# from ang^2/ps to nm^2/s (can be changed)
#scale_factor = 10000

# keep Dloc as 1/ps
scale_factor = 1

np.savetxt("msd_out.txt", Dloc/scale_factor)

for k in range (0,len(bin_coord)-1):
    print("In the interval %.2f-%.2f"%(bin_coord[k],bin_coord[k+1]))
    print("The computed local diffusivity is %.5E +/- %.5E 1/ps \n" %(np.mean(Dloc[1:,k],axis=0)/scale_factor,np.std(Dloc[1:,k],axis=0)/scale_factor))

####
# Plot part
#for k in range (0,len(bin_coord)-1):
#
#    Dplot[k]=np.sum(Dloc, axis=0)
#    
#    plt.plot((bin_coord[k]+bin_coord[k+1])/2, (Dplot[k]/scale_factor), 'b--',
#            linewidth=2.0, label='%.1f-%.1f' %(bin_coord[k], bin_coord[k+1]))
#    plt.legend(loc='lower right')
#    plt.xlabel('coord')
##    plt.ylabel('MSD ($\AA^2/ps$)')
#    plt.ylabel('Dloc ($cm^2/s$)')
#    plt.show()
