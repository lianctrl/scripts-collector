#!/usr/bin/env python
"""
Author:Luca Sagresti

Date: V1.0 11/24/19
      V2.0 03/19/20

Details:Script computing 1D local MSD along a speicfied axis
	for reference or all the sodium atoms. It is also possible to
	tune the window of time and the binning specifics of the
	chosen coordinate

"""
######################################
#               USAGE
# python script.py file.pdb file.xtc start_pos end_pos spacing window lag yes/no axis
#
#    SCRIPT CURRENT VERSION IS V2.0
######################################

import sys; from sys import argv

import numpy as np

import matplotlib.pyplot as plt

import matplotlib as mpl

import MDAnalysis as mda

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


#error if argv less than 8 arguments
if len(sys.argv) != 10:
    sys.exit("\nYou did not provide the correct number\
of arguments in the command line.\nUSAGE IS:\n\
python script.py file.pdb file.xtc start_pos end_pos spacing window lag yes/no axis")

##########################################################

# File reading

# N.B. IT IS IMPORTANT TO IMPORT THE TRAJECTORY (XTC) UNWRAPPED WITH
# NO PBC JUMPS OTHERWISE LOCAL_MSD COULD BE WRONG COMPUTED

# Path to the pdb file
PDB = argv[1]

# Path to the pbc nojump xtc file (use gmx trjconv -pbc nojump)
XTC = argv[2]

u = mda.Universe(PDB, XTC)

####
#Space discretization

start_pos = float(argv[3]) # in Angstrom

end_pos   = float(argv[4]) # in Angstrom

if (start_pos<end_pos):

    print("Script written to assign start pos the greater value: \
swap between end and start pos. Run continue...")

    start_pos, end_pos = end_pos, start_pos

spacing = float(argv[5]) # in Angstrom

if (spacing>start_pos or spacing<0):

    sys.exit("Wrong spacing inserted")

number_points = np.int(np.abs((end_pos-start_pos)/spacing))

bin_coord = np.linspace(start_pos,end_pos,number_points+1)

####
#Time discretization

# Time between frames (ps)
delta_t = u.trajectory[1].time-u.trajectory[0].time

# Window of time recording MSD
T   = np.float(argv[6]) #300 # ps

#how many picoseconds you want to skip (crucial param)
lag = np.float(argv[7]) # ps

#check if user provided a credible lag time
if (lag < delta_t):
    sys.exit("\n lag time less than time between frames!!")

shift = np.int(lag/delta_t)

N = np.int(T/(delta_t*shift))

# time array used to plot msd
t = np.linspace(0,T,np.int(T/delta_t))

#check on the group selected
if (argv[8]=="yes"):
    print ("Assuming you want to study only the reference ion")
    NA_atoms = u.select_atoms("resname NAR")

elif (argv[8]=="no"):
    print("Assuming you want to study the non-restrained NA ions")
    NA_atoms = u.select_atoms("resname NA")

elif (argv[8]=="both"):
    print("Assuming you want to study the non-restrained NA ions")
    NA_atoms = u.select_atoms("resname NA or resname NAR")

else:
    sys.exit("\n Unknown command for selection group")

#check on the axis selected
if(argv[9]=="x"):
    print("x axis chosen")
    axis = 0

elif(argv[9]=="y"):
    print("y axis chosen")
    axis = 1

elif(argv[9]=="z"):
    print("z axis chosen")
    axis = 2

else:
    sys.exit("\n Unknown axis inserted")

# Reading data from trajectories

oneD_coordinate_NA = np.zeros((len(NA_atoms),len(u.trajectory)))

box_lenghts = u.dimensions

oneD_dimension = box_lenghts[axis]

m=0

for frm in u.trajectory[:]:

    oneD_coordinate_NA[:,m] = NA_atoms.positions[:,axis]

    m+=1

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

    for l in range (0,len(NA_atoms)):

        for i in range (0,len(u.trajectory)):

            #Tricky part to ensure unwrapped trajectory being correctly read

            N_Box = np.floor(oneD_coordinate_NA[l][i]/oneD_dimension)

            lim_inf = oneD_dimension* N_Box + bin_coord[j+1]

            lim_sup = oneD_dimension* N_Box + bin_coord[j]

            # check if position is inside the current bin
            if (lim_inf<=oneD_coordinate_NA[l][i]<lim_sup and i<(len(u.trajectory)- N * shift)):

                D_temp[:]+=compute_msd(oneD_coordinate_NA[l], i, N)

                count+=1

#       insert some statistic here to give an idea of the goodness of the algorithm

#    print (count)

    Dloc[:,j] = D_temp[:]/float(count)

# from ang^2/ps to nm^2/s (can be changed)
scale_factor = 10000

np.savetxt("msd_out.txt", Dloc/scale_factor)

for k in range (0,len(bin_coord)-1):
    print("In the interval %.2f-%.2f"%(bin_coord[k],bin_coord[k+1]))
    print("The computed local diffusivity is %.5E +/- %.5E cm^2/s \n" %(np.mean(Dloc[1:,k],axis=0)/scale_factor,np.std(Dloc[1:,k],axis=0)/scale_factor))

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
