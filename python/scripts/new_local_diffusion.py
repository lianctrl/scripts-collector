#!/usr/bin/env python
"""
Author: Luca Sagresti

Date:   V1.0 03/20/20

Details:Script computing 1D local MSD along a speicfied axis
        for reference or all the sodium atoms. It is also possible to
        tune the lag time and the binning specifics of the
        chosen coordinate

"""
######################################
#               USAGE
# python script.py file.pdb file.xtc start_pos end_pos spacing lag yes/no axis
#
#    SCRIPT CURRENT VERSION IS V1.0
######################################

import sys; from sys import argv

import numpy as np

import MDAnalysis as mda

#error if argv less than 9 arguments
if len(sys.argv) != 9:
    sys.exit("\nYou did not provide the correct number\
of arguments in the command line.\nUSAGE IS:\n\
python script.py file.pdb file.xtc start_pos end_pos spacing lag yes/no axis")

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

#how many picoseconds you want to skip
lag = np.float(argv[6]) # ps

#check if user provided a credible lag time
if (lag < delta_t):
    sys.exit("\n lag time less than time between frames!!")

shift = np.int(lag/delta_t)

#check on the group selected
if (argv[7]=="yes"):
    print ("Assuming you want to study only the reference ion")
    NA_atoms = u.select_atoms("resname NAR")

elif (argv[7]=="no"):
    print("Assuming you want to study the non-restrained NA ions")
    NA_atoms = u.select_atoms("resname NA")

elif (argv[7]=="both"):
    print("Assuming you want to study the non-restrained NA ions")
    NA_atoms = u.select_atoms("resname NA or resname NAR")

else:
    sys.exit("\n Unknown command for selection group")

#check on the axis selected
if(argv[8]=="x"):
    print("x axis chosen")
    axis = 0

elif(argv[8]=="y"):
    print("y axis chosen")
    axis = 1

elif(argv[8]=="z"):
    print("z axis chosen")
    axis = 2

else:
    sys.exit("\n Unknown axis inserted")

# Reading data from trajectories

oneD_coordinate_NA = np.zeros((len(NA_atoms),len(u.trajectory)))

box_lenghts = u.dimensions

oneD_dimension = box_lenghts[axis]

j=0

for frm in u.trajectory[:]:

    oneD_coordinate_NA[:,j] = NA_atoms.positions[:,axis]

    j+=1

# main loop on the atoms, then on the single atom positions and finally on the binning coordinate

count = np.zeros((len(NA_atoms),len(bin_coord)-1))

Ddyn = np.zeros((len(NA_atoms),len(bin_coord)-1))

for m in range (0,len(NA_atoms)):

    for n in range(0,len(u.trajectory)-shift):

        temp = oneD_coordinate_NA[m][n+shift]-oneD_coordinate_NA[m][n]

        for i in range (0,len(bin_coord)-1):

            N_Box = np.floor(oneD_coordinate_NA[m][n]/oneD_dimension)

            lim_inf = oneD_dimension* N_Box + bin_coord[i+1]

            lim_sup = oneD_dimension* N_Box + bin_coord[i]

            if (lim_inf<=oneD_coordinate_NA[m][n]<lim_sup):

                count[m][i] += 1.0

                Ddyn[m][i]  += ((temp * temp) / (2.0 * lag))

tot_count = np.sum(count,axis=0)

temp_D = np.sum(Ddyn,axis=0)

Dloc = temp_D/tot_count

scale_factor = 10000 # from A^2/ps to cm^2/s

print("\n All right, no problem soldier!\n")

np.savetxt("Dloc_out.txt", Dloc/scale_factor)

for k in range(0,len(bin_coord)-1):
    print ("%.5E cm^2/s in the interval %2.3f-%2.3f" % (Dloc[k]/scale_factor,bin_coord[k],bin_coord[k+1]))
