#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
V1.0 04/07/2020
V2.0 04/09/2020 reshaping all the atoms in one array
                choice of the resoultion method: Higham or inverse method

Authors: L. Sagresti and G. Brancato and L. Peri
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
if len(sys.argv) != 10:
    sys.exit("\nYou did not provide the correct number\
of arguments in the command line.\nUSAGE IS:\n\
python3 script.py file.pdb file.xtc start_pos end_pos spacing lag yes/no axis method")

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

if (start_pos>end_pos):

    print("Script written to assign start pos the lower value: \
swap between end and start pos. Run continue...")

    start_pos, end_pos = end_pos, start_pos

spacing = float(argv[5]) # in Angstrom

if (spacing>end_pos or spacing<0):

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
    print("Assuming you want to study the non-restrained NA ions and NAR")
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

m=0
k=0

# load trajectory
for frm in u.trajectory[:]:

    if ( ( m % int(len(u.trajectory)/10)  ) == 0 ) :

        print(" %i / 10 of file loaded- %s" %(k, datetime.datetime.now()))

        k+=1

    oneD_coordinate_NA[:,m] = NA_atoms.positions[:,axis]

    m+=1

oneD_coordinate_NA = np.reshape(oneD_coordinate_NA,len(u.trajectory)*len(NA_atoms))

#  Hummer diffusion

binning_space = bin_coord[1] - bin_coord[0]

print('Computing diffuson with Hummer Method for %d bins' %(len(bin_coord)-1))

#Function declaration
def exp_R_Hummer(nshift, binning, traj):

    exp_R = np.zeros((len(binning)-1, len(binning)-1))

    binning_space = binning[1] - binning[0]

    ix = -1 * np.ones(len(traj[:-nshift]))

    n_Box = np.floor(traj[:-nshift]/oneD_dimension)

    lim_inf = np.zeros(len(traj[:-nshift]))

    lim_sup = np.zeros(len(traj[:-nshift]))

    for i in range(0,len(binning)-1):

        lim_inf = oneD_dimension * n_Box + binning[i]

        lim_sup = oneD_dimension * n_Box + binning[i+1]

        #find indexes of trajectory where positions are between the binning
        #spaces and put in the index of the correct binning
        ix = np.where(np.logical_and(traj[:-nshift] < lim_sup[:] , traj[:-nshift] >= lim_inf[:]), i, ix)

    ix = np.where(ix==-1, len(binning)-2, ix)

    ix_delta = np.zeros(len(traj[:-nshift]))

    xj = ix_delta

    print('Computing exponential of the Hummer Matrix')

    for n in range(0,len(traj[:-nshift])):

        ix_delta[n] = traj[n] - binning[np.int(ix[n])] - (n_Box[n] * oneD_dimension) + 0.5 * binning_space

        xj[n] = traj[n+nshift] - ix_delta[n]

    jx = -1 * np.ones(len(traj[:-nshift]))

    lim_inf = np.zeros(len(traj[:-nshift]))

    lim_sup = np.zeros(len(traj[:-nshift]))

    for i in range(0,len(binning)-1):

        lim_inf = oneD_dimension * n_Box + binning[i]

        lim_sup = oneD_dimension * n_Box + binning[i+1]

        #find indexes of trajectory where positions are between the binning
        #spaces and put in the index of the correct binning
        jx = np.where(np.logical_and(xj[:] < lim_sup[:] , xj[:] >= lim_inf[:]), i, jx)

    jx = np.where(jx==-1, len(binning)-2, jx)

    print('Initialising exp(tR) at t = %d ps' %(nshift * delta_t))

    for n in range(0,len(traj[:-nshift])):

        # sum one for every occurrence in that binning interval previously passed through ix and jx
        exp_R[np.int(jx[n])][np.int(ix[n])]+=1

    Ri = np.zeros(len(binning)-1)

    for i in range(0,len(exp_R)):

        for j in range(0,len(exp_R)):

            Ri[i] = Ri[i] + exp_R[j][i]

    for i in range(0,len(exp_R)):

       for j in range(0,len(exp_R)):

           if Ri[j]!=0:
               #normalize the rates, to have them < 1
               exp_R[i][j] = np.float(np.float(exp_R[i][j]) / np.float(Ri[j]))

    return exp_R

# Start the Main

exp_R = np.zeros((len(bin_coord)-1, len(bin_coord)-1))

print("Started computing the exponential matrix- %s" %(datetime.datetime.now()))

exp_R = exp_R_Hummer(shift, bin_coord, oneD_coordinate_NA)

p,_   = np.histogram(oneD_coordinate_NA, bins=bin_coord, density=True)

delta_Q = binning_space

D_Hum   = np.zeros(len(bin_coord)-2)

'''
#Here a straightforward problem:
#matrix exp_R_sym (exp_(tR~)) it is not symmetric
#we can enforce it by hand, however eigenvalues could be < 0,
#and the logarithm cannot be taken, so the inversion trick at
#the end of the comment section is needed


exp_R_sym=np.zeros([len(bin_coord)-1,len(bin_coord)-1])

for i in range(len(exp_R)):

    for j in range(len(exp_R)):

        exp_R_sym[i][j] = np.sqrt(1/p[j])*exp_R[i][j]*np.sqrt(p[i])

#exp_R_sym it is not symmetric so we enforce it to be

for i in range(len(exp_R)):

    for j in range(i+1):

        exp_R_sym[i][j] = exp_R_sym[j][i] = float((exp_R_sym[i][j]+exp_R_sym[j][i])/2)


#w, v=np.linalg.eig(exp_R_sym)
#print((w))
'''

method = argv[9]

if (method=="hummer"):

	print("Using Hummer method to evaluate the matrix R")

	print('Computing time derivtive O(dt^4)')

	dt_exp_R = (-1 * exp_R_Hummer(shift+2, bin_coord, oneD_coordinate_NA)+ 8 * exp_R_Hummer(shift+1, bin_coord, oneD_coordinate_NA)- 8 * exp_R_Hummer(shift-1, bin_coord, oneD_coordinate_NA) + exp_R_Hummer(shift-2, bin_coord, oneD_coordinate_NA))/(12 * delta_t)

	print('Inverting exp_(tR) to get R=(d/dt exp(tR))*(exp(-tR)')

	R = np.dot(dt_exp_R, np.linalg.inv(exp_R))

	for i in range(0,len(D_Hum)):

		if (p[i+1]==0):

			D_Hum[i]=0

		else:

			D_Hum[i]  = (delta_Q**2) * (R[i+1][i] * np.sqrt(p[i]/p[i+1]))#+R[i][i+1]*np.sqrt(p[i+1]/p[i]))/2

#			D_Hum2[i] = (delta_Q**2) * (R[i][i+1] * np.sqrt(p[i+1]/p[i]))

			#Manually enforcing detailed balance

elif(method=="higham"):

	print("Using Higham algorithm for the taylor expansion of the matrix logarithm")

	R = (1/(shift * delta_t)) * logm(exp_R)

	for i in range(0,len(D_Hum)):

		if (p[i+1]==0):

			D_Hum[i]=0

		else:

			D_Hum[i]  = (delta_Q**2) * (R[i+1][i] * np.sqrt(p[i]/p[i+1]))#+R1[i][i+1]*np.sqrt(p1[i+1]/p1[i]))/2

#			D_Hum2[i] = (delta_Q**2) * (R1[i][i+1] * np.sqrt(p[i+1]/p[i]))

			#Manually enforcing detailed balance

else:
	sys.exit("\n Unknown method inserted, only hummer or higham are available!")

print("Resumes of the parameters used:")
print("Bin width used is %.4f Ang, lagtime used is %1.2f ps"%(binning_space,lag))
print("Everything ok, saving data in progress...")

#binning centers
bin_centers = (bin_coord[1:] + bin_coord[:-1]) / 2

scale_factor = 10000 #from A^2/ps to cm^2/s

#plot_data = np.array([bin_centers,D_Hum/scale_factor])

np.savetxt("Hummer_Diffusion.dat", D_Hum/scale_factor)

#not possible to save like that since are colums of different length
#np.savetxt("plot_local_diffusion.dat", plot_data.T)
