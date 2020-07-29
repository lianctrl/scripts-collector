#!/usr/bin/env python                                                                      


import sys; from sys import argv

import numpy as np

import MDAnalysis as mda
######################################
#               USAGE 
# python passage.py file.pdb file.xtc output.dat
#
#
######################################

PDB = argv[1]

XTC = argv[2]

Output_file = argv [3]

u = mda.Universe(PDB, XTC)

NA_atom=u.select_atoms("resname NAR and resid 1")

f=open(Output_file,"w")


for frm in u.trajectory[:]:

    z_coordinate_NAR=NA_atom.positions[:,2]

    f.write(str(z_coordinate_NAR[0])+"\n")


f.close()

