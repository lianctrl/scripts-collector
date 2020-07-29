#!/usr/bin/env python3

import numpy as np

delr=0.0005 #nm

rcut=1.3 #nm

rnull=0.05 #nm

k=0.811 #k*C6=C4 for LJ 12-6-4

nbins=np.int((rcut+1)/delr)+1

for i in range(0,nbins):

    r=delr*i

    if r <= rnull:

        print (r,0.0,0.0,0.0,0.0,0.0,0.0)

    else:

        print (r,1.0/r,1/r**2,(-1/r**4)*(1/r**2+k),(-2/r**5)*(3/r**2+2*k),1/r**12,12/r**13)

# to run: python3 gen_table.py >> table.xvg
