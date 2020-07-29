#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
V1.0 04/07/2020

Author: Luca Sagresti
"""

import numpy as np
import matplotlib.pyplot as plt; import matplotlib as mpl
from scipy.interpolate import CubicSpline
import sys; from sys import argv

#error if argv less than 3 arguments
if len(sys.argv) != 3:
    sys.exit("\nYou did not provide the correct number\
of arguments in the command line.\nUSAGE IS:\n\
python script.py 0025M.dat 1M.dat")


fig=plt.figure()
plt.title('Results from Hummer Method')
x1 = (binning1[1:]+binning1[:-1])/2 #bin centers
#ax2=fig.add_subplot(122)
x1 = (x1[1:]+x1[:-1])/2
ax2.plot(x, cs_md_diff(x), label='From MD')
ax2.plot(x, D_Hum, label='From Hummer')
ax2.plot(x1, D1, label='From scipy')
ax2.legend()
plt.grid()
plt.show()

cs_Humm=CubicSpline(x1, D1)

# detailed balance check

err=np.zeros((len(binning1)-1, len(binning1)-1))

for i in range(len(err)):
    err[i][i]=R1[i][i]
    for j in range(len(err)):
        if j!=i:
            err[i][i]+=R[j][i]
    for j in range(i+1, len(err)):
        err[i][j]=err[j][i]=R1[i][j]-R1[j][i]*float(p1[i]/p1[j])


plt.figure()
plt.imshow(err, cmap='jet')
#plt.imshow(R1)
plt.colorbar()
plt.show()
