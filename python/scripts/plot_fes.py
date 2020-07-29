#!/usr/bin/env python

import sys; from sys import argv

import numpy as np

from scipy.interpolate import CubicSpline

import matplotlib.pyplot as plt

import matplotlib as mpl

if len(sys.argv) != 5:
    sys.exit("\nYou did not provide the correct number\
of arguments in the command line.\n\
USAGE IS: python script.py 0025Mres.dat 1Mres.dat 0025M_PMF.dat 1M_PMF.dat ")

res_0025M = np.genfromtxt(argv[1])

res_1M = np.genfromtxt(argv[2])

data_0025M = np.genfromtxt(argv[3])

data_1M = np.genfromtxt(argv[4])

cs1 = CubicSpline(data_0025M[:,0], data_0025M[:,1])

cs2 = CubicSpline(data_1M[:,0], data_1M[:,1])

# Starting to plot
fig = plt.figure()

xs = np.arange(0.0, 90.0, 0.1)

ys1 = cs1(xs) #- abs(np.amin(cs(xs)))

ys2 = cs2(xs)


x_stch1 = np.array([64.15 for i in range(10000)])

y_stch1 = np.linspace(0.0,20000,10000)

x_ench1 = np.array([62.94 for i in range(10000)])

y_ench1 = np.linspace(0.0,20000,10000)

x_stch2 = np.array([55.68 for i in range(10000)])

y_stch2 = np.linspace(0.0,20000,10000)

x_ench2 = np.array([54.47 for i in range(10000)])

y_ench2 = np.linspace(0.0,20000,10000)

x_up   = np.array([75 for i in range(10000)])

y_up   = np.linspace(0.0,20000,10000)

x_down = np.array([15 for i in range(10000)])

y_down = np.linspace(0.0,20000,10000)

#KbT1 = 0.0083144621 * 298 #KJoule/mol
#KbT2 = 1.380649e-23 * 298 #Joule

ax1 = fig.add_subplot(121) 
ax1.plot(res_0025M[:,0],res_0025M[:,1], label="0025M 38 events")  # res distr from MD 0025M
ax1.plot(res_1M[:,0],res_1M[:,1], label="1M 18 events")  # res distr from MD 1M
ax1.plot(x_up,y_up, 'r-',label="enter CNT")  
ax1.plot(x_down,y_down,'r-', label="exit CNT")
ax1.plot(x_stch1,y_stch1,'g+',markersize=0.5, label="start ch layer")  
ax1.plot(x_ench1,y_ench1,'g+',markersize=0.5, label="end ch layer")
ax1.plot(x_stch2,y_stch2,'g+',markersize=0.5, label="start ch layer")  
ax1.plot(x_ench2,y_ench2,'g+',markersize=0.5, label="end ch layer")

ax1.set_title('Res times-Qtot=-0.84')
ax1.legend(loc='best')
ax1.set_ylim([0.0,2e4])

ax1.set_xlabel('Z coordinate')
ax1.set_ylabel('Counts')
plt.grid()

ax2 = fig.add_subplot(122) 
ax2.plot(xs,ys1, label="Spline 0025M 38 events")  # Free energy from MD 0025M
ax2.plot(xs,ys2, label="Spline 1M 18 events")  # Free energy from MD 1M
ax2.plot(x_up,y_up, 'r-',label="enter CNT")  
ax2.plot(x_down,y_down,'r-', label="exit CNT")
ax2.plot(x_stch1,y_stch1,'go',markersize=2.5, label="start ch layer")  
ax2.plot(x_ench1,y_ench1,'go',markersize=2.5, label="end ch layer")
ax2.plot(x_stch2,y_stch2,'go',markersize=2.5, label="start ch layer")  
ax2.plot(x_ench2,y_ench2,'go',markersize=2.5, label="end ch layer")

ax2.set_title('PMF-Qtot=-0.84')
ax2.legend(loc='best')
ax2.set_ylim([0.0,30.0])

ax2.set_xlabel('Z coordinate')
ax2.set_ylabel('PMF [KJoule/mol]')
plt.grid()
plt.show()
