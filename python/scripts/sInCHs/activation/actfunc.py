import numpy as np

from . import actprotocol
from . import actode

from scipy.integrate import odeint

def Act_func(x, teta, S0):

    # conductance parameters

    #EK      = 0.0    # mV
    gK_max  = 33.2     # nS

    # Membrane capacitance
    Cm      = 1.0    # microF cm^-2

    #define activation sequence protocol

    Vmax = 60.0 #mV

    increment = 10.0 #mV

    Vtest = actprotocol.Act_Protocol(Vmax, increment)

    # Time discretiztion

    tend = 1.05

    ttest_i = 0.50 #time at which you start record the current

    ttest_f = 0.55 #time at which you end record the current

    Npoints = len(x)

    Points_per_sec = np.int(Npoints/tend)


    # prepare empty arrays
    Open_states = np.zeros((Npoints,len(Vtest)))

    max_conductance = np.zeros(len(Vtest))


    for i in range (0,len(Vtest)):
        gamma = np.append(teta, Vtest[i])

        f = lambda S,t: actode.ode_Act(S, t, gamma)

        r = odeint(f, S0, x)

        Open_states[:,i] = r[:,10]

        max_conductance[i] = gK_max * np.amax(r[np.int(Points_per_sec*ttest_i):np.int(Points_per_sec*ttest_f),10])


    g_gmax = (max_conductance)/(np.amax(max_conductance))

    return (g_gmax)

