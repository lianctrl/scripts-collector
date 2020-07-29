import numpy as np

from . import inacprotocol
from . import inacode

from scipy.integrate import odeint

def Inac_func(x, teta, S0):

    # conductance parameters

    EK      = -90.0    # mV
    gK_max  = 33.2     # nS

    # Membrane capacitance
    Cm      = 1.0    # microF cm^-2


    Vtesting =  60.0 #mV


    #define activation sequence protocol

    Vmax = 60.0 #mV

    increment = 10.0 #mV

    Vtest = inacprotocol.Inact_Protocol(Vmax, increment)

    # Time discretiztion

    tend = 3.00

    ttest = 2.00 #time at which you start record the current

    Npoints = len(x)

    Points_per_sec = np.int(Npoints/tend)


    # prepare empty arrays
    Open_states = np.zeros((Npoints,len(Vtest)))

    max_conductance = np.zeros(len(Vtest))

    max_currents = np.zeros(len(Vtest))

    for i in range (0,len(Vtest)):
        gamma = np.append(teta, Vtest[i])

        f = lambda S,t: inacode.ode_Inac(S, t, gamma)

        r = odeint(f, S0, x)

        Open_states[:,i] = r[:,10]

        max_conductance[i] = gK_max * np.amax(r[np.int(Points_per_sec*ttest):,10]-r[np.int(Points_per_sec*ttest)-1,10])

        max_currents[i] = max_conductance[i] * (Vtesting - EK)

    I_Imax = (max_currents)/(np.amax(max_currents))

    return (I_Imax)

