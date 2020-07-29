import numpy as np

from . import recode
from . import recprotocol

from scipy.integrate import odeint

def Rec_func(x, teta, S0):

    # conductance parameters

    EK      = -90.0    # mV
    gK_max  = 33.2      # nS

    # Membrane capacitance
    Cm      = 1.0    # microF cm^-2

    #define activation sequence protocol

    tmax = 0.300 #s

    increment = 0.015 #s

    t_pulse = recprotocol.Rec_Protocol(tmax, increment)

    tini_prep = 1.00 #s

    tend_prep = 2.00 #s

    Vtesting = 60.0 #mV

    # Time discretiztion

    tend = 4.00 #s

    Npoints = len(x)

    Points_per_sec = np.int(Npoints/tend)

    # prepare empty arrays
    Open_states = np.zeros((Npoints,len(t_pulse)))

    max_conductance = np.zeros(len(t_pulse))

    max_currents = np.zeros(len(t_pulse))

    max_conductance_prep = np.zeros(len(t_pulse))

    max_currents_prep = np.zeros(len(t_pulse))

    for i in range (0,len(t_pulse)):

        gamma = np.append(teta, t_pulse[i])

        f = lambda S,t: recode.ode_Rec(S, t, gamma)

        r = odeint(f, S0, x)

        Open_states[:,i] = r[:,10]

        max_conductance[i] = gK_max * np.amax(r[np.int(Points_per_sec*(tend_prep+t_pulse[i])):,10])

        max_conductance_prep[i] = gK_max * np.amax(r[np.int(Points_per_sec*tini_prep):np.int(Points_per_sec*tend_prep),10])

        # Compute the current proportional to the open channel conductance and potential applied

        max_currents[i] = max_conductance[i] * (Vtesting - EK)

        max_currents_prep[i] = max_conductance_prep[i] * (Vtesting - EK) # nS * mV = pA

    I_Imax = np.true_divide(max_currents,max_currents_prep)

    return (I_Imax)


