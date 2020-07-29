import numpy as np

def ode_Act (S, t, p):

    C0=S[0]
    C1=S[1]
    C2=S[2]
    C3=S[3]
    C4=S[4]
    I0=S[5]
    I1=S[6]
    I2=S[7]
    I3=S[8]
    I4=S[9]
    O =S[10]


    #constants

    T = 291.0 #K or 18 degree celsius
    e =  1.602176634 * (10**-19.0) # C
    K_B = 1.380649 * (10**-23.0) # J*K^-1

    exp_factor = (e/(K_B * T)) * (10**-3)

    #Voltage sequences

    V = 0.0 #mV

    if 0 <= t < 0.5:
        V = -90.0 # mV

    if 0.5 <= t < 0.55:
        V = p[11] # Vtest

    if 0.55 <= t <= 1.05:
        V = -50.0 # mV

    if t > 1.05 :
        V = -50.0 #mV


    #voltage dependent rate constants

    alpha = p[0] * np.exp(p[1] * (V * exp_factor))
    beta = p[2] * np.exp(-1.0 * p[3] * (V * exp_factor))
    k_CO = p[4] * np.exp(p[5] * (V * exp_factor))
    k_OC = p[6] * np.exp(-1.0 * p[7] * (V * exp_factor))


    k_CI = p[8]

    k_IC = p[9]

    f = p[10]

    # ODEs

    dC0dt = beta * C1 + (k_IC/(f**4.0)) * I0 - (k_CI*(f**4.0) + 4.0 * alpha) * C0
    dC1dt = 4.0 * alpha * C0 + 2.0 * beta * C2 + (k_IC/(f**3.0)) * I1 - (k_CI*(f**3.0) + beta + 3.0 * alpha) * C1
    dC2dt = 3.0 * alpha * C1 + 3.0 * beta * C3 + (k_IC/(f**2.0)) * I2 - (k_CI*(f**2.0) + 2.0 * beta + 2.0 * alpha) * C2
    dC3dt = 2.0 * alpha * C2 + 4.0 * beta * C4 + (k_IC/f) * I3 - (k_CI*f + 3.0 * beta + 1.0 * alpha) * C3
    dC4dt = 1.0 * alpha * C3 + k_OC * O + k_IC * I4 - (k_CI + k_CO + 4.0 * beta) * C4

    dI0dt = beta * f * I1 + (k_CI*(f**4.0)) * C0 - (k_IC/(f**4.0) + 4.0 * (alpha/f)) * I0
    dI1dt = 4.0 * (alpha/f) * I0 + 2.0 * beta * f * I2 + (k_CI*(f**3.0)) * C1 - (k_IC/(f**3.0) + beta * f + 3.0 * (alpha/f)) * I1
    dI2dt = 3.0 * (alpha/f) * I1 + 3.0 * beta * f * I3 + (k_CI*(f**2.0)) * C2 - (k_IC/(f**2.0) + 2.0 * beta * f + 2.0 * (alpha/f)) * I2
    dI3dt = 2.0 * (alpha/f) * I2 + 4.0 * beta * f * I4 + (k_CI*f) * C3 - (k_IC/f + 3.0 * beta * f + 1.0 * (alpha/f)) * I3
    dI4dt = 1.0 * (alpha/f) * I3 + k_CI * C4 - (k_IC + 4.0 * beta * f) * I4

    dOdt = k_CO * C4 - k_OC * O

    return (dC0dt, dC1dt, dC2dt, dC3dt, dC4dt, dI0dt, dI1dt, dI2dt, dI3dt, dI4dt, dOdt)

