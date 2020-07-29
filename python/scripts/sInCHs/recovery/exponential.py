import numpy as np

def Simple_exp_eq(x,tau, tau_0 = 15.0):
    return (1.0 - np.exp(-(x - tau_0 )/ tau))

