import numpy as np

def Boltz_eq(x, V1_2, k):
    return (1.0/(1.0 + np.exp(-(x-V1_2)/k)))

