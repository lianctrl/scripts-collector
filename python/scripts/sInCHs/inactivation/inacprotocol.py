import numpy as np

def Inact_Protocol(max_V, DeltaV):
    Vhold = -90.0 #mV

    Vtest = np.linspace(Vhold,max_V,np.abs(np.int((max_V-Vhold)/DeltaV))+1)

    return (Vtest)
