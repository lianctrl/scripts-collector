import numpy as np

def Rec_Protocol(max_t, Deltat):
    min_t = 0.015 #s

    t_pulse = np.linspace(min_t,max_t,np.abs(np.int((max_t-min_t)/Deltat))+1)

    return (t_pulse) 

