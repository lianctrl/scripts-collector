
import numpy as np
#from scipy.integrate import odeint
from scipy.optimize import curve_fit
from scipy.optimize import minimize
#from scipy.optimize import least_squares

from sInCHs.activation import actfunc
from sInCHs.activation import boltzmann
from sInCHs.activation import actprotocol
from sInCHs.inactivation import inacfunc
from sInCHs.inactivation import invboltzmann
from sInCHs.inactivation import inacprotocol
from sInCHs.recovery import recfunc
from sInCHs.recovery import exponential
from sInCHs.recovery import recprotocol

from SALib.sample import saltelli
from SALib.analyze import sobol

import pandas as pd


###############################################################
##################### Import Experimental Data ################
folder="~/PHD/NOTEBOOK/Data_exp_Kv/"

# Activation dataset

dataset_act = pd.read_csv(folder+"Act_datasets.csv", skiprows=1)

act_data = dataset_act.to_numpy()

act_data_WT = act_data[:,1]

act_data_S390N = act_data[:,3]

# Inactivation dataset

dataset_inact = pd.read_csv(folder+"Inact_datasets.csv", skiprows=1)

inact_data = dataset_inact.to_numpy()

inact_data_WT = inact_data[:,1]

inact_data_S390N = inact_data[:,3]

# Recovery dataset

dataset_rec = pd.read_csv(folder+"Rec_datasets.csv", skiprows=1)

rec_data = dataset_rec.to_numpy()

rec_data_WT = rec_data[:,1]

rec_data_S390N = rec_data[:,3]

######################### End part of Import ########################
#####################################################################

######################## Cost Function Building #####################


def residual_Tot(p, S0):

    # Time discretiztion
    tini = 0.0

    tend = 4.00

#    ttest = 0.56 #time at which you start record the current

    Npoints = 200000

    Points_per_sec = np.int(Npoints/tend)

    # time array
    t = np.linspace(tini,tend,Npoints)

####################### START OF THE 3 VOLTAGE PROTOCOLS ###################

    # DEFINE ACTIVATION SEQUENCE PROTOCOL

    Vmax = 60.0 #mV

    increment = 10.0 #mV

    Vtest_act = actprotocol.Act_Protocol(Vmax, increment)

    #experimental data for hemiactivation and slope

    exp_Vhemi_act = -22.70 #mV

    exp_k_act = 13.00 #mV

    #evaluation of experimental points trough the declared Boltzmann fit

    exp_data_act = boltzmann.Boltz_eq(Vtest_act,exp_Vhemi_act,exp_k_act)

    # here data from perfect fit
    #sq_err_act = np.sum(np.subtract(exp_data_act,actfunc.Act_func(t,p,S0))**2)

    #here experimental data imported with Pandas in initial section
    sq_err_act = np.sum(np.subtract(act_data_WT, actfunc.Act_func(t,p,S0))**2)

    act_cost_func = (1.0/len(Vtest_act)) * sq_err_act


    # DEFINE INACTIVATION SEQUENCE PROTOCOL

#    Vmax = 60.0 #mV

#    increment = 10.0 #mV

    Vtest_inact = inacprotocol.Inact_Protocol(Vmax, increment)

    #experimental data for hemiactivation and slope

    exp_Vhemi_inact = -49.60 #mV

    exp_k_inact = 5.10 #mV

    #evaluation of experimental points trough the declared Boltzmann fit

    exp_data_inact = invboltzmann.Boltz_eq_in(Vtest_inact,exp_Vhemi_inact,exp_k_inact)

    # here data from perfect fit
    #sq_err_inact = np.sum(np.subtract(exp_data_inact,inacfunc.Inac_func(t,p,S0))**2)

    #here experimental data imported with Pandas in initial section
    sq_err_inact = np.sum(np.subtract(inact_data_WT,inacfunc.Inac_func(t,p,S0))**2)

    inact_cost_func = (1.0/len(Vtest_inact)) * sq_err_inact


    # DEFINE RECOVERY SEQUENCE PROTOCOL

    tmax = 300.0 #ms

    increment = 15.0 #ms

    t_pulse = recprotocol.Rec_Protocol(tmax, increment)

    #experimental data for recovery time

    exp_tau = 47.0 #ms

    #evaluation of experimental points trough the declared exponential fit

    exp_data_rec = exponential.Simple_exp_eq(t_pulse,exp_tau)

    # here data from perfect fit
    #sq_err_rec = np.sum(np.subtract(exp_data_rec,recfunc.Rec_func(t,p,S0))**2)

    #here experimental data imported with Pandas in initial section
    sq_err_rec = np.sum(np.subtract(rec_data_WT,recfunc.Rec_func(t,p,S0))**2)

    rec_cost_func = (1.0/len(t_pulse)) * sq_err_rec


    # sum of the three protocols cost function


    return (act_cost_func + inact_cost_func + rec_cost_func)

##################### SENSITIVITY ANALYSIS #######################################

C0_0 = 439.0
C1_0 = 258.8
C2_0 = 57.2
C3_0 = 5.6
C4_0 = 0.2
I0_0 = 12.8
I1_0 = 55.3
I2_0 = 89.4
I3_0 = 64.2
I4_0 = 17.2
O_0  = 0.1

# Pack up the initial conditions:

C0 = [C0_0, C1_0, C2_0, C3_0, C4_0, I0_0, I1_0, I2_0, I3_0, I4_0, O_0]


problem ={
    'num_vars' : 11,
    'names' : ['alpha_0','alpha_1','beta_0','beta_1','KCO_0','KCO_1','KOC_0','KOC_1','KCI','KIC','f'],
    'bounds' : [[100.0,10000.0],
                [0.0,1.0],
                [0.0,400.0],
                [0.0,10.0],
                [200.0,1000.0],
                [0.0,1.0],
                [0.0,300.0],
                [0.0,1.0],
                [0.0,400.0],
                [0.0,200.0],
                [0.0,0.8]]
}

param_values=saltelli.sample(problem,500)

Y = np.zeros([param_values.shape[0]])

for i, X in enumerate(param_values):

    Y[i] = residual_Tot(X,C0)

# Substitute possible Nan with the mean of non-Nan values
Y[np.isnan(Y)] = np.nanmean(Y)

Si = sobol.analyze(problem, Y, print_to_console=True)
