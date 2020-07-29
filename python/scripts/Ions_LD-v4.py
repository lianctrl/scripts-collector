#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 14:22:58 2020

@author: lorenzo
"""

import numpy as np
import matplotlib.pyplot as plt; import matplotlib as mpl
from scipy.interpolate import CubicSpline
import datetime
from scipy.optimize import curve_fit

# In[1]: Tuning Parameters

skipped_index=1000

coordination_bin_width=0.05 # width of a bin for the histogram of coordination number- you need to have a
                           # statistics of at least one for each bin!  Original=0.05
#Diffusion computing 
                          
lagtime  = 0.072        # Lag time in ps, Md timestep = 0.004 ps, #original =0.016
md_dt = 0.004        # time step in ps
ilagtime = int(lagtime / md_dt)   # Lag time converted into index position

binning_numbers = [50]   #binning numbers for diffusion (original=50)

#MFTP computing

t_nbins = 30    # number of bins to divide the MFPT into for BOTH MD and LD

#LD

Sym_type='short'        # 'long' to perform long LD
                       #'short' to perform short LD
                       #'both' to perform short both
                       
ld_dt                   = md_dt/2.0      # md_dt=0.004
                       
ld_t_end_long           = 100000.0     # long ld end (picoseconds)

ld_t_end_short         = 1000.0 # picoseconds (MD ends ant around 500 ps)
ld_ntrials             = 2000

#%%
print("\nreading input coordination file- should take 2m10s for a 100ns, 25000000 lines")

data=np.loadtxt("/Users/lorenzo/Google Drive/LORENZO/Files/gromos_spc-Ca-data.dat")

npoints = len(data)
print("\nsuccessfully read %d data" %(len(data)))

dynamics = data


# In[2]:



fig = plt.figure()
ax0 = fig.add_subplot(111)

xperquestoplot=np.linspace(0,len(data)*0.000004,len(data)) #so it is in ns
ax0.plot(xperquestoplot,data)

ax0.set_ylabel('Coordination Number')
ax0.set_xlabel('time [ns]')
ax0.set_title('Coordination Number obtained in MD')
plt.show()


md_x_inf = np.min(dynamics)
md_x_sup = np.max(dynamics)
md_x_bins = np.arange(md_x_inf, md_x_sup, coordination_bin_width)
#md_x_bins_max_precision = np.arange(md_x_inf, md_x_sup, 0.01)
print("Calculating FES for %d bins of width %3.2f" %(len(md_x_bins),coordination_bin_width))
print("Range of coordination number obtained in MD:  [ %.3f  :  %.3f ]" % (md_x_inf, md_x_sup))

hist, _ = np.histogram(dynamics, bins=md_x_bins, density=True)

plt.figure()
plt.hist(dynamics, md_x_bins)
plt.xlim([5,10])
plt.xlabel('Coordination Number')
plt.ylabel('Statistics')
plt.title('Statistics of coordination number')
plt.yscale('log')
#plt.legend(loc='best')
plt.show()



# In[3]:

# Compute Free Energy from MD trajectory and get spline function

x_free_ener = (md_x_bins[1:] + md_x_bins[:-1]) / 2

tmp = -1.0 * np.log(hist)
tmp = tmp + abs(np.amin(tmp))
free_ener = tmp


# now free_ener is in units of KbT. Therefore:
KbT = 0.0083144621 * 300   #KJoule/mol
free_ener = free_ener * KbT
# now free_ener is in KJoule/mol

# compute spline of free energy
cs = CubicSpline(x_free_ener, free_ener)


# In[4]: FE from  Metadynamics


tmp=np.loadtxt("/Users/lorenzo/Google Drive/LORENZO/Files/gromos_spc-Ca-fes_from_hills.dat")

x_free_energy  = tmp[:,0]
free_energy = tmp[:,1]

free_energy=free_energy+abs(np.amin(free_energy))   # this free_energy is already in KJoule/mol
cs_free_energy = CubicSpline(x_free_energy, free_energy)

# In [5]: # Plot Free Energy A(x) from MD and compare with Metadynamics

fig = plt.figure()

xs = np.arange(md_x_inf, md_x_sup, 0.01)
ys = cs(xs) + abs(np.amin(cs(xs)))
ys2 = cs_free_energy(xs) + abs(np.amin(cs_free_energy(xs)))

ax0 = fig.add_subplot(121) 
ax0.plot(x_free_ener, free_ener, label="Boltzmann statistics of MD runs")  # Free energy from MD, discrete points
ax0.plot(x_free_energy, free_energy, label="Metadynamics")  # Free energy from MD, discrete points

ax1 = fig.add_subplot(122) # instead of plt.subplot(2, 2, 1)
ax1.plot(xs, ys, label="Boltzmann statistics of MD runs")  # Free energy from MD, discrete points
ax1.plot(xs, ys2, label="Metadynamics")  # Free energy from MD, discrete points

ax0.set_title('Free Energy')
ax1.set_title('Spline')
ax0.legend(loc='best')
ax1.legend(loc='best')
ax0.set_ylim([-np.max(free_ener)*0.1,np.max(free_ener)*1.5])
ax1.set_ylim([-np.max(free_ener)*0.1,np.max(free_ener)*1.5])

ax0.set_xlabel('Coordination Number')
ax1.set_xlabel('Coordination Number')
ax0.set_ylabel('Free Energy [KJoule/mol]')

plt.show()

# In[6] : Compute  and plot (Coordination dependent) Diffusion Coefficent

dynamics = data[skipped_index:]   # 1D array of MD simulation

#binning_numbers = [ 50 ]
data_number = len(dynamics[:-ilagtime])

md_msd  = np.zeros(((len(binning_numbers)),(np.max(binning_numbers)-1)))
md_nmsd = np.zeros(((len(binning_numbers)),(np.max(binning_numbers)-1)))
md_diff = np.zeros(((len(binning_numbers)),(np.max(binning_numbers)-1)))
md_x_diff = np.zeros(((len(binning_numbers)),(np.max(binning_numbers)-1)))



for j in range(len(binning_numbers)):
    md_x_bins = np.linspace(md_x_inf,md_x_sup,binning_numbers[j])
    for i in range((len(md_x_bins)-1)):
        md_x_diff[j][i] = (md_x_bins[1+i] + md_x_bins[i]) / 2



k=0 # counter inutile
print("will compute Diffusion of %d data points with these numbers of bins: %s" %(data_number,binning_numbers))
print("it should take 7 minutes with 3 bins numbers")
print("\n")
for i in range(data_number):
    x_ini = dynamics[i]
    x_end = dynamics[i+ilagtime]
    msd = (x_end - x_ini)**2

    if ( ( i % int(data_number/10)  ) == 0 ) :
            print(" %d / 10 of file analysed- %s" %(k, datetime.datetime.now()))
            k=k+1
    for j in range(len(binning_numbers)):
        nbins = binning_numbers[j]
        md_binstep = (md_x_sup - md_x_inf) / float(nbins-1)
          
        idx = int((x_ini - md_x_inf) / md_binstep)
        
        if (idx == (nbins-1)):
            idx = nbins-2
            
        md_msd[j][idx] = md_msd[j][idx] + msd
        md_nmsd[j][idx] = md_nmsd[j][idx] + 1
    
for j in range(len(binning_numbers)):
    for i in range(nbins-1):
        if md_nmsd[j][i] > 0 :
            md_msd[j][i] = md_msd[j][i] / md_nmsd[j][i]
        else :
            md_msd[j][i] = 0.0
        md_diff[j][i] = md_msd[j][i] / (2.0 * lagtime)



fig1 = plt.figure()
ax1 = fig1.add_subplot(211)
ax2 = fig1.add_subplot(212)

for j in range(len(binning_numbers)):

    nmax = binning_numbers[j]
    if md_x_diff[j][nmax-2] == 0 :
        nmax = nmax - 1
    ax1.plot(md_x_diff[j][:nmax-1], md_diff[j][:nmax-1], label=str(nmax)+' bins', linewidth=1.5)
    ax1.set_xlim(np.min(md_x_diff[j][:])-1.0,np.max(md_x_diff[j][:])+0.25)
    ax1.legend()
    ax1.legend(loc='best')
    ax1.set_ylabel('diffusion coefficent [1/ps]')
    plt.title('Diffusion')
for j in range(len(binning_numbers)):
    nmax = binning_numbers[j]
    ax2.semilogy(md_x_diff[j][:nmax-1], md_nmsd[j][:nmax-1], label=str(nmax)+' bins', linewidth=1.5)

ax2.set_xlim(np.min(md_x_diff[j][:])-1.0,np.max(md_x_diff[j][:])+0.25)
ax2.set_xlabel('coordination number')
ax2.set_ylabel('Number of events')
ax2.legend()
plt.show()


#Sanity Check
plt.figure()
plt.plot((md_x_diff[j][:nmax-1]), (md_diff[j][:nmax-1])/md_diff.max(), label='Diffusivity')
plt.plot(xs, ys/ys.max(),label='Free Energy')
plt.legend()
plt.show()


discretization_chosen=0  # corresponds to the 4th line in previous 2 plots

# Plot final discretized diffusion chosen, D(x), and corresponding spline function
print("final diffusion chosen is the one with %d bins" %(binning_numbers[discretization_chosen]))
md_x_diff_final=md_x_diff[discretization_chosen][:binning_numbers[discretization_chosen]-1]
md_diff_final=md_diff[discretization_chosen][:binning_numbers[discretization_chosen]-1]

LD_max_allowed_value=xs[-1]

fig = plt.figure()
axes = fig.add_subplot(111)

cs_md_diff = CubicSpline(md_x_diff_final,md_diff_final)

xs = np.arange(md_x_inf,md_x_sup, 0.01)
ys_md_diff = cs_md_diff(xs)
for i in range(len(xs)):
    
    if ys_md_diff[i] < 0 :
        print("\n\nAttention! value %d of diffusion constant is negative and will be set to diffusion average" %i)
        ys_md_diff[i] =0
        
        plt.axvline(x=xs[i-1], color="red")
        plt.annotate('negative diffusion: '+str(format(xs[i-1], ".1f")),
             xy=(xs[i-1], np.max(ys_md_diff)*0.9), xycoords='data',
             xytext=(-140, +30), textcoords='offset points', fontsize=16, color='red',
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))


for i in range(len(xs)):
    if ys_md_diff[i] == 0:
        ys_md_diff[i]=np.average(ys_md_diff)


axes.plot(md_x_diff_final, md_diff_final, label="Diff")
axes.plot(xs, ys_md_diff, label="Spline Diff")


axes.set_xlim((np.min(xs)-(np.max(xs)-np.min(xs))*0.1),np.max(xs)+(np.max(xs)-np.min(xs))*0.1)
axes.set_ylim(0,np.max(ys_md_diff)*1.2)
axes.set_xlabel('coordination number')
axes.set_ylabel('diffusion constant [1/ps]')
axes.legend()
plt.show()




# In[7]:

# ## define 3 minima between which we will check the dynamics

x_ini   = 7.9
x_delta = 0.2 # discretization of order parameter for MFPT
x_left = 7.0
x_right = 8.75




#names = [ "left minimum", "initial point", "right minimum" ]
names = [ "Left", "Start", "Right" ]
xpoints = [ x_left, x_ini, x_right ]


#print(x_free_energy)

def Closest_value( array_x, value):
    return np.argmin(abs(array_x- value))


def Maxdistance( array, index_start, index_end ):  # calculates the maximum POSITIVE distance (array[j]-array[index_start]) with index_start<j<=index_end
    
    if index_start > index_end :
        index0=index_end
        index_end=index_start
        index_start=index0
    
    tmp=array [index_start : index_end]        
    return np.max(tmp-array[index_start])
    
        
fig = plt.figure()
ax2 = fig.add_subplot(111)  
ax2.plot(x_free_ener, free_ener, label="Boltzmann statistics of MD runs")  # Free energy from MD, discrete points
ax2.plot(x_free_energy, free_energy, label="Metadynamics")  # Free energy from MD, discrete points
ax2.plot(xs, ys, label="Spline of MD runs")  # Free energy from MD, discrete points
ax2.legend()
ax2.legend(loc='best')
ax2.set_xlabel('coordination number')
ax2.set_ylabel('free energy [KJoule/mol]')
ax2.set_title('Location of the minima in the free energy')
maxplot=np.max(free_energy)
maxplot2=np.max(free_ener)
if maxplot2 > maxplot :
    maxplot=maxplot2
ax2.set_ylim(0-maxplot*0.2,maxplot*1.2)

where_to_put_y = np.zeros(3)
where_to_put_y2 = np.zeros(3)
where_to_put_x = np.zeros(3)
index_free_energy = np.zeros(3,dtype='int')
index_free_ener = np.zeros(3,dtype='int')
where_to_put_x[0]=x_left; 
where_to_put_x[1]=x_ini; 
where_to_put_x[2]=x_right


for i in range(3):
    index_free_energy[i]=Closest_value( x_free_energy, where_to_put_x[i])
    where_to_put_y[i]=free_energy[Closest_value( x_free_energy, where_to_put_x[i])]
    index_free_ener[i]=Closest_value( x_free_ener, where_to_put_x[i])
    where_to_put_y2[i]=free_ener[Closest_value( x_free_ener, where_to_put_x[i])]
    
print("                                   A(x_left)-A(x_ini)         A(x_right)-A(x_ini)")

for i in range(3):
        
        plt.annotate(str(where_to_put_x[i])+', '+str(format(where_to_put_y[i], ".1f")),
             xy=(where_to_put_x[i], where_to_put_y[i]), xycoords='data',
             xytext=(+20, +30), textcoords='offset points', fontsize=16, color='#ff7f0e',
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

        plt.annotate(str(where_to_put_x[i])+', '+str(format(where_to_put_y2[i], ".1f")),
             xy=(where_to_put_x[i], where_to_put_y2[i]), xycoords='data',
             xytext=(+20, -30), textcoords='offset points', fontsize=16, color='#1f77b4',
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
        
print("\nMD- point to point                  %.2f KJoule/mol           %.2f KJoule/mol"
      %(where_to_put_y[0]-where_to_put_y[1],where_to_put_y[2]-where_to_put_y[1]) )
print("MD- max delta in range              %.2f KJoule/mol           %.2f KJoule/mol"
      %(Maxdistance(free_ener, index_free_ener[1], index_free_ener[0]),
        Maxdistance(free_ener, index_free_ener[1], index_free_ener[2])))
print("\nMetadynamics- point to point        %.2f KJoule/mol           %.2f KJoule/mol"
      %(where_to_put_y2[0]-where_to_put_y2[1],where_to_put_y2[2]-where_to_put_y2[1]) )
print("Metadynamics- max delta in range    %.2f KJoule/mol           %.2f KJoule/mol"
      %(Maxdistance(free_energy, index_free_energy[1], index_free_energy[0]),
        Maxdistance(free_energy, index_free_energy[1], index_free_energy[2])))

plt.show()


'''
fileout = open('MFPT-MD.txt','w') 
fileout.write("                                   A(x_left)-A(x_ini)         A(x_right)-A(x_ini)\n")
fileout.write("\nMD- point to point                  %.2f KJoule/mol           %.2f KJoule/mol\n"
      %(where_to_put_y[0]-where_to_put_y[1],where_to_put_y[2]-where_to_put_y[1]) )
fileout.write("MD- max delta in range              %.2f KJoule/mol           %.2f KJoule/mol\n"
      %(Maxdistance(free_ener, index_free_ener[1], index_free_ener[0]),
        Maxdistance(free_ener, index_free_ener[1], index_free_ener[2])))
fileout.write("\nMetadynamics- point to point        %.2f KJoule/mol           %.2f KJoule/mol\n"
      %(where_to_put_y2[0]-where_to_put_y2[1],where_to_put_y2[2]-where_to_put_y2[1]) )
fileout.write("Metadynamics- max delta in range    %.2f KJoule/mol           %.2f KJoule/mol\n"
      %(Maxdistance(free_energy, index_free_energy[1], index_free_energy[0]),
        Maxdistance(free_energy, index_free_energy[1], index_free_energy[2])))

fileout.close()
'''




# In[8]: Compute Mean First Passage Time from central minimum x_ini to left x_left and right x_right minima, from MD





def fitfunction(t, tau, A):
    # tau is in ps only if t is
    initial_size=A # self-correlation with self is 1
    y= initial_size*np.exp(-t/tau)
    return y


#Compute MFPT from x_ini to x_left and x_right FROM MD


xl = x_ini - x_delta    # will consider only starting points in the range   x_ini-x_delta,x_ini+x_delta
xr = x_ini + x_delta

max_time=0.0
mfpt_l = []
mfpt_r = []
idx = 0

for i in range(skipped_index,len(data)):
    if i < idx:
        continue
    if data[i] > xl and data[i] < xr:
        
        for j in range(i+1,len(data)):
            if data[j] > x_right:
                mfpt_r.append(j - i)
                idx = j
                break
            if data[j] < x_left:
                mfpt_l.append(j - i)
                idx = j
                break

mfpt_l = [x * md_dt for x in mfpt_l]
mfpt_r = [x * md_dt for x in mfpt_r]

#%%
fig = plt.figure()
axes = fig.add_subplot(121)
t_nbins_fit=20

max_time=np.max([np.max(mfpt_l), np.max(mfpt_r)])
t_bins = np.linspace(0.1,max_time,t_nbins_fit)

hist, _ = np.histogram(mfpt_l, bins=t_bins, density=True)
init_vals = [1.0, hist[0]] #[interval,100,1/100] 
best_valsR, covarR = curve_fit(fitfunction, (t_bins[1:]+t_bins[:-1])/2, hist, p0=init_vals)

print("max time obtained for left minimum: %8.2f ps \n" %(np.max(mfpt_l)))
print("max time obtained for right minimum: %8.2f ps \n" %(np.max(mfpt_r)))

axes.set_yscale('log')
axes.bar(t_bins[:-1], hist, alpha = 0.5, width=t_bins[0]-t_bins[1] ) #, label = "$P_{MD}(x)$ " + str(interval) )
axes.plot((t_bins[1:]+t_bins[:-1])/2, fitfunction((t_bins[1:]+t_bins[:-1])/2, *best_valsR),'r')
axes.set_title('MFPT- Center (%.1f) to Left (%.1f)' %(x_ini,x_left))
axes.set_ylabel('probability')
axes.set_xlabel('time [ps]')
ax2 = fig.add_subplot(122)

hist, _ = np.histogram(mfpt_r, bins=t_bins, density=True)
init_vals = [1.0, hist[0]] #[interval,100,1/100]
    
best_valsL, covarL = curve_fit(fitfunction, (t_bins[1:]+t_bins[:-1])/2, hist, p0=init_vals) # , bounds=((interval-epsilon,-np.inf,-np.inf), (interval+epsilon,np.inf,np.inf)) )
ax2.set_yscale('log')
ax2.bar(t_bins[:-1], hist, alpha = 0.5,  width=t_bins[0]-t_bins[1])  #, label = "$P_{MD}(x)$ " + str(interval) )
ax2.plot((t_bins[1:]+t_bins[:-1])/2, fitfunction((t_bins[1:]+t_bins[:-1])/2, *best_valsL),'r')
ax2.set_title('MFPT- Center (%.1f) to Right (%.1f)' %(x_ini,x_right))
ax2.set_xlabel('time [ps]')
plt.show()
#fig.savefig('MFPT-MD.jpg')


print("We had %d events from Center to left,  mean time: %.3f ps      tau: %.2f +/- %.2f ps" %(len(mfpt_l), np.mean(mfpt_l),best_valsL[0],np.sqrt(covarL[0,0])))
print("We had %d events from Center to right, mean time: %.3f ps      tau: %.2f +/- %.2f ps" %(len(mfpt_r), np.mean(mfpt_r),best_valsR[0],np.sqrt(covarR[0,0])))


'''
fileout = open('MFPT-MD.txt','w') 
fileout.write("MD- We had %d events from Center to left,  mean time: %.3f ps      tau: %.2f +/- %.2f ps\n" %(len(mfpt_l), np.mean(mfpt_l),best_valsL[0],np.sqrt(covarL[0,0])))
fileout.write("MD- We had %d events from Center to right, mean time: %.3f ps      tau: %.2f +/- %.2f ps\n" %(len(mfpt_r), np.mean(mfpt_r),best_valsR[0],np.sqrt(covarR[0,0])))
fileout.close()
'''

    # In[]
# Experiments w/ Diffusion Zones


i_lag_zones=4
lagtime_zones=i_lag_zones*md_dt

data_number = len(dynamics[:-i_lag_zones])

#Zones=[dynamics.min(),6.9, 7.3,7.7, 8.3, 8.6, 8.9, dynamics.max()]
#Zones=np.array([dynamics.min(),6.9, 7.1,7.3,7.4, 7.5,7.7,7.8,7.9,8,8.1,8.2, 8.4, 8.6,8.75, 8.9, dynamics.max()])
#Zones=np.array([dynamics.min(),6.9,7, 7.1,7.25,7.3,7.4, 7.5,7.55,7.6,7.7,7.8,7.87,7.9, 8,8.1,8.2,8.3, 8.4,8.5, 8.6,8.75, 8.9, dynamics.max()])
Zones=np.array([dynamics.min(),6.9,7, 7.1,7.25,7.3,7.4, 7.5,7.55,7.6,7.7,7.75,7.8,7.82,7.87,7.9,7.93,7.97, 8,8.1,8.2,8.3, 8.4,8.5, 8.6,8.75, 8.9, dynamics.max()])

md_msd  = np.zeros(len(Zones)-1)
md_nmsd = np.zeros(len(Zones)-1)
md_diff = np.zeros(len(Zones)-1)
md_x_diff = np.zeros(len(Zones)-1)


k=0 # counter inutile
print("will compute Diffusion of %d data points with these numbers of zones: %s" %(data_number,len(Zones)))
print("\n")
for i in range(data_number):
    x_ini = dynamics[i]
    x_end = dynamics[i+i_lag_zones]
    msd = (x_end - x_ini)**2

    if ( ( i % int(data_number/10)  ) == 0 ) :
            print(" %d / 10 of file analysed w/ %d zones - %s" %(k,len(Zones), datetime.datetime.now()))
            k=k+1
    
        
          
    for j in range(len(Zones)-1):
        if x_ini< Zones[(j+1)%len(Zones)] and x_ini>= Zones[j]:
            idx=j
            break
        elif j==len(Zones)-2:
            idx=j
            break
            
    md_msd[idx] = md_msd[idx] + msd
    #print('%d %f'%(idx, msd))
    md_nmsd[idx] = md_nmsd[idx] + 1
    
#print(md_msd)

for i in range(len(Zones)-1):
    if md_nmsd[i] > 0 :
        md_msd[i] = md_msd[i] / md_nmsd[i]
    else :
        md_msd[i] = 0.0       
    md_diff[i] = md_msd[i] / (2.0 * lagtime_zones)



fig1 = plt.figure()
ax1 = fig1.add_subplot(211)
ax2 = fig1.add_subplot(212)

x=((Zones[1:]+Zones[:-1])*0.5)

ax1.plot(x, md_diff)
#ax1.set_xlim(np.min(md_x_diff[j][:])-1.0,np.max(md_x_diff[j][:])+0.25)
#ax1.legend(loc='best')
ax1.set_ylabel('diffusion coefficent [1/ps]')
plt.title('Diffusion w/ %d zones' %len(Zones))


ax2.semilogy(x, md_nmsd)
#ax2.set_xlim(np.min(md_x_diff[j][:])-1.0,np.max(md_x_diff[j][:])+0.25)
ax2.set_xlabel('coordination number')
ax2.set_ylabel('Number of events')
#ax2.legend()
plt.show()


#Sanity Check
plt.figure()
plt.plot(x, md_diff/md_diff.max(), label='Diffusivity')
plt.plot(xs, ys/ys.max(),label='Free Energy')
plt.legend()
plt.show()


cs_zones=CubicSpline(x, md_diff)

plt.figure()

plt.plot(x, md_diff, label='Diffusivity')
plt.plot(x, cs_zones(x),'r', label='Spline' )
plt.show()

# In[]
#  Hummer diffusion

from scipy.linalg import logm

binning=np.linspace(dynamics.min(), dynamics.max(), 30)  #original 25
binning_space=binning[0]-binning[1]

print('Computing diffuson with Hummer Method for %d bins' %(len(binning)-1))

n_lag=4

lag_H=n_lag*md_dt

KbT = 0.0083144621 * 300   #KJoule/mol

def exp_R_Hummer(n_lag, binning):
    
    exp_R=np.zeros((len(binning)-1, len(binning)-1))
    binning_space=binning[0]-binning[1]
    
    ix=-1*np.ones(len(dynamics[:-n_lag]))
    for i in range(len(binning)-1):
        ix=np.where(np.logical_and(dynamics[:-n_lag]< binning[(i+1)] , dynamics[:-n_lag]>= binning[i]), i, ix)
    #print(min(ix))
    ix=np.where(ix==-1, len(binning)-2, ix)
    
    ix_delta=np.zeros(len(dynamics[:-n_lag]))
    xj=ix_delta

    print('Computing exponential Hummer Matrix')

    for n in range(len(dynamics[:-n_lag])):
        ix_delta[n]=dynamics[n]- binning[int(ix[n])]+ 0.5*binning_space
        xj[n]=dynamics[n+n_lag]-ix_delta[n]
    
    jx=-1*np.ones(len(dynamics[:-n_lag]))
    for i in range(len(binning)-1):
        jx=np.where(np.logical_and(xj< binning[(i+1)] , xj>= binning[i]), i, jx)
    jx=np.where(jx==-1, len(binning)-2, jx)
 
    print('Initialising exp(tR) at t = %d' %(n_lag))
    for n in range(len(dynamics[:-n_lag])):
        exp_R[int(jx[n])][int(ix[n])]+=1
    
    '''     
    #Traduzione plain di come lo fa bracncato nello script c, NON usare perchè non ottimizzato   
    for n in range(len(dynamics[:-n_lag])):
        x_0=dynamics[n]
        for i in range(len(binning)):
            if x_0< binning[(i+1)%len(binning)] and x_0>= binning[i]:
                ix=i
                break
            elif i==len(binning)-1:
                ix=i
                break
   
        ix_delta=dynamics[n]- binning[ix]+ 0.5*binning_space
        xj=dynamics[n+n_lag]-ix_delta
       
        jx=0
    
        for j in range(len(binning)):
            if xj< binning[(j+1)%len(binning)] and xj>= binning[j]:
                jx=j
                break
            elif j==len(binning)-1:
                jx=j
                break
        
        if jx==0 or ix==0:
            print('Problem with Hammer at %d' %n)
        else:
            exp_R[jx][ix]+=1
        '''
    Ri=np.zeros(len(binning)-1)

    for i in range(len(exp_R)):
        for j in range(len(exp_R)):
            Ri[i] = Ri[i] + exp_R[j][i]
        
    for i in range(len(exp_R)):
       for j in range(len(exp_R)):
           if Ri[j]!=0:
               exp_R[i][j]= float(float(exp_R[i][j] )/ float(Ri[j]))
        
    return exp_R

exp_R=np.zeros((len(binning)-1, len(binning)-1))
exp_R=exp_R_Hummer(n_lag, binning)
p,_=np.histogram(dynamics[skipped_index:], bins=binning, density=True)
'''
#Così è come la calcoma Hummer nell'articolo, ma ci sono problemi:
    #la matrice exp_R_sym (ezp_(tR~)) NON è simmetrica come invece dovrebbe
    #se la simmetrizziamo noi a mano gli autovalori (w) NON sono strettamente positivi
    #quindi non ne posso prendere il logaritmo (in teoria DEVONO essere positivi perchè R~ è simmetrica e ne prendo l'esponenziale)

exp_R_sym=np.zeros([len(binning)-1,len(binning)-1])
for i in range(len(exp_R)):
    for j in range(len(exp_R)):
        exp_R_sym[i][j]=np.sqrt(1/p[j])*exp_R[i][j]*np.sqrt(p[i])

#La matrice DOVREBBE essere simmetrica, ma NON lo è
#Quindi la SIMMETRIZZIAMO noi

for i in range(len(exp_R)):
    for j in range(i+1):
        exp_R_sym[i][j]=exp_R_sym[j][i]=float((exp_R_sym[i][j]+exp_R_sym[j][i])/2)


w, v=np.linalg.eig(exp_R_sym)
print((w))
'''

print('Computing time derivtive O(h^4)')
dt_exp_R=(-1*exp_R_Hummer(n_lag+2, binning)+8*exp_R_Hummer(n_lag+1, binning)-8*exp_R_Hummer(n_lag-1, binning)+exp_R_Hummer(n_lag-2, binning))/(12*md_dt)
print('Inverting exp_(tR) to get R=(d/dt exp(tR))*(exp(-tR)')
R=np.dot(dt_exp_R, np.linalg.inv(exp_R))



print('Computing Free Energy and Diffusivity with Hummer method for %d bins' %(len(binning-1)))
Delta_Q=binning[1]-binning[0]
F_Hum=-KbT*np.log(p/Delta_Q)
D_Hum=np.zeros(len(binning)-2)

for i in range(len(D_Hum)):
    D_Hum[i]=(Delta_Q**2)*(R[i+1][i]*np.sqrt(p[i]/p[i+1]))#+R[i][i+1]*np.sqrt(p[i+1]/p[i]))/2
#%%
    
binning1=np.linspace(dynamics.min(), dynamics.max(), 30)  #27 con n_lag=3 riproduce la msd, 30 valore migliore
n_lag1=3
D1=np.zeros(len(binning1)-2)
R1=(1/(n_lag1*md_dt))*logm(exp_R_Hummer(n_lag1, binning1))   

Delta_Q1=binning1[1]-binning1[0]
p1,_=np.histogram(dynamics[skipped_index:], bins=binning1, density=True)


for i in range(len(D1)):
    D1[i]=(Delta_Q1**2)*(R1[i+1][i]*np.sqrt(p1[i]/p1[i+1]))#+R1[i][i+1]*np.sqrt(p1[i+1]/p1[i]))/2
    #D2[i]=(Delta_Q1**2)*(R1[i][i+1]*np.sqrt(p1[i+1]/p1[i]))
    #Manually enforcing detailed balance



fig=plt.figure()
plt.title('Results from Hummer Method')
ax1=fig.add_subplot(121)
x=(binning[1:]+binning[:-1])/2 #bin centers
x1=(binning1[1:]+binning1[:-1])/2 #bin centers
ax1.plot(x, cs_free_energy(x), label='From Metadynamics')
ax1.plot(x, F_Hum, label='From Hummer')
ax1.legend()
ax2=fig.add_subplot(122)
x=(x[1:]+x[:-1])/2 #middle points between bin centers
x1=(x1[1:]+x1[:-1])/2
ax2.plot(x, cs_md_diff(x), label='From MD')
ax2.plot(x, D_Hum, label='From Hummer')
ax2.plot(x1, D1, label='From scipy')
ax2.legend()
plt.show()

cs_Humm=CubicSpline(x1, D1)


'''
np.save('MD_spline', cs_md_diff)
np.save('Zones_spline', cs_zones)
np.save('Humm_spline', cs_Humm)
np.save('FE_spline',cs_free_energy)
'''

# In[]: betailed balance check

err=np.zeros((len(binning1)-1, len(binning1)-1))

for i in range(len(err)):
    err[i][i]=R1[i][i]
    for j in range(len(err)):
        if j!=i:
            err[i][i]+=R[j][i]
    for j in range(i+1, len(err)):
        err[i][j]=err[j][i]=R1[i][j]-R1[j][i]*float(p1[i]/p1[j])
        

plt.figure()
plt.imshow(err, cmap='jet')
#plt.imshow(R1)
plt.colorbar()
plt.show()

# In[9]:  LD (Langevin Dynamics) options


""" both short and long LDs """

ld_t_ini                = 0.0
#ld_dt                   = 0.004       # same as the one of md
ld_sqrt_dt              = 0.063245553203367

""" long LD"""
#ld_t_end_long           = 10000.0     # long ld end (picoseconds)   # more than 1.2 microsec hangs the ram
ld_y_ini                = x_ini
#ld_fpt_y_left          = x_left      # Stop ld at y = ld_fpt_y
#ld_fpt_y_right         = x_right     # Stop ld at y = ld_fpt_y

ld_max_allowed_y_long  = 9.5  # max value admitted before stopping ld
ld_min_allowed_y_long  = 5.0          # min value admitted before stopping ld


""" short LDs"""
#ld_t_end_short         = 500.0 # picoseconds
#ld_ntrials             = 100
ld_max_allowed_y_short = x_right+0.2  # max value admitted before stopping ld
ld_min_allowed_y_short = x_left-0.2   # min value admitted before stopping ld





# In[10]: LD propagator function- python version


# langevin dynamics propagator requires 3 parameters: ld_D, ld_D1, ld_U1,
#    calculated from MD (D,D1) or metadynamics (U1). 
    
#    Since we do not YET have converged metadynamics runs, we use the boltzmann
#    statistics of the MD for now
    
    # before, we were setting it to 1 (wrong)
    #x_bins = np.linspace(ld_fpt_y_left,ld_fpt_y_right,nbins)
    #g = np.linspace(1.0,1.0,nbins)
    #cs_free_energy = CubicSpline(x_bins, g)
    
# in original propagator function
#ld_U1 = cs_free_energy(ld_yold,1)     # this way, the potential is taken from Metadynamics (as Giuseppe did)
# in propagator function, now (by Tommy)
#ld_U1 = cs(ld_yold,1)        # this way, the potential is taken from Boltzmann statistics on MD runs








def Propagate(ld_y_ini, ld_steps, ld_dt, ld_yn, dW, ld_fpt_y_left, ld_fpt_y_right, action_for_extremes):
    """ Initialize some variables """
    ld_yn[0] = ld_y_ini
    ld_fpt_t_l = 0.0
    ld_fpt_t_r = 0.0
    #ld_sqrt_dt = np.sqrt(ld_dt)
    percentage_analysed = 0

    
    
    for i in range(1,ld_steps):
        ld_yold = ld_yn[i-1]

        # Problem specific parameters

        ld_D  = cs_md_diff(ld_yold)             #  [1/ps]
        ld_D1 = cs_md_diff(ld_yold,1)           #  [1/ps]
        ld_U1 = cs_free_energy(ld_yold,1)/KbT  # this way, the potential is taken from Metadynamics (as Giuseppe did)
        #ld_U1 = cs(ld_yold,1)/KbT               #  [num]     this way, the potential is taken from Boltzmann statistics on MD runs

        # End of Problem specific parameters

        
        
        
        ld_a  = (ld_D1 - ld_D * ld_U1)  
        ld_b  = np.sqrt(2.0*ld_D)

        # Euler
        ld_ynew = ld_yold + ld_a * ld_dt + ld_b * dW[i]
        
        
        if ( ld_steps == int((ld_t_end_long - ld_t_ini) / ld_dt) ): # do only for long LD
            if ( ( i % int(ld_steps/10)  ) == 0 ) :
                percentage_analysed=percentage_analysed+1
                print(" %d / 10 of LD dynamics ran- %s" %(percentage_analysed, datetime.datetime.now()))

        ld_yn[i] = ld_ynew


        
        if ld_ynew < ld_fpt_y_left:
            if action_for_extremes == "reflect" :
                #ld_ynew = ld_ynew + 2*abs(ld_fpt_y_left-ld_y_new)
                ld_ynew = 2*ld_fpt_y_left-ld_ynew
            elif action_for_extremes == "count" :
                ld_fpt_t_l = i * ld_dt
                #if ld_fpt_y_left == ld_min_allowed_y_long :
                #print("attention- particle from LD escaped from range [%.2f,%.2f] at time %.3f ns- last value %.2f" %(ld_fpt_y_left,ld_fpt_y_right,i*ld_dt/1000,ld_ynew))
                break
            else :
                print("error! action_for_extremes not defined")
                break

        elif ld_ynew > ld_fpt_y_right:
            if action_for_extremes == "reflect" :
                #ld_ynew = ld_ynew - 2*abs(ld_fpt_y_right-ld_y_new)
                ld_ynew = 2*ld_fpt_y_right-ld_ynew
            elif action_for_extremes == "count" :
                ld_fpt_t_r = i * ld_dt
                #if ld_fpt_y_right == ld_max_allowed_y_long :
                #print("attention- particle from LD escaped from range [%.2f,%.2f] at time %.3f ns- last value %.2f" %(ld_fpt_y_left,ld_fpt_y_right,i*ld_dt/1000,ld_ynew))
                break
            else :
                print("error! action_for_extremes not defined")
                break
    return ld_fpt_t_l, ld_fpt_t_r 


def Propagate_NEW(ld_y_ini, ld_steps, ld_dt, ld_NEW, dW, ld_fpt_y_left, ld_fpt_y_right, action_for_extremes):
    ld_yn[0] = ld_y_ini
    ld_fpt_t_l = 0.0
    ld_fpt_t_r = 0.0
    #ld_sqrt_dt = np.sqrt(ld_dt)
    percentage_analysed = 0

    chosen_spline_diff=cs_md_diff   #cs_md_diff     for local MSD w/ even binning
                                 #cs_zones       for local MSD w/ zones
                                 #cs_Humm        for Hummer method
    
    for i in range(1,ld_steps):
        ld_yold = ld_yn[i-1]

        # Problem specific parameters

        ld_D  = chosen_spline_diff(ld_yold)             #  [1/ps]
        ld_D1 = chosen_spline_diff(ld_yold,1)           #  [1/ps]
        ld_U1 = cs_free_energy(ld_yold,1)/KbT  # this way, the potential is taken from Metadynamics (as Giuseppe did)
        #ld_U1 = cs(ld_yold,1)/KbT               #  [num]     this way, the potential is taken from Boltzmann statistics on MD runs

        # End of Problem specific parameters

        
        
        
        ld_a  = (ld_D1 - ld_D * ld_U1)  
        ld_b  = np.sqrt(2.0*ld_D)

        # Euler
        ld_ynew = ld_yold + ld_a * ld_dt + ld_b * dW[i]
        
        
        if ( ld_steps == int((ld_t_end_long - ld_t_ini) / ld_dt) ): # do only for long LD
            if ( ( i % int(ld_steps/10)  ) == 0 ) :
                percentage_analysed=percentage_analysed+1
                print(" %d / 10 of LD dynamics ran- %s" %(percentage_analysed, datetime.datetime.now()))

        ld_yn[i] = ld_ynew


        
        if ld_ynew < ld_fpt_y_left:
            if action_for_extremes == "reflect" :
                #ld_ynew = ld_ynew + 2*abs(ld_fpt_y_left-ld_y_new)
                ld_ynew = 2*ld_fpt_y_left-ld_ynew
            elif action_for_extremes == "count" :
                ld_fpt_t_l = i * ld_dt
                #if ld_fpt_y_left == ld_min_allowed_y_long :
                #print("attention- particle from LD escaped from range [%.2f,%.2f] at time %.3f ns- last value %.2f" %(ld_fpt_y_left,ld_fpt_y_right,i*ld_dt/1000,ld_ynew))
                break
            else :
                print("error! action_for_extremes not defined")
                break

        elif ld_ynew > ld_fpt_y_right:
            if action_for_extremes == "reflect" :
                #ld_ynew = ld_ynew - 2*abs(ld_fpt_y_right-ld_y_new)
                ld_ynew = 2*ld_fpt_y_right-ld_ynew
            elif action_for_extremes == "count" :
                ld_fpt_t_r = i * ld_dt
                #if ld_fpt_y_right == ld_max_allowed_y_long :
                #print("attention- particle from LD escaped from range [%.2f,%.2f] at time %.3f ns- last value %.2f" %(ld_fpt_y_left,ld_fpt_y_right,i*ld_dt/1000,ld_ynew))
                break
            else :
                print("error! action_for_extremes not defined")
                break
    return ld_fpt_t_l, ld_fpt_t_r 

# In[11]: perform Long LD from one minimum to adjacent ones

if Sym_type == 'long' or Sym_type == 'both':

    print("Performing long LD- %d ns- should take 30 minutes for a 100ns run" %(ld_t_end_long/1000))


    ld_steps  = int((ld_t_end_long - ld_t_ini) / ld_dt)
    ld_yn    = np.zeros(ld_steps)

    dW = np.random.normal(loc = 0.0, scale=ld_sqrt_dt,size=ld_steps)    


    if LD_max_allowed_value < ld_max_allowed_y_long :
        print("\n\nATTENTION: max value allowed for long LD (%.2f) exceeds the range in which D is defined- setting it to %.2f\n\n" %(ld_max_allowed_y_long,LD_max_allowed_value))
        ld_max_allowed_y_long=LD_max_allowed_value


    fpt_t_l, fpt_t_r = Propagate(ld_y_ini, ld_steps, ld_dt, ld_yn, dW, ld_min_allowed_y_long, ld_max_allowed_y_long, "reflect")


    fig = plt.figure()
    ax0 = fig.add_subplot(111) 

    xperquestoplot=np.linspace(0,float(ld_t_end_long/1000),float(ld_t_end_long/ld_dt)) #so it is in ns
    ax0.plot(xperquestoplot,ld_yn)
    ax0.set_ylabel('Coordination Number')
    ax0.set_xlabel('time [ns]')
    ax0.set_title('Coordination Number obtained in LD long dynamics')
    #plt.legend(loc='best')
    plt.show()


# In[12]: From Long LD recover A(x)

if Sym_type == 'long' or Sym_type == 'both':
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    t_nbins = binning_numbers[discretization_chosen] # same as the one used for Diffusion constant
    t_bins = np.linspace(ld_min_allowed_y_long, ld_max_allowed_y_long, t_nbins)

    ld_yn=ld_yn[~np.isnan(ld_yn)]  # let's take over all these NaN, shall we?
    hist, _ = np.histogram(ld_yn, bins=t_bins, density=True)
    tmp = -1.0 * np.log(hist)
    tmp = tmp + abs(np.amin(tmp))
    tmp = tmp * KbT

    ax1.plot(x_free_ener, free_ener, label="Obtained from MD run")  # Free energy from MD, discrete points
    ax1.plot(t_bins[:t_nbins-1]+(t_bins[1]-t_bins[0])/2, tmp, label="Recovered from LD run")
    ax1.set_ylabel('Free Energy [KJoule/mol]')
    ax1.set_xlabel('coordination number')


    ypoints = [ cs(xpoints[0]), cs(xpoints[1]), cs(xpoints[2]) ]
    ax1.scatter(xpoints, ypoints, s=100, color='red')
    for i, txt in enumerate(names):
        ax1.annotate(txt, xy=(xpoints[i], ypoints[i]), xycoords='data',
                     xytext=(+0, +30), textcoords='offset points', fontsize=16, color='red')

    plt.legend(loc='best')
    plt.show()


# In[21]:  From Long LD recover D(x)

if Sym_type == 'long' or Sym_type == 'both':

    fig = plt.figure()
    axes = fig.add_subplot(111)

    tothist=np.sum(hist)
    axes.bar(t_bins[:-1]+(t_bins[1]-t_bins[0])/2,hist/tothist, width=(t_bins[1]-t_bins[0]))
    axes.plot(t_bins[:-1]+(t_bins[1]-t_bins[0])/2,hist/tothist, '--k')
    tempspline=CubicSpline(t_bins[:-1]+(t_bins[1]-t_bins[0])/2, hist/tothist)
    ypoints = [ tempspline(xpoints[0]), tempspline(xpoints[1]), tempspline(xpoints[2]) ]
    axes.scatter(xpoints, ypoints, s=100, color='red')
    for i, txt in enumerate(names):
        axes.annotate(txt, xy=(xpoints[i], ypoints[i]), xycoords='data',
                      xytext=(+15, +10), textcoords='offset points', fontsize=16, color='red')

    axes.set_title('histogram of the coordination numbers encountered in Langevin Dynamics')
    #plt.legend(loc='best')
    plt.show()


    print("using the same discretization as before- %d bins" %(binning_numbers[discretization_chosen]))
    nbins=binning_numbers[discretization_chosen]

    ld_x_inf = np.min(ld_yn)
    ld_x_sup = np.max(ld_yn) + 0.01
    ld_x_bins = np.linspace(ld_x_inf,ld_x_sup,nbins)
    ld_binstep = (ld_x_sup - ld_x_inf) / float(nbins-1)

    print("Analysing previous Langevin Dynamics, that was between %.2f and %.2f, to find diffusion constant" %(ld_x_inf,ld_x_sup))


    ld_msd  = np.zeros(nbins-1) 
    ld_nmsd = np.zeros(nbins-1)
    ld_diff = np.zeros(nbins-1) 

    ld_x_diff = (ld_x_bins[1:] + ld_x_bins[:-1]) / 2


    for i in range(len(ld_yn[:-ilagtime])):
        x_ini = ld_yn[i]
        x_end = ld_yn[i+ilagtime]
        msd = (x_end - x_ini)**2
    
        idx = int((x_ini - ld_x_inf) / ld_binstep)
       
        ld_msd[idx] = ld_msd[idx] + msd
        ld_nmsd[idx] = ld_nmsd[idx] + 1
    
    for i in range(nbins-1):
        ld_msd[i] = ld_msd[i] / ld_nmsd[i]
        ld_diff[i] = ld_msd[i] / (2.0 * lagtime)
    
    fig = plt.figure()
    axes = fig.add_subplot(111)

    cs_md_diff = CubicSpline(md_x_diff_final,md_diff_final)
    axes.plot(ld_x_diff, ld_diff, label="Langevin Dynamics")
    axes.plot(md_x_diff_final, md_diff_final, label="Molecular Dynamics")

    tempspline=CubicSpline(ld_x_diff, ld_diff)
    ypoints = [ tempspline(xpoints[0]), tempspline(xpoints[1]), tempspline(xpoints[2]) ]
    axes.scatter(xpoints, ypoints, s=100, color='red')
    for i, txt in enumerate(names):
        axes.annotate(txt, xy=(xpoints[i], ypoints[i]), xycoords='data',
                      xytext=(+0, +30), textcoords='offset points', fontsize=16, color='red')

    axes.set_xlim([t_bins[0],t_bins[-1]])
    axes.set_xlabel('coordination number')
    axes.set_ylabel('diffusion constant [1/ps]')
    axes.legend()
    plt.show()




# In[23]: Compute Mean First Passage Time from central minimum x_ini to left x_left and right x_right minima, from Long LD

if Sym_type == 'long' or Sym_type == 'both':
    
    xl = x_ini - x_delta    # will consider only starting points in the range   x_ini-x_delta,x_ini+x_delta
    xr = x_ini + x_delta

    max_time=0.0
    mfpt_l_long_ld = []
    mfpt_r_long_ld = []

    idx = 0

    print('Computing MFPT from long LD, %d data points, could take a while' %len(ld_yn))

    k=0
    for i in range(1,len(ld_yn)):
        if ( ( (i+1) // int(len(ld_yn)/10)  ) >= k ) :
            print(" %d / 10 of long LD analysed- %s" %(k, datetime.datetime.now()))
            k=k+1
    
        if i < idx:
            continue
        if ld_yn[i] > xl and ld_yn[i] < xr:
        
            j_r=np.where(ld_yn[i:]>x_right, ld_yn[i:]-x_right, np.inf).argmin()
            j_l=np.where(ld_yn[i:]<x_left, ld_yn[i:]-x_left, np.inf).argmin()
        
            if j_r==0:
                j_r=np.inf
            
            if j_l==0:
                j_l=np.inf
        
            if j_r <j_l:
                mfpt_r_long_ld.append(j_r)
                idx = j_r + i
            elif j_l<j_r:
                mfpt_l_long_ld.append(j_l)
                idx = j_l + i
        

    mfpt_l_long_ld = [x * md_dt for x in mfpt_l_long_ld]
    mfpt_r_long_ld = [x * md_dt for x in mfpt_r_long_ld]


    max_time = 0.0
    if len(mfpt_l_long_ld) > 0 and len(mfpt_r_long_ld) > 0:
        max_time=np.max([np.max(mfpt_l_long_ld),np.max(mfpt_r_long_ld)])
        print("max time for a transition: %.2f" %max_time)
        t_bins = np.linspace(0.1,max_time,t_nbins)
    elif len(mfpt_r_long_ld) == 0:
        if  len(mfpt_r_long_ld) != 0:
            max_time=np.max(mfpt_r_long_ld)
        else: print('NO DATA')
    elif len(mfpt_l_long_ld) == 0:
        if  len(mfpt_l_long_ld) != 0:
            max_time=np.max(mfpt_l_long_ld)
        else: print('NO DATA')



    fig = plt.figure()
    axes = fig.add_subplot(121)
    if len(mfpt_l_long_ld) == 0 :
        print("attention! left minimum was never reached")
    else :
        if np.max(mfpt_l_long_ld) == np.min(mfpt_l_long_ld) :
            print("attention! we had transition to left minimum %d times, but for one time: %.3f ps" %(len(mfpt_l_long_ld),np.max(mfpt_l_long_ld)))
        else :
            hist, _ = np.histogram(mfpt_l_long_ld, bins=t_bins, density=True)
            axes.bar((t_bins[1:]+t_bins[:-1])/2, hist, width=(t_bins[1:]+t_bins[:-1]))
            axes.set_title('MFPT- Center (%.1f) to Left (%.1f)' %(x_ini,x_left))
            axes.set_ylabel('probability')
            axes.set_xlabel('time [ps]')
            init_vals = [1.0, hist[0]] #[interval,100,1/100]
            best_valsL_long_ld, covarL_long_ld = curve_fit(fitfunction, (t_bins[1:]+t_bins[:-1])/2, hist, p0=init_vals) # , bounds=((interval-epsilon,-np.inf,-np.inf), (interval+epsilon,np.inf,np.inf)) )
            axes.plot((t_bins[1:]+t_bins[:-1])/2, fitfunction((t_bins[1:]+t_bins[:-1])/2, *best_valsL_long_ld))
            print("We had %d events from Center to left,  mean time: %.3f ps      tau: %.2f +/- %.2f ps" %(len(mfpt_l_long_ld), np.mean(mfpt_l_long_ld),best_valsL_long_ld[0],np.sqrt(covarL_long_ld[0,0])))

    if len(mfpt_r_long_ld) == 0 :
        print("attention! right minimum was never reached")
    else :
        if np.max(mfpt_r_long_ld) == np.min(mfpt_r_long_ld) :
            print("attention! we had transition to right minimum %d times, but for one time: %.3f ps" %(len(mfpt_r_long_ld),np.max(mfpt_r_long_ld)))
        else :
            ax2 = fig.add_subplot(122)
            hist, _ = np.histogram(mfpt_r_long_ld, bins=t_bins, density=True)
            ax2.bar((t_bins[1:]+t_bins[:-1])/2, hist, width=(t_bins[1:]+t_bins[:-1]))
            ax2.set_title('MFPT- Center (%.1f) to Right (%.1f)' %(x_ini,x_right))
            ax2.set_xlabel('time [ps]')
            init_vals = [1.0, hist[0]] #[interval,100,1/100]
            best_valsR_long_ld, covarR_long_ld = curve_fit(fitfunction, (t_bins[1:]+t_bins[:-1])/2, hist, p0=init_vals) # , bounds=((interval-epsilon,-np.inf,-np.inf), (interval+epsilon,np.inf,np.inf)) )
            axes.plot((t_bins[1:]+t_bins[:-1])/2, fitfunction((t_bins[1:]+t_bins[:-1])/2, *best_valsL_long_ld))
            print("We had %d events from Center to right, mean time: %.3f ps      tau: %.2f +/- %.2f ps" %(len(mfpt_r_long_ld), np.mean(mfpt_r_long_ld),best_valsR_long_ld[0],np.sqrt(covarR_long_ld[0,0])))

    plt.show()


    plt.figure()
    plt.plot((t_bins[1:]+t_bins[:-1])/2, hist)
    plt.show()

'''
    fileout = open('MFPT-LD_long.txt','w') 
    fileout.write("Long LD- We had %d events from Center to left,  mean time: %.3f ps      tau: %.2f +/- %.2f ps\n" %(len(mfpt_l_long_ld), np.mean(mfpt_l_long_ld),best_valsL_long_ld[0],np.sqrt(covarL_long_ld[0,0])))
    fileout.write("Long LD- We had %d events from Center to right, mean time: %.3f ps      tau: %.2f +/- %.2f ps\n" %(len(mfpt_r_long_ld), np.mean(mfpt_r_long_ld),best_valsR_long_ld[0],np.sqrt(covarR_long_ld[0,0])))
    fileout.close()
'''

# In[]: many short LD

ld_steps  = int((ld_t_end_short - ld_t_ini) / ld_dt)
ld_yn    = np.zeros(ld_steps)

max_time=0.0
mfpt_l_short_ld = []
mfpt_r_short_ld = []
idx = 0

if Sym_type == 'short' or Sym_type == 'both':

    print('\n Performing %d short LD' %ld_ntrials)
    print('Sholuld run about 10 trials a minute with %d pionts \n' %ld_t_end_short )
    
    for itrial in range(1,ld_ntrials):
        """ Get langevin dynamics """
        fpt_t_l = 0
        fpt_t_r = 0
        dW = np.random.normal(loc = 0.0, scale=ld_sqrt_dt,size=ld_steps)
        fpt_t_l, fpt_t_r = Propagate_NEW(ld_y_ini, ld_steps, ld_dt , ld_yn, dW, ld_min_allowed_y_short, ld_max_allowed_y_short, "count")
                       

        if fpt_t_l != 0:
            mfpt_l_short_ld.append(fpt_t_l)
        
        if fpt_t_r != 0:
            mfpt_r_short_ld.append(fpt_t_r)
        
        if ((itrial+1) % int(ld_ntrials/20) == 0): #change to 20
            print("trial %d / %d done!" %(itrial+1,ld_ntrials))
            print(datetime.datetime.now())
        


    max_time = 0.0
    if len(mfpt_l_short_ld) > 0 and len(mfpt_r_short_ld) > 0 :
        max_time=np.max([np.max(mfpt_l_short_ld), np.max(mfpt_r_short_ld)])
    elif len(mfpt_r_short_ld) == 0:
        if  len(mfpt_r_short_ld) != 0:
            max_time=np.max(mfpt_r_short_ld)
        else: print('NO DATA')
    elif len(mfpt_l_short_ld) == 0:
        if  len(mfpt_l_short_ld) != 0:
            max_time=np.max(mfpt_l_short_ld)
        else: print('NO DATA')
     
        print("max time for a transition: %.2f" %max_time)

    t_bins = np.linspace(0.1,max_time,t_nbins)

#%%
    fig = plt.figure()
    axes = fig.add_subplot(121)
    if len(mfpt_l_short_ld) == 0 :
        print("attention! left minimum was never reached")
    else :
        if np.max(mfpt_l_short_ld) == np.min(mfpt_l_short_ld) :
            print("attention! we had transition to left minimum %d times, but for one time: %.3f ps" %(len(mfpt_l_short_ld),np.max(mfpt_l_short_ld)))
        else :
            t_nbins_fit=25
            t_bins_fit=np.linspace(0.1,max_time,t_nbins_fit)
            hist, _ = np.histogram(mfpt_l_short_ld, bins=t_bins_fit, density=True)
            axes.set_yscale('log')
            axes.bar((t_bins_fit[1:]+t_bins_fit[:-1])/2, hist, width=t_nbins_fit)
            axes.set_title('MFPT- Center (%.1f) to Left (%.1f)' %(x_ini,x_left))
            axes.set_ylabel('probability')
            axes.set_xlabel('time [ps]')
            init_vals = [1.0, hist[0]] #[interval,100,1/100]
            best_valsL_short_ld, covarL_short_ld = curve_fit(fitfunction, (t_bins_fit[1:]+t_bins_fit[:-1])/2, hist, p0=init_vals) 
            axes.plot((t_bins_fit[1:]+t_bins_fit[:-1])/2, fitfunction((t_bins_fit[1:]+t_bins_fit[:-1])/2, *best_valsR),'y', label='MD fit')
            axes.plot((t_bins_fit[1:]+t_bins_fit[:-1])/2, fitfunction((t_bins_fit[1:]+t_bins_fit[:-1])/2, *best_valsL_short_ld), 'r', label='LD fit')
            print("We had %d events from Center to left,  mean time: %.3f ps      tau: %.2f +/- %.2f ps" %(len(mfpt_l_short_ld), np.mean(mfpt_l_short_ld),best_valsL_short_ld[0],np.sqrt(covarL_short_ld[0,0])))

    if len(mfpt_r_short_ld) == 0 :
            print("attention! right minimum was never reached")
    else :
        if np.max(mfpt_r_short_ld) == np.min(mfpt_r_short_ld) :
            print("attention! we had transition to right minimum %d times, but for one time: %.3f ps" %(len(mfpt_r_short_ld),np.max(mfpt_r_short_ld)))
        else :
            ax2 = fig.add_subplot(122)
            
            hist, _ = np.histogram(mfpt_r_short_ld, bins=t_bins, density=True)
            ax2.set_yscale('log')
            ax2.bar((t_bins[1:]+t_bins[:-1])/2, hist, width=t_nbins)
            ax2.set_title('MFPT- Center (%.1f) to Right (%.1f)' %(x_ini,x_right))
            ax2.set_xlabel('time [ps]')
            init_vals = [1.0, hist[0]] #[interval,100,1/100]
            best_valsR_short_ld, covarR_short_ld = curve_fit(fitfunction, (t_bins[1:]+t_bins[:-1])/2, hist, p0=init_vals) # , bounds=((interval-epsilon,-np.inf,-np.inf), (interval+epsilon,np.inf,np.inf)) )
            ax2.plot((t_bins[1:]+t_bins[:-1])/2, fitfunction((t_bins[1:]+t_bins[:-1])/2, *best_valsR_short_ld), 'r', label='LD fit')
            ax2.plot((t_bins[1:]+t_bins[:-1])/2, fitfunction((t_bins[1:]+t_bins[:-1])/2, *best_valsR),'y', label='MD fit')
            print("We had %d events from Center to right, mean time: %.3f ps      tau: %.2f +/- %.2f ps" %(len(mfpt_r_short_ld), np.mean(mfpt_r_short_ld),best_valsR_short_ld[0],np.sqrt(covarR_short_ld[0,0])))
    plt.legend()
    plt.show()
'''
    fig.savefig('MFPT-LD_short_local_MSD.jpg')

    fileout = open('MFPT-LD_short_local_MSD.txt','w') 
    fileout.write("Short LDs- We had %d events from Center to left,  mean time: %.3f ps      tau: %.2f +/- %.2f ps\n" %(len(mfpt_l_short_ld), np.mean(mfpt_l_short_ld),best_valsL_short_ld[0],np.sqrt(covarL_short_ld[0,0])))
    fileout.write("Short LDs- We had %d events from Center to right, mean time: %.3f ps      tau: %.2f +/- %.2f ps\n" %(len(mfpt_r_short_ld), np.mean(mfpt_r_short_ld),best_valsR_short_ld[0],np.sqrt(covarR_short_ld[0,0])))
    fileout.close()
    print(datetime.datetime.now())
'''

# In[]:
#Attempt to use Fokker Plank
import scipy
'''
def gaussian(x, mu, sig):
    return (1/sig*np.sqrt(2*np.pi))*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

FP_dt=md_dt*1.0     #md_dt=0.004 ps
dx_FP=0.001

FP_time_length= 200.0 #ps
FP_steps=int(FP_time_length/FP_dt)

x_FP=np.arange(dynamics.min(), dynamics.max(), dx_FP)
KbT = 0.0083144621 * 300   #KJoule/mol


def tridiag_FP(x, Sign ,periodic_bc=False):
    chosen_spline=cs_zones   #cs_md_diff     for local MSD w/ even binning
                                 #cs_zones       for local MSD w/ zones
                                 #cs_Humm        for Hummer method
    N=len(x)
    x_mid=np.arange(dynamics.min()-dx_FP/2, dynamics.max()+dx_FP/2, dx_FP)
    #CFL numbers
    alpha=(chosen_spline(x_mid,1)+chosen_spline(x_mid)*cs_free_energy(x_mid,2)/KbT)*(FP_dt/dx_FP)
    delta=(chosen_spline(x_mid))*(FP_dt/dx_FP**2)
    rho=(chosen_spline(x,1)*cs_free_energy(x,1)+chosen_spline(x)*cs_free_energy(x,2))*(FP_dt/KbT)
    
    print(alpha.max(), delta.max())
    
    if Sign>0:
            a=(delta[0:-1]-alpha[0:-1]/2)/2
            b=(delta[1:]+alpha[1:]/2)/2
            c=1-((delta[1:]+delta[:-1])+rho)/2
    elif Sign <0:
            a=-(delta[0:-1]-alpha[0:-1]/2)/2
            b=-(delta[1:]+alpha[1:]/2)/2
            c=1+((delta[1:]+delta[:-1])+rho)/2
    
    
    diagonals = [a[0:-1], c, b[1:]]   # 3 diagonals, CN on reaction
    A=scipy.sparse.diags(diagonals, [-1,0,1], format='csc').toarray()
    if periodic_bc== False:
        return A
    elif periodic_bc ==True:
        A[0][N-1]=a[0]
        A[N-1][0]=b[-1]
        return A

Matrix_FP=tridiag_FP(x_FP, 1)
Matrix_inv_FP=np.linalg.inv(tridiag_FP(x_FP, -1))
Matrix_CN_FP=np.dot(Matrix_inv_FP, Matrix_FP)
Matrix_CN_FP=np.where(Matrix_CN_FP>10**-8, Matrix_CN_FP, 0)
#initial conditions

x_ini   = 7.9
x_delta = 0.2 # discretization of order parameter for MFPT
x_left = 7.0
x_right = 8.75

mask_L=x_FP<x_left
mask_R=x_FP>x_right

mfpt_FP_l=[]
mfpt_FP_r=[]
y_FP=gaussian(x_FP, x_ini, x_delta/10.0)
y_FP=y_FP/y_FP.sum()
p_left=y_FP[mask_L].sum()
p_right=y_FP[mask_R].sum()
past=tmp=np.zeros(len(y_FP))
for i in range(FP_steps):
    
    y_FP=Matrix_CN_FP.dot(y_FP)
    y_FP=y_FP/y_FP.sum() #keep normalization
        
    
    if(i%10000==0):
        plt.figure()
        plt.plot(x_FP, y_FP)
        plt.show()
    #computing MFPT
    if y_FP[mask_L].sum()>p_left: 
        mfpt_FP_l.append(i*(y_FP[mask_L].sum()-p_left)*FP_dt/dx_FP) #so it's in [ps]
        p_left=y_FP[mask_L].sum()
    if y_FP[mask_R].sum()>p_right: 
        mfpt_FP_r.append(i*(y_FP[mask_R].sum()-p_right)*FP_dt/dx_FP) #so it's in [ps]
        p_right=y_FP[mask_R].sum()   


print(np.sum(mfpt_FP_l), np.sum(mfpt_FP_r))   

plt.figure()
plt.plot(x_FP, np.log(y_FP))
plt.show()
'''
