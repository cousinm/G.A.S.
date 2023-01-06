#!/usr/bin/env python
# Extrapolate stellar lifetimes for metallicity bin Z = 0.0001
# Lifetimes are given by Portinari+98 for a large set of stellar mass and for 5 Metallicity Bins
#
# import library
import glob                      # to list files in a path
import os                        # to build path
import numpy as npy              # array 
import matplotlib.pyplot as plt  # plot
import math                      # for mathematical operations
from scipy.interpolate import InterpolatedUnivariateSpline as interp # for interpolation
#
# Path Portinari+98
original_data_path = 'sp/data/originals/Padova_Tracks_Portinari'
retreated_data_path = 'sp/data/retreated'
#
# Define metallicity list
MetBins = [0.0004, 0.004, 0.008, 0.02, 0.05]
# Define metallicity ref
MetBinsRef = ['0004','004','008','02','05']
#
# Load original data
filename = os.path.join(original_data_path, 'stellar_lifetimes.dat')
# Initial mass
initialMass = npy.loadtxt(filename, unpack=True, usecols=[0])
# Stellar life times
tmp = npy.loadtxt(filename, unpack=True, usecols=[1, 2, 3, 4, 5])
stellarLifeTimes = {}
for iz, Z in enumerate(MetBinsRef):
    stellarLifeTimes[Z] = tmp[iz] * 1.e-9  # Convert in Gyr
#
# Extrapolate (as a function of metalicity) for all stellar initial masses at Z = 0.0001
sltZ0001 = []
for m, item in enumerate(initialMass):
    slt = [stellarLifeTimes[Z][m] for Z in MetBinsRef]
    s = interp(npy.log10(MetBins), slt, k=1)
    val = float(s(math.log10(0.0001)))
    sltZ0001.append(val)
stellarLifeTimes['0001'] = npy.array(sltZ0001)
#
# Add the new Met bin
MetBins = npy.insert(MetBins, 0, 0.0001)
# Add the new Met bin ref
MetBinsRef = npy.insert(MetBinsRef, 0, '0001')
#
# Create interpolation spline (order 2) for each metalicity
interpol = {}
for Z, slt in stellarLifeTimes.items():
    s = interp(initialMass, slt, k=2)
    interpol[Z] = s
# 
# Oversample data
dm = 0.1
Mmin = 0.8
Mmax = 100.
NMass = int((Mmax - Mmin)/dm) +1
# Create the new initial mass list
initialMass = npy.linspace(Mmin, Mmax, NMass)
# Oversample
for iz, Z in enumerate(MetBinsRef):
    stellarLifeTimes[Z] = interpol[Z](initialMass)
#
# Plot
plot_path = 'sp/plots'
plt.figure(figsize=(6, 5))
color = ['SkyBlue', 'LightSeaGreen', 'YellowGreen', 'Gold', 'DarkOrange', 'FireBrick']
#
for iz, Z in enumerate(MetBinsRef):
    plt.plot(initialMass, stellarLifeTimes[Z], c=color[iz], label=f'Z = 0.{Z}') 
plt.xlim([0.8, 100])                     # xrange
plt.xlabel('$m_{\star}$ [$M_{\odot}$]')  # label for x axis
plt.ylabel(r'$\tau_{\star}$ [Gyr]')      # label for y axis
plt.xscale('log')                        # use log scale for x axis
plt.yscale('log')                        # use log scale for y axis
plt.grid()
plt.legend(loc='best',prop={'size':10})        # display legend    
# Save the plot in pdf format
filename = os.path.join(plot_path, 'stellar_lifetimes.pdf')
plt.savefig(filename) 
#
# Save data
for z, item in enumerate(MetBinsRef):
    # Reformat data
    D = []
    D.append(initialMass)
    for Z in MetBinsRef:
        D.append(stellarLifeTimes[Z])
    filename = os.path.join(retreated_data_path, 'extended_Stellar_Lifetimes.dat')
    hdr =  'Stellar lifetimes (in Gyr) \n'
    hdr += '1) Mass, 2) Z = 0.0001, 3) Z = 0.0004, 4) Z = 0.004, 5) Z = 0.008, 6) Z = 0.02, 7) Z = 0.05'
    npy.savetxt(filename, npy.column_stack(D), header=hdr, fmt='%15.6e'*7)    
