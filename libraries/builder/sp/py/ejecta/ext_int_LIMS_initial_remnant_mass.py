#!/usr/bin/env python
#
# Interpolate final mass for all initial stellar mass at Z = 0.0004
# Extrapolate final mass for all initial stellar mass at Z = 0.05
#
# Import library
import os                        									 # to build path
import numpy as npy              									 # array
import math                      									 # for mathematical operations
import matplotlib.pyplot as plt  									 # plot
from scipy.interpolate import InterpolatedUnivariateSpline as interp # for interpolation
#
# path Karakas+10
original_data_path = 'sp/data/originals/Karakas'
retreated_data_path = 'sp/data/retreated'
#
# Create MetBins
MetBins = [0.0001, 0.004, 0.008, 0.02]
MetBinsRef = ['0001', '004', '008', '02']
#
# Load the original initial_final_mass table
filename = os.path.join(original_data_path, 'initial_final_mass_LIMS.dat')
#
initialMass = npy.loadtxt(filename, unpack=True, usecols=[0])
tmp = npy.loadtxt(filename, unpack=True, usecols=[1, 2, 3, 4])
finalMass = {}
for iz, Z in enumerate(MetBinsRef):
    finalMass[Z] = tmp[iz]
del tmp            # delete tmp array
#
# Interpolate (as a function of the metalicity)
# the final mass for all initial stellar mass at Z = 0.0004
fM0004 = []
for im, M in enumerate(initialMass):
    fM = [finalMass[Z][im] for Z in MetBinsRef]
    intrpl = interp(npy.log10(MetBins), fM, k=1)
    fM = float(intrpl(math.log10(0.0004)))
    if fM < M:
        fM0004.append(fM)
    else:
        raise Exception('Final mass > initial mass')
finalMass['0004'] = npy.array(fM0004)
#
# Add the new metallicty bin
MetBins = npy.insert(MetBins, 1, 0.0004)
MetBinsRef = npy.insert(MetBinsRef, 1, '0004')
#
# Extrapolate (as a function of the metallicity)
# the final mass for all initial stellar mass at Z = 0.05
fM05 = []
for im, M in enumerate(initialMass):
    fM = [finalMass[Z][im] for Z in MetBinsRef]
    intrpl = interp(npy.log10(MetBins), fM, k=1)
    fM = float(intrpl(math.log10(0.05)))
    if fM < M:
        fM05.append(fM)
    else:
        raise Exception('Final mass > initial mass')
finalMass['05'] = npy.array(fM05)
#
# Add the new metallicty bin
MetBins = npy.insert(MetBins,5,0.05)
MetBinsRef = npy.insert(MetBinsRef,5,'05')
#
# plot
plot_path = 'sp/plots'
plt.figure(figsize = (5.5,5))
# color reference
color = ['SkyBlue', 'LightSeaGreen', 'YellowGreen',
         'Gold','DarkOrange','FireBrick']
#
for iz, Z in enumerate(MetBinsRef):
    legend = f'Z = 0.{Z}'
    plt.plot(initialMass, finalMass[Z], c=color[iz], label=legend)
    if Z in ['0004', '05']: # Z = 0.0004, 0.05,
                            # all remnant masses are extrapolated
        plt.scatter(initialMass, finalMass[Z], s=40.,
                    marker='o', facecolor='none', edgecolor=color[iz])
    else:
        plt.scatter(initialMass, finalMass[Z], s=40., c=color[iz],
                    marker='o', edgecolor=color[iz])
#
plt.ylim([0.5,1.1])                               # yrange
plt.xlabel('Initial mass ($M_{\odot}$)')       # label for x axis
plt.ylabel('Remnant mass ($M_{\odot}$)')       # label for y axis	
plt.grid(True,which="both",ls="-")             # display a grid with all sub log tics
plt.legend(loc='best',prop={'size':10})        # display legend                            
filename = os.path.join(plot_path,'initial_remnant_mass_LIMS.pdf')
plt.savefig(filename)             	           # save the figure
#
# Save data
# Create filename
filename = os.path.join(retreated_data_path, 'initial_remnant_mass_LIMS.dat')
D = []
D.append(initialMass)
for Z in MetBinsRef:
    D.append(finalMass[Z])
hdr = 'Initial .vs. remnant mass for a set of metallicity bins\n'
hdr += 'All remnant stellar masses at Z = 0.0004 are interpolated\n'
hdr += 'All remnant stellar masses at Z = 0.05 are extrapolated\n'
hdr += '1) Minit, 2) Z=0.0001 3) Z=0.0004 (interp), 4) Z=0.004, 5) Z=0.008, 6) Z=0.02, 7) Z=0.05 (extrap)'
npy.savetxt(filename, npy.column_stack(D), header=hdr, fmt='%15.6e'*7)
