#!/usr/bin/env python
#
# Interpolate the final mass at Z = 0.05 for stellar mass m = 30.
# Extrapolate the final mass for Z = 0.0001 for all initial stellar mass
#
# import library
import os                        # to build path
import math
from glob import glob
import numpy as npy              # array
import matplotlib.pyplot as plt  # plot
from scipy.interpolate import InterpolatedUnivariateSpline as interp
#
# Path to original data from Portinari+98
original_data_path = 'sp/data/originals/Padova_Tracks_Portinari'
retreated_data_path = 'sp/data/retreated'
#
# Read the mass table
filename = os.path.join(original_data_path, 'mass.dat')
Mass = npy.loadtxt(filename, unpack=True)
#
# Read the original initial/final mass data file
filename = os.path.join(original_data_path, 'initial_final_mass_MS.dat')
finalMass = npy.loadtxt(filename, unpack=True, usecols=[1, 2, 3, 4, 5])
# finalMass are structured as follow: finalMass[Z][m]
# Z from 0.0004 to Z = 0.05
#
# Create MetBins
# The metalicity bin 0.0001 is not given in the study of Portinar+98,
# We have to interpolate data according to other metalicity bins
MetBins = [0.0004, 0.004, 0.008, 0.02, 0.05]
MetBinsRef = ['0004', '004', '008', '02', '05']
#
# Interpolate the final mass associated to m = 30. at Z = 0.05:
# finalMass[4][3]
initialMassList = [Mass[i] for i in [0, 1, 2, 4, 5, 6, 7]]
finalMassList = [finalMass[4][i] for i in [0, 1, 2, 4, 5, 6, 7]]
finalMassInterpFromMass = interp(initialMassList, finalMassList, k=1)
ZList = [MetBins[z] for z in [0, 1, 2, 3]]
finalMassList = [finalMass[z][3] for z in [0, 1, 2, 3]]
finalMassInterpFromMet = interp(ZList, finalMassList, k=1)
fM30M = finalMassInterpFromMass(30.0)
fM30Z = finalMassInterpFromMet(0.05)
fM30 = (fM30M + fM30Z) / 2.
#
# Replace the "0" value given in input
if fM30 < 30.0:
    finalMass[4][3] = fM30
else:
    msg = ('Final mass > initial mass (Z = %.4f, m = %.3f, mf = %.3f)'
           % (0.05, 30.0, fM30))
    raise Exception(msg)

# Interpolate final mass for all masses at Z = 0.0001
finalMassZ0001 = []
for im, m in enumerate(Mass):
    initialMassList = [finalMass[iz][im] for iz, _ in enumerate(MetBins)]
    finalMassInterp = interp(MetBins, initialMassList, k=1)
    fMZ0001 = float(finalMassInterp(0.0001))
    if fMZ0001 < m:
        finalMassZ0001.append(fMZ0001)
    else:
        fm = 0.999*m
        msg = ('Final mass > initial mass (Z = %.4f, m = %.3f, mf = %.3f),\
                reset to m = Z = %.4f'
               % (0.0001, m, fMZ0001, fm))
        finalMassZ0001.append(fm)
        print(msg)

# Add the new metallicity bin
MetBins = npy.insert(MetBins, 0, 0.0001)
MetBinsRef = npy.insert(MetBinsRef, 0, '0001')

# Insert this new MetBins
newfinalMass = [finalMassZ0001]
for z, item in enumerate(finalMass):
    newfinalMass.append(finalMass[z])
#
# plot
plot_path = 'sp/plots'
#
plt.figure(figsize=(5.5, 5))
# color reference
color = ['SkyBlue', 'LightSeaGreen', 'YellowGreen',
         'Gold', 'DarkOrange', 'FireBrick']
#
for z, item in enumerate(MetBins):
    legend = 'Z = 0.'+MetBinsRef[z]
    plt.plot(Mass, newfinalMass[z], c=color[z], label=legend)
    if z == 0:  # Z = 0.0001, all dat are extrapolated
        plt.scatter(Mass, newfinalMass[z], s=40.,
                    marker='s', facecolor='none', edgecolor=color[z])
    elif z == 5:  # Z = 0.05, final mass for M = 30 is interpolated
        x = [Mass[m] for m in [0, 1, 2, 4, 5, 6, 7]]
        y = [newfinalMass[z][m] for m in [0, 1, 2, 4, 5, 6, 7]]
        plt.scatter(x, y, s=40., c=color[z], marker='o', edgecolor=color[z])
        x = [Mass[m] for m in [3]]
        y = [newfinalMass[z][m] for m in [3]]
        plt.scatter(x, y, s=40., marker='s',
                    facecolor='none', edgecolor=color[z])
    else:
        plt.scatter(Mass, newfinalMass[z], s=40., c=color[z],
                    marker='o', edgecolor=color[z])
#
plt.ylim([0, 80])                               # yrange
plt.xlabel(r'Initial mass ($M_{\odot}$)')       # label for x axis
plt.ylabel(r'Final mass ($M_{\odot}$)')         # label for y axis
plt.grid(True, which="both", ls="-")            # display a grid
plt.legend(loc=2, prop={'size': 10})            # display legend
filename = os.path.join(plot_path, 'initial_final_mass_MS.pdf')
plt.savefig(filename)             	            # save the figure
#
# Save data
# Create filename
filename = os.path.join(retreated_data_path, 'initial_final_mass_MS.dat')
D = []
D.append(Mass)
for z, item in enumerate(newfinalMass):
    D.append(newfinalMass[z])
hdr = 'Initial .vs. final mass for a set of metallicity bins\n'
hdr += 'The final mass at Z = 0.05 for m = 30. is interpolated\n'
hdr += 'Final mass at Z = 0.0001 for all mass are extrapolated\n'
hdr += '1) Minit, 2) Z=0.0001 (extrap), 3) Z=0.0004,\
        4) Z=0.004, 5) Z=0.008, 6) Z=0.02, 7) Z=0.05'
npy.savetxt(filename, npy.column_stack(D), header=hdr, fmt='%15.6e'*7)
