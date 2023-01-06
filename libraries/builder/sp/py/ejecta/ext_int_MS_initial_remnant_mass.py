#!/usr/bin/env python
#
# Interpolate the remnant mass at Z = 0.004 for stellar mass m = 60.
# Interpolate the remnant mass at Z = 0.008 for stellar mass m = 100.
# Interpolate the remnant mass at Z = 0.05 for stellar mass m = 30.
# Extrapolate the final mass for Z = 0.0001 for all initial stellar mass
#
# Import library
import os                        # to build path
import numpy as npy              # array
import matplotlib.pyplot as plt  # plot
from scipy.interpolate import InterpolatedUnivariateSpline as interp
#
# Path to original data from Portinari+98
data_path = 'sp/data/originals/Padova_Tracks_Portinari'
retreated_data_path = 'sp/data/retreated'
#
# read the mass table
Mass = npy.loadtxt(os.path.join(data_path, 'mass.dat'), unpack=True)
#
# Read the original initial_final_mass data file
filename = os.path.join(data_path, 'initial_remnant_mass_MS.dat')
remnantMass = npy.loadtxt(filename, unpack=True, usecols=[1, 2, 3, 4, 5])
# remnantMass are structured as follow: remnantMass[Z][m]
# Z from 0.0004 to Z = 0.05
#
# Create MetBins
# The metalicity bin 0.0001 is not given in the study of Portinar+98,
# We have to interpolate these data.
MetBins = [0.0004, 0.004, 0.008, 0.02, 0.05]
MetBinsRef = ['0004', '004', '008', '02', '05']
#
# Interpolate the remnant mass associated to m = 60. at Z = 0.004:
# remnantMass[1][5]
initialMassList = [Mass[i] for i in [0, 1, 2, 3, 4, 6, 7]]
remnantMassList = [remnantMass[1][i] for i in [0, 1, 2, 3, 4, 6, 7]]
remnantMassInterpFromMass = interp(initialMassList, remnantMassList, k=1)
ZList = [MetBins[z] for z in [0, 2, 3, 4]]
remnantMassList = [remnantMass[z][5] for z in [0, 2, 3, 4]]
remnantMassInterpFromMet = interp(ZList, remnantMassList, k=1)
rM60M = remnantMassInterpFromMass(60.0)
rM60Z = remnantMassInterpFromMet(0.004)
rM60 = (0.3*rM60M + 0.7*rM60Z)
#
# Replace the "0" value given in input
if rM60 < 60.0:
    remnantMass[1][5] = float(rM60)
else:
    msg = ('Final mass > initial mass (Z = %.4f, m = %.3f, mf = %.3f)'
           % (0.004, 60.0, rM60))
    raise Exception(msg)
#
# Interpolate the remnant mass associated to m = 100. at Z = 0.008:
# remnantMass[2][6]
initialMassList = [Mass[i] for i in [0, 1, 2, 3, 4, 5, 7]]
remnantMassList = [remnantMass[2][i] for i in [0, 1, 2, 3, 4, 5, 7]]
remnantMassInterpFromMass = interp(initialMassList, remnantMassList, k=1)
ZList = [MetBins[z] for z in [0, 1, 3, 4]]
remnantMassList = [remnantMass[z][6] for z in [0, 1, 3, 4]]
remnantMassInterpFromMet = interp(ZList, remnantMassList, k=1)
rM100M = remnantMassInterpFromMass(100.0)
rM100Z = remnantMassInterpFromMet(0.008)
rM100 = rM100M  # (rM100M + rM100Z) / 2.
#
# Replace the "0" value given in input
if rM100 < 100.0:
    remnantMass[2][6] = float(rM100)
else:
    msg = ('Final mass > initial mass (Z = %.4f, m = %.3f, mf = %.3f)'
           % (0.008, 100.0, rM100))
    raise Exception(msg)
#
# Interpolate the remnant mass associated to m = 30. at Z = 0.05:
# remnantMass[4][3]
initialMassList = [Mass[i] for i in [0, 1, 2, 4, 5, 6, 7]]
remnantMassList = [remnantMass[4][i] for i in [0, 1, 2, 4, 5, 6, 7]]
remnantMassInterpFromMass = interp(initialMassList, remnantMassList, k=1)
ZList = [MetBins[z] for z in [0, 1, 2, 3]]
remnantMassList = [remnantMass[z][3] for z in [0, 2, 1, 3]]
remnantMassInterpFromMet = interp(ZList, remnantMassList, k=1)
rM30M = remnantMassInterpFromMass(30.0)
rM30Z = remnantMassInterpFromMet(0.05)
rM30 = (rM30M + rM30Z) / 2.
#
# Replace the "0" value given in input
if rM30 < 30.0:
    remnantMass[4][3] = float(rM30)
else:
    msg = ('Final mass > initial mass (Z = %.4f, m = %.3f, mf = %.3f)'
           % (0.05, 30.0, rM30))
    raise Exception(msg)
#
# Interpolate remnant mass for all masses at Z = 0.0001
remnantMassZ0001 = []
for im, m in enumerate(Mass):
    initialMassList = [remnantMass[iz][im] for iz, _ in enumerate(MetBins)]
    remnantMassInterp = interp(MetBins, initialMassList, k=1)
    fMZ0001 = float(remnantMassInterp(0.0001))
    if fMZ0001 < m:
        remnantMassZ0001.append(fMZ0001)
    else:
        fm = 0.999*m
        msg = ('Remnant mass > initial mass (Z = %.4f, m = %.3f, mf = %.3f),\
                reset to m = Z = %.4f'
               % (0.0001, m, fMZ0001, fm))
        remnantMassZ0001.append(fm)
        print(msg)

# Add the new metallicity bin
MetBins = npy.insert(MetBins, 0, 0.0001)
MetBinsRef = npy.insert(MetBinsRef, 0, '0001')

# Insert this new MetBins
newRemnantMass = [remnantMassZ0001]
for z, item in enumerate(remnantMass):
    newRemnantMass.append(remnantMass[z])
#
# plot
plot_path = 'sp/plots'
#
# Check with final mass
# final mass correspond to the stellar mass at the end of the wind phase
# We must have final mass > remnant mass
# Load final mass
filename = os.path.join(retreated_data_path, 'initial_final_mass_MS.dat')
finalMass = npy.loadtxt(filename, unpack=True, usecols=[1, 2, 3, 4, 5, 6])
# Run onto metalliity bins
for iz, z in enumerate(MetBins):
    # run onto mass
    for im, m in enumerate(finalMass):
        fM = finalMass[iz][im]
        rM = newRemnantMass[iz][im]
        if  fM < rM:
            msg = ('Final mass < remnant mass (Z = %.4f, fm = %.3f, rm = %.3f)'
                % (z, fM, rM))
            raise Exception(msg)
#
# Save data
filename = os.path.join(retreated_data_path, 'initial_remnant_mass_MS.dat')
D = []
D.append(Mass)
for z, item in enumerate(newRemnantMass):
    D.append(newRemnantMass[z])
    
hdr = 'Initial .vs. remnant mass for a set of metallicity bins\n'
hdr += 'The remnant mass at Z = 0.004 for am = 60. is interpolated\n'
hdr += 'The remnant mass at Z = 0.008 for m = 100. is interpolated\n'
hdr += 'The remnant mass at Z = 0.05 for m = 30. is interpolated\n'
hdr += 'Remnant mass at Z = 0.0001 for all mass are extrapolated\n'
hdr += '1) Minit, 2) Z=0.0001 (extrap), 3) Z=0.0004,\
        4) Z=0.004, 5) Z=0.008, 6) Z=0.02, 7) Z=0.05'
npy.savetxt(filename,npy.column_stack(D),header=hdr,fmt='%15.6e'*7)
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
    plt.plot(Mass, newRemnantMass[z], c=color[z], label=legend)
    if z == 0:  # Z = 0.0001, all data are extrapolated
        plt.scatter(Mass, newRemnantMass[z], s=40.,
                    marker='s', facecolor='none', edgecolor=color[z])
    elif z == 2: # Z = 0.004, remnant mass at m = 60 is interpolated
        x = [Mass[m] for m in [0, 1, 2, 3, 4, 6, 7]]
        y = [newRemnantMass[z][m] for m in [0, 1, 2, 3, 4, 6, 7]]
        plt.scatter(x, y, s=40., c=color[z], marker='o', edgecolor=color[z])
        plt.scatter(Mass[5], newRemnantMass[z][5], s=40., marker='o',
                    facecolor='none', edgecolor=color[z])
    elif z == 3: # Z = 0.008, remnant mass at m = 100 is interpolated
        x = [Mass[m] for m in [0, 1, 2, 3, 4, 5, 7]]
        y = [newRemnantMass[z][m] for m in [0, 1, 2, 3, 4, 5, 7]]
        plt.scatter(x, y, s=40., c=color[z], marker='o', edgecolor=color[z])
        plt.scatter(Mass[6], newRemnantMass[z][6], s=40., marker='o',
                    facecolor='none', edgecolor=color[z])
    elif z == 5: # Z = 0.05, remnant mass at m = 30 is interpolated
        x = [Mass[m] for m in [0, 1, 2, 4, 5, 6, 7]]
        y = [newRemnantMass[z][m] for m in [0, 1, 2, 4, 5, 6, 7]]
        plt.scatter(x, y, s=40., c=color[z], marker='o', edgecolor=color[z])
        plt.scatter(Mass[3], newRemnantMass[z][3], s=40., marker='o',
                    facecolor='none', edgecolor=color[z])
    else:
        plt.scatter(Mass, newRemnantMass[z], s=40., c=color[z], 
                    marker='o', edgecolor=color[z])

plt.ylim([0, 80])                               # yrange
plt.xlabel(r'Initial mass ($M_{\odot}$)')       # label for x axis
plt.ylabel(r'Remnant mass ($M_{\odot}$)')       # label for y axis
plt.grid(True, which="both", ls="-")            # display a grid
plt.legend(loc=2, prop={'size': 10})            # display legend
filename = os.path.join(plot_path, 'initial_remnant_mass_MS.pdf')
plt.savefig(filename)             	            # save the figure
