#!/usr/bin/env python
#
# import library
from os import path
from numpy import linspace, loadtxt, savetxt, column_stack, log10
import matplotlib.pyplot as plt                                      # plot
import matplotlib.gridspec as gridspec                               # for multiplot grid
from scipy.interpolate import InterpolatedUnivariateSpline as interp # for interpolation
#
original_data_path = 'cooling/data/originals'
retreated_data_path = 'cooling/data/retreated'
#
# Create Log10 Temp bins
lT_min = 1.0 
lT_max = 9.0
Nbins = 500
lTbins = linspace(lT_min, lT_max, Nbins)
#
# Rijcke
# Define metalicity log10 bins
lMetBinsRijcke = [-4.0, -2.0, -1.0, -0.5, 0.0, 0.5]
dataFilesRijcke = {-4.0: 'Rijcke_2013_Fe_H_-4_0_z_0.dat',
                   -2.0: 'Rijcke_2013_Fe_H_-2_0_z_0.dat',
                   -1.0: 'Rijcke_2013_Fe_H_-1_0_z_0.dat',
                   -0.5: 'Rijcke_2013_Fe_H_-0_5_z_0.dat',
                    0.0: 'Rijcke_2013_Fe_H_0_0_z_0.dat',
                    0.5: 'Rijcke_2013_Fe_H_0_5_z_0.dat'
                  }
lLambdaRijcke = {}
#
# Gnat & Sternberg
metBinsGnatSternberg = [1.e-3, 1.e-2, 1.e-1, 1., 2.0]             # 1 = Zsun
dataFilesGnatSternberg = 'cooling_curves_Gnat_Sternberg.dat'
lTGnatSternberg = log10(loadtxt(path.join(original_data_path, dataFilesGnatSternberg),
                                unpack=True, usecols=[0]))
lLambda = log10(loadtxt(path.join(original_data_path, dataFilesGnatSternberg),
                                     unpack=True, usecols=[1, 2, 3 ,4, 5]))
# Reformat data in a dict
lLambdaGnatSternberg = {}
for iz, Z in enumerate(metBinsGnatSternberg):
    lLambdaGnatSternberg[float(log10([Z]))] = lLambda[iz]
#
fig = plt.figure(figsize = (6., 5.))
gs = gridspec.GridSpec(1 ,1)
ax  = plt.subplot(gs[0])
gs.update(left=0.12, right=0.85, bottom=0.1, top=0.96)
#
# Load original Rijcke data and re-interpolate on lTbins
for lZ, filename in dataFilesRijcke.items():
    # Load
    lT = log10(loadtxt(path.join(original_data_path, filename), unpack=True, usecols=[0]))
    # Rijcke data have done as nh² Lambda (erg/s/cm3) with a nH = 0.01 cm-3
    # +4 in log scale to compensate
    lLambda	= loadtxt(path.join(original_data_path, filename), unpack=True, usecols=[1]) + 4.0
    # Interpolation through full temperature scale
    interpol = interp(lT[::-1], lLambda[::-1], k=1)
    lLambdaRijcke[lZ] = interpol(lTbins)
    #
    plt.plot(lT, lLambda, lw=2.0, ls='-',label='Rijcke original')
    plt.plot(lTbins, lLambdaRijcke[lZ], lw=1.0, ls=':',label='Rijcke interpolated')
    #
    # Add Gnat & Sternberg data if availlable
    if lZ in lLambdaGnatSternberg:
        plt.plot(lTGnatSternberg, lLambdaGnatSternberg[lZ], lw=2.0, ls='-',label='Gnat & Sternberg')
#
ax.set_xlabel(r'log $T$ [K]')                                  # label for x axis
ax.set_ylabel(r'log $\Lambda$ $[erg\cdot s^{-1}\cdot cm^3]$')  # label for x axis
plt.ylim([-27.0, -20.0])
plt.xlim([1.0, 9.0])
#
# Save the plot in pdf format
# path figure
plot_path = 'cooling/plots'
plt.savefig(path.join(plot_path,'Cooling_efficiency_Rijcke.pdf'))
plt.close(fig)
#
# Save reformated data
filename  = path.join(retreated_data_path, 'cooling_curves_Rijcke.dat')
# reformat data
D = [lTbins]
Ncols = 1
hdr = 'Rijcke et al cooling efficiency\n'
hdr += '1) log Temp [K]\n'
hdr += '2) log Lambda [Fe/H] = -4.0 [erg/s.cm^3]\n'
hdr += '3) log Lambda [Fe/H] = -2.0 [erg/s.cm^3]\n'
hdr += '4) log Lambda [Fe/H] = -1.0 [erg/s.cm^3]\n'
hdr += '5) log Lambda [Fe/H] = -0.5 [erg/s.cm^3]\n'
hdr += '6) log Lambda [Fe/H] = +0.0 [erg/s.cm^3]\n'
hdr += '7) log Lambda [Fe/H] = +0.5 [erg/s.cm^3]\n'
for lLambda in lLambdaRijcke.values():
    Ncols += 1
    D.append(lLambda)
    savetxt(filename, column_stack(D), header=hdr, fmt='%15.6e'*Ncols)
