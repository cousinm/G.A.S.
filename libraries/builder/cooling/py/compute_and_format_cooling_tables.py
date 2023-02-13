#!/usr/bin/env python
#
# import library
from os import path
from numpy import linspace, loadtxt, savetxt, column_stack, log10, array, meshgrid, amax, amin
from math import log10, floor, ceil
import matplotlib.pyplot as plt         # plot
import matplotlib.cm as cm              # colormap
import matplotlib.gridspec as gridspec  # for multiplot grid
import matplotlib.colors as colors      # colors
import matplotlib.colorbar as cb        # colorbar
from scipy.interpolate import interp2d  # for interpolation
#
retreated_data_path = 'cooling/data/retreated'
#
# Load data  Rijcke 2013
author = 'Rijcke 2013'
filename = 'cooling_curves_Rijcke.dat'
lTempBins_orig = loadtxt(path.join(retreated_data_path, filename),
                         unpack=True,usecols=[0])
lLambda_orig = loadtxt(path.join(retreated_data_path, filename),
                       unpack=True, usecols=[1, 2, 3, 4, 5, 6])
#
# Metallicity table
# Rijcke 2013 metalicities are given in solar scale
# Convert in metal mass fraction according to Zsun = 0.02
Zsun = 0.02
lZsun = log10(Zsun)
lMetBins_orig = array([-4.0, -2.0, -1.0, -0.5, 0.0, 0.5]) + lZsun
#
# Interpolation as a 2D grid
interp = interp2d(lTempBins_orig, lMetBins_orig, lLambda_orig)
#
# Create Bin table
# Temperature
lTemp_min = 2.0  # 100K 
lTemp_max = 8.0  # 10^8K
nlTemp = 201
lTempBins = linspace(lTemp_min, lTemp_max, nlTemp)
dlTemp = lTempBins[1] - lTempBins[0]
print(f'dlTemp={dlTemp}')
# Metallicity
nlMet = 51
lMetBins = linspace(-5., -1., nlMet) # Metal mass fraction 
dlMet = lMetBins[1] - lMetBins[0]
lMetMin = min(lMetBins)
lMetMax = max(lMetBins)
print(f'dlMet={dlMet}')
#
# Generate Cooling efficiency grid
lLambda = {}
lLambda_min = 1.e10
lLambda_max = -1.e10
for iz, _ in enumerate(lMetBins):
    lLambda[iz] = interp(lTempBins, lMetBins[iz])
    min = amin(lLambda[iz])
    if min < lLambda_min:
        lLambda_min = min
    max = amax(lLambda[iz])
    if max > lLambda_max:
        lLambda_max = max
lLambda_min = floor(lLambda_min)
lLambda_max = ceil(lLambda_max)
#
# Compute instability parameter: 2.0 - dlLambda / dlT
dlLambda_dlT = {}
dlLambda_dlT_min = 1.e10
dlLambda_dlT_max = -1.e10
for iz, _ in enumerate(lMetBins):
    dlLdlT = []
    for it, lT in enumerate(lTempBins):
        if ((it > 0) and (it < nlTemp - 1)):
            # Center derivate
            dlnLambda = lLambda[iz][it + 1] - lLambda[iz][it - 1]
            val = dlnLambda/(2.*dlTemp)
        if (it == 0):
            # Upper derivative
            dlnLambda = lLambda[iz][it + 1] - lLambda[iz][it]
            val = dlnLambda/dlTemp
        if (it == nlTemp - 1):
            dlnLambda = lLambda[iz][it] - lLambda[iz][it - 1]
            val = dlnLambda/dlTemp
        dlLdlT.append(2.0 - val)
    dlLambda_dlT[iz] = array(dlLdlT)
    #
    # Evaluate min, max values
    min = amin(dlLambda_dlT[iz])
    if min < dlLambda_dlT_min:
        dlLambda_dlT_min = min
    max = amax(dlLambda_dlT[iz])
    if max > dlLambda_dlT_max:
        dlLambda_dlT_max = max
dlLambda_dlT_min = floor(dlLambda_dlT_min)
dlLambda_dlT_max = ceil(dlLambda_dlT_max)
#
# Save Data
# Cooling_efficiency
filename = path.join(retreated_data_path, 'cooling_efficiency.in')
# Header
lTempBinsStr = ' '.join([f'{lt:.2f}' for lt in lTempBins])
lMetBinsStr = ' '.join([f'{lz:.2f}' for lz in lMetBins])
hdr = ('New G.A.S. format\n'
       'In this file, each line of the header is starting by #\n'
       '   data are saved after the keyword ----\n'
       'Cooling efficiency [erg/s.cm^3] associated to different metallicities\n' 
       '  computed from Rijcke 2013\n'
       'This data file is build as follow:\n'
       'nTempBins nMetBins\n'
       'dlogTemp dlogMet\n'
       'List of Temperatures (in log) [K]             : TempBins\n'
       'List of metalicity bins (in log) [mass ratio] : MetBins\n'
       'A blank line starting with #\n'
       'And for each following data line (i: Temp index): logLambda(i,1:nMetBins) in erg/s.cm^3\n'
       '----\n'
       f'{nlTemp} {nlMet}\n'
       f'{dlTemp:.2f} {dlMet:.2f}\n'
       f'{lTempBinsStr}\n'
       f'{lMetBinsStr}\n')
# reformat data
D = []
Ncols = 0
for iz, _ in enumerate(lMetBins):
    Ncols += 1
    D.append(lLambda[iz])
savetxt(filename, column_stack(D), header=hdr, fmt='%15.6e'*Ncols)
#
# Plot
#
fig = plt.figure(figsize = (6.,5.))
gs = gridspec.GridSpec(1,1)
ax  = plt.subplot(gs[0])
gs.update(left=0.12, right=0.85, bottom=0.1, top=0.96)
#
# Define the colormaps
cmap = cm.hot
#
bound_z = linspace(lLambda_min, lLambda_max, 9)
norm_z  = colors.Normalize(vmin=amin(bound_z),vmax=amax(bound_z))
ax2     = fig.add_axes([0.86, 0.1, 0.01, 0.86]) # [left,bottom,width,height]
cb_z    = cb.ColorbarBase(ax2, cmap=cmap, spacing='proportional',
                          ticks=bound_z, boundaries=bound_z,
                          format='%3.1f', label=r'log $\Lambda$ [$erg/s\cdot cm^3$]')
# create a ScalarMappable and initialize a data structure
cb_sm   = cm.ScalarMappable(cmap=cmap, norm=norm_z)
cb_sm.set_array([])
#
ax.imshow(array(D), vmin=lLambda_min, vmax=lLambda_max, \
          extent=[lTemp_min, lTemp_max, lMetMin-lZsun, lMetMax-lZsun],
          origin='lower', interpolation='nearest', cmap=cmap,
          aspect='auto', zorder = 0)
#
ax.set_xlabel(r'log $T$ [K]', fontsize=14)    # label for x axis
ax.set_ylabel(r'log $Z$ [$Z_{\odot}$]', fontsize=14)    # label for x axis
#
# Save the plot in pdf format
path_plot = 'cooling/plots'
plt.savefig(path.join(path_plot, 'cooling_efficiency_grid.pdf'))
plt.close(fig)
#
# Thermal instabillities
filename  = path.join(retreated_data_path, 'thermal_instability.in')
# Header
hdr = ('New G.A.S. format\n'
       'In this file, each line of the header is starting by #\n'
       '   data are saved after the keyword ----\n'
       'Thermal_instability parameter (2. - dlnLambda/dlnT) associated to different metallicities\n'
       '  computed from Rijcke 2013\n'
       'This data file is build as follow:\n'
       'nTempBins nMetBins\n'
       'dlogTemp dlogMet\n'
       'List of Temperatures (in log) [K]             : TempBins\n'
       'List of metalicity bins (in log) [mass ratio] : MetBins\n'
       'A blank line starting with #\n'
       'And for each following data line (i: Temp index): 2. - dlnLambda/dlnT(i,1:nMetBins)\n'
       '----\n'
       f'{nlTemp} {nlMet}\n'
       f'{dlTemp:.2f} {dlMet:.2f}\n'
       f'{lTempBinsStr}\n'
       f'{lMetBinsStr}\n')
# Reformat data
D = []
Ncols = 0
for iz, _ in enumerate(lMetBins):
    Ncols += 1
    D.append(dlLambda_dlT[iz])
savetxt(filename, column_stack(D), header=hdr, fmt='%15.6e'*Ncols)
#
fig = plt.figure(figsize = (6.,5.))
gs = gridspec.GridSpec(1,1)
ax  = plt.subplot(gs[0])
gs.update(left=0.12, right=0.85, bottom=0.1, top=0.96)
# define the colormaps
cmap = cm.hot
#
bound_z = linspace(dlLambda_dlT_min, dlLambda_dlT_max, 10)
norm_z  = colors.Normalize(vmin=amin(bound_z),vmax=amax(bound_z))
ax2     = fig.add_axes([0.86, 0.1, 0.01, 0.86]) # [left,bottom,width,height]
cb_z    = cb.ColorbarBase(ax2, cmap=cmap, spacing='proportional',
                          ticks=bound_z, boundaries=bound_z,
                          format='%3.1f', label=r'$2 -dln\Lambda/dlnT$')
# create a ScalarMappable and initialize a data structure
cb_sm   = cm.ScalarMappable(cmap=cmap, norm=norm_z)
cb_sm.set_array([])
#
# Contours
# Grid for contours
lTgrid, lZgrid = meshgrid(lTempBins, lMetBins-lZsun)
ax.contour(lTgrid, lZgrid, array(D), array([-2.5, 0.]), ls='-',
           origin='lower', colors='royalblue', linewidths=2.0, zorder = 1)
#
ax.imshow(array(D), vmin=dlLambda_dlT_min, vmax=dlLambda_dlT_max, \
          extent=[lTemp_min, lTemp_max, lMetMin-lZsun, lMetMax-lZsun], origin='lower',
          interpolation='nearest', cmap=cmap, aspect='auto', zorder = 0)
#
ax.set_xlabel(r'log $T$ [K]', fontsize=14)    # label for x axis
ax.set_ylabel(r'log $Z$ [$Z_{\odot}$]', fontsize=14)    # label for x axis
#
# Save the plot in pdf format
# path figure
plt.savefig(path.join(path_plot, 'thermal_instability_grid.pdf'))
plt.close(fig)
