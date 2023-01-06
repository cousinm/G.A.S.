#!/usr/bin/env python
# Interpolate abundances for metalicity bin Z = 0.0004
# Extrapolate abundances for metalicity bin Z = 0.05
#
# import library
import os                        # to build path
import numpy as npy              # array 
import matplotlib.pyplot as plt  # plot
import math                      # for mathematical operations
from scipy.interpolate import InterpolatedUnivariateSpline as interp # for interpolation
#
# path Karakas
original_data_path = 'sp/data/originals/Karakas'
#
# filename
filename = os.path.join(original_data_path, 'initial_abundance.dat')
# Read the metallicity table
Metalicity = npy.loadtxt(filename, unpack=True, usecols=[0])
# Convert in log scale for interpolation
logMetalicity = npy.log10(Metalicity)
#
# Create element list
Elements = ['H1', 'He4', 'C12', 'N14', 'O16', 'Fe56']
Elements_Tex = ['$^{1}H$', '$^{4}He$', '$^{12}C$', '$^{14}N$', '$^{16}O$', '$^{56}Fe$']
#
# Read abundance table
X0 = npy.loadtxt(filename, unpack=True, usecols=[1, 2, 3, 4, 5, 6])
# Convert in log scale for interpolation
logX0 = npy.log10(X0)
#
# Test mass fraction and metallicity values
for iz, Z in enumerate(Metalicity):
    # Compute metal mass fraction Zmax = 1 - X0_H1 - X0_He4 and compared with Z
    Zmax = 1. - X0[0][iz] - X0[1][iz]
    # Sum mass fraction associated to metal elements take into account (i.e C12, N14, O16, Fe56)
    Zmin = sum([X0[e][iz] for e in [2, 3, 4, 5]])
    # Test
    if (Zmax < Z or Z < Zmin):
        X0_H1_He4 = X0[0][iz] + X0[1][iz]
        dX = X0_H1_He4 - (1. - Z)
        X0[0][iz] = X0[0][iz]*(1. -dX/X0_H1_He4)
        X0[1][iz] = X0[1][iz]*(1. -dX/X0_H1_He4)
        # recompute Zmax
        Zmax = 1. - X0[0][iz] - X0[1][iz]
#
# Interpolate/Extrapolate abundances
# Use log scale
NewX0 = []
for ielt, elt in enumerate(Elements):
    tmp = logX0[ielt]
    interpol = interp(logMetalicity, tmp, k=2)
    val = interpol(math.log10(0.0004))
    tmp = npy.insert(tmp, 1, val)
    val = interpol(math.log10(0.05))
    tmp = npy.append(tmp, val)
    NewX0.append(10.**tmp)
#
# Add metallicity Bins
Metalicity = npy.insert(Metalicity, 1, 0.0004)
Metalicity = npy.insert(Metalicity, 5, 0.05)
#
# Test mass fraction and metallicity values
for iz, Z in enumerate(Metalicity):
    # Compute metal mass fraction Zmax = 1 - X0_H1 - X0_He4 and compared with Z
    Zmax = 1. - NewX0[0][iz] - NewX0[1][iz]
    # Sum mass fraction associated to metal elements take into account (i.e C12, N14, O16, Fe56)
    Zmin = sum([NewX0[e][iz] for e in [2, 3, 4, 5]])
    # Test
    if (Zmax < Z or Z < Zmin):
        X0_H1_He4 = NewX0[0][iz] + NewX0[1][iz]
        dX = X0_H1_He4 - (1. - Z)
        NewX0[0][iz] = NewX0[0][iz]*(1. -dX/X0_H1_He4)
        NewX0[1][iz] = NewX0[1][iz]*(1. -dX/X0_H1_He4)
        # recompute Zmax
        Zmax = 1. - NewX0[0][iz] - NewX0[1][iz]
#
# Define primordial abundance
Metalicity = npy.insert(Metalicity, 0, 0.)
X0_primordial = [0.78, 0.22, 0., 0. ,0., 0. ]
for i, X0 in enumerate(NewX0):
    NewX0[i] = npy.insert(X0, 0, X0_primordial[i])
#
# plots
plot_path = 'sp/plots'
fig, ax = plt.subplots(figsize=(5.5, 5))
color = ['Crimson', 'Purple', 'SlateGray', 'YellowGreen',
         'DodgerBlue', 'SaddleBrown']
#
for ielt, elt in enumerate(Elements):
    legend = Elements_Tex[ielt]
    x = [math.log10(Metalicity[i]) for i in [1, 3, 4, 5]]
    y = [math.log10(NewX0[ielt][i]) for i in [1, 3, 4, 5]]
    plt.scatter(x, y, s=40., c=color[ielt], marker='o',
                edgecolor=color[ielt], label=legend)
    x = [math.log10(Metalicity[i]) for i in [2, 6]]
    y = [math.log10(NewX0[ielt][i]) for i in [2, 6]]
    plt.scatter(x, y, s=40., c=color[ielt], marker='o',
                facecolor='none', edgecolor=color[ielt])
xmin = -5
xmax = -1
plt.xlim([xmin, xmax])                                               # xrange
ax.arrow(xmin, math.log10(X0_primordial[0]), 0.2, 0., head_width=0.05,
         head_length=0.1, fc=color[0], ec=color[0])
ax.arrow(xmin, math.log10(X0_primordial[1]), 0.2, 0., head_width=0.05,
         head_length=0.1,fc=color[1], ec=color[1])
plt.xlabel('Metallicity $log_{10}(Z)$ (mass fraction)')              # label for x axis
plt.ylabel('Element abundances $log_{10}(\chi_i)$ (mass fraction)')  # label for y axis
plt.legend(loc='best')                                               # display legend
plt.grid(True,which="both",ls="-")                                   # display a grid with all sub log tics                               
filename = os.path.join(plot_path, 'element_abundances.pdf')
plt.savefig(filename)             	                                 # save the figure
#
# Save data
retreated_data_path = 'sp/data/retreated'
filename = os.path.join(retreated_data_path, 'element_abundances.dat')
# create the header
hdr = 'Element abundances (retreated from Karakas+2010)'
hdr += 'for different metallicity bins\n'
hdr += 'X0_i is expressed in mass fraction\n'
hdr += 'Values for Z = 0.0004 and Z = 0.05 are interpolated and extrapolated\n'
hdr += '1) Z (mass fraction) |  2) X0_H1 |  3) X0_He4 |  4) X0_C12 |  5) X0_N14 |  6) X0_O16  | 7) X0_Fe56'
# Print data
D = []
D.append(Metalicity)
for X0 in NewX0:
    D.append(X0)
npy.savetxt(filename, npy.column_stack(D), header=hdr, fmt='%17.8e'*7)
