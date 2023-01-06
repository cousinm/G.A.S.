#!/usr/bin/env python
# join (LIMS + QMS + MS) tot ejecta associated to a single stellar population
#
import os                        # to build path
import numpy as npy              # array
from scipy.interpolate import InterpolatedUnivariateSpline as interp # for interpolation
import matplotlib.pyplot as plt  # plot
import matplotlib.gridspec as gridspec
#
Elements = ['H1', 'He4', 'C12', 'N14', 'O16', 'Fe56']
fullElements = ['H1', 'He4', 'C12', 'N14', 'O16', 'Fe56', 'Total']
ZElements = ['C12', 'N14', 'O16', 'Fe56']
#
# dDefine MetBins
MetBins = [0.0000, 0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]
# Define metallicity bin references
MetBinsRef = ['0000', '0001', '0004', '004', '008', '02', '05']

# Create the data structure
final_Ejct = {}
#
retreated_data_path = 'sp/data/retreated'
#
# Load LIMS
for Z in MetBinsRef:
    #
    filename = os.path.join(retreated_data_path, f'selected_Final_Ejecta_LIMS_Z{Z}.dat')
    # Read Ejecta
    tmp = npy.loadtxt(filename, unpack=True, usecols=([1, 2, 3, 4, 5, 6, 7]))
    final_Ejct[Z] = {}
    for ielt, elt in enumerate(fullElements):
        final_Ejct[Z][elt] = tmp[ielt]
    # Z=0.0000 data are specific initial mass
    if Z in ['0000']:
        initialMassZ0000 = npy.loadtxt(filename, unpack=True, usecols=([0]))
    else:
        initialMass = npy.loadtxt(filename, unpack=True, usecols=([0]))
#
# Add MS
for Z in MetBinsRef:
    #
    filename = os.path.join(retreated_data_path, f'selected_Final_Ejecta_MS_Z{Z}.dat')
    tmp = npy.loadtxt(filename, unpack=True, usecols=([0]))
    if Z in ['0000']:
        initialMassZ0000 = npy.concatenate([initialMassZ0000, tmp])
    if Z in ['0001']:
        initialMass = npy.concatenate([initialMass, tmp])
    # Read Ejecta
    tmp = npy.loadtxt(filename, unpack=True, usecols=([1, 2, 3, 4, 5, 6, 7]))
    for ielt, elt in enumerate(fullElements):
        final_Ejct[Z][elt] = npy.concatenate([final_Ejct[Z][elt], tmp[ielt]])
#
# Create interpolation spline (order 2) for each elements + total
# Ejecta are interpolated using log scale
interpol = {}
for iz, Z in enumerate(MetBinsRef):
    interpol[Z] = {}
    for elt, ej in final_Ejct[Z].items():
        # Reset null values to 1.e-20
        ej = npy.where(ej <= 0., 1.e-20, ej)
        if Z in ['0000']:
            Mass = initialMassZ0000
        else:
            Mass = initialMass
        s = interp(Mass, npy.log10(ej), k=2)
        interpol[Z][elt] = s
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
    for elt, ej in final_Ejct[Z].items():
        ej = 10.**(interpol[Z][elt](initialMass))
        ej = npy.where(ej <= 1.e-20, 0., ej)
        final_Ejct[Z][elt] = ej
#
# Check:
# - H1 production = 0.
# - Total ejected mass < initial mass
# Apply correction factor if necessary
# Compute: 
#  - The metal ejected mass 1. - H1 - He4
#
# Load initial abundance table
filename = os.path.join(retreated_data_path, 'element_abundances.dat')
# Data are stored following 1) elements 2) metallicity bins, XO[elt][Z]
tmp = npy.loadtxt(filename, unpack=True, usecols=([1, 2, 3, 4, 5, 6])).T
X0 = {}
for iz, Z in enumerate(MetBinsRef):
    X0[Z] = {}
    for ielt, elt in enumerate(Elements):
        X0[Z][elt] = tmp[iz][ielt]
#
for iz, Z in enumerate(MetBinsRef):
    metals = []
    for im, M in enumerate(initialMass):
        # Total ejecta
        Mej = final_Ejct[Z]['Total'][im]
        if Mej > M:
            # Correction
            final_Ejct[Z]['Total'][im] = 0.99*M
            Mej = final_Ejct[Z]['Total'][im]
        #
        # H1 production
        H1_ej = final_Ejct[Z]['H1'][im]
        H1_orig = Mej * X0[Z]['H1']
        if (H1_ej > H1_orig):
            dm = H1_ej - H1_orig
            # Apply correction
            final_Ejct[Z]['H1'][im] -= dm
            # Apply correction on Total
            final_Ejct[Z]['Total'][im] -= dm
            Mej = final_Ejct[Z]['Total'][im]
        #
        # Ejected mass for all tracking elements
        sumElts = sum([final_Ejct[Z][elt][im] for elt in Elements])
        if (Mej > 0.) and (sumElts >= Mej):
            dm = sumElts - Mej
            for elt in Elements:
                # Apply correction
                w = final_Ejct[Z][elt][im] / sumElts
                final_Ejct[Z][elt][im] -= dm * w
        #
        # Ejected metals
        MZej = Mej - final_Ejct[Z]['H1'][im] - final_Ejct[Z]['He4'][im]
        # Tests
        sumZElts = sum([final_Ejct[Z][elt][im] for elt in ZElements])
        if (Mej > 0.) and (sumZElts >= MZej):
            dm = sumZElts - MZej
            for elt in ZElements:
                w = final_Ejct[Z][elt][im] / sumZElts
                final_Ejct[Z][elt][im] -= dm * w
        #
        metals.append(MZej)
    # Set metal ejected mass
    final_Ejct[Z]['Metals'] = npy.array(metals)
#
# New Elements with Metals
fullElements = ['H1', 'He4', 'C12', 'N14', 'O16', 'Fe56', 'Metals', 'Total']
#
# Save data in files
for Z in MetBinsRef:
    D = []
    D.append(npy.array(initialMass))
    for elt in fullElements:
        D.append(final_Ejct[Z][elt])
    filename = os.path.join(retreated_data_path, f'selected_Final_Ejecta_Z{Z}.dat')
    hdr = 'Final ejected mass (during the TP-AGB phase)\n'
    hdr += '  computed for a set of metallicity and a range of low/intermediate mass by Karakas+2010'
    hdr += 'Are listed here only main ISM elements H1, He4, C12, N14, O16 and Fe56\n'
    hdr += 'Metal eand Total ejected mass are also added\n'
    hdr += f'Z=0.{Z}\n'
    hdr += '1) Mass, 2) Ej_H1, 3) Ej_He4, 4) Ej_C12, 5) Ej_N14, 6) Ej_O16, 7) Ej_Fe56, 8) Ej_Z 9) Ej_tot'
    npy.savetxt(filename, npy.column_stack(D), header=hdr, fmt='%17.8e'*9)
#
# Plot
plot_path = 'sp/plots'
plt.figure(figsize = (16, 9))
gs = gridspec.GridSpec(2, 4)
gs.update(wspace=0.08, hspace=0.05,
          top=0.95, bottom=0.05,
          left=0.05, right=0.95)
# Legend
# tot  : Black
# Z    : DimGray
# H1   : Crimson
# He4  : Purple
# C12  : SlateGray
# N14  : YellowGreen
# O16  : DodgerBlue
# Fe56 : SaddleBrown
color = ['Crimson', 'Purple', 'SlateGray',
         'YellowGreen', 'DodgerBlue',
         'SaddleBrown','DimGray','Black']
# Legends
Elements_Tex = ['$^{1}H$', '$^{4}He$', '$^{12}C$',
                '$^{14}N$', '$^{16}O$', '$^{56}Fe$',
                '$Metals$','$Total$']
#
for iz, Z in enumerate(MetBinsRef):
    ax = plt.subplot(gs[iz])
    for ielt, elt in enumerate(fullElements):
        # Plot interpolation spline
        ej_ = final_Ejct[Z][elt]
        ej_ = npy.where(ej_ < 1.e-20, 0., ej_)
        ax.plot(initialMass, ej_, ls='-', lw=1.5,
                color=color[ielt], label=Elements_Tex[ielt])
    #
    ax.set_xlabel('Stellar initial mass ($M_{\odot}$)')      # label for x axis
    ax.set_ylabel('Final ejected Mass ($M_{\odot}$)')        # label for y axis
    ax.set_xlim([0.8, 100.])
    ax.set_ylim([5.e-6, 100.])
    ax.set_yscale('log')           # use log scale for y axis
    ax.set_xscale('log')           # use log scale for y axis
    ax.grid(True,ls="-")           # display a grid with all sub log tics, which="both"
    # turn on/off x and/or y axis tics and labels
    if iz < 4:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_ylim([1.e-6, 100.])
    if (iz == 0) or (iz == 4):
        ax.yaxis.tick_left()
        ax.yaxis.set_label_position('left')
    if (iz == 3) or (iz == 6):
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')
    if (iz == 1) or (iz == 2) or (iz==5):
        ax.set_yticklabels([])
        ax.set_ylabel('')
    if Z in ['05']:
        plt.legend(ncol=3, loc='best', prop={'size':10})
    plt.text(0.03, 0.93, f'Z = 0.{Z}', size=12.,
             transform=ax.transAxes,
             bbox=dict(facecolor='white', edgecolor='white'))
#
# Save the plot in pdf format
filename = os.path.join(plot_path,'selected_Final_Ejecta.pdf')
plt.savefig(filename)
