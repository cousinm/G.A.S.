#!/usr/bin/env python
#
# Karakas data are been divided in different files corresponding to metallicity and mass
# This routine allows to create new metallicity bins
# interpolate final mass for all initial stellar mass at Z = 0.0004
# extrapolate final mass for all initial stellar mass at Z = 0.05
#
# import library
from copy import deepcopy
import os                        # to build path
import numpy as npy              # array
import math                      # for math operations
from scipy.interpolate import InterpolatedUnivariateSpline as interp # for interpolation
#
# define MetBins
MetBins = [0.0001, 0.004, 0.008, 0.02]
# define metallicity bin references
MetBinsRef = ['0001', '004', '008', '02']
#
# path Karakas+10
original_data_path = 'sp/data/originals/Karakas'
retreated_data_path = 'sp/data/retreated'
#
# load the initial mass table
initialMass = npy.loadtxt(os.path.join(original_data_path, 'mass.dat'))
#
Elements = ['H1', 'He4', 'C12', 'N14', 'O16', 'Fe56']
#
# Read original files
# Complete with total ejected mass, deduced from the initial/remnant masses
tot_Ejct = {}
#
# Read remnant mass
filename = os.path.join(retreated_data_path, 'initial_remnant_mass_LIMS.dat')
# Data are stored following 1) metallicity 2) mass bins
# Data for Z = 0.0004 and Z = 0.05 previously extrapolated are not red
remnantMass = npy.loadtxt(filename, unpack=True, usecols=([1, 3, 4, 5]))
tot_Ejct = {}  # tot_Ejct[Z][m][elt]
for iz, Z in enumerate(MetBinsRef):
    # Run onto various metalicity reference files
    filename = f'selected_Final_Ejecta_LIMS_Z{Z}.dat'
    filename = os.path.join(original_data_path, filename)
    tmp = npy.loadtxt(filename, unpack=True, usecols=[1, 2, 3, 4, 5, 6])
    tot_Ejct[Z] = {}
    for ielt, elt in enumerate(Elements):
        tot_Ejct[Z][elt] = tmp[ielt]
    tot_Ejct[Z]['Total'] = initialMass - remnantMass[iz]
#
# Add Total to Elements
Elements.append('Total')
#
# Interpolate (as a function of the metallicity)
# the ejected mass for all initial stellar mass at Z = 0.0004
tot_Ejct['0004'] = {}
for elt in Elements:
    ej_ = []
    for im, M in enumerate(initialMass):
        ej = npy.array([tot_Ejct[z][elt][im] for z in MetBinsRef])
        ej = npy.where(ej <= 0., 1.e-20, ej)
        intrpl = interp(npy.log10(MetBins), npy.log10(ej), k=1)
        ej = 10.**(intrpl(math.log10(0.0004)))
        ej = npy.where(ej <= 1.e-15, 0., ej)
        ej_.append(float(ej))
    tot_Ejct['0004'][elt] = npy.array(ej_)
#
# Extrapolate (as a function of the metallicity)
# the ejected mass for all initial stellar mass at Z = 0.05
tot_Ejct['05'] = {}
for elt in Elements:
    ej_ = []
    for im, M in enumerate(initialMass):
        ej = npy.array([tot_Ejct[z][elt][im] for z in MetBinsRef])
        ej = npy.where(ej <= 0., 1.e-20, ej)
        intrpl = interp(npy.log10(MetBins), npy.log10(ej), k=1)
        ej = 10.**(intrpl(math.log10(0.05)))
        ej = npy.where(ej <= 1.e-15, 0., ej)
        ej_.append(float(ej))
    tot_Ejct['05'][elt] = npy.array(ej_)
#
# Add the new metallicty bin
MetBins = npy.insert(MetBins, 1, 0.0004)
MetBinsRef = npy.insert(MetBinsRef, 1, '0004')
#
# Add the new metallicty bin
MetBins = npy.insert(MetBins, 5, 0.05)
MetBinsRef = npy.insert(MetBinsRef, 5, '05')
#
# Save data in files
for Z in MetBinsRef:
    D = []
    D.append(npy.array(initialMass))
    for elt in Elements:
        D.append(tot_Ejct[Z][elt])
    filename = os.path.join(retreated_data_path, f'selected_Final_Ejecta_LIMS_Z{Z}.dat')
    hdr = 'Final ejected mass (during the TP-AGB phase)\n'
    hdr += '  computed for a set of metallicity and a range of low/intermediate mass by Karakas+2010'
    hdr += 'Are listed here only main ISM elements H1, He4, C12, N14, O16 and Fe56\n'
    hdr += 'Total-ejected mass is also added\n'
    hdr += f'Z=0.{Z}\n'
    hdr += '1) Mass, 2) Ej_H1, 3) Ej_He4, 4) Ej_C12, 5) Ej_N14, 6) Ej_O16, 7) Ej_Fe56, 8) Ej_tot'
    npy.savetxt(filename, npy.column_stack(D), header=hdr, fmt='%17.8e'*8)
