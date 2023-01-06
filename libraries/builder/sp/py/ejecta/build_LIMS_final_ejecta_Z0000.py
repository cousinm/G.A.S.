#!/usr/bin/env python
# join (LIMS1 + LIMS2) tot ejecta associated to a single stellar population
# for Z = 0.000
#
# import library
import os                                                            # to build path
import numpy as npy                                                  # array 
#
# Data path
retreated_data_path = 'sp/data/retreated'
#
Elements = ['H1', 'He4', 'C12', 'N14', 'O16','Fe56', 'Total']
#
# Read LIMS1
filename = os.path.join(retreated_data_path, 'selected_Final_Ejecta_LIMS1_Z0000.dat')
# Load data
initialMassLIMS1 = npy.loadtxt(filename, unpack=True, usecols=[0])
tmp = npy.loadtxt(filename, unpack=True, usecols=[1, 2, 3, 4, 5, 6, 7])
tot_Ejct = {}
for ielt, elt in enumerate(Elements):
    tot_Ejct[elt] = tmp[ielt]
#
# Read LIMS2
filename = os.path.join(retreated_data_path, 'selected_Final_Ejecta_LIMS2_Z0000.dat')
# Load data
initialMassLIMS2 = npy.loadtxt(filename, unpack=True, usecols=[0])
tmp = npy.loadtxt(filename, unpack=True, usecols=[1, 2, 3, 4, 5, 6, 7])
for ielt, elt in enumerate(Elements):
    tot_Ejct[elt] = npy.concatenate([tot_Ejct[elt], tmp[ielt]])
#
# Create joined mass table
initialMass = npy.append(initialMassLIMS1, initialMassLIMS2)
#
# Reformat data
D = []
D.append(initialMass)
for elt in Elements:
    D.append(tot_Ejct[elt])
filename = os.path.join(retreated_data_path, 'selected_Final_Ejecta_LIMS_Z0000.dat')
hdr = 'Final stellar ejecta (AGB)\n'
hdr += 'Computed for a range of low and intermediate mass stars\n'
hdr += ' by Campbell+2008 and Gil-Pons+2013\n'
hdr += 'Initial stellar mass from 0.85 to 9 Msun\n'
hdr += 'Are listed here only main ISM elements H1, He4, C12, N14, O16 and Fe56\n'
hdr += 'Total-ejected mass is also added\n'
hdr += 'Z=0.0000\n'
hdr += '1) Mass, 2) Ej_H1, 3) Ej_He4, 4) Ej_C12, 5) Ej_N14, 6) Ej_O16, 7) Ej_Fe56, 8) Ej_tot'
npy.savetxt(filename, npy.column_stack(D), header=hdr, fmt='%17.8e'*8)
