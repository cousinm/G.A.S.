#!/usr/bin/env python
# Reformat Campbell LIMS ejecta
# import library
import glob                      # to list files in a path
import os                        # to build path
import numpy as npy              # array
#
# path Campbell
original_data_path = 'sp/data/originals/Gil-Pons_Campbell'
retreated_data_path = 'sp/data/retreated'
#
# Load original data file
filename = os.path.join(original_data_path, 'mass[Campbell].dat')
initialMass = npy.loadtxt(filename, unpack=True, usecols=[0])
#
filename = os.path.join(original_data_path, 'initial_remnant_mass[Campbell].dat')
remnantMass = npy.loadtxt(filename, unpack=True, usecols=[1])
#
filename = os.path.join(original_data_path, 'yields_LIMS_metal_free[Campbell].dat')
# WARNING in input data yields are given in mass fraction
# of each species in the total ejecta.
# To compute total ejected mass, we have to take into the total ejected mass
# given by (M0 - Mrem)
# Ej_frac[m][elt] ejected mass fraction associated to a specie
AElts = npy.loadtxt(filename, unpack=True, usecols=[2])
Ej_frac = npy.loadtxt(filename, unpack=True, usecols=[4, 5, 6, 7]).T
# Define elements table
Elements = ['H1', 'He4', 'C12', 'N14', 'O16','Fe56', 'Total']
# Table de correspondance
A2Elt = {1: 'H1', 4: 'He4', 12: 'C12', 14: 'N14', 16: 'O16', 56: 'Fe56'}
#
# Convert ejecta fration in ejected mass
tot_Ejct = {}
tot_Ejct['Total'] = npy.zeros(initialMass.size)
total = 0.
for ia, A in enumerate(AElts):
    ej = []
    for im, M in enumerate(initialMass):
        Menv = initialMass[im]-remnantMass[im]
        tot_Ejct['Total'][im] += Ej_frac[ia][im] * Menv
        if A in list(A2Elt.keys()):
            ej.append(Ej_frac[ia][im] * Menv)
    if ej:
        elt = A2Elt.get(A)
        tot_Ejct[elt] = npy.array(ej)

# Reformat data
D = []
D.append(initialMass)
for elt in Elements:
    D.append(tot_Ejct[elt])
filename = os.path.join(retreated_data_path, 'selected_Final_Ejecta_LIMS1_Z0000.dat')
hdr = 'Final stellar ejecta (AGB)\n'
hdr += 'Computed for a range of low and intermediate stars by Campbell+2008\n'
hdr += ' with initial stellar mass 0.85, 1.00, 2.00 and 3.00 Msun\n'
hdr += 'Are listed here only main ISM elements H1, He4, C12, N14, O16 and Fe56\n'
hdr += 'Total-ejected mass is added\n'
hdr += 'Z=0.000\n'
hdr += '1) Mass, 2) Ej_H1, 3) Ej_He4, 4) Ej_C12, 5) Ej_N14, 6) Ej_O16, 7) Ej_Fe56, 8) Ej_tot'
npy.savetxt(filename, npy.column_stack(D), header=hdr, fmt='%17.8e'*8)
