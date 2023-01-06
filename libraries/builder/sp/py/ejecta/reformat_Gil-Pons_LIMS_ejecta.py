#!/usr/bin/env python
# Reformat Gil-Pons LIMS ejecta
#
import os                        # to build path
import numpy as npy              # array
#
# Path Gil-Pons
original_data_path = 'sp/data/originals/Gil-Pons_Campbell'
retreated_data_path = 'sp/data/retreated'
# 
Elements = ['H1', 'He4', 'C12', 'N14', 'O16']
#
# Load original data file
filename = os.path.join(original_data_path, 'ejecta_LIMS_Z00001[Gil-Pons].dat')
# Ej_tmp[elt][m] ejected mass associated to a specie
tmp = npy.loadtxt(filename, unpack=True, usecols=[2, 3, 4, 5, 6])
tot_Ejct = {}
for ielt, elt in enumerate(Elements):
    tot_Ejct[elt] = tmp[ielt]
#
initialMass = npy.loadtxt(filename, unpack=True, usecols=[0])
envelopMass = npy.loadtxt(filename, unpack=True, usecols=[1])
#
# Add Fe56 and total ejected mass
total = []
for im, M in enumerate(initialMass):
    total_im = 0.
    for elt in Elements:
        ej = tot_Ejct[elt][im]
        total_im += ej
    total.append(total_im)
tot_Ejct['Fe56'] = npy.zeros(initialMass.size)
tot_Ejct['Total'] = npy.array(total)
#
# Add Fe56 and Total to Elements
Elements.append('Fe56')
Elements.append('Total')
#
# Reformat data
D = []
D.append(initialMass)
for elt in Elements:
    D.append(tot_Ejct[elt])
#
filename = os.path.join(retreated_data_path, 'selected_Final_Ejecta_LIMS2_Z0000.dat')
hdr = 'Final stellar ejecta (AGB)\n'
hdr += 'Computed for a range of low and intermediate stars by Gil-Pons+2013\n'
hdr += ' with initial stellar mass 4.0, 5.0, 6.0, 7.0, 8.0 and 9.0 Msun\n'
hdr += 'Are listed here only main ISM elements H1, He4, C12, N14, O16 and Fe56\n'
hdr += 'Total-ejected mass is added\n'
hdr += 'Z=0.00001\n'
hdr += '1) Mass, 2) Ej_H1, 3) Ej_He4, 4) Ej_C12, 5) Ej_N14, 6) Ej_O16, 7) Ej_Fe56, 8) Ej_tot'
npy.savetxt(filename, npy.column_stack(D), header=hdr,fmt='%17.8e'*8)
#
D = []
D.append(initialMass)
filename = os.path.join(retreated_data_path, 'mass[Gil-Pons].dat')
hdr = 'Stellar initial mass Gil-Pons+13'
npy.savetxt(filename, npy.column_stack(D), header=hdr, fmt='%17.8e')
#
D = []
D.append(initialMass)
D.append(initialMass - envelopMass)
filename = os.path.join(retreated_data_path, 'initial_remnant_mass[Gil-Pons].dat')
hdr = 'Stellar initial mass and remnant associated mass Gil-Pons+13'
npy.savetxt(filename, npy.column_stack(D), header=hdr, fmt='%17.8e'*2)
