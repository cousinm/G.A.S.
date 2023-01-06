#!/usr/bin/env python
#
import os                              # to build path
import numpy as npy                    # array
from scipy.interpolate import InterpolatedUnivariateSpline as interp # for interpolation
#
original_data_path = 'sp/data/originals/Heger'
retreated_data_path = 'sp/data/retreated'
#
# Load original data file
filename = os.path.join(original_data_path,'ejected_mass.dat')
data = npy.loadtxt(filename, unpack=True, usecols=[0, 1, 3, 5, 6], comments='#') 
Cut = npy.loadtxt(filename, unpack=True, usecols=[2], comments='#', dtype=str) 
Mass = data[0]
Energy = data[1]
Mix = data[2]
Aelt = data[3]
Ej = data[4]
#
# Elements table
Elements = ['H1', 'He4', 'C12', 'N14', 'O16', 'Fe56']
Zelements = Elements[2:]
#
# Correspondance table
A2Elts = {1: 'H1', 4: 'He4', 12: 'C12', 14: 'N14', 16: 'O16', 56: 'Fe56'}
#
# Create total ejecta table
tot_Ejct = {}
# Create elements
for val in A2Elts.values():
    tot_Ejct[val] = []
# Add total
tot_Ejct['Total'] = []
#
# Run on datas and extract usefull parameters
# Start with the first mass
initialMass = [Mass[0]]
remnantMass = []
total = 0.
for i, m in enumerate(Mass):
    if (Energy[i] == 0.9 and
        Cut[i] == 'S4' and
            Mix[i] == 0.00251):
        ej = Ej[i]
        if ej < 1.e-20:
            ej = 0.
        # Sum total ejected mass
        total += ej
        elt = A2Elts.get(Aelt[i], None)
        if elt is not None:
            tot_Ejct[elt].append(ej)
        if m not in initialMass:
            initialMass.append(m)
        if Mass[i+1] not in initialMass:
            # Dump total ejected mass for previous mass
            tot_Ejct['Total'].append(total)
            # Deduce remnant mass from initial mass and total ejected mass
            remnantMass.append(m - total)
            total = 0.

# Dumps latest total ejected mass
tot_Ejct['Total'].append(total)
# Deduce remnant mass from initial mass and total ejected mass
remnantMass.append(m - total)
total = 0.
#
# Add Total to Element lits
Elements.append('Total')
#
# Select Portinari closest initial mass and limit to this set of masses
# Undersample
# Load Portinari+98 initial mass list
filename = os.path.join(retreated_data_path, 'initial_final_mass_MS.dat')
MSMass = npy.loadtxt(filename, unpack=True, usecols=[0], comments='#') 
# Limit to M = 100
MSMass = MSMass[:-1]
# Select closest mass
selection = []
for m in MSMass:
   selection.append(npy.argmin(npy.abs(npy.array(initialMass) - m)))
#
# Keep only selected mass
initialMass = npy.array(initialMass)[selection]
remnantMass = npy.array(remnantMass)[selection]
for elt, ej in tot_Ejct.items():
    tot_Ejct[elt] = npy.array(ej)[selection]
#
# Save data
D = []
D.append(initialMass)
for elt in Elements:
    D.append(tot_Ejct[elt])
filename = os.path.join(retreated_data_path, 'selected_Final_Ejecta_MS_Z0000.dat')
hdr = 'Final ejecta (SNII)\n'
hdr += 'Computed for a range of massive stars by Heger+10\n'
hdr += ' with initial stellar mass from 10 to 100 Msun\n'
hdr += 'Are listed here only main ISM elements H1, He4, C12, N14, O16 and Fe56\n'
hdr += 'Total-ejected mass is also added\n'
hdr += 'Ej_i < 1.e-10 Msun are set to 0\n'
hdr += 'Z=0.000\n'
hdr += '1) Mass, 2) Ej_H1, 3) Ej_He4, 4) Ej_C12, 5) Ej_N14, 6) Ej_O16, 7) Ej_Fe56, 8) Ej_tot'
npy.savetxt(filename, npy.column_stack(D), header=hdr, fmt='%17.8e'*8)
#
# Save remnant mass as function of the initial mass
D = []
D.append(initialMass)
D.append(remnantMass)
filename = os.path.join(retreated_data_path, 'initial_remnant_mass_MS_Z0000.dat')
hdr = 'Initial and Remnant mass\n'
hdr += 'Computed for a range of massive stars by Heger+10\n'
hdr += ' with initial stellar mass from 10 to 100 Msun\n'
hdr += '1) Minit, 2) Mrem'
npy.savetxt(filename, npy.column_stack(D), header=hdr, fmt='%17.8e'*2)
