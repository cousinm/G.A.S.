#!/usr/bin/env python
#
# Data are extracted from Portinari et al 1998
# Interpolate the total ejecta for all elements at Z = 0.004 for stellar mass m = 60.
# Interpolate the total ejecta for all elements at Z = 0.008 for stellar mass m = 100.
# Interpolate the total ejecta for all elements at Z = 0.05 for stellar mass m = 30.
# Extrapolate the total ejecta for all elements for Z = 0.0001 for all initial stellar mass
#
# import library
import os                                                            # to build path
import numpy as npy                                                  # array
import math                                                          # for math operations
from scipy.interpolate import InterpolatedUnivariateSpline as interp # for interpolation
#
# Path Portinari+98
original_data_path = 'sp/data/originals/Padova_Tracks_Portinari'
retreated_data_path = 'sp/data/retreated'
#
# Read the mass table
initialMass = npy.loadtxt(os.path.join(original_data_path, 'mass.dat'), unpack=True)
#
# Define elements list
Elements = ['H1', 'He4', 'C12', 'N14', 'O16', 'Fe56']
#
# Define initial metallicity list
MetBins = [0.0004, 0.004, 0.008, 0.02, 0.05]
# Define metallicity ref
MetBinsRef = ['0004','004','008','02','05']
#
# Create a dict for ejecta indexed by metallicity and elements
# Complete with total ejected mass, deduced from the initial/remnant masses
tot_Ejct = {}
# Read remnant mass
filename = os.path.join(retreated_data_path, 'initial_remnant_mass_MS.dat')
# Data are stored following 1) metallicity 2) mass bins
# Remnant mass for Z=0.0001, previously extrapolated is not used
remnantMass = npy.loadtxt(filename, unpack=True, usecols=([2, 3, 4, 5, 6]))
#
for iz, Z in enumerate(MetBinsRef):
    # Read original data
    filename = f'selected_Total_Ejecta_MS_Z{Z}.dat'
    # Load data
    tmp = npy.loadtxt(os.path.join(original_data_path, filename),
                      unpack=True, usecols=[1, 2, 3, 4, 5, 6])
    dct = {}
    for ielt, elt in enumerate(Elements):
        dct[elt] = tmp[ielt]
    dct['Total'] = initialMass - remnantMass[iz]
    tot_Ejct[Z] = dct
#
# Add Total to Elements
Elements.append('Total')
#
# Define guess values for empty data, those values will be then corrected
# according to mass conservation laws
#
# Interpolate the total ejected mass for all elements at Z = 0.004
# for stellar mass m = 60: tot_Ejct['004'][:][5]
for elt, ej in tot_Ejct['004'].items():
    m = [initialMass[i] for i in [0, 1, 2, 3, 4, 6, 7]]
    e = [math.log10(ej[i]) for i in [0, 1, 2, 3, 4, 6, 7]]
    # Interpolate the total ejecta for a given element for stellar mass m = 60
    intrpl = interp(m, e, k=2)
    ej60 = 10.**(intrpl(60.0))
    # Replace the data in tot_Ejct
    tot_Ejct['004'][elt][5] = ej60
#
# Interpolate the total ejected mass for all elements at Z = 0.008
# for stellar mass m = 100: tot_Ejct['008'][:][6]
for elt, ej in tot_Ejct['008'].items():
    m = [initialMass[i] for i in [0, 1, 2, 3, 4, 5, 7]]
    e = [math.log10(ej[i]) for i in [0, 1, 2, 3, 4, 5, 7]]
    # Interpolate the total ejecta for a given element for stellar mass m = 60
    intrpl = interp(m, e, k=1)
    ej100 = 10.**(intrpl(100.0))
    # Replace the data in tot_Ejct
    tot_Ejct['008'][elt][6] = max(0., ej100)
#
# Interpolate the total ejected mass for all elements at Z = 0.05
# for stellar mass m = 30: tot_Ejct['05'][:][3]
for elt, ej in tot_Ejct['05'].items():
    # Create the Ejecta table for interpolation as function of the stellar mass (at a given metallicity) 
    m = [initialMass[i] for i in [0, 1, 2, 4, 5, 6, 7]]
    e = [math.log10(ej[i]) for i in [0, 1, 2, 4, 5, 6, 7]]
    # Interpolate the total ejecta for a given element for stellar mass m = 60
    intrpl = interp(m, e, k=2)
    ej30 = 10.**(intrpl(30.0))
    # Replace the data in tot_Ejct
    tot_Ejct['05'][elt][3] = max(0., ej30)
#
# Extrapolate ejecta (as function of the metallicity)
# for all elements and all masses at Z = 0.0001
# To do that we complete the Portinari+98 data set with total ejecta
# computed by Heger+10 for metal free stars
# Load data
filename = f'selected_Final_Ejecta_MS_Z0000.dat'
# Extract initial mass
massZ0000 = npy.loadtxt(os.path.join(retreated_data_path, filename),
                        unpack=True, usecols=[0])
# Extract ejecta for all elements and masses
tmp = npy.loadtxt(os.path.join(retreated_data_path, filename),
                  unpack=True, usecols=[1, 2, 3, 4, 5, 6, 7])
# Interpolate Heger metal-free total ejecta
dct = {}
for ielt, elt in enumerate(Elements):
    # Interpolation is performed with log scale
    # Set null values to 1.e-20
    tmp[ielt] = npy.where(tmp[ielt] <= 0., 1.e-20, tmp[ielt])
    intrpl = interp(massZ0000, npy.log10(tmp[ielt]), k=1)
    ej = 10.**(intrpl(initialMass))
    # Reset null values
    ej = npy.where(ej <= 1.e-20, 0., ej)
    dct[elt] = ej
tot_Ejct['0000'] = dct
#
# Add the new Met bin
MetBins = npy.insert(MetBins, 0, 1.e-50)
MetBinsRef = npy.insert(MetBinsRef, 0, '0000')
#
tot_Ejct['0001'] = {}
for elt in Elements:
    ej_ = []
    for im, _ in enumerate(initialMass):
        ej = npy.array([tot_Ejct[z][elt][im] for z in MetBinsRef])
        ej = npy.where(ej <= 0., 1.e-20, ej)
        intrpl = interp(npy.log10(MetBins), npy.log10(ej), k=1)
        ej = 10.**(intrpl(math.log10(0.0001)))
        ej = npy.where(ej <= 1.e-20, 0., ej)
        ej_.append(float(ej))
    tot_Ejct['0001'][elt] = npy.array(ej_)
#
# Add the new Met bin
MetBins = npy.insert(MetBins, 1, 0.0001)
MetBinsRef = npy.insert(MetBinsRef, 1, '0001')
# Remove metal-free data
MetBins = npy.delete(MetBins, 0)
MetBinsRef = npy.delete(MetBinsRef, 0)
tot_Ejct.pop('0000')
#
# Save data
for iz, Z in enumerate(MetBinsRef):
# Reformat data
    D = []
    D.append(initialMass)
    for elt in Elements:
        D.append(tot_Ejct[Z][elt])
    filename = os.path.join(retreated_data_path, f'selected_Final_Ejecta_MS_Z{Z}.dat')
    hdr = 'Total stellar ejecta (SNII + over + wind)\n'
    hdr += 'Computed for a range of massive stars by Potinari+98 with Padova Tracks for a set of metallicity\n'
    hdr += '  and a range of initial stellar mass from 12 to 100 Msun\n'
    hdr += 'Are listed here only main ISM elements H1, He4, C12, N14, O16 and Fe56\n'
    hdr += 'Total-ejected mass is also added\n'
    hdr += f'Z=0.{Z}\n'
    hdr += '1) Mass, 2) Ej_H1, 3) Ej_He4, 4) Ej_C12, 5) Ej_N14, 6) Ej_O16, 7) Ej_Fe56, 8) Ej_tot'
    npy.savetxt(filename, npy.column_stack(D), header=hdr, fmt='%17.8e'*8)
