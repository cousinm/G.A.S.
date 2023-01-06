#!/usr/bin/env python
# stellar_lifetimes.py
# Allow to load stellar lifetimes and compute interpolation splines
#
import os                                                            # to build path
import numpy as npy                                                  # array
from scipy.interpolate import InterpolatedUnivariateSpline as interp # for interpolation

def load():
    #
    # Define metallicity list
    MetBinsRef = ['0001', '0004', '004', '008', '02', '05']
    #
    data_path = 'sp/data/retreated'
    #
    filename = os.path.join(data_path, 'extended_Stellar_Lifetimes.dat')
    # Load data
    initialMass = npy.loadtxt(filename, unpack=True, usecols=[0])
    # Load stellar lifetimes
    tmp = npy.loadtxt(filename, unpack=True, usecols=[1, 2, 3, 4, 5, 6])
    slt = {}
    for iz, Z in enumerate(MetBinsRef):
        slt[Z] = tmp[iz]
    #
    # Stellar lifetime at Z=0.0000 is assumed to be similar to Z=0.0001
    slt['0000'] = slt['0001']
    #
    # Add the new metallicity bin
    MetBinsRef = npy.insert(MetBinsRef, 0, '0000')
    #
    # Create interpolation splines for stellar lifetimes
    slt2m = {}
    m2slt = {}
    for Z, s in slt.items():
        m2slt[Z] = interp(initialMass, s, k=1)
        slt2m[Z] = interp(s[::-1], initialMass[::-1], k=1)
    # return
    return slt2m, m2slt
