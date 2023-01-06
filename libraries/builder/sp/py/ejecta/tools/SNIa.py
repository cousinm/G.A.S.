#!/usr/bin/env python
#
import os                  # to build path
import numpy as npy        # array
import math
import tools.IMF as imf    # initial mass funcito
#
def load():
    #
    # Define Elements
    Elements = ['H1', 'He4', 'C12', 'N14', 'O16', 'Fe56', 'Metals', 'Total']
    #
    # Load SNIa ejecta
    data_path = 'sp/data/originals/Iwamoto'
    filename = os.path.join(data_path, 'SNIa_synthesized_mass.dat')
    tmp = npy.loadtxt(filename, unpack=True, usecols=[1])
    SNIa_ejct = {}
    for ielt, elt in enumerate(Elements):
        SNIa_ejct[elt] = tmp[ielt]
    #
    return SNIa_ejct
#
def f(mu):
    g = 2.0
    return 2.**(1. + g) * (1. + g) * mu**g
#
def norm():
    return 0.05
#
def integrator(m, imfRef):
    #
    # Define parameters for SNIa binnary systems
    MB_min = 3.0       # in Msun
    MB_max = 16.0      # in Msun
    # if MB_max = 16, then the secondary component m2 can be more massive than MB_max/2 = 8.0
    M_max = MB_max/2.  # in Msun
    M_min = 1.0        # in Msun
    # define the distribution function
    # Binnary SNIa systems
    # In the binnary system @ t m2 = m
    m2 = m
    dm = 0.5 # Msun
    # Init to null value
    i = 0.
    #
    if (m2 > M_min) and (m2 < M_max):
        # There are binnary system which can generate SNIa
        # Compute integration borders
        MB_inf = max(2.*m2, MB_min)
        MB_sup = 0.5 * MB_max + m2

        if MB_inf < MB_sup:
            # create mass integration table
            Nsteps = int(math.floor((MB_sup - MB_inf) / dm)) +1
            MB_scale = npy.linspace(MB_inf, MB_sup, Nsteps)
            for im in range(len(MB_scale)-1):
                MB_d  = MB_scale[im]
                int_d = imf.IMF(MB_d, imfRef) * f(m2 / MB_d) / MB_d**2.
                #
                MB_u  = MB_scale[im + 1]
                int_u = imf.IMF(MB_u, imfRef) * f(m2 / MB_u) / MB_u**2.
                dMB = MB_u - MB_d
                i += 0.5 * dMB * (int_u + int_d)
        if MB_inf == MB_sup:
            i += imf.IMF(MB_sup, imfRef) * f(0.5) / MB_sup**2.
    #
    return i
