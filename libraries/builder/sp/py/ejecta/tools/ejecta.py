#!/usr/bin/env python
#
import tools.IMF as imf    # initial mass funciton
import tools.slt as slt    # stellar life-times
import tools.SNIa as snIa  # SNIa
#
from scipy.interpolate import InterpolatedUnivariateSpline as interp # for interpolation
from os import path
from numpy import loadtxt
#
# Define metallicity list
MetBins = [0.0000, 0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]
MetBinsRef = ['0000', '0001', '0004', '004', '008', '02', '05']
#
# Define Elements
Elements = ['H1', 'He4', 'C12', 'N14', 'O16', 'Fe56', 'Metals', 'Total']
#
# Load stellar lifestimes
slt2m, m2slt = slt.load()
#
# Load SNIa ejected mass
SNIa_ejct = snIa.load()
#
# Load ejected mass
global m2ej
m2ej = None
#
# Define mass range
mMin = 0.8
mMinSNII = 8.0
mMax = 100.
#
# Create time table
dt = 1.e-4      # Gyr
ddt = dt/4.


def load():
    #
    data_path = 'sp/data/retreated'
    #
    # Load data
    final_Ejct = {}
    for Z in MetBinsRef:
        final_Ejct[Z] = {}
        filename = path.join(data_path, f'selected_Final_Ejecta_Z{Z}.dat')
        # Initial mass
        initialMass = loadtxt(filename, unpack=True, usecols=[0])
        # Ejecta
        tmp = loadtxt(filename, unpack=True, usecols=[1, 2, 3, 4, 5, 6, 7, 8])
        for ielt, elt in enumerate(Elements):
            final_Ejct[Z][elt] = tmp[ielt]
    #
    # Create interpolation splines
    m2ej = {}
    for Z in MetBinsRef:
        m2ej[Z] = {}
        for elt in Elements:
            m2ej[Z][elt] = interp(initialMass, final_Ejct[Z][elt], k=1)
    # return     
    return m2ej


def dMdt(t, m, Z, ddt):
    r = 0.
    if m <= mMax and m > mMin:
        td = t - ddt
        tu = t + ddt
        md = slt2m[Z](td)
        mu = slt2m[Z](tu)
        r = abs(mu - md) / abs(tu - td)
    return r


def compute(t, Z, elt, imfRef):
    """
    Compute instantaneous final ejecta rate
    In output are given the full average ejecta rate
    and the associated specific SNIa average ejecta rate
    """
    global m2ej
    if m2ej is None:
        m2ej = load()
    #
    ejr = 0.
    ejr_SNIa = 0.
    #
    m = slt2m[Z](t)
    dM_dt = dMdt(t, m, Z, ddt)
    # Instantaneous final ejecta rate
    ejr = m2ej[Z][elt](m) * imf.IMF(m, imfRef) * dM_dt / m
    #
    # SNII
    SNII = 0.
    if (m>=8. and m<=mMax):
        SNII = imf.IMF(m, imfRef) * dM_dt / m
    #
    SNIa = dM_dt * snIa.integrator(m, imfRef)
    ejr_SNIa = snIa.norm() * SNIa_ejct[elt] * SNIa
    ejr = (1. - snIa.norm()) * ejr + ejr_SNIa
    #
    # SNIa + SNII
    SNR = SNII + snIa.norm() * SNIa
    return ejr, ejr_SNIa, SNR


def compute_average(t0, t1, Z, elt, imfRef):
    """
    Compute average ejecta rate produced between t0 and t1
    """
    # Compute timescale
    timescale = [t0]
    t = t0
    while t < t1:
        t = min(t1, t + dt)
        timescale.append(t)
    #
    ejct = []
    ejctSNIa = []
    SNRate = []
    av_ejct = 0.
    av_ejctSNIa = 0.
    av_Mej = 0.
    av_SNR = 0.
    #
    for it, t in enumerate(timescale):
        #
        ejr, ejr_SNIa, SNR = compute(t, Z, elt, imfRef)
        SNRate.append(SNR)
        ejct.append(ejr)
        ejctSNIa.append(ejr_SNIa)
        #
        # Trapez integration
        if it > 0:
            cdt = timescale[it] - timescale[it - 1]
            av_ejct += cdt * (ejct[it - 1] + ejct[it])/2.
            av_Mej += cdt * av_ejct
            av_ejctSNIa += cdt * (ejctSNIa[it - 1] + ejctSNIa[it])/2.
            av_SNR += cdt * (SNRate[it - 1] + SNRate[it])/2.
    #
    # Normalization
    av_ejct /= (t1 - t0)
    av_ejctSNIa /= (t1 - t0)
    av_Mej /= (t1 - t0)
    av_SNR /= (t1 - t0)
    #
    return timescale, av_ejct, av_ejctSNIa, av_Mej, av_SNR, ejct, ejctSNIa
