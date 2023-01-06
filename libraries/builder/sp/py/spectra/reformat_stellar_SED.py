
# import library
import matplotlib.pyplot as plt         # plot
import matplotlib.gridspec as gridspec  # for multiplot grid
from os import path, environ
from numpy import loadtxt, insert, array, log10, trapz, savetxt, column_stack, ravel, where
from math import isnan
from scipy.interpolate import InterpolatedUnivariateSpline as interp # for interpolation
#
# Data path
path_BC03 = environ.get('bc03', None)
if path_BC03 is None:
    raise Exception("Please set the 'bc03' environnement variable")
#
path_spectra = path.join(path_BC03, 'spectra')
path_model = path.join(path_BC03, 'models/Padova1994_updated')
#
retreated_data_path = 'sp/data/retreated'
#
# Path to dust (to get complete wavelenght table)
path_dust = 'dust/data/retreated'
#
# Load full wavelenght table
filename = 'wavelenghts.in'
waves = loadtxt(path.join(path_dust, filename),
                unpack=True, usecols=[0])  # micron
#
# Create MetBins
MetBins = [0.0000, 0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]
MetBinsRef = ['0000', '0001', '0004', '004', '008', '02', '05']
MetBinLabs = ['22', '22', '32', '42', '52', '62', '72']
#
IMF_list = {'Chabrier+03': 'chab', 'Kroupa+93': 'kroup', 'Salpeter+55': 'salp'}
#
# Loop over IMFs
for imf in IMF_list:
    #
    new_SED = {}
    original_SED = {}
    SSP_luminosity = {}
    # Loop over metalicity
    #
    for iz, Z in enumerate(MetBinsRef):
        #
        # Load original spectra
        # Define base filename
        filename = f'bc2003_lr_BaSeL_Z{Z}_ssp.spec'
        #
        # Original wavelenghts
        original_waves = []
        with open(path.join(path_spectra, imf, filename)) as file:
            lines = file.readlines()
            #
            for il, line in enumerate(lines):
                if il == 2:
                    # Extract records for references in .4color file
                    records = array([int(r) - 1 for r in line.split()[2:]])
                if il == 3:
                    # Extract ssp ages [Gyr]
                    ages = array([float(a)/1.e9 for a in line.split()[2:]])
                    original_SED[Z] = [[] for _ in ages]
                if il > 5:
                    sline = line.split()
                    w = float(sline[0])   # Wavelenght [amgstrom]
                    # Convert wavelenght in microns
                    original_waves.append(w/1.e4)  # [microns]
                    #
                    for ia, _ in enumerate(ages):
                        # wavelengths is given in first column, so ages start at 1: ia + 1
                        # Normalisation, spectrum intensity is initialy given in Lsun/Angstom convert in Lsun
                        original_SED[Z][ia].append(float(sline[ia + 1]) * w)  # [lamb * I_lamb]
        #
        # Load associated normalization data
        filename = f'bc2003_lr_BaSeL_m{MetBinLabs[iz]}_{IMF_list[imf]}_ssp.4color'
        data = loadtxt(path.join(path_model, imf, filename), unpack=True)
        #
        # One additional spectrum is produced at 0 Gyr without any associated normalization
        # in the .4color file
        # Then data are filtered to records
        check_ages = insert(10.**(data[0])/1.e9, 0, 0.)[records]
        #
        # Get stellar mass normalisation, columns 6, 7, so id 5, 6
        # The SSP models are normalized to a total mass of 1Msun in stars at age t = 0
        # Then data are filtered to records
        living_remnant_mass_original = insert(data[5] + data[6], 0, 1.)[records]
        #
        # Build new SED with normalisation and interpolation over full wavelenght range
        new_SED[Z] = {}
        # Individual bolometric luminosity
        SSP_luminosity[Z] = []
        #
        for ia, age in enumerate(ages):
            #
            # Normalization according to the original scaling (each spectrum is then given for M*_liv + M_rem = 1 Msun)
            new_SED[Z][ia] = original_SED[Z][ia] / living_remnant_mass_original[ia]
            #
            # Interpolation over the complete wavelenght range (dust + sp)
            # Interpolation has to be done in log-log
            lw = log10(original_waves)
            ls = where(new_SED[Z][ia] <= 0., 1.e-50, new_SED[Z][ia])
            ls = log10(ls)
            intrpl = interp(lw, ls, k=1)
            iSED = 10.**(intrpl(log10(waves)))
            new_SED[Z][ia] = iSED
            #
            # # Plot  SED
            # fig = plt.figure(figsize = (7, 6))
            # gs = gridspec.GridSpec(1, 1)
            # gs.update(left=0.12, right=0.9, bottom=0.06, top=0.94)
            # #
            # ax = plt.subplot(gs[0])
            # #
            # ax.set_xscale('log')
            # plt.xlim([8.e-3, 2.5e2])
            # plt.xlabel('wavelenghts $\lambda$ [$\mu m$]')  # label for x axis
            # #
            # plt.grid(True,ls=":",which='both') # display a grid with all sub log tics, which="both"
            # #
            # plt.ylim([-6., 3.])
            # plt.ylabel(r'log $\lambda I_{\lambda}$')       # label for y axis
            # #
            # text = ' iage = %3.3i ,age = %7.5f' %(ia, age)
            # plt.text(3.e-2, -0.5, text,size=14., bbox=dict(facecolor='white', edgecolor='none'))
            # plt.plot(original_waves, log10(original_SED[Z][ia]), '-', color='blue', lw=2.0)
            # plt.plot(waves, log10(new_SED[Z][ia]), '--', color='green', lw=1.0)
            # #
            # plt.show()
            # plt.close()
            #
            # Compute elementary Luminosity of a given population
            Lum = trapz(new_SED[Z][ia] / waves, waves)
            if isnan(Lum):
                msg = f'Luminosity for Z={Z} and iage={ia} is NAN'
                raise ValueError(msg)
            SSP_luminosity[Z].append(Lum)
    #
    # Create data file
    filename  = f'sp_SED_[BC03]_{imf}.in'
    # Reformat data
    ageBins = ' '.join([f'{a}' for a in ages])
    MetBinsStr = ' '.join([f'{z}' for z in MetBins])
    waveBins = ' '.join([f'{w}' for w in waves])
    tmp = ravel([[f'{l}' for l in lumz] for lumz in SSP_luminosity.values()])
    LumBins = ' '.join([f'{t}' for t in tmp])
    #
    D = []
    Ncols = 0
    for iz, Z in enumerate(MetBinsRef):
        Ncols += 1
        D.append([])
        for ia, _ in enumerate(ages):
            for l, _ in enumerate(waves):
                D[iz].append(new_SED[Z][ia][l])
    mMin = 0.8
    mMax = 100.0
    #
    hdr = ('G.A.S. input file\n'
           'In this file, each header line is starting by #\n'
           '   data are saved after the ---- limit\n'
           'Are listed here for a set of initial metallicities:\n'
           'SSP individual SED as a function of the SSP age:\n'
           '     - for an instantaneous bursts (SSP) of 1 Msun\n'
           f'     - following a {imf} IMF (Mmin = {mMin}, Mmax = {mMax} Msun)\n'
           'At a given age, the normalization is done for M*_liv + M_rem = 1 Msun\n'
           'The SED luminosity is given in solar unit (Lsun) with Lsun = 3.826 x 10^33 ergs/s\n'
           'This data file has build as follow:\n'
           'nAgeBins nMetBins nlambda\n'
           'List of stellar ages [Gyr]             : AgeBins\n'
           'List of metalicity bins [mass ratio]   : MetBins\n'
           'List of wavelenghts [microns]          : waves\n'
           'List of luminosity for each age bin    : LumBins(1:nAgeBins, 1:nMetBins) in Lsun\n'
           'A blank line\n'
           'And for each following data line (i: age index, j: lambda index): SB(i,j,1:nMetBins) in Lsun\n'
           '----\n'
           f'{len(ages)} {len(MetBinsRef)} {len(waves)}\n'
           f'{ageBins}\n'
           f'{MetBinsStr}\n'
           f'{waveBins}\n'
           f'{LumBins}\n')
    savetxt(path.join(retreated_data_path, filename), column_stack(D), header=hdr, fmt='%15.6e'*Ncols)
