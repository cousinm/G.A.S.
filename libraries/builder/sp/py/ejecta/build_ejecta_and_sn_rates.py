#!/usr/bin/env python
# Compute ejecta rates associated to a single stellar population
#
# import library
import matplotlib.pyplot as plt         # plot
import matplotlib.gridspec as gridspec  # for multiplot grid
#
from os import path, environ
from tools.IMF import IMF_List
from tools.ejecta import MetBins, MetBinsRef, Elements, mMin, mMax, compute_average
from numpy import savetxt, column_stack, array
#
# Define pathes
path_BC03 = environ.get('bc03', None)
if path_BC03 is None:
    raise Exception("Please set the 'bc03' environnement variable")
#
retreated_data_path = 'sp/data/retreated'
#
# Load original BC03 filtered timescales
filename = 'filtered_original_BC03_timescale_Gyr.dat'
with open(path.join(path_BC03, filename), 'r') as file:
    lines = file.readlines()
    timescale = [float(t) for t in lines[2].split(',')]
# Add a additionnal step at 15.25 Gyr to allow a computation until 15 Gyr (included)
timescale.append(15.25)
#
for imfRef in IMF_List:
    print('> IMF: ' + imfRef)
    #
    # Create data table
    # Full
    ejctRate = {}
    ejctRateSNIa = {}
    #
    # BC03 average age bins
    BC03_ejctRate = {}               # Full ejecta
    BC03_ejctRateSNIa = {}           # Ejection + metal creation from SNIa
    BC03_living_remnant_mass = {}    # living + remnant stellar mass
    BC03_SNR = {}                    # Average SNIa + SNII rates
    #
    for Z in MetBinsRef:
        print(f"  > Z = 0.{Z}")
        ejctRate[Z] = {}
        ejctRateSNIa[Z] = {}
        BC03_ejctRate[Z] = {}
        BC03_ejctRateSNIa[Z] = {}
        BC03_living_remnant_mass[Z] = []
        BC03_SNR[Z] = []
        #
        for elt in Elements:
            print(f"    > elt = {elt}")
            ejctRate[Z][elt] = []
            ejctRateSNIa[Z][elt] = []
            BC03_ejctRate[Z][elt] = []
            BC03_ejctRateSNIa[Z][elt] = []
            #
            fullTimescale = []
            #
            for it, t in enumerate(timescale[:-1]):
                #
                t0 = t
                t1 = timescale[it + 1]
                ts, av_ejct, av_ejctSNIa, av_Mej, av_SNR, ejr, ejrSNIa = compute_average(t0, t1, Z, elt, imfRef)
                #
                fullTimescale += ts[:-1]
                ejctRate[Z][elt] += ejr[:-1]
                ejctRateSNIa[Z][elt] += ejrSNIa[:-1]
                BC03_ejctRate[Z][elt].append(av_ejct)
                BC03_ejctRateSNIa[Z][elt].append(av_ejctSNIa)
                #
                if elt in ['Total']:
                    BC03_SNR[Z].append(av_SNR)
                    if it == 0:
                        Mrem = 1. - av_Mej
                    else:
                        Mrem = BC03_living_remnant_mass[Z][it - 1] - av_Mej
                    BC03_living_remnant_mass[Z].append(Mrem)
    #
    # Save data
    #
    # Ordered elements
    orderedElments = ['Total', 'Metals', 'H1', 'He4', 'C12', 'N14', 'O16', 'Fe56']
    #
    # Living and remnant stellar mass
    filename = path.join(retreated_data_path, f'BC03_living_remnant_mass_{imfRef}.dat')
    D =[timescale[:-1]]
    for Z in MetBinsRef:
        D.append(BC03_living_remnant_mass[Z])
    hdr = 'Living + remnant stellar mass as function of the BC03 selected timesteps\n'
    hdr += '  computed for a set of metallicity\n'
    hdr += f' and according to {imfRef} IMF\n'
    hdr += '1) TimeStep, 2) Z=0.0000, 3) Z=0.0001, 4) Z=0.0004, 5) Z=0.004, 6) Z=0.008, 7) Z=0.02, 8) Z=0.05'
    savetxt(filename, column_stack(D), header=hdr, fmt='%17.8e'*8)
    #
    # Average ejecta rate as function of the SSP age for the different elements and metalicity
    filename = path.join(retreated_data_path, f'sp_MLR_[BC03]_{imfRef}.in')
    MetBinsStr = ' '.join([f'{z}' for z in MetBins])
    hdr = ('G.A.S. input file\n'
        'In this file, each header line is starting by #\n'
        '   data are saved after the ---- limit\n'
        'Are listed here, for a set of initial metallicities:\n'
        'Total, Metals(Z), H1, He4, C12, N14, O16 and Fe56 average mass loss rates (mlr) [Msun/Gyr]\n'
        '   as a function of the SSP age\n'
        '     - for an instantaneous bursts (SSP) of 1 Msun\n'
        f'     - following a {imfRef} IMF (Mmin = {mMin}, Mmax = {mMax} Msun)\n'
        'At a given age, the normalization is done for M*_liv + M_rem = 1 Msun\n'
        'This data file has build as follow:\n'
        'nAgeBins, nMetBins and nElts\n'
        'List of metalicity bins: MetBins\n'
        'A blank line\n'
        'And for each following data line : AgeBin(i), Total_mlr(i,1:nMetBins), Z_mlr(i,1:nMetBins), H1_mlr(i,1:nMetBins), He4_mlr(i,1:nMetBins) ...\n'
        '  ... C12_mlr(i,1:nMetBins), N14_mlr(i,1:nMetBins), O16_mlr(i,1:nMetBins), Fe56_mlr(i,1:nMetBins)\n'
        '----\n'
        f'{len(timescale[:-1])} {len(MetBinsRef)} {len(orderedElments) -2}\n'
        f'{MetBinsStr}\n'
        )
    D = []
    D.append(timescale[:-1])
    Ncols = 1
    for elt in orderedElments:
        for Z in MetBinsRef:
            Ncols += 1
            # Normalisation to the living + remnant current stellar mass
            D.append(array(BC03_ejctRate[Z][elt])/array(BC03_living_remnant_mass[Z]))
    # Write the file
    savetxt(filename, column_stack(D), header=hdr, fmt='%17.8e'*Ncols)
    #
    # Save average SN rates (SNIa + SNII)
    filename = path.join(retreated_data_path, f'sp_SNR_[BC03]_{imfRef}.in')
    hdr = ('G.A.S. input file\n'
        'In this file, each header line is starting by #\n'
        '   data are saved after the ---- limit\n'
        'Are listed here for a set of initial metallicities:\n'
        'SN Rates (SNIa + SNII) [nb/Gyr] as a function of the SSP age:\n'
        '     - for an instantaneous bursts (SSP) of 1 Msun\n'
        f'     - following a {imfRef} IMF (Mmin = {mMin}, Mmax = {mMax} Msun)\n'
        'At a given age, the normalization is done for M*_liv + M_rem = 1 Msun\n'
        'This data file has build as follow:\n'
        'nAgeBins, nMetBins\n'
        'List of metalicity bins: MetBins\n'
        'A blank line\n'
        'And for each following data line : AgeBin(i), SNR(i,1:nMetBins)\n'
        '----\n'
        f'{len(timescale[:-1])} {len(MetBinsRef)}\n'
        f'{MetBinsStr}\n'
        )
    D = []
    D.append(timescale[:-1])
    Ncols = 1
    for elt in orderedElments:
        for Z in MetBinsRef:
            Ncols += 1
            # Normalisation to the living + remnant current stellar mass
            D.append(array(BC03_SNR[Z])/array(BC03_living_remnant_mass[Z]))
    # Write the file
    savetxt(filename, column_stack(D), header=hdr, fmt='%17.8e'*Ncols)
    #
    # Plot
    plot_path = 'sp/plots'
    plt.figure(figsize = (16, 9))
    gs = gridspec.GridSpec(2, 4)
    gs.update(wspace=0.08, hspace=0.05,
            top=0.95, bottom=0.05,
            left=0.05, right=0.95)
    # Legend
    # tot  : Black
    # Z    : DimGray
    # H1   : Crimson
    # He4  : Purple
    # C12  : SlateGray
    # N14  : YellowGreen
    # O16  : DodgerBlue
    # Fe56 : SaddleBrown
    color = ['Crimson', 'Purple', 'SlateGray',
            'YellowGreen', 'DodgerBlue',
            'SaddleBrown','DimGray','Black']
    # Legends
    Elements_Tex = ['$^{1}H$', '$^{4}He$', '$^{12}C$',
                    '$^{14}N$', '$^{16}O$', '$^{56}Fe$',
                    '$Metals$','$Total$']
    #
    for iz, Z in enumerate(MetBinsRef):
        ax = plt.subplot(gs[iz])
        for ielt, elt in enumerate(Elements):
            ax.plot(timescale[:-1], BC03_ejctRate[Z][elt], ls='-', color=color[ielt])
        #
        ax.set_xlabel('Times ($Gyr$)')                                # label for x axis
        ax.set_ylabel('Instantaneous ejecta rate ($M_{\odot}/Gyr$)')  # label for y axis
        # ax.set_xlim([0.8, 100.])
        ax.set_ylim([1.e-6, 30.])
        ax.set_yscale('log')           # use log scale for y axis
        ax.set_xscale('log')           # use log scale for y axis
        ax.grid(True,ls="-")           # display a grid with all sub log tics, which="both"
        # turn on/off x and/or y axis tics and labels
        if iz < 4:
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position('top')
            # ax.set_ylim([1.e-6, 100.])
        if (iz == 0) or (iz == 4):
            ax.yaxis.tick_left()
            ax.yaxis.set_label_position('left')
        if (iz == 3) or (iz == 6):
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position('right')
        if (iz == 1) or (iz == 2) or (iz==5):
            ax.set_yticklabels([])
            ax.set_ylabel('')
        if Z in ['05']:
            plt.legend(ncol=3, loc='best', prop={'size':10})
        plt.text(0.7, 0.93, f'Z = 0.{Z}', size=12.,
                transform=ax.transAxes,
                bbox=dict(facecolor='white', edgecolor='white'))
    #
    # Save the plot in pdf format
    filename = path.join(plot_path, f'selected_Final_Instantaneous_Ejecta_Rates_{imfRef}.pdf')
    plt.savefig(filename)
    plt.close()
    #
    # Plot evolution of the remnant mass
    plt.figure(figsize=(5.5, 5))
    # color reference
    color = ['LightBlue', 'SkyBlue', 'LightSeaGreen', 'YellowGreen',
            'Gold', 'DarkOrange', 'FireBrick']
    #
    for iz, Z in enumerate(MetBinsRef):
        legend = f'Z = 0.{Z}'
        plt.plot(timescale[:-1], BC03_living_remnant_mass[Z], c=color[iz], label=legend)
    #
    # plt.ylim([0.4, 1.05])                              # yrange
    plt.xscale('log')                                    # use log scale for x axis
    plt.xlabel(r'SSP Age [Gyr]')                         # label for x axis
    plt.ylabel(r'Living + remnant mass ($M_{\odot}$)')   # label for y axis
    plt.grid(True, which="both", ls="-")                 # display a grid
    plt.legend(loc='best', prop={'size': 10})            # display legend
    #
    filename = path.join(plot_path, f'BC03_living_remnant_mass_{imfRef}.pdf')
    plt.savefig(filename)
    plt.close()
    #
    # Plot evolution of the SN Rates (SNIa + SNII)
    plt.figure(figsize=(5.5, 5))
    #
    for iz, Z in enumerate(MetBinsRef):
        legend = f'Z = 0.{Z}'
        plt.plot(timescale[:-1], BC03_SNR[Z], c=color[iz], label=legend)
    #
    #plt.ylim([0.4, 1.05])                         # yrange
    plt.xscale('log')                              # use log scale for x axis
    plt.yscale('log')                              # use log scale for x axis
    plt.xlabel(r'SSP Age [Gyr]')                   # label for x axis
    plt.ylabel(r'SN Rates (SNIa + SNII) [/Gyr]')   # label for y axis
    plt.grid(True, which="both", ls="-")           # display a grid
    plt.legend(loc='best', prop={'size': 10})      # display legend
    filename = path.join(plot_path, f'BC03_SNR_{imfRef}.pdf')
    plt.savefig(filename)
    plt.close()
