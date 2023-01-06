from os import  path
from numpy import savetxt
#
path_spectra = '/home/morgane/Numerique/bc03/spectra'
#
# Load original spectra
# Define base filename
filename = f'Chabrier+03/bc2003_lr_BaSeL_Z05_ssp.spec'
#
# Wavelenghts
wavelenghts = []
with open(path.join(path_spectra, filename)) as file:
    lines = file.readlines()
    #
    for il, line in enumerate(lines):
        if il > 5:
            sline = line.split()
            wavelenghts.append(float(sline[0])/1.e4)   # wavelenght [microns]
#
# Save
filename = path.join(path_spectra, 'BC03_wavelenghts.dat')
hdr = 'BC03 wavelenght [micron]\n'
hdr += f'{len(wavelenghts)}'
savetxt(filename, wavelenghts, header=hdr, fmt='%17.8e'*1)
