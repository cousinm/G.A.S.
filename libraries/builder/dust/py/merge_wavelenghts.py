#!/usr/bin/env python
# reformat dustem data for eGalICS input
#
# import library
from numpy import loadtxt, savetxt
from os import path
from scipy.interpolate import InterpolatedUnivariateSpline as interp # for interpolation

# Path to original data
original_data_path = 'dust/data/originals'
#
# Path to retreated data
retreated_data_path = 'dust/data/retreated'
#
# Read original DUSTEM wavelenghts
filename = 'LAMBDA.DAT'
ir_waves = loadtxt(path.join(original_data_path, filename),
                   unpack=True, usecols=[0])  # micron
#
min_ir_waves = min(ir_waves)
max_ir_waves = max(ir_waves)
print('> (N, min_ir_waves, max_ir_waves) = '
      f'({len(ir_waves)}, {min_ir_waves}, {max_ir_waves})')
#
# Read original BC03 wavelenghts
filename = 'BC03_wavelenghts.dat'
uv_opt_nir_waves = loadtxt(path.join(original_data_path, filename),
                           unpack=True,usecols=[0])  # micron
min_uv_opt_nir_waves = min(uv_opt_nir_waves)
max_uv_opt_nir_waves = max(uv_opt_nir_waves)
print ('> (N, min_uv_opt_nir_waves, max_uv_opt_nir_waves) = '
       f'({len(uv_opt_nir_waves)}, {min_uv_opt_nir_waves}, {max_uv_opt_nir_waves})')
#
# Until 10 microns BC03 wavelenghts are used
BC03_wave_limit = 10. # microns
lbc03 = 0
waves = []
BC03_indexes_in_waves = []
while (uv_opt_nir_waves[lbc03] <= BC03_wave_limit):
    waves.append(uv_opt_nir_waves[lbc03])
    BC03_indexes_in_waves.append(lbc03)
    lbc03 = lbc03 + 1
#
# Complete full wavelenght table
l   = lbc03
lir = 0
# Looking for the first ir wave next BC03_wave_limit
while (ir_waves[lir] < BC03_wave_limit):
    lir = lir + 1
#
while ((lir < len(ir_waves)) and (lbc03 < len(uv_opt_nir_waves))):
        wave_bc03 = uv_opt_nir_waves[lbc03]
        wave_ir   = ir_waves[lir]
        if (wave_ir < wave_bc03):
            # Use ir_wave
            waves.append(wave_ir)
            l = l + 1
            lir = lir + 1
        else:
            # The next BC03 wavelenght is too close
            # we must take into account BC03 wavelenght
            waves.append(wave_bc03)
            BC03_indexes_in_waves.append(lbc03)
            lbc03 = lbc03 + 1
            l = l + 1
            dw = wave_ir - waves[-1]
            if dw < (ir_waves[lir + 1] - wave_ir)/2.:
                lir = lir + 1
#
# Save full wavelenghts
filename = 'wavelenghts.in'
hdr = 'Complete wavelenght table [micron]\n'
hdr += f'{len(waves)}'
savetxt(path.join(retreated_data_path, filename), waves, header=hdr, fmt='%15.6e')
#
# save BC03_in_waves
filename = 'BC03_indexes_in_waves.dat'
hdr = 'Indexes of the original BC03 wavelengt table included in the full wavelentgh table\n'
hdr += f'{len(BC03_indexes_in_waves)}'
savetxt(path.join(retreated_data_path, filename), BC03_indexes_in_waves,
        header=hdr, fmt='%4i')
