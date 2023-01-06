#!/bin/bash
echo '> Configure virtual env'
pipenv install
pipenv shell

echo '> Erase previous plots and retreated data'
rm ./sp/plots/*
rm ./sp/data/retreated/*

echo '> Create initial abundance'
echo '  > ext_int_element_abundance'
python ./sp/py/ejecta/ext_int_element_abundance.py

echo '> Extrapolate/Interpolate final/remnant MS'
echo '  > ext_int_MS_initial_final_mass'
python ./sp/py/ejecta/ext_int_MS_initial_final_mass.py
echo '  > ext_int_MS_initial_remnant_mass'
python ./sp/py/ejecta/ext_int_MS_initial_remnant_mass.py

echo '> Extrapolate/Interpolate final/remnant LIMS'
echo '  > ext_int_LIMS_initial_remnant_mass'
python ./sp/py/ejecta/ext_int_LIMS_initial_remnant_mass.py

echo '> Reformat original data'
echo '  > reformat_Campbell_LIMS_ejecta'
python ./sp/py/ejecta/reformat_Campbell_LIMS_ejecta.py
echo '  > reformat_Gil-Pons_LIMS_ejecta'
python ./sp/py/ejecta/reformat_Gil-Pons_LIMS_ejecta.py
echo '  > reformat_Heger_final_ejecta'
python ./sp/py/ejecta/reformat_Heger_final_ejecta.py
echo '> Build LIMS ejecta Z=0'
echo '  > build_LIMS_final_ejecta_Z0000'
python ./sp/py/ejecta/build_LIMS_final_ejecta_Z0000.py

echo '> Extrapolate/Interpolate ejecta LIMS'
echo '  > ext_int_LIMS_final_ejecta'
python ./sp/py/ejecta/ext_int_LIMS_final_ejecta.py

echo '> Extrapolate/Interpolate ejecta MS'
echo '  > ext_int_MS_final_ejecta'
python ./sp/py/ejecta/ext_int_MS_final_ejecta.py

echo '> Build final ejecta'
echo '  > build_complete_ejecta'
python ./sp/py/ejecta/build_complete_ejecta.py

echo '> Buid stellar lifetimes'
echo '  > ext_int_stellar_lifetimes'
python ./sp/py/ejecta/ext_int_stellar_lifetimes.py
