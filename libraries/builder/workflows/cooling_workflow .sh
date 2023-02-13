#!/bin/bash
echo '> Configure virtual env'
pipenv install
pipenv shell

echo '> Erase previous plots and retreated data'
rm ./cooling/plots/*
rm ./cooling/data/retreated/*

echo '> Reformat original data'
echo '  > reformat_Rijcke_cooling_curves'
python ./cooling/py/reformat_Rijcke_cooling_curves.py

echo '> Compute cooling efficiency and thermal instability factor'
echo '  > compute_and_format_cooling_tables'
python ./cooling/py/compute_and_format_cooling_tables.py