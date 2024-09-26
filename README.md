# OLEDs Kinetics Model with Transient-PL and Roll-Off Fits

Developed by Rebecca Hanscam

## Overview
This code builds flexible kinetic models to simulate or fit to photoluminescence and electroluminescence experiments for OLEDs. Using the kinetic equations involving reactions between electron-hole pair charges and bright and dark triplets, this code locates optimal sets of rate constants for a given model using these species. This code is designed to facilitate searches over thousands of different rate parameters for many different models.

## Dependencies:
- Python
- Scipy
- Numpy
- f2py3

## Code structure
- `input_analyze.py` : define desired model, controls optimization variables, parameters
- `kinetics_fit.py` : builds oled class, runs optimizer, wrapper for steady-state code
- `kinetics.F` : runs steady-state for kinetic equations

## Fit to Experimental Data
Provide experimental data in text files labeled `devicename_pl.txt` or `devicename_el.txt`. The PL data should be given in two columns, one for time (in microseconds) and one for PL (in counts). The roll-off data should be given in two columns, one for current density (mA/cm^2) and one for EQE. If only a prediction of the PL or roll-off curves is desired, provide similar data files including the time range (for PL) or current density range (for roll-off) with place-holder PL or EQE data.

## How to run
1. Update `input_example.py` with correct data files and initial rates

2. Update `kinetics.F` to reflect the kinetic equations in the desired model

3. Compile fortran library for python to call

	f2py3 -c -m kinetics kinetics.F

4. Run input file

	python input_analyze.py

See `examples` directory for example input files and basic usage options:
- PL analysis or fit
- Roll-off analysis or fit
- Scan of PL rate parameters
Note: All examples use made up experimental PL and roll-off data.
