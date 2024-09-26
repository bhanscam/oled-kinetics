Example of Parameter Search for Roll-Off Fit

**(Using made up roll-off data)**

## How to run
1. Update `kinetics.F` to reflect the kinetic equations in the desired model

2. Compile fortran library for python to call

	f2py3 -c -m kinetics kinetics.F

3. Update `input_rolloff.py` with correct data files and rates included in the model

4. Run input file

	python input_rolloff.py
