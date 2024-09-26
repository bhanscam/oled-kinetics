# Example of Parameter Search for PL Fit

**(Using made up PL data)

## How to run
1. Update `kinetics.F` to reflect the kinetic equations in the desired model

2. Compile fortran library for python to call

	`f2py3 -c -m kinetics kinetics.F`

3. Update `input_pl_scan.py` with correct data files and rates included in the model

4. Make changes to `scan.py` based on number of rates in the model the defined parameter space

5. Build list of rates of inputs

	`python scan.py`

6. Update `submit_scan.sh` with the total number of jobs to run using the output of `scan.py`

7. Run input file (many many times)

	`sbatch submit_scan.sh`

8. Analyze the run logs of the scan once complete. Output is a list of optimized rates and their associated rmse value.

	`bash analyze_scan.sh N_nodes out_file_name`

