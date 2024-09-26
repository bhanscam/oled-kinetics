from kinetics_fit import *
import sys

# Example input file to generate PL curves for a set of three OLEDs

'''
Possible Reactions
a)           -> ex       K_abs  : Exciton absorption
b)        ex ->          K_rad  : Exciton Radiative Decay
c)        ex ->          K_ET   : Exciton Nonradiative Decay
d)        dx ->          K_DT   : Dark Exciton Nonradiative Decay
e)        ch ->          K_CT   : Charge Leakage
f)        ex -> dx       K_BD   : Exciton conversion, bright to dark
g)        dx -> ex       K_DB   : Exciton conversion, dark to bright
h)        ex -> el + ho  K_sepB : Charge Separation, bright
i)        dx -> el + ho  K_sepD : Charge Separation, dark
j)   el + ho -> ex       K_recB : Charge Recombination, bright
k)   el + ho -> dx       K_recD : Charge Recombination, dark
l)   ex + ex -> ex       K_TTA  : Triplet-Triplet Annihilation
     ex + dx -> ex
     ex + dx -> dx
     dx + dx -> dx
m)   ex + ch -> ch       K_TCA  : Triplet-Charge Annihilation
     dx + ch -> ch
'''

# create list of input parameters based on scan info
node = int(sys.argv[1])
core = int(sys.argv[2])
total_cores = int(sys.argv[3])
job = ( ( node - 1 ) * total_cores ) + core
scan_params = np.loadtxt("scan_params.txt")
params_list = scan_params[int(job-1)]

# list experimental data files, in same order as device_rates lists
# Note: input file should only contain the relevant section of data
#       and can be normalized or raw
devices_exp = ["exciplex_ag","exciplex_agal","exciplex_al"]

# Build rate dictionary of PL rates
# Set any rates unused in the model to zero
# Set any device specific rates to zero
PLrates = {
	"a_Kabs"  : 0.0,
	"b_Krad"  : 0.0,
	"c_K_ET"  : 0.0,
	"d_K_DT " : 0.0,
	"e_K_CT " : 0.0,
	"f_K_BD " : 0.0,
	"g_K_DB " : 0.0,
	"h_KsepB" : params_list[0],
	"i_KsepD" : 0.0,
	"j_KrecB" : params_list[1],
	"k_KrecD" : 0.0,
	"l_K_TTA" : 0.0, 
	"m_K_TCA" : 0.0,
}
# Host specific PL rates
exciplex_PLrates = {
}

simple_PLrates = {
}

# Coupling layer specific rates
Ag_PLrates = {
	"c_K_ET"  : params_list[2],
}
     
AgAl_PLrates = {
	"c_K_ET"  : params_list[3],
}

Al_PLrates = { 
	"c_K_ET"  : params_list[4],
}

layer_PLrates = [Ag_PLrates, AgAl_PLrates, Al_PLrates]
host_PLrates = [exciplex_PLrates]

# initial relative concentrations of each species
Concs = {
	"ex": 0.005,
	"ch": 0.0,
	"dx": 0.005,
}

# Build OLED object and run transient-PL fit/analysis
''' fit=pl         : read in experimental times and [ex] from input file
    optimize=False : generates model PL curve for given set of rate parameters 
    optimize=True  : runs a Nelder-Meade optimization over the given rate
                     parameters to fit to the provided experimental PL data
    plot=True      : prints to file model PL data and times for plotting 
'''
my_pl_model = oled( Concs, devices_exp, fit='pl' )
my_pl_model.pl_fit( PLrates, host_PLrates, layer_PLrates, optimize=True, plot=False )
