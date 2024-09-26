from kinetics_fit import *

# Example input file to generate PL and Roll-off curves for a set
# of three OLEDs

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

# list experimental data files, in same order as device_rates lists
# Note: input file should only contain the relevant section of data
#       and can be normalized or raw
devices_exp = ["exciplex_ag","exciplex_agal","exciplex_al"]

# Build rate dictionary of fixed PL rates
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
	"h_KsepB" : 1.33E+06,
	"i_KsepD" : 0.0,
	"j_KrecB" : 2.43E+09,
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
	"c_K_ET"  : 1.35e+06,
}
     
AgAl_PLrates = {
	"c_K_ET"  : 1.06e+06,
}

Al_PLrates = { 
	"c_K_ET"  : 7.00e+05,
}

layer_PLrates = [Ag_PLrates, AgAl_PLrates, Al_PLrates]
host_PLrates = [exciplex_PLrates]

# Build dictionary for rates optimized with roll-off
# Provide initial guesses for relative magnitude
ELrates = {
}

# Host specific EL rates
exciplex_ELrates = {
	"l_K_TTA" : 6.04E+10, 
	"m_K_TCA" : 5.90E+08,
}

simple_ELrates = {
}

# Coupling layer specific rates
Ag_ELrates = {
	"e_K_CT " : 7.12E+08,
}
     
AgAl_ELrates = {
	"e_K_CT " : 7.59E+08,
}

Al_ELrates = {
	"e_K_CT " : 1.20E+09,
}

layer_ELrates = [Ag_ELrates, AgAl_ELrates, Al_ELrates]
host_ELrates = [exciplex_ELrates]

# initial relative concentrations of each species
Concs = {
	"ex": 0.005,
	"ch": 0.0,
	"dx": 0.005,
}

# Build OLED object and run roll-off fit/analysis
''' fit=el           : read in experimental J (current density) and EQE values
                       from input file
    optimize=False   : generates model roll-off curve for given set of rate parameters 
    optimize=True    : runs a Nelder-Meade optimization over the given rate
                       parameters to fit to the provided experimental roll-off data
    plot=True        : prints to file model roll-off data and J values for plotting
'''
my_el_model = oled( Concs, devices_exp, fit='el' )
my_el_model.rolloff_fit( PLrates, host_PLrates, layer_PLrates, ELrates, host_ELrates, layer_ELrates, optimize=False, plot=True )
