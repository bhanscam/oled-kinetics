import numpy as np
import kinetics
from scipy.optimize import minimize
from collections import OrderedDict
import math

np.set_printoptions(precision=10)

def sn( num ):
	'''
	Print a given number in scientific notation
	'''
	return "{:.3e}".format( num )

def read_exp_files( devices_exp, fit ):
	'''
	Read experimental data files into arrays
	Input: list of file names, without "_pl.txt" or "_el.txt" extension
	'''
	if fit == 'pl':
		kind = "_pl.txt"
	if fit == 'el':
		kind = "_el.txt"

	# if PL then Jvals=time and eqe_exp=exciton concentration
	Jvals = []
	eqe_exp = []
	for filename in devices_exp:
		Jtemp = []
		EQEtemp = []
		with open(filename+kind, 'r') as file:
			for line in file:
				line = line.strip()
				words = line.split()
				Jtemp.append(float(words[0]))
				EQEtemp.append(float(words[1]))
		Jvals.append(Jtemp)
		eqe_exp.append(EQEtemp)
	
	return Jvals, eqe_exp

class oled:

	def __init__( self, Concs_dict, device_names, fit ):

		# add arrays to oled object, convert things to lists for fortran
		self.Concs = [Concs_dict["ex"],Concs_dict["ch"],Concs_dict["dx"]]

		self.device_names = device_names
		self.Nspecies = len(self.Concs)

		# use device names to make list of device layers and their indices
		layer_names = [device.split('_')[1] for device in device_names]
		host_names = [device.split('_')[0] for device in device_names]
		self.layer_names = layer_names
		self.host_names = host_names

		layer_device_indices = []
		for layer in layer_names:
			if layer == 'ag':
				layer_device_indices.append(0)
			elif layer == 'agal':
				layer_device_indices.append(1)
			elif layer == 'al':
				layer_device_indices.append(2)
			else:
				print("coupling layer label not recognized")
				exit()
		self.layer_device_indices = layer_device_indices

		# organize device info
		host_device_indices = []
		rec_sep_ratios = []	# use if constraining k_rec/k_sep ratio
		for host in host_names:
			if host == 'exciplex':
				host_device_indices.append(0)
				rec_sep_ratios.append(5.16e+01)
			elif 'exciplex' in host_names and host == 'simple':
				host_device_indices.append(1)
				rec_sep_ratios.append(7.79e+02)
			elif 'exciplex' not in host_names and host == 'simple':
				host_device_indices.append(0)
				rec_sep_ratios.append(7.79e+02)
			else:
				print("host label not recognized")
				exit()
		self.host_device_indices = host_device_indices
		self.rec_sep_ratios = rec_sep_ratios

		# read in experimental data files
		if fit == 'el':
			Jvals, eqe_exp = read_exp_files( device_names, fit='el' )
			for device in range(len(Jvals)):
	
				# normalize exp data as percentage loss
				EQEdev = eqe_exp[device]
				EQEdev = [(float(i)/max(EQEdev))*100 for i in EQEdev]
				eqe_exp[device] = np.array(EQEdev)
			self.eqe_exp = eqe_exp
			self.Jvals = Jvals

		if fit == 'pl':
			times, pl_exp = read_exp_files( device_names, fit='pl' )
			for device in range(len(times)):
				# convert to microseconds
				times_dev = times[device]
				times_dev = [t*1e-6 for t in times_dev]
				times[device] = times_dev

			self.times = times
			self.pl_exp = pl_exp

		return

	def pl_rmse( self, eqe_model, device_ind, plot=False ):
		'''
		Calculate RMSE between model eqe and experimental eqe
		'''
		# normalize both exp curve and model
		# shift eqe model data so there are no negative numbers
		if min(eqe_model) < 0:
			eqe_model = [(float(i)-min(eqe_model))/(eqe_model[0]-min(eqe_model)) for i in eqe_model]
		else:
			eqe_model = [float(i)/eqe_model[0] for i in eqe_model]
		eqe_model = np.array(eqe_model)

		# calculate the rmse between exp and model
		rmse = np.sqrt(((eqe_model - self.pl_exp[device_ind]) ** 2).mean())

		if math.isnan( rmse ):
			print( "species populations got too big to handle" )
			exit()

		if plot:
			np.savetxt(self.device_names[device_ind]+"_pl_model.txt",eqe_model)
			np.savetxt(self.device_names[device_ind]+"_pl_exp.txt",self.pl_exp[device_ind])
			np.savetxt(self.device_names[device_ind]+"_pl_times.txt",self.times[device_ind])

			print(self.device_names[device_ind]," rmse is ", rmse)

		return rmse


	def pl_objective( self, ELrates_list, plot=False ):
		'''
		Objective function for Nelder-Mead that runs steady state to calculate eqe and compare to experimental data
		input: list of optimizable rates {PL1,PL2,PL3,...,EL1,EL2,...,ag-EL1,ag-EL2,...,agal-EL1,agal-EL2,...,al-EL1,al-EL2,...}
		output: rmse between model roll-off curve and experimental
		'''

		rmse = 0

		# enforce that the rates be positive numbers because Nelder-Mead can't handle bounds on the parameters
		if any(rate < 0 for rate in ELrates_list):
			rmse = 1e3
		else:
			#for device_ind in range(len(self.device_ELrates_dicts)):
			for device_ind in range(len(self.device_names)):

				# update EL rates in full rate dictionary
				i = 0
				for key in self.ELkeys:
					self.rates_dict[key] = ELrates_list[i]
					i += 1
				j = 0
				for key in self.device_ELkeys:
					self.rates_dict[key] = ELrates_list[i+(self.host_device_indices[device_ind]*len(self.device_ELkeys))+j]
					j += 1
				k = 0
				for key in self.layer_ELkeys:
					self.rates_dict[key] = ELrates_list[i+(len(self.device_ELkeys)*len(self.device_ELrates_dicts))+(self.layer_device_indices[device_ind]*len(self.layer_ELkeys))+k]
					k += 1
				self.rates_dict = OrderedDict(sorted(self.rates_dict.items()))
				rates_list = [val for key,val in self.rates_dict.items()]

				# fortran is dumb so embed rates and times in lists of determined length
				fortran_rates = [0*i for i in range(30)]
				for i in range(len(rates_list)):
					fortran_rates[i] = rates_list[i]
				Ntimes = len(self.times[device_ind])
				fortran_times = [0*i for i in range(1000)]
				for i in range(Ntimes):
					fortran_times[i] = self.times[device_ind][i]

				# run steady-state for given rates and calculate rmse
				pl = kinetics.steady_state_pl( Ntimes, fortran_times, self.Concs, fortran_rates )
				pl = pl[:Ntimes]
				rmse += self.pl_rmse( pl, device_ind, plot=plot )

		# print rmse and optimizable rates
		self.rmse_final = rmse
		if rmse != 1e3:
			print_rmse = "rmse = " + str(rmse)
			i = 0
			for key in self.ELkeys:
				print_rmse += " " + key + ":" + str(sn(ELrates_list[i]))
				i += 1
			for device_ind in range(len(self.times)):
				for j in range(len(self.device_ELkeys)):
					print_rmse += " " + self.device_names[device_ind] + ":" + self.device_ELkeys[j] + ": " + str(sn(ELrates_list[i+(self.host_device_indices[device_ind]*len(self.device_ELkeys))+j]))
				for k in range(len(self.layer_ELkeys)):
					print_rmse += " " + self.device_names[device_ind] + ":" + self.layer_ELkeys[k] + ": " + str(sn(ELrates_list[i+(len(self.device_ELkeys)*len(self.device_ELrates_dicts))+(self.layer_device_indices[device_ind]*len(self.layer_ELkeys))+k]))
			print(print_rmse)
		return rmse

	def pl_fit( self, ELrates_dict, device_ELrates_dicts, layer_ELrates_dicts, optimize=False, plot=False ):
		'''
		Wrapper for PL rate optimization
		Note: here "EL" labels are used to indicate optimizable rates
		'''
		# alphabetize and store each optimizable dictionary in oled object
		self.ELrates_dict = OrderedDict(sorted(ELrates_dict.items()))
		for dev_ind in range(len(device_ELrates_dicts)):
			device_ELrates_dicts[dev_ind] = OrderedDict(sorted(device_ELrates_dicts[dev_ind].items()))
		self.device_ELrates_dicts = device_ELrates_dicts
		for lay_ind in range(len(layer_ELrates_dicts)):
			layer_ELrates_dicts[lay_ind] = OrderedDict(sorted(layer_ELrates_dicts[lay_ind].items()))
		self.layer_ELrates_dicts = layer_ELrates_dicts

		# build a dictionary of all the nonzero rates and sort alphabetically
		rates_dict = ELrates_dict.copy()
		rates_dict = {key:val for key,val in rates_dict.items() if val!=0 }
		rates_dict = OrderedDict(sorted(rates_dict.items()))
		self.ELkeys = list(rates_dict.keys())
		rates_dict.update(device_ELrates_dicts[0])
		rates_dict.update(layer_ELrates_dicts[0])
		rates_dict = {key:val for key,val in rates_dict.items() if val!=0 }
		self.rates_dict = OrderedDict(sorted(rates_dict.items()))
		self.device_ELkeys = list(self.device_ELrates_dicts[0].keys())
		self.layer_ELkeys = list(self.layer_ELrates_dicts[0].keys())
		print("full dictionary of rates included in the model: ",self.rates_dict)
		print("optimizable rates: ",self.ELkeys)
		print("host specific optimizable rates: ",self.device_ELkeys)
		print("layer specific optimizable rates: ",self.layer_ELkeys)
		print()

		# turn dictionary of optimizable rates into list for optimizer
		# [r1, r2, ..., Ag_r1, Ag_r2, ..., AgAl_r1, AgAl_r2, ..., Al_r1, Al_r2, ...]
		ELrates_dict = OrderedDict(sorted(ELrates_dict.items()))
		ELrates_list = [val for key,val in ELrates_dict.items() if val!=0]
		for device_dict in device_ELrates_dicts:
			device_dict = OrderedDict(sorted(device_dict.items()))	# sort device specific rates
			[ELrates_list.append(val) for key,val in device_dict.items()]
		for layer_dict in layer_ELrates_dicts:
			layer_dict = OrderedDict(sorted(layer_dict.items()))	# sort layer specific rates
			[ELrates_list.append(val) for key,val in layer_dict.items()]

		if optimize:
			# optimize rates
			result = minimize( self.pl_objective, ELrates_list, method='Nelder-Mead', options={'disp': True, 'maxiter': 1e4} )
			print("Nelder-Mead done!\nfinal optimized rates = ",result.x) 
			self.ELrates = result.x		

			## run steady-state on final conditions
			final_rmse = self.pl_objective( result.x, plot=plot )
		else:	
			## run steady-state on final conditions
			final_rmse = self.pl_objective( ELrates_list, plot=plot )

		return final_rmse


	def rolloff_rmse( self, eqe_model, device_ind, plot=False ):
		'''
		Calculate RMSE between model eqe and experimental eqe
		'''
		# normalize and rescale eqe model for direct comparison to exp data
		eqe_model_fit = [(float(i)/max(eqe_model))*100 for i in eqe_model]
		eqe_model_fit = np.array(eqe_model_fit)

		# calculate the rmse between exp and model
		rmse = np.sqrt((abs(eqe_model_fit - self.eqe_exp[device_ind]) ** 2).mean())

		if math.isnan( rmse ):
			print( "species populations got too big to handle" )
			exit()

		if plot:
			np.savetxt(self.device_names[device_ind]+"_eqe_model.txt",eqe_model)
			np.savetxt(self.device_names[device_ind]+"_eqe_exp.txt",self.eqe_exp[device_ind])
			np.savetxt(self.device_names[device_ind]+"_eqe_Jvals.txt",self.Jvals[device_ind])

			print(self.device_names[device_ind]," rmse is ", rmse)

		return rmse

	def rolloff_objective( self, ELrates_list, plot=False ):
		'''
		Objective function for Nelder-Mead that runs steady state to calculate eqe and compare to experimental data
		input: list of optimizable rates {EL1,EL2,...,ag-EL1,ag-EL2,...,agal-EL1,agal-EL2,...,al-EL1,al-EL2,...}
		output: rmse between model roll-off curve and experimental
		'''

		# enforce that the rates be positive numbers because Nelder-Mead can't handle bounds on the parameters
		if any(rate < 0 for rate in ELrates_list):
			rmse = 1e3
		else:
			rmse = 0
			# loop over each device: Ag, AgAl, Al
			for device_ind in range(len(self.device_names)):

				# update EL rates in full rate dictionary
				i = 0
				for key in self.ELkeys:
					self.rates_dict[key] = ELrates_list[i]
					i += 1
				j = 0
				for key in self.host_ELkeys:
					self.rates_dict[key] = ELrates_list[i+(self.host_device_indices[device_ind]*len(self.host_ELkeys))+j]
					j += 1
				k = 0
				for key in self.layer_ELkeys:
					self.rates_dict[key] = ELrates_list[i+(len(self.host_ELkeys)*len(self.host_ELrates_dicts))+(self.layer_device_indices[device_ind]*len(self.layer_ELkeys))+k]
					k += 1

				# use a constrained k_rec/k_sep ratio
				if self.set_rec_sep and ("k_KrecD" in self.rates_dict):
					self.rates_dict['i_KsepD'] = self.rates_dict['k_KrecD'] / self.rec_sep_ratios[device_ind]

				# use the correct device specific PL rates
				host_ind = self.host_device_indices[device_ind]
				for key in self.host_PLkeys:
					self.rates_dict[key] = self.host_PLrates_dicts[host_ind][key]

				lay_ind = self.layer_device_indices[device_ind]
				for key in self.layer_PLkeys:
					self.rates_dict[key] = self.layer_PLrates_dicts[lay_ind][key]

				self.rates_dict = OrderedDict(sorted(self.rates_dict.items()))
				rates_list = [val for key,val in self.rates_dict.items()]

				# fortran is annoying so embed rates and current densities in lists of determined length
				fortran_rates = [0*i for i in range(30)]
				for i in range(len(rates_list)):
					fortran_rates[i] = rates_list[i]
				Njvals = len(self.Jvals[device_ind])
				fortran_Jvals = [0*i for i in range(1000)]
				for i in range(Njvals):
					fortran_Jvals[i] = self.Jvals[device_ind][i]
		
				# run steady-state for given rates and calculate rmse
				eqe,Js,fail = kinetics.steady_state_el( Njvals, fortran_Jvals, self.Concs, fortran_rates )
				if fail != 0:
					print("Did not reach steady-state within 1 sec time evolution")
					exit()
				eqe = eqe[:Njvals]	# remove extra zeros
				Js = Js[:Njvals]	# remove extra zeros
				rmse += self.rolloff_rmse( eqe, device_ind, plot=plot )

		# print rmse and optimizable rates
		self.rmse_final = rmse
		if rmse != 1e3:
			print_rmse = "rmse = " + str(rmse)
			i = 0
			for key in self.ELkeys:
				print_rmse += " " + key + ":" + str(sn(ELrates_list[i]))
				i += 1
			for device_ind in range(len(self.device_names)):
				for j in range(len(self.host_ELkeys)):
					print_rmse += " " + self.device_names[device_ind] + ":" + self.host_ELkeys[j] + ": " + str(sn(ELrates_list[i+(self.host_device_indices[device_ind]*len(self.host_ELkeys))+j]))
					if self.set_rec_sep and self.host_ELkeys[j]=='k_KrecD':
						print_rmse += " " + self.device_names[device_ind] + ":i_KsepD: " + str(sn(self.rates_dict['i_KsepD']))
				for k in range(len(self.layer_ELkeys)):
					print_rmse += " " + self.device_names[device_ind] + ":" + self.layer_ELkeys[k] + ": " + str(sn(ELrates_list[i+(len(self.host_ELkeys)*len(self.host_ELrates_dicts))+(self.layer_device_indices[device_ind]*len(self.layer_ELkeys))+k]))
					if self.set_rec_sep and self.layer_ELkeys[j]=='k_KrecD':
						print_rmse += " " + self.device_names[device_ind] + ":i_KsepD: " + str(sn(self.rates_dict['i_KsepD']))
			print(print_rmse)
		return rmse

	def rolloff_fit( self, PLrates_dict, host_PLrates_dicts, layer_PLrates_dicts, ELrates_dict, host_ELrates_dicts, layer_ELrates_dicts, optimize=False, plot=False, set_rec_sep=False ):
		'''
		Wrapper for roll-off rate optimization
		'''
		self.PLrates_dict = PLrates_dict
		self.set_rec_sep = set_rec_sep

		# alphabetize and store each optimizable dictionary in oled object
		self.ELrates_dict = OrderedDict(sorted(ELrates_dict.items()))
		for dev_ind in range(len(self.device_names)):
			lay_ind = self.layer_device_indices[dev_ind]
			host_ind = self.host_device_indices[dev_ind]
			layer_PLrates_dicts[lay_ind] = OrderedDict(sorted(layer_PLrates_dicts[lay_ind].items()))
			host_PLrates_dicts[host_ind] = OrderedDict(sorted(host_PLrates_dicts[host_ind].items()))
			layer_ELrates_dicts[lay_ind] = OrderedDict(sorted(layer_ELrates_dicts[lay_ind].items()))
			host_ELrates_dicts[host_ind] = OrderedDict(sorted(host_ELrates_dicts[host_ind].items()))

		self.layer_PLrates_dicts = layer_PLrates_dicts
		self.host_PLrates_dicts = host_PLrates_dicts
		self.layer_ELrates_dicts = layer_ELrates_dicts
		self.host_ELrates_dicts = host_ELrates_dicts
		self.layer_PLkeys = list(self.layer_PLrates_dicts[0].keys())
		self.host_PLkeys = list(self.host_PLrates_dicts[0].keys())

		# build a dictionary of all te nonzero rates and sort alphabetically
		rates_dict = self.PLrates_dict.copy()
		rates_dict.update(ELrates_dict)
		rates_dict.update(self.host_PLrates_dicts[0])
		rates_dict.update(self.layer_PLrates_dicts[0])
		rates_dict.update(layer_ELrates_dicts[0])
		rates_dict.update(host_ELrates_dicts[0])
		rates_dict = {key:val for key,val in rates_dict.items() if val!=0 }
		self.rates_dict = OrderedDict(sorted(rates_dict.items()))
		self.layer_ELkeys = list(self.layer_ELrates_dicts[0].keys())
		self.host_ELkeys = list(self.host_ELrates_dicts[0].keys())
		self.ELkeys = list(self.ELrates_dict.keys())
		print("full dictionary of rates included in the model: ",self.rates_dict)
		print("optimizable rates: ",self.ELkeys)
		print("host specific optimizable rates: ",self.host_ELkeys)
		print("layer specific optimizable rates: ",self.layer_ELkeys)
		print()

		# turn dictionary of optimizable rates into list for optimizer
		# [r1, r2, ..., Ag_r1, Ag_r2, ..., AgAl_r1, AgAl_r2, ..., Al_r1, Al_r2, ...]
		ELrates_dict = OrderedDict(sorted(ELrates_dict.items()))
		ELrates_list = [val for key,val in ELrates_dict.items() if val!=0]
		for host_dict in host_ELrates_dicts:
			host_dict = OrderedDict(sorted(host_dict.items()))	# sort device specific rates
			[ELrates_list.append(val) for key,val in host_dict.items()]
		for layer_dict in layer_ELrates_dicts:
			layer_dict = OrderedDict(sorted(layer_dict.items()))	# sort device specific rates
			[ELrates_list.append(val) for key,val in layer_dict.items()]

		if optimize:
			# optimize rates
			result = minimize( self.rolloff_objective, ELrates_list, method='Nelder-Mead', options={'disp': True,'maxiter': 1000} )
			print("Nelder-Mead done!\nfinal optimized rates = ",result.x) 
	
			# save and plot results of fit
			nonzero_PLrates = {key:val for key,val in self.PLrates_dict.items() if val!=0}
			print("PL fixed rates: ",nonzero_PLrates)
			print("PL fixed host specific rates: ", self.host_PLrates_dicts)
			print("PL fixed layer specific rates: ", self.layer_PLrates_dicts)
			eqe_final = self.rolloff_objective( result.x, plot=plot )
		else:
			# save and plot results of fit
			nonzero_PLrates = {key:val for key,val in self.PLrates_dict.items() if val!=0}
			print("PL fixed rates: ",nonzero_PLrates)
			print("PL fixed host specific rates: ", self.host_PLrates_dicts)
			print("PL fixed layer specific rates: ", self.layer_PLrates_dicts)
			eqe_final = self.rolloff_objective( ELrates_list, plot=plot )

		return 

