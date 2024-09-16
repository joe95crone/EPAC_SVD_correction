# EPAC Correction Script
import sys
import subprocess
import re
import numpy as np
import time
from string import Template
from itertools import chain

# READ_ME
# run as python3 correction.py <lattice> <beamline> <energy (GeV)> <correction case> <SINGLE|MULTI|BUNCH_JITTER|MULTI_BUNCH> <TIK|TRUNC>
# classes for response matrix determination, error assignment + running, correction running
# class to take in all of the required data
# main should be able to integrate all classes & make whole code run just from command line inputs

class Lattice_Search:
	# class variables
	
	# constructor
	def __init__(self, name):
		self.name = name
	
	# defines a temporary lattice via a simple elegant run
	# uses save lattice with output_seq = 2, just needs a save_lattice command nothing else!
	# output_seq =2 produces lattice in one long string L0001 to search for
	def lattice_generator(lattice, beamline, energy):
		subprocess.call(["./lat_gen.sh",str(lattice),str(beamline),str(energy)]) 
		
	
	# read in temp.lte & take lattice
	# run lattice generator from within this function
	def define_beamline(self, lattice, beamline, energy):
		removal = ['L0001','LINE',str(beamline.upper())]
		elements = []
		self.lattice_generator(lattice, beamline, energy)
		with open('temp.lte','r') as text_file:
			lines = text_file.readlines()
			for line in lines:
				if "L0001:" in line:
					#'\w+' means "a word character (a-z etc.) repeated one or more times
					 elements.append(re.findall(r"[\w']+",line))
				elif "&" in line:
					elements.append(re.findall(r"[\w]+",line))
				elif ")" in line:
					elements.append(re.findall(r"[\w]+",line))
			elements = list(chain.from_iterable(elements))
			elements = list(filter(lambda x: x not in removal, elements))
			return elements	
	
	# read in lattice and produce list of magnets
	# only for specified magnet types
	# these are the magnets the errors are placed onto
	def define_magnets(self, lattice, beamline, energy):
		elements = self.define_beamline(self, lattice, beamline, energy)
		mag_elements = []
		with open('temp.lte','r') as text_file:
			lines = text_file.readlines()
			for line in lines:
				for k in range(len(elements)):
					if "{0}: KQUAD".format(elements[k]) in line:
						mag_elements.append(elements[k])
					elif "{0}: QUAD".format(elements[k]) in line:
						mag_elements.append(elements[k])
					elif "{0}: EKICKER".format(elements[k]) in line:
						mag_elements.append(elements[k])
					elif "{0}: SBEN".format(elements[k]) in line:
						mag_elements.append(elements[k])
					elif "{0}: RBEN".format(elements[k]) in line:
						mag_elements.append(elements[k])
			return mag_elements
			
	def define_correctors(self, lattice, beamline, energy, cor_case):
		elements = self.define_beamline(self, lattice, beamline, energy)
		cor_elements = [[],[],[]]
		with open('temp.lte','r') as text_file:
			lines = text_file.readlines()
			for line in lines:
				for k in range(len(elements)):
					if "{0}: KQUAD".format(elements[k]) in line:
						cor_elements[0].append(elements[k])
					elif "{0}: EKICKER".format(elements[k]) in line:
						cor_elements[1].append([elements[k]])
					elif "{0}: QUAD".format(elements[k]) in line:
						cor_elements[2].append([elements[k]])
			# massive switch case to determine 
			if cor_case == "PMQ":
				if energy == 1000:
					n = 4
					cor_elements = [cor_elements[0][i:i + n] for i in range(0, len(cor_elements[0]), n)]
				elif energy == 5000:
					n = 2
					cor_elements = [cor_elements[0][i:i + n] for i in range(0, len(cor_elements[0]), n)]
			elif cor_case == "DED" and beamline != "PMQarray":
				cor_elements = cor_elements[1]
			elif cor_case == "EMQ" and beamline != "PMQarray":
				cor_elements = cor_elements[2]
			elif cor_case == "PMQ_DED" and beamline != "PMQarray":
				if energy == 1000:
					n = 4
					cor_elements = [cor_elements[0][i:i + n] for i in range(0, len(cor_elements[0]), n)] + cor_elements[1]
				elif energy == 5000:
					n = 2
					cor_elements = [cor_elements[0][i:i + n] for i in range(0, len(cor_elements[0]), n)] + cor_elements[1]
			elif cor_case == "PMQ_EMQ" and beamline != "PMQarray":
				if energy == 1000:
					n = 4
					cor_elements = [cor_elements[0][i:i + n] for i in range(0, len(cor_elements[0]), n)] + cor_elements[2]
				elif energy == 5000:
					n = 2
					cor_elements = [cor_elements[0][i:i + n] for i in range(0, len(cor_elements[0]), n)] + cor_elements[2]
			elif cor_case == "DED_EMQ" and beamline != "PMQarray":
				cor_elements = cor_elements[1] + cor_elements[2]
			elif cor_case == "PMQ_DED_EMQ" and beamline != "PMQarray":
				if energy == 1000:
					n = 4
					cor_elements = [cor_elements[0][i:i + n] for i in range(0, len(cor_elements[0]), n)] + cor_elements[1] + cor_elements[2]
				elif energy == 5000:
					n = 2
					cor_elements = [cor_elements[0][i:i + n] for i in range(0, len(cor_elements[0]), n)] + cor_elements[1] + cor_elements[2]	
			else:
				print("Invalid correction case selection")
				exit()
			
			return cor_elements
	
	# function to get the no. BPMs in the lattice	
	def no_BPMs(self, lattice, beamline, energy):
		elements = self.define_beamline(self, lattice, beamline, energy)
		nBPM = 0
		for i in range(len(elements)):
			if elements[i] == "BPM":
				nBPM += 1
		return nBPM
	
	# function to get the corrector limits
	def define_corrector_limits(correctors, case, energy):
		PMQlim = 2e-3
		DEDlim = 3e-3 # @ 1 GeV
		EMQlim = 3e-3 # @ 1 GeV
		lim = []
		nDED = 0
		for i in range(len(correctors)):
			if len(correctors[i]) > 1:
				lim.append(PMQlim)
			elif len(correctors[i]) == 1 and "DED" in case and nDED < 2:
				lim.append(DEDlim*(1000/energy))
				nDED = nDED +1
			else:
				lim.append(EMQlim*(1000/energy))
		return lim
			
# writes SDDS parameter files of passed data, assuming data is of correct format! 
# 2D array following column of data specified in the header
def write_params(dat):
	param = {
    	'header': 'SDDS1\n&parameter name=InputFile, description="Response Bump Correction", type=string, &end\n&column name=ElementName, type=string,  &end\
	\n&column name=ElementParameter, type=string,  &end\n&column name=ParameterValue, type=double,  &end\
	\n&column name=ParameterMode, type=string,  &end\n&data mode=ascii, &end\n! page number 1\npipe',
    	'page_no': '\t'+' '+str(dat.shape[0]),
    	'data': '\n'.join([' '.join(np.append([str(dat[i,j]) for j in range(dat.shape[1])],str("differential"))) for i in range(dat.shape[0])])
	}

	with open('param/template.txt', 'r') as infile:
    		src = Template(infile.read())
    		result = src.substitute(param)
	
	with open('param/mod.param', 'w') as outfile:
		outfile.write(result)
	
	

class Response:
	# Class Variable
	Rrootpath = 'bpm_files/response/{0}_{1}_{2}_{3}_bpm.dat'
	Rsavepath = 'response_matrices/rmatrix_{0}_{1}GeV_{2}_{3}.dat'
	invRsavepath = 'response_matrices/inv_rmatrix_{0}_{1}GeV_{2}_{3}.dat'
	
	standard_path = '../FBPIC_Beam/Saved_Beams/OF2_EPAC_KDE_1GeV.sdds'
	#hb_path = '../FBPIC_Beam/2D_SCAN/EPAC_2D_SCAN_{0}.sdds'
	hb_path = '../FBPIC_Beam/3D_SCAN_HB/EPAC_3D_1GeV_HB_{0}.sdds'
	# constructor	
	def __init__(self, name):
		self.name = name	
	
	# run the sims to produce elegant response matrix
	# void function! only executes other code	
	def response_run(self, lattice, beamline, energy, seed, err_vals, correctors, kick_vals, HV_FLAG, beam_dist=0):
		# switch for beam distribution
		if beam_dist == 0:
			path = self.standard_path
		elif beam_dist != 0:
			path = self.hb_path.format(beam_dist)
		
		# switch for kick based on energy
		kick_en = 1
		if energy == 1000:
			kick_en = 1
		elif energy == 5000:
			kick_en = 2
		for i in range(len(correctors)):
			response_kick = []
			for j in range(len(correctors[i])):
				if HV_FLAG == "X":
					response_kick.append([correctors[i][j], "DX" if "QUAD" in correctors[i][j] else "HKICK", str(kick_vals[0] if "QUAD" in correctors[i][j] else kick_vals[kick_en])])
				elif HV_FLAG == "Y":
					response_kick.append([correctors[i][j], "DY" if "QUAD" in correctors[i][j] else "VKICK", str(kick_vals[0] if "QUAD" in correctors[i][j] else kick_vals[kick_en])])	
			response_kick = np.array(response_kick)
			write_params(response_kick)
			subprocess.call(["./response.sh",str(lattice),str(beamline),str(energy),str(seed),str(err_vals[0]),str(err_vals[1]),str(err_vals[2]),str(err_vals[3]),str(err_vals[4]),str('param/mod.param'),str(i+1),str(HV_FLAG),str(path),str(beam_dist)])	
	
	
	# produces weighted response matrix
	# row weights are BPM weights, column weights are corrector weights
	def r_matrix(self, correctors, BPM_weight, COR_weight, nBPM, HV_FLAG, seed, beam_dist=0):
		Rmatrix = np.empty((len(correctors),nBPM))
		for k in range(len(correctors)):
			# column weighting (correctors)
			Rmatrix[k] = np.sqrt(BPM_weight[k])*np.genfromtxt(self.Rrootpath.format(seed,k+1,HV_FLAG,beam_dist),skip_header=1,dtype=float)[:,1 if HV_FLAG == "X" else 2]
		# row weighting (BPMs)
		return (Rmatrix.T * np.sqrt(COR_weight)).T
		
	# produces + saves response matrix	
	def save_r_matrix(self, beamline, energy, HV_FLAG, correctors, BPM_weight, COR_weight, nBPM, seed, beam_dist=0):
		np.savetxt(self.Rsavepath.format(beamline,energy/10**3,HV_FLAG,beam_dist),self.r_matrix(self,correctors,BPM_weight,COR_weight,nBPM,HV_FLAG,seed,beam_dist))
		return self.r_matrix(self,correctors,BPM_weight,COR_weight,nBPM,HV_FLAG,seed,beam_dist)
		
	# SVD inverse of a matrix + cutting
	def SVD_inv_r_matrix(self, matrix, cutoff, svd_method):	
		# choice of tikhonov regularization or truncated SVD
		if svd_method == 'TRUNC': 
			invertRmatrix = np.linalg.pinv(matrix,rcond=(cutoff))
			return invertRmatrix
		elif svd_method == 'TIK':
			# SVD of matrix
			U, S, Vh = np.linalg.svd(matrix, full_matrices=False)
			# Tikhonov regularization parameter
			alphaTik = np.mean(S)/150
			# SVD components for inverse of matrix
			V = Vh.T
			Uh = U.T
			STik = np.diag(S/(S**2 + alphaTik))
			invertRmatrix = np.linalg.multi_dot([V, STik, Uh])
			return invertRmatrix
		else:
			print("Invalid SVD method selection <TIK|TRUNC>")
		
	# produces + saves inverse response matrix	
	def save_inv_r_matrix(self, beamline, energy, HV_FLAG, correctors, BPM_weight, COR_weight, nBPM, seed, cutoff, svd_method, beam_dist=0):
		np.savetxt(self.invRsavepath.format(beamline,energy/10**3,HV_FLAG,beam_dist),self.SVD_inv_r_matrix(self,self.r_matrix(self,correctors,BPM_weight,COR_weight,nBPM,HV_FLAG,seed,beam_dist),cutoff,svd_method))
		return self.SVD_inv_r_matrix(self,self.r_matrix(self,correctors,BPM_weight,COR_weight,nBPM,HV_FLAG,seed,beam_dist),cutoff,svd_method)

class Error: 
	# class variables
	standard_path = '../FBPIC_Beam/Saved_Beams/OF2_EPAC_KDE_1GeV.sdds'
	#hb_path = '../FBPIC_Beam/2D_SCAN/EPAC_2D_SCAN_{0}.sdds'
	hb_path = '../FBPIC_Beam/3D_SCAN_HB/EPAC_3D_1GeV_HB_{0}.sdds'
	# constructor
	def __init__(self, name):
		self.name = name
	
	#runs error run with passed seed for any beamline
	# void function
	def error_run(self, lattice, beamline, energy, ErrVals, seed, trial, beam_dist=0):
		if beam_dist == 0:
			path = self.standard_path
		elif beam_dist != 0:
			path = self.hb_path.format(beam_dist)
		# running of error run		
		subprocess.call(["./error.sh",str(lattice),str(beamline),str(energy),str(seed),str(ErrVals[0]),str(ErrVals[1]),str(ErrVals[2]),str(ErrVals[3]),str(ErrVals[4]),str(seed),str(trial),str(path),str(beam_dist)]) 
		
class Correction:
	# class variables
	ErrBPMpath = 'bpm_files/error/{0}_{1}_{2}_bpm.dat'	
	corsavepath = 'correction_run/{0}_{1}_{2}_{3}_{4}_{5}_corrector.dat'
	CorBPMpath = 'bpm_files/correct/{0}_{1}_{2}_bpm.dat'
	
	standard_path = '../FBPIC_Beam/Saved_Beams/OF2_EPAC_KDE_1GeV.sdds'
	#hb_path = '../FBPIC_Beam/2D_SCAN/EPAC_2D_SCAN_{0}.sdds'
	hb_path = '../FBPIC_Beam/3D_SCAN_HB/EPAC_3D_1GeV_HB_{0}.sdds'
	# constructor
	def __init__(self, name):
		self.name = name
	
	# function to get the corrector values
	def corrector_vals(self, energy, seed, trial, correctors, KickVals, invRxmatrix, invRymatrix, beam_dist=0):
		BPM = []
		cor_kick = []
		BPMpath = ' '
		if trial > 0:
			BPMpath = self.CorBPMpath
			# get mean BPM values
			BPM.append([np.genfromtxt(BPMpath.format(seed,trial-1,beam_dist),skip_header=1,dtype=float)[:,1], np.genfromtxt(BPMpath.format(seed,trial-1,beam_dist),skip_header=1,dtype=float)[:,2]])
		else:
			BPMpath = self.ErrBPMpath
			# get mean BPM values
			BPM.append([np.genfromtxt(BPMpath.format(seed,trial,beam_dist),skip_header=1,dtype=float)[:,1], np.genfromtxt(BPMpath.format(seed,trial,beam_dist),skip_header=1,dtype=float)[:,2]])
		# get mean BPM values
		BPM = np.array(BPM).squeeze()
		# get list of initial misalignment kick values
		# switch for kick based on energy
		kick_en = 1
		if energy == 1000:
			kick_en = 1
		elif energy == 5000:
			kick_en = 2
		for i in range(len(correctors)):
			cor_kick.append(-KickVals[0 if "QUAD" in correctors[i] else kick_en])
		cor_kick = np.array(cor_kick)
		correction = np.array([cor_kick*np.dot(BPM[0],invRxmatrix), cor_kick*np.dot(BPM[1],invRymatrix)])
		return correction
	
	# function to save the corrector values + return them
	def save_corrector_vals(self, lattice, beamline, energy, cor_case, seed, trial, correctors, KickVals, invRxmatrix, invRymatrix, beam_dist=0):
		np.savetxt(self.corsavepath.format(beamline,energy/10**3,cor_case,seed,trial,beam_dist),self.corrector_vals(self,energy,seed,trial,correctors,KickVals,invRxmatrix,invRymatrix,beam_dist))
		return self.corrector_vals(self, energy, seed, trial, correctors, KickVals, invRxmatrix, invRymatrix, beam_dist)
		
	def corrector_run(self, lattice, beamline, energy, seed, err_vals, mag_cor, cor_vals, trial, beam_dist=0):
		if beam_dist == 0:
			path = self.standard_path
		elif beam_dist != 0:
			path = self.hb_path.format(beam_dist)
		cor_vals_param = []
		cor_type = [["DX", "DY"], ["HKICK", "VKICK"]]
		cor_mag_flat = list(chain.from_iterable(mag_cor))
		for i in range(cor_vals.shape[0]):
			for j in range(len(mag_cor)):
				for k in range(len(mag_cor[j])):
					cor_vals_param.append([mag_cor[j][k], cor_type[0 if len(mag_cor[j]) > 1 else 1][i], str(cor_vals[i][j])])			
		cor_vals_param = np.array(cor_vals_param)	
		write_params(cor_vals_param)
		subprocess.call(["./correct.sh",str(lattice),str(beamline),str(energy),str(seed),str(err_vals[0]),str(err_vals[1]),str(err_vals[2]),str(err_vals[3]),str(err_vals[4]),str('param/mod.param'),str(trial-1),str(path),str(beam_dist)])	
	
	# untested!!!!
	def corrector_run_beamvar(self, lattice, beamline, energy, seed, err_vals, mag_cor, cor_vals, trial, beam_dist=0):
		if beam_dist == 0:
			print("Fail - not a modified bunch")
		else:
			path = self.hb_path.format(beam_dist)
			subprocess.call(["./correct.sh",str(lattice),str(beamline),str(energy),str(seed),str(err_vals[0]),str(err_vals[1]),str(err_vals[2]),str(err_vals[3]),str(err_vals[4]),str('param/mod.param'),str(trial-1),str(path),str(beam_dist)])
		
if __name__ == "__main__":
	#-------				
	# SETUP
	#-------
	
	# clear previous run
	subprocess.call(["./remove_correction.sh"])
	
	# RUN SETTINGS
	# BEAMLINE SETTINGS
	# take in the lattice, beamline, energy (convert to MeV) & wether we do a single run or multi-seed from command line
	lattice = str(sys.argv[1])
	beamline = str(sys.argv[2])
	energy = int(sys.argv[3])*10**3
	
	# CORRECTION SETTINGS
	#takes in correction case <PMQ>|<DED>|<EMQ>|<PMQ_DED>|<PMQ_EMQ>|<DED_EMQ>|<PMQ_DED_EMQ>
	cor_case = str(sys.argv[4])
	
	# SVD METHOD
	# either <TIK|TRUNC> for tikhonov regularization or truncated SVD
	svd_method = str(sys.argv[6])
	
	# LATTICE SETUP
	lattice_elements = Lattice_Search.define_beamline(Lattice_Search, lattice, beamline, energy) 
	mag_elements = Lattice_Search.define_magnets(Lattice_Search, lattice, beamline, energy) 
	mag_cor = Lattice_Search.define_correctors(Lattice_Search, lattice, beamline, energy, cor_case)
	nBPM = Lattice_Search.no_BPMs(Lattice_Search, lattice, beamline, energy)
	
	# LIMITS
	limits = Lattice_Search.define_corrector_limits(mag_cor, cor_case, energy)
	
	# WEIGHTING
	BPM_weight = np.ones(nBPM)
	#BPM_weight = np.linspace(0.1,1,nBPM)
	
	COR_weight = np.ones(len(mag_cor))
	#COR_weight = np.array([0, 1, 1])
	
	# RUN SWITCH
	# haven't included the option to run a single errored LWFA source bunch because I'm not sure when it would be used!
	# SINGLE performs correction on a single bunch (by default this is the standard bunch distribution)
	# BUNCH_JITTER performs correction on the single bunch (default standard) then runs the errored bunches through this to get a feel for jitter
	# MULTI runs multiple correction runs with different seeds
	# MULTI_BUNCH runs multiple correction runs with different seeds and different bunches 
	# <SINGLE|BUNCH_JITTER|MULTI|MULTI_BUNCH>
	trial_case = str(sys.argv[5]) 
	
	#-----------
	# CONSTANTS
	#-----------
	
	# ERROR
	# error tolerances (transverse misalignments only currently)
	ErrXYinit = 100e-6
	ErrXPYPinit = 10e-6 # should be 10 microrad as standard
	ErrXYpmq = 50e-6
	ErrFSEpmq = 0.01
	ErrXYemq = 100e-6
	error_vals = [ErrXYinit, ErrXPYPinit, ErrXYpmq, ErrFSEpmq, ErrXYemq]
	
	# RESPONSE
	# apply the response_run function to generate all of the bpm positions for response matrix
	# kick values
	DXYpmq = 1e-3
	DKICK1GeV = 1e-3 
	DKICK5GeV = 0.2e-3
	KickVals = [DXYpmq, DKICK1GeV, DKICK5GeV]
	# invert response matrix with cut-off
	# cutoff = 0.4
	if beamline == 'PMQarray':		
		cutoff = 0.3 # believed to get ~2 single values for PMQ array#
	else:
		cutoff = 0.05 # believed to get ~2 single values for AF 5m
	
	# switch cases here i.e. loop starts here + need a failure else
	if trial_case == "SINGLE":
		tic = time.perf_counter()
		seed_magic = 42
		# first error run
		Error.error_run(Error, lattice, beamline, energy, error_vals, seed_magic, 0)
		# response matrix trajectory bumps in X and Y
		Response.response_run(Response, lattice, beamline, energy, seed_magic, error_vals, mag_cor, KickVals, "X")
		Response.response_run(Response, lattice, beamline, energy, seed_magic, error_vals, mag_cor, KickVals, "Y")
		# get the response matrix + save it
		# note the use of the self class here requiring the class name to be passed too! 
		RxMat = Response.save_r_matrix(Response, beamline, energy, "X", mag_cor, BPM_weight, COR_weight, nBPM, seed_magic)
		RyMat = Response.save_r_matrix(Response, beamline, energy, "Y", mag_cor, BPM_weight, COR_weight, nBPM, seed_magic)	
		# production of inverse response matrices
		invRxMat = Response.save_inv_r_matrix(Response, beamline, energy, "X", mag_cor, BPM_weight, COR_weight, nBPM, seed_magic, cutoff, svd_method)
		invRyMat = Response.save_inv_r_matrix(Response, beamline, energy, "Y", mag_cor, BPM_weight, COR_weight, nBPM, seed_magic, cutoff, svd_method)
		# corrector values + corrector run
		ntrial = 20
		for trial in range(0,ntrial):
			new_correct = Correction.save_corrector_vals(Correction, lattice, beamline, energy, cor_case, seed_magic, trial, mag_cor, KickVals, invRxMat, invRyMat)
			if np.all(new_correct == 0):
				break
			else:
				# corrector values are the previous correction + the new correction!
				if trial == 0:
					cor_vals = new_correct
				else:
					cor_vals = cor_vals + new_correct
					for j in range(cor_vals.shape[0]):
						for i in range(cor_vals.shape[1]):
							if abs(cor_vals[j,i]) > limits[i] and cor_vals[j,i] > 0:
								cor_vals[j,i] = limits[i]
							elif abs(cor_vals[j,i]) > limits[i] and cor_vals[j,i] < 0:
								cor_vals[j,i] = -limits[i]
				# run with these corrector values
				Correction.corrector_run(Correction, lattice, beamline, energy, seed_magic, error_vals, mag_cor, cor_vals, trial+1)
		print(cor_vals)
		Correction.corrector_run(Correction, lattice, beamline, energy, seed_magic, error_vals, mag_cor, cor_vals, 0)
		toc = time.perf_counter()
		print("sim time: ",toc-tic)
	if trial_case == "BUNCH_JITTER":
		tic = time.perf_counter()
		seed_magic = 12
		# first error run
		Error.error_run(Error, lattice, beamline, energy, error_vals, seed_magic, 0)
		# response matrix trajectory bumps in X and Y
		Response.response_run(Response, lattice, beamline, energy, seed_magic, error_vals, mag_cor, KickVals, "X")
		Response.response_run(Response, lattice, beamline, energy, seed_magic, error_vals, mag_cor, KickVals, "Y")
		# get the response matrix + save it
		# note the use of the self class here requiring the class name to be passed too! 
		RxMat = Response.save_r_matrix(Response, beamline, energy, "X", mag_cor, BPM_weight, COR_weight, nBPM, seed_magic)
		RyMat = Response.save_r_matrix(Response, beamline, energy, "Y", mag_cor, BPM_weight, COR_weight, nBPM, seed_magic)	
		# production of inverse response matrices
		invRxMat = Response.save_inv_r_matrix(Response, beamline, energy, "X", mag_cor, BPM_weight, COR_weight, nBPM, seed_magic, cutoff, svd_method)
		invRyMat = Response.save_inv_r_matrix(Response, beamline, energy, "Y", mag_cor, BPM_weight, COR_weight, nBPM, seed_magic, cutoff, svd_method)
		# corrector values + corrector run
		ntrial = 20
		for trial in range(0,ntrial):
			new_correct = Correction.save_corrector_vals(Correction, lattice, beamline, energy, cor_case, seed_magic, trial, mag_cor, KickVals, invRxMat, invRyMat)
			if np.all(new_correct == 0):
				break
			else:
				# corrector values are the previous correction + the new correction!
				if trial == 0:
					cor_vals = new_correct
				else:
					cor_vals = cor_vals + new_correct
					for j in range(cor_vals.shape[0]):
						for i in range(cor_vals.shape[1]):
							if abs(cor_vals[j,i]) > limits[i] and cor_vals[j,i] > 0:
								cor_vals[j,i] = limits[i]
							elif abs(cor_vals[j,i]) > limits[i] and cor_vals[j,i] < 0:
								cor_vals[j,i] = -limits[i]
				# run with these corrector values
				Correction.corrector_run(Correction, lattice, beamline, energy, seed_magic, error_vals, mag_cor, cor_vals, trial+1)
		
		#Correction.corrector_run(Correction, lattice, beamline, energy, seed_magic, error_vals, mag_cor, cor_vals, ntrial)
		nbunch = 30
		for bun in range(1,nbunch+1):
			Correction.corrector_run_beamvar(Correction, lattice, beamline, energy, seed_magic, error_vals, mag_cor, cor_vals, 0, bun)
		toc = time.perf_counter()
		print(cor_vals)
		print("sim time: ",toc-tic)
	elif trial_case == "MULTI":
		nseed = 50
		maxrand = 1e8
		seed_vals = []
		sim_time = []
		for i in range(nseed):
			tic = time.perf_counter()
			seed = np.random.randint(maxrand)
			# first error run
			Error.error_run(Error, lattice, beamline, energy, error_vals, seed, 0)
			# response matrix trajectory bumps in X and Y
			Response.response_run(Response, lattice, beamline, energy, seed, error_vals, mag_cor, KickVals, "X")
			Response.response_run(Response, lattice, beamline, energy, seed, error_vals, mag_cor, KickVals, "Y")
			# get the response matrix + save it
			# note the use of the self class here requiring the class name to be passed too! 
			RxMat = Response.save_r_matrix(Response, beamline, energy, "X", mag_cor, BPM_weight, COR_weight, nBPM, seed)
			RyMat = Response.save_r_matrix(Response, beamline, energy, "Y", mag_cor, BPM_weight, COR_weight, nBPM, seed)	
			# production of inverse response matrices
			invRxMat = Response.save_inv_r_matrix(Response, beamline, energy, "X", mag_cor, BPM_weight, COR_weight, nBPM, seed, cutoff, svd_method)
			invRyMat = Response.save_inv_r_matrix(Response, beamline, energy, "Y", mag_cor, BPM_weight, COR_weight, nBPM, seed, cutoff, svd_method)
			# corrector values + corrector run
			ntrial = 5
			for trial in range(0,ntrial):
				new_correct = Correction.save_corrector_vals(Correction, lattice, beamline, energy, cor_case, seed, trial, mag_cor, KickVals, invRxMat, invRyMat)
				if np.all(new_correct == 0):
					break
				else:
					# corrector values are the previous correction + the new correction!
					if trial == 0:
						cor_vals = new_correct
					else:
						cor_vals = cor_vals + new_correct
						for j in range(cor_vals.shape[0]):
							for i in range(cor_vals.shape[1]):
								if abs(cor_vals[j,i]) > limits[i] and cor_vals[j,i] > 0:
									cor_vals[j,i] = limits[i]
								elif abs(cor_vals[j,i]) > limits[i] and cor_vals[j,i] < 0:
									cor_vals[j,i] = -limits[i]
					# run with these corrector values
					Correction.corrector_run(Correction, lattice, beamline, energy, seed, error_vals, mag_cor, cor_vals, trial+1)
			#print(cor_vals)
			toc = time.perf_counter()
			sim_time.append(toc-tic)
			seed_vals.append(seed)
		print(sum(sim_time))
		# write out times
		time_vals = open("correction_run/sim_time.dat","w")
		for i in range(len(sim_time)):
			time_vals.write(str(sim_time[i])+"\n")
	elif trial_case == "MULTI_BUNCH":
		nseed = 20
		nbunch = 11
		maxrand = 1e8
		seed_vals = []
		sim_time = []
		for bun in range(1,nbunch+1):
			for i in range(nseed):
				tic = time.perf_counter()
				seed = np.random.randint(maxrand)
				# first error run
				Error.error_run(Error, lattice, beamline, energy, error_vals, seed, 0, bun)
				# response matrix trajectory bumps in X and Y
				Response.response_run(Response, lattice, beamline, energy, seed, error_vals, mag_cor, KickVals, "X", bun)
				Response.response_run(Response, lattice, beamline, energy, seed, error_vals, mag_cor, KickVals, "Y", bun)
				# get the response matrix + save it
				# note the use of the self class here requiring the class name to be passed too! 
				RxMat = Response.save_r_matrix(Response, beamline, energy, "X", mag_cor, BPM_weight, COR_weight, nBPM, seed, bun)
				RyMat = Response.save_r_matrix(Response, beamline, energy, "Y", mag_cor, BPM_weight, COR_weight, nBPM, seed, bun)	
				# production of inverse response matrices
				invRxMat = Response.save_inv_r_matrix(Response, beamline, energy, "X", mag_cor, BPM_weight, COR_weight, nBPM, seed, cutoff, svd_method, bun)
				invRyMat = Response.save_inv_r_matrix(Response, beamline, energy, "Y", mag_cor, BPM_weight, COR_weight, nBPM, seed, cutoff, svd_method, bun)
				# corrector values + corrector run
				ntrial = 10
				for trial in range(0,ntrial):
					new_correct = Correction.save_corrector_vals(Correction, lattice, beamline, energy, cor_case, seed, trial, mag_cor, KickVals, invRxMat, invRyMat, bun)
					if np.all(new_correct == 0):
						break
					else:
						# corrector values are the previous correction + the new correction!
						if trial == 0:
							cor_vals = new_correct
						else:
							cor_vals = cor_vals + new_correct
							for j in range(cor_vals.shape[0]):
								for i in range(cor_vals.shape[1]):
									if abs(cor_vals[j,i]) > limits[i] and cor_vals[j,i] > 0:
										cor_vals[j,i] = limits[i]
									elif abs(cor_vals[j,i]) > limits[i] and cor_vals[j,i] < 0:
										cor_vals[j,i] = -limits[i]
						# run with these corrector values
						Correction.corrector_run(Correction, lattice, beamline, energy, seed, error_vals, mag_cor, cor_vals, trial+1, bun)
				# print(cor_vals)
				toc = time.perf_counter()
				sim_time.append(toc-tic)
				seed_vals.append(seed)
		print(sum(sim_time))
		# write out times
		time_vals = open("correction_run/sim_time.dat","w")
		for i in range(len(sim_time)):
			time_vals.write(str(sim_time[i])+"\n")
	else:
		print("Invalid trial case selection")
		exit()	
