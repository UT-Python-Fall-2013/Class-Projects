###########################################
# Random Utility Models with Phylogenies ##
# Python3/Cython/C						 ##
# Nate Pope								 ##
# Requires: 							 ##
#		-GSL libraries					 ##
#		-Ceygen							 ##
#		-NumPy  						 ##
#		-matplotlib						 ##
#		-Cython							 ##
#		-time			 				 ##
###########################################
import numpy as npy
import poyla as po
import time
import matplotlib.pyplot as plt

class RUMP:
	"""Phylogenetic random utility models"""
	def __init__(self, Y, X, Omega, fixef, nMCMC = 1000, nChains = 3, thin = 1, burnin = 100):
		"""Y is matrix of observations, X is matrix of covariates, Omega is inverse of phylogenetic covariance matrix, 
		fixef is array with indicies of fixed effects. If N observations of Q plants under P covariates:
		Y is N x Q, X is N x P, Omega is P_fixef x P_fixef"""
		## Model settings
		self.nChains = nChains
		self.nMCMC = nMCMC
		
		## Model dimensionality
		Q = Y.shape[1]
		N = Y.shape[0]
		P = X.shape[1]
		
		## MCMC sampling of model
		self.runtime = npy.zeros([nChains])
		beta = npy.empty([nMCMC * thin * nChains, Q-1, P])
		lamb = npy.empty([nMCMC * thin * nChains, Q-1])
		for i in range(nChains):
			startTime = time.time()
			betaSamples, lambSamples = po.MCMCsamp(nMCMC*thin, Y, X, Omega, fixef)
			endTime = time.time()
			for q in range(Q-1):
				lamb[(i*nMCMC*thin):(i*nMCMC*thin+nMCMC*thin),q] = npy.array(lambSamples)[q,:]
				for p in range(P):
					beta[(i*nMCMC*thin):(i*nMCMC*thin+nMCMC*thin),q,p] = npy.array(betaSamples)[p,q,:]
			self.runtime[i] = endTime - startTime
			print("Chain", i+1, "sampled in", self.runtime[i], "seconds.")	## progress report
		
		## create MCMC objects for each parameter
		self.sampleList = {}
		for q in range(Q-1):
			lambName = "lambda_"+str(q)
			self.sampleList[lambName] = MCMC(lamb[:,q], lambName, nMCMC, nChains)
			for p in range(P):
				betaName = "beta_"+str(q)+"_"+str(p)
				self.sampleList[betaName] = MCMC(beta[:,q,p], betaName, nMCMC, nChains)
	
	def summary(self):
		"""add summary statistics, numb obs, etc"""
		print("Convergence Diagnostics:\n")
		self.convergence()
		print("Quantiles:\n")
		self.quantiles()
		
	def samples(self, name = False):
		"""returns samples for all parameters"""
		keys = self.sampleList.keys()
		if name != False:
			return(self.sampleList[name].samples)
		else:
			formatList = ["f4"]*len(keys)
			samples = npy.empty(self.nMCMC*self.nChains, dtype = {'names':list(keys), 'formats':formatList})
			for name in keys:
				samples[name] = self.sampleList[name].samples.flatten()
			return(samples)
		
	def save(self, filename):
		"""save output to .csv file"""
		npy.savetxt(filename, self.samples(), delimiter=",", header = ",".join(self.samples().dtype.names) )
	
	def quantiles(self):
		"""Gives convergence criteria for all parameters"""
		print( "{0:<10}{1:<10}{2:<10}{3:<10}{4:<10}{5:<10}{6:<10}{7:<10}".format("Param.","Mean", "Std.Dev.", "0.025", "0.250", "0.500", "0.750", "0.975") )
		for i in self.sampleList.keys():
			self.sampleList[i].Quantile(title = False)
	
	def convergence(self):
		"""Gives convergence criteria for all parameters"""
		print( "{0:<10}{1:<10}{2:<10}{3:<10}{4:<10}{5:<10}".format("Param.","Samples","Chains","Within","Between","Rhat") )
		for i in self.sampleList.keys():
			self.sampleList[i].Rhat(title = False)
	
	def plot(self, name, type = "d"):
		"""Simply wraps class MCMC methods, to plot a given parameter. Valid plot types are 'd' for density and 't' for trace"""
		if type == "d":
			self.sampleList[name].Density()
		elif type == "t":
			self.sampleList[name].Trace()
		else:
			print("Invalid plot type. Valid types 'd' for density and 't' for trace.")

# constructor function
class MCMC:
	def __init__(self, input, name, nMCMC, nChains):
		"""Takes flattened vector input and wraps into array; where columns are chains, rows are samples."""
		self.samples = npy.reshape(input, [nMCMC, nChains], "F")
		self.name = name
		self.nMCMC = nMCMC
		self.nChains = nChains
	
	def Rhat(self, title = True):
		"""Calculates Gelman-Rubin convergence criterion"""
		## calculate within-, between- chain variance, use to calculate scale reduction factor Rhat
		mu_all = npy.mean(self.samples)
		mu_chain = npy.mean(self.samples, axis = 0)
		var_chain = npy.empty([self.nChains])
		for i in range(self.nChains):
			var_chain[i] = npy.sum((self.samples[:,i] - mu_chain[i])**2)/(self.nMCMC - 1)
		B = npy.sum((mu_chain - var_chain)**2) * self.nMCMC / (self.nChains - 1)
		W = npy.sum(var_chain)/self.nChains
		#W = npy.mean(npy.std(self.samples, axis = 1, ddof = 1)**2)
		#B = (self.nMCMC*npy.mean((npy.mean(self.samples) - npy.mean(self.samples, axis = 1))**2)) ## not quite right
		Var = (1-1/self.nMCMC)*W + B/self.nMCMC
		Rhat = npy.sqrt(Var / W)
		## print labels if desired
		if title == True:
			print( "{0:<10}{1:<10}{2:<10}{3:<10}{4:<10}{5:<10}".format("Param.","Samples","Chains","Within","Between","Rhat") )
		## print convergence diagnostic Rhat
		print( "{0:<10}{1:<10.4g}{2:<10.4g}{3:<10.4g}{4:<10.4g}{5:<10.4g}".format(self.name,self.nMCMC,self.nChains,W,B,Rhat) )
	
	def Quantile(self, title = True):
		"""Calculates quantiles, mean, standard deviation."""
		## calculate mean, sd, quantiles
		values = npy.percentile(self.samples, [2.5, 25, 50, 75, 97.5] )
		mean = npy.mean(self.samples)
		sd = npy.std(self.samples)
		## print labels if desired
		if title == True:
			print( "{0:<10}{1:<10}{2:<10}{3:<10}{4:<10}{5:<10}{6:<10}{7:<10}".format("Param.","Mean", "Std.Dev.", "0.025", "0.250", "0.500", "0.750", "0.975") )
		## print quantiles
		print( "{0:<10}{1:<10.4g}{2:<10.4g}{3:<10.4g}{4:<10.4g}{5:<10.4g}{6:<10.4g}{7:<10.4g}".format(self.name, mean, sd, values[0], values[1], values[2], values[3], values[4]))
	
	def Density(self):
		"""Plots density."""
		plt.hist(self.samples.flatten())
		plt.title("Samples of "+self.name)
		plt.ylabel("# Samples")
		plt.xlabel("Value")
		plt.show()
	
	def Trace(self):
		"""Plots trace."""
		for i in range(self.nChains):
			plt.plot(npy.arange(self.nMCMC), self.samples[:,i])
		plt.title("Trace of "+self.name)
		plt.ylabel("Value")
		plt.xlabel("Iteration")
		plt.show()

#################
## Example Use ##
#################
# load in example data
Y = npy.genfromtxt("sampleData/Y.csv", delimiter="\t", dtype = float)
X = npy.genfromtxt("sampleData/X.csv", delimiter="\t")
Omega = npy.genfromtxt("sampleData/Omega.csv", delimiter="\t")
Omega = npy.linalg.inv(Omega)	#for memoryview version
fixef = npy.array([0])

# create object, sample
model_fit = RUMP(Y, X, Omega, fixef, 1000, 3)

# extract raw samples
model_fit.samples()

# save to file
model_fit.save("test.csv")

# check convergence
model_fit.convergence()

# get quantiles
model_fit.quantiles()

# extract summary
model_fit.summary()

# plot trace
model_fit.plot("lambda_1", "t")

# plot density
model_fit.plot("lambda_1", "d")
