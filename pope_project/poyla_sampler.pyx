# cython: profile=True

#####################################################
#	MCMC functions-Poyla-Gamma multinomial model 	#
# 	Nate Pope										#
#####################################################

##################
## Dependencies ##
##################
cimport cython
import sys
## import numpy
cimport numpy as np
import numpy as np
## import linear algebra for memoryviews
cimport ceygen.elemwise as ce
cimport ceygen.core as cc
cimport ceygen.llt as clt
cimport ceygen.lu as clu
cimport ceygen.reductions as cr 
import ceygen.elemwise as ce
import ceygen.core as cc
import ceygen.llt as clt
import ceygen.lu as clu
import ceygen.reductions as cr 

##################################
## import basic transformations ##
##################################
cdef extern from "math.h": #nogil:
	double log "log"(double)
	double exp "exp"(double)
	double sqrt "sqrt"(double)
	double abs "fabs"(double)	## note that 'abs' coverts to integer, eg. abs(0.5) = 0
cdef double pi = 3.14159265

##############################################
## Import random number generators from GSL ##
##############################################
cdef extern from "gsl/gsl_rng.h":#nogil:
	ctypedef struct gsl_rng_type:
		pass
	ctypedef struct gsl_rng:
		pass
	gsl_rng_type *gsl_rng_mt19937
	gsl_rng *gsl_rng_alloc(gsl_rng_type * T)
  
cdef gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937)

cdef extern from "gsl/gsl_randist.h":#nogil:
	double unif "gsl_rng_uniform"(gsl_rng * r)
	double unif_interval "gsl_ran_flat"(gsl_rng * r,double,double)	## syntax; (seed, lower, upper)
	double gamma "gsl_ran_gamma"(gsl_rng * r,double,double)			## syntax; (seed, shape, scale) ... scale is 1/rate
	double gaussian "gsl_ran_gaussian"(gsl_rng * r,double)			## syntax; (seed, standard deviation)
	double exponential "gsl_ran_exponential"(gsl_rng * r,double)	## syntax; (seed, mean) ... mean is 1/rate

cdef extern from "gsl/gsl_cdf.h":
	double pnorm "gsl_cdf_ugaussian_P"(double)

############################################
## Poyla-Gamma random variate generation ###
############################################
cdef double t = 0.64		## truncation parameter. may influence efficiency ... ?
## Poyla-Gamma sampler
cpdef double[:] PG(double[:] ROWSUM, double[:] INPUT): #nogil:
	"""z > 0"""
	### ridiculous list of declared variables, see if you can shorten
	cdef double z = 0			## holder for i-th input
	cdef double N = 0			## holder for i-th rowsum
	cdef double K = 0			## holder for transformation of z
	cdef double mu = 0			## inverse of z, aka mean of inverse-Gaussian distribution
	cdef double X = 0			## holder for random draw of inverse-Gaussian/truncated exponential
	cdef long n = 0				## counter for squeezing iterations at end
	cdef double S = 0			## 0-degree piecewise coefficient
	cdef double Y = 0			## threshold for squeezing algorithm at end 
	cdef double alpha = 0		## threshold for accept/reject of chi-squared variate
	cdef double E = 0			## first draw from exponential used to construct chi-squared variate
	cdef double E2 = 0			## second draw from exponential used to construct chi-squared variate
	cdef double b = 0			## used to construct threshold for z
	cdef double a = 0			## used to construct threshold for z
	cdef double xb = 0			## transformation of b, maybe can combine with b?
	cdef double xa = 0			## transformation of a, maybe can combine with a?
	cdef double x0 = 0			## transformation of K, but can't combine with k
	cdef double qdivp = 0		## used to construct threshold for z, could combine with x0?
	cdef double mass_texp = 0	## transformation of qdivp, could combine with qdivp?
	cdef double PoGa = 0		## holder for aggregate Poyla-Gamma, eg. sum of X's, maybe combine with INPUT[ob]?
	cdef long nvar = 0			## counter for number of Poyla gamma draws
	cdef long ob = 0			## counter for number of observations
	### at this point we would add loop for all observations. So that:
	for ob in xrange(INPUT.shape[0]):
		PoGa = 0	## have to reset PoGa for each obs
		nvar = 0	## have to reset nvar for each obs
		N = ROWSUM[ob]
		z = abs(INPUT[ob]) * 0.5
		K = 0.125* pi*pi + 0.5* z*z	# this was error
		X, Y = 0, 0
		S = 1
		b = sqrt(1.0 / t) * (t * z - 1)
		a = sqrt(1.0 / t) * (t * z + 1) * -1.0
		x0 = log(K) + K * t
		xb = x0 - z + log(pnorm(b))
		xa = x0 + z + log(pnorm(a))
		qdivp = 4 / pi * ( exp(xb) + exp(xa) )
		mass_texp = 1.0 / (1.0 + qdivp)
		#print(mass_texp)#debug
		while nvar < N:
			if unif(r) < mass_texp: ## aka p/(p+q)
				#print(mass_texp)#debug
				## generate proposal from truncated exponential
				X = t + exponential(r, 1) / K
				#print("exp", X/4)#debug
			else:
				## generate proposal from truncated inverse Gaussian
				X = t+1 	## intialize X so that greater than t
				if z < 1/t:
					## when z is small, e.g. 1/z > t, then TIG() is well approximated by chi-squared distribution and accept/reject step is used.
					alpha = 0 ## set to zero so that condition not initially met
					#E = 0 ## set E, E2 so that sampler redraws with each new obs (doesn't reuse same E/E2)
					#E2 = 1
					while unif(r) > alpha:	## rejection criteria
						E = exponential(r, 1)
						E2 = exponential(r, 1)
						while (E*E) > (2*E2 / t):	# corrected this, conditional was other way
							E = exponential(r, 1)
							E2 = exponential(r, 1)
							#print(E,E2)#debug
						X = 1 + t*E
						X = t / (X*X)
						alpha = exp( -0.5 * z*z * X)
					#print("chi", X/4)#debug
				else:
					mu = 1/z
					## when z is not small, generate proposal from TIG(.) using rejection sampling.
					while X > t: ## rejection criteria
						Y = gaussian(r, 1)
						Y = Y*Y
						X = mu + 0.5*mu*mu*Y - 0.5*mu*sqrt(4*mu*Y + (mu*Y)*(mu*Y))
						if unif(r) > (mu/(mu+X)): 
							X = mu*mu/X
					#print("IG", X/4)#debug
			## having generated a proposal X, we decide whether to accept X
			## using series sampling method of Devroye (2009)
			S = a_n(0.0,X)
			Y = unif(r)*S
			n = 0
			while True:	## infinite loop
				n += 1
				#print(n, X, S, Y)#debug
				if n%2 == 1:
					S = S - a_n(n,X)
					if Y <= S:
						PoGa += (X * 0.25)
						#print("Sq", X)#debug
						nvar += 1
						break
				else:
					S = S + a_n(n,X)
					if Y > S: 
						break
				#if n > 1000:#debug
				#	sys.exit("Stuck in Poyla-Gamma sampler, exiting.")
		INPUT[ob] = PoGa  ## sum of independent Poyla-Gamma(1,x) draws
	##@ END LOOP OVER OBSERVATIONS
	return(INPUT)	# return vector of inputs, now Poyla-Gamma draws

cpdef inline double a_n(double n, double x):
	"""x > 0"""
	cdef double K = (n + 0.5) * pi
	if 0 < x <= t:
		return( exp( -1.5 * (log(0.5*pi) + log(x)) + log(K) - 2.0 * (n+0.5)*(n+0.5) / x ) ) # copied this from jesse's code, seems to give correct result
	#	return( pi*(n+0.5) * (2/(pi*x))**1.5 * exp(-(2*(n + 0.5)**2)/x) )
	if x > t:
		return( K * exp( -0.5 * K*K * x) )
	#	return( pi*(n+0.5) * exp(-0.5 * ((n+0.5)**2 * pi**2 * x) ) * x )
		
###############################################################################
##  Gibbs sampler using Polya-Gamma mixture for multinomial model			 ##
###############################################################################

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def MCMCsamp(long nMCMC, double[:,:] Y, double[:,:] X, double[:,:] Omega, long[:] fixef, double shape = 0.001, double rate = 0.001):
	"""Omega should be inverse"""
	##@ C-STYLE TYPE DECLARATIONS
	## model dimensionality
	cdef int Q = Y.shape[1]										## number of choices
	cdef int N = X.shape[0]										## number of observations
	cdef int P = X.shape[1]										## number of regression coefficients
	## indices
	cdef long n, i, q, f, j										## integers for loops
	cdef long [:,:] 	oQ = np.zeros([Q, Q-1], dtype = int)	## indices for each value of q, that are all Q except q
	cdef long [:] 		re = np.delete(np.arange(P), fixef)		## index for random effects
	## model parameters
	cdef double [:,:] 	v = np.zeros([P, P])					## prior variance-covariance matrix for beta
	cdef double [:,:] 	Z = np.zeros([N, Q-1])					## linear predictor
	cdef double [:,:] 	lam = np.zeros([Q-1, nMCMC])			## phylogenetic signal (scales Omega)
	cdef double [:,:,:] Lambda = np.zeros([N, N, Q-1])			## mixing weights for linear predictor
	cdef double [:,:,:] beta = np.zeros([P, Q, nMCMC])			## regression coefficients
	cdef double [:,:]	betaQ = np.zeros([P, Q-1])				## regression coefficients for all Q choices, excluding choice q 
	cdef double [:]		eta = np.zeros([P - len(fixef)])		## random effect subset of beta
	cdef double [:]		etaT = eta.T							## tranpose of eta
	cdef double [:] 	C = np.zeros([N])						## denominator of simplex
	cdef double [:,:] 	D = np.zeros([N,Q-1])					## exponent of linear predictor for all Q choices, excluding choice q
	cdef double [:] 	B = np.zeros([P])						## location parameters for conditional distribution of beta
	cdef double [:] 	T = np.zeros([P])						## random standard normal draws, to feed into beta
	cdef double [:,:] 	V = np.zeros([P, N])					## holder for Xt_Lam
	cdef double [:,:] 	L = np.zeros([P, P])					## cholesky decomposition of V
	cdef double [:,:] 	Xt = np.transpose(X)					## transpose of model matrix X
	cdef double [:,:] 	Xt_Lam = np.zeros([P,N])				## transpose of model matrix X %*% mixing weights Lambda
	cdef double [:,:] 	Xt_Lam_X = np.zeros([P,P])				## transpose of model matrix X %*% mixing weights Lambda %*% model matrix X
	cdef double [:,:] 	I = np.eye(N)							## identity matrix to initialize Lambda
	cdef double [:] 	m = np.zeros([N])						## 'fitted' values for linear predictor, eg. X %*% beta
	cdef double [:] 	p = np.zeros([N])						## probabilities used to generate random variates of linear predictor Z
	cdef double [:]		kappa = np.zeros([N])					## kappa is the observed freq. for a given event minus the total freq. for that event
	cdef double [:]		rowSumsY = np.zeros([N])				## row sums of Y -- e.g. the total freq. for each event ... should be 1 for this model
	
	cdef long reMi = re[0]										## lower index for random effects
	cdef long reMa = re[re.shape[0]-1]+1						## upper index for random effects

	## decomposition of Omega [assumes Omega has already been inverted]
	cdef double[:] eigO = np.linalg.eig(Omega)[0]		## eigenvalues of Omega
	cdef double[:,:] O = np.eye(Omega.shape[0])			## O is a diagonal matrix for the eigen values -> maybe can combine with G somehow
	cdef double[:,:] G = np.eye(Omega.shape[0])			## G is a diagonal matrix ...
	for i in xrange(Omega.shape[0]):
		G[i,i] = G[i,i]*eigO[i]							## ... of the eigenvalues
	cdef double[:,:] Pl = np.linalg.eig(Omega)[1]		## eigenvectors of Omega
	cdef double[:,:] Plt = Pl.T							## transpose of Pl
	cdef double[:,:] PlO = np.zeros([Omega.shape[0],Omega.shape[0]])		## for Pl %*% O
	cdef double sumL = np.log(np.sum(eigO))				## sum of log-eigenvalues

	## prior parameters for lam
	#cdef double shape = shape									## shape parameter
	#cdef double rate = rate									## rate parameter
	#cdef double shapeUpdate = 0								## shape parameter for conditional distribution of lam
	#cdef double rateUpdate = 0									## rate parameter for conditional distribution of lam

	## for slice sampler
	cdef double left = 0							## lower boundary for uniform prior on lam
	cdef double right = 1							## upper boundary for uniform prior on lam
	cdef double AUX = 0								## auxiliary variable, proportional to conditional probability of lam
	cdef double offset = 0							## auxiliary variable
	cdef double fA = 0								## conditional distribution of lam evaluated at lower limit
	cdef double fB = 0								## conditional distribution of lam evaluated at upper limit
	cdef double fAUX = 0							## conditional distribution of lam evaluated at AUX
	
	##@ GENERATE Q LEAVE-ONE-OUT INDICES
	for i in xrange(Q):
		for j in xrange(Q):
			if j < i: oQ[i, j] = j
			if j > i: oQ[i, j-1] = j

	##@ CALCULATE ROWSUMS OF Y
	cr.rowwise_sum(Y, out=rowSumsY)
	
	##@ PRIOR PARAMETERS
	## fix variance to high value for fixed effects
	for f in xrange(fixef.shape[0]):
		v[fixef[f], fixef[f]] = 1 / 100			## parameterized as precision (inverse of variance) 
	v[reMi:reMa,reMi:reMa] = Omega
	
	##@ INITIALIZE LAM, LAMBDA
	for q in xrange(Q-1):
		lam[q,0] = unif_interval(r, left, right)	## initialize lam
		## initialize Lambda to identity matrix
		Lambda[:,:,q] = I## this may cause problems later on. CHeck that is OK with kappa below
		
	#with nogil:
	##@ BEGIN GIBBS SAMPLER
	for i in xrange(1, nMCMC):					## loop over each MCMC iteration
		
		for q in xrange(Q-1):								## loop over each category			
			#@ STEP 1. UPDATE LAM WITH SLICE SAMPLER ## this doesn't take much time
			eta = beta[reMi:reMa, q, i-1]				## eta is random effect subset of beta
			AUX = condL(lam[q,i-1], G, eta, etaT, Pl, Plt, sumL, 0)## evaluate "proportional probability" at current value of lam (offset is 0)
			offset = unif_interval(r, 0, AUX)			## generate random number from interval 0-AUX
			fA = condL(left, G, eta, etaT, Pl, Plt, sumL, offset)	## evaluate "proportional probability" minus offset at lower limit
			fB = condL(right, G, eta, etaT, Pl, Plt, sumL, offset)	## evaluate "proportional probability" minus offset at upper limit
			fAUX = AUX - offset							## evaluate "proportional probability" minus offset at current value of lam
			# the following assumes that the distribution is monotonic with regards to the space fAUX--fA or fAUX--fB ... may or may not be true. can't prove that function condL is monotonic over all l/eta/Omega
			if fAUX * fA < 0:													## if fAUX and fA have different signs ...
				aRoot = zeroin(left, lam[q,i-1], fA, fAUX, G, eta, etaT, Pl, Plt, sumL, offset)	## set the lower boundary to where the density equals the offset
				#print(aRoot)#debug
			else:																## otherwise,
				aRoot = left													## set to the predefined boundary (e.g. 0)
			if fAUX * fB < 0:													## if fAUX and fB have different signs ...
				bRoot = zeroin(lam[q,i-1], right, fAUX, fB, G, eta, etaT, Pl, Plt, sumL, offset)  ## set the upper boundary to where the density equals the offset
				#print(bRoot)#debug
			else:																## otherwise,
				bRoot = right													## set to the predefined boundary (e.g. 1)
			lam[q,i] = unif_interval(r, aRoot, bRoot)							## random variate from interval aRoot--bRoot, is new draw of lam
			for f in xrange(eigO.shape[0]):
				O[f,f] = eigO[f] ** lam[q,i]									## O is holder for scaled diagonalization of Omega
			PlO = cc.dot_mm(Pl, O)							## using eta as holder for Pl %*% O
			v[reMi:reMa, reMi:reMa] = cc.dot_mm(PlO, Plt)	## update v to reflect new lam ... v[re,re] = Pl %*% G^lam %*% Plt
			
			#@ STEP 2. UPDATE C
			betaQ = beta[:, oQ[q,0]:(oQ[q,Q-2]+1), i-1]			## betaQ contains all beta except for current category (q) ## should this be i-1??
			D = cc.dot_mm(X, betaQ)								## D = X %*% betaQ 		# why doesn't out argument work??
			for n in xrange(N):
				for f in xrange(Q-1):
					D[n,f] = exp(D[n,f])						## D = exp(X %*% betaQ)
			cr.rowwise_sum(D, out = C)							## C = rowsum(exp(X %*% betaQ))
			for n in xrange(N):
				C[n] = log(C[n])								## C = log(rowsum(exp(X %*% betaQ))
			
			##@ STEP 3. UPDATE BETA
			## build covariance matrix for conditional distribution of beta
			cc.dot_mm(Xt, Lambda[:,:,q], out = Xt_Lam)		## Xt_Lam = X.T %*% Lambda
			cc.dot_mm(Xt_Lam, X, out = Xt_Lam_X)			## Xt_Lam_X = X.T %*% Lambda %*% X
			ce.add_mm(Xt_Lam_X, v, out = Xt_Lam_X)			## Xt_Lam_X = X.T %*% Lambda %*% X + v
			clu.iinv(Xt_Lam_X) 								## Xt_Lam_X = (X.T %*% Lambda %*% X + v)^-1 ... thus variance-covariance for new draw of beta
			## build location vector for conditional distribution of beta
			ce.add_vs(Y[:,q], -0.5, out = kappa)			## kappa = Y - 0.5
			ce.multiply_vv(rowSumsY, kappa, out = kappa)	## kappa = Y*rowSumsY - rowSumsY*0.5 ... will always be negative?
			cc.dot_mv(Lambda[:,:,q], C, out = m)			## m = Lambda %*% C, returns [N], co-opting m as holder ...
			#ce.multiply_vs(m, -1, out = m)					## m = -m
			ce.add_vv(kappa, m, out = m)					## m = (Y-rowSumsY/2) - Lambda %*% C ... returns [N] .. 
			cc.dot_mv(Xt, m, out = T) 						## T = Xt %*% m ... returns [P], co-opting T as holder ...
			cc.dot_mv(Xt_Lam_X, T, out = B)					## XtLamX %*% Xt %*% m ... returns [P] ... vector for beta
			## draw new beta from conditional distribution ... beta ~ MVN(B,V)
			for f in xrange(P):
				T[f] = gaussian(r, 1.0)							## T ~ N(0,1) ... thus random normal variate for each regression coefficient
			clt.cholesky(Xt_Lam_X, out = L)						## L is cholesky decomp. so that Xt_Lam_X = L %*% L.T
			cc.dot_mv(L, T, out = Xt_Lam_X[:,0])				## Xt_Lam_X[:,0] = L %*% T ... use Xt_Lam_X as 'workspace' for T
			ce.add_vv(B, Xt_Lam_X[:,0], out = beta[:, q, i])	## beta = B + L %*% T ... thus new draw for beta
			
			##@ STEP 4. UPDATE LAMBDA
			## compute new value of linear predictor
			cc.dot_mv(X, beta[:,q,i], out = m)				## m = X %*% beta
			ce.multiply_vs(C, -1, out = C)					## C = -C
			ce.add_vv(m, C, out = m)						## m = X %*% beta - C
			## draw Lambda from conditional distribution ...
			m = PG(rowSumsY, m)
			for n in xrange(N):
				Lambda[n,n,q] = m[n]						## assign new Lambdas
				##@ end loop over obs	
			##@ END LOOP OVER Q
	##@ END GIBBS SAMPLER
	##@ END NO-GIL
	return(beta, lam)
##@ END FUNCTION

## function that describes conditional distribution of lam ... could probably be optimized better ... maybe with pointers
cdef double condL(double l, double[:,:] G, double[:] eta, double[:] etaT, double[:,:] Pl, double[:,:] Plt, double sumL, double offset):	
	cdef long i	## iterator
	## maybe want to make slice sampler work on log scale, for underflow reasons with high-dimensional eta??
	for i in xrange(G.shape[0]):
		G[i,i] = G[i,i] ** l		## l does the actually scaling
	cdef double ev = -0.5 * cc.dot_vv(cc.dot_vm(cc.dot_vm(cc.dot_vm(etaT, Pl), G), Plt), eta)	## exponential term
	cdef double sl = l * sumL																	## determinant term
	return(exp(sl + ev) - offset)	## check this is correct??

cdef double zeroin(double ax, double bx, double fa, double fb, double[:,:] G, double[:] eta, double[:] etaT, double[:,:] Pl, double[:,:] Plt, double sumL, double offset):
	"""Use Brent's zero-finding algorithm to find the root of conditional distribution of lam.
		ax : left border of the range wherein to seek the root
		bx : right border of the range wherein to seek the root
		fa : the function evaluated at the lower boundary
		fb : the function evaluated at the upper boundary
	All other terms are for 'condL' function"""
	##@ C-STYLE TYPE DECLARATIONS
	cdef double a, b, c, fc		## Abscissae, f(c) */
	cdef double prev_step		 ## Holder for result of previous step
	cdef double tol_act			 ## Actual tolerance
	cdef double p				 ## Interpolation step is calcu-
	cdef double q				 ## lated in the form p/q; division operations is delayed until the last moment	*/
	cdef double new_step		 ## Step at this iteration
	cdef double EPSILON = 2**-53 ## machine precision eps
	cdef long MaxIter = 200
	cdef double Tol = 0.0001220703
	## interpolation threshholds
	cdef double t1
	cdef double cb
	cdef double t2

	##@ ROOT-FINDING ALGORITHM
	a = ax
	b = bx
	c = a
	fc = fa	
	## First test if we have found a root at an endpoint
	if fa == 0.0:
		return a
	if fb ==  0.0:
		return b
	## Otherwise, find root
	for i in xrange(MaxIter):		## Main iteration loop
		prev_step = b-a	## Distance from the last but one to the last approximation
		if abs(fc) < abs(fb):				
			## Swap data for b to be the best approximation
			a = b
			b = c
			c = a
			fa = fb
			fb = fc
			fc = fa
		tol_act = 2*EPSILON*abs(b) + Tol/2
		new_step = (c-b)/2
		if abs(new_step) <= tol_act or fb == 0.0:
			return b			## Acceptable approx. is found
		## Decide if the interpolation can be tried
		if abs(prev_step) >= tol_act and abs(fa) > abs(fb): ## If prev_step was large enough and was in true direction, interpolation may be tried
			cb = c - b
			if a == c: ## If we have only two distinct points linear interpolation can only be applied */
				t1 = fb/fa
				p = cb * t1
				q = 1.0 - t1
			else:		## Quadric inverse interpolation
				q = fa/fc
				t1 = fb/fc
				t2 = fb/fa
				p = t2 * ( cb * q * (q-t1) - (b-a) * (t1-1.0) )
				q = (q-1.0) * (t1-1.0) * (t2-1.0)
			if p > 0.0:
				q = -q			## p was calculated with the opposite sign; make p positive
			else:
				p = -p			## and assign possible minus to q
			if p < (0.75*cb*q-abs(tol_act*q)/2) and p < abs(prev_step*q/2):	## If b+p/q falls in [b,c] and isn't too large it is accepted
					new_step = p/q
		## If p/q is too large then the bisection procedure can reduce [b,c] range to greater extent
		if abs(new_step) < tol_act:	## Adjust the step to be not less than tolerance
			if new_step > 0.0:
				new_step = tol_act
			else:
				new_step = -tol_act
		a = b
		fa = fb				## Save the previous approx.
		b += new_step
		fb = condL(b, G, eta, etaT, Pl, Plt, sumL, offset)		## Do step to a new approxim.
		if (fb > 0 and fc > 0) or (fb < 0 and fc < 0):
			c = a
			fc = fa ## Adjust c for it to have a sign opposite to that of b
	##@ FAILED -- MAX ITERATIONS REACHED
	## print("Brent's algorithm failed to converge on solution")
	return b
##@ END BRENTS ALGORITHM
