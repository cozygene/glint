"""
Coded by Omer Weissbrod and Reut yedidim
"""
from utils import common, tools
from module import Module
import numpy as np
import scipy.linalg as la
import scipy.optimize as optimize
import scipy.stats as stats
import time
import sys
import scipy.linalg.blas as blas
import logging

AVAILABLE_KINSHIPS = ["refactor"]

def negLLevalLong(logdelta, s, Uy, UX, logdetXX, reml=True):
    Sd = s + np.exp(logdelta)
    UyS = Uy / Sd
    yKy = UyS.T.dot(Uy) 
    logdetK = np.log(Sd).sum()
    null_ll, beta_0, null_F = lleval(Uy, UX, Sd, yKy, logdetK, logdetXX, reml=reml)
    return -null_ll
        
def findLogDelta(U, s, phe, covars, numIntervals=100, ldeltamin=-5, ldeltamax=5, reml=True):
    """
    s - eigenvalues 
    U - eigenvectors
    phe - phenotype vector
    covars - covariates matrix
    """
    logging.info("computing log delta...")
    #Prepare required matrices
    Uy = np.dot(U.T, phe).flatten()
    UX = np.dot(U.T, covars)
    XX = covars.T.dot(covars)
    [Sxx,Uxx] = la.eigh(XX)
    logdetXX  = np.log(Sxx).sum()
    
    nllgrid = np.ones(numIntervals+1) * np.inf
    ldeltagrid = np.arange(numIntervals+1) / (numIntervals*1.0) * (ldeltamax-ldeltamin) + ldeltamin 
    for i in xrange(nllgrid.shape[0]):
        nllgrid[i] = negLLevalLong(ldeltagrid[i], s, Uy, UX, logdetXX, reml=reml)

    # find minimum
    nll_min = nllgrid.min()
    ldeltaopt_glob = ldeltagrid[nllgrid.argmin()]
    nllmin = nll_min.copy()
    
    # more accurate search around the minimum of the grid search
    for i in xrange(1, nllgrid.shape[0]-1):
        if (nllgrid[i]<nllgrid[i-1] and nllgrid[i]<nllgrid[i+1]):
            ###print (ldeltagrid[i-1],nllgrid[i-1]), (ldeltagrid[i],nllgrid[i]), (ldeltagrid[i+1],nllgrid[i+1])
            try:
                ldeltaopt,nllopt,iter,funcalls = optimize.brent(negLLevalLong,(s, Uy, UX, logdetXX, reml),(ldeltagrid[i-1],ldeltagrid[i],ldeltagrid[i+1]),full_output=True)
                if nllopt<nllmin:           
                    nllmin=nllopt
                    ldeltaopt_glob=ldeltaopt
            except: pass
    
    logging.debug("log delta is %s", ldeltaopt_glob)
    return ldeltaopt_glob
    

def lleval(Uy, UX, Sd, yKy, logdetK, logdetXX, reml=True):
    N = Uy.shape[0]
    D = UX.shape[1]

    UXS = UX / np.lib.stride_tricks.as_strided(Sd, (Sd.size, D), (Sd.itemsize,0))   
    XKy = UXS.T.dot(Uy) 
    
    XKX = UXS.T.dot(UX) 
    [SxKx,UxKx]= la.eigh(XKX)   
    i_pos = SxKx>1E-10
    beta = np.dot(UxKx[:,i_pos], (np.dot(UxKx[:,i_pos].T, XKy) / SxKx[i_pos]))
    r2 = yKy-XKy.dot(beta)

    if reml:
        logdetXKX = np.log(SxKx).sum()      
        sigma2 = r2 / (N - D)
        ll =  -0.5 * (logdetK + (N-D)*np.log(2.0*np.pi*sigma2) + r2/sigma2 + logdetXKX - logdetXX)
        detXKX = SxKx.prod()
        if (D==2):
            F = (beta[0]**2) * detXKX / (XKX[1,1] * sigma2)
        elif (D>2):
            F = ((UxKx[0,i_pos]**2) / SxKx[i_pos]).sum()
            F = beta[0]**2 / (F*sigma2)
        elif (D==1):
            F = (XKy[0]**2) / (XKX[0,0] * sigma2)
        else: F=0       
    else:
        sigma2 = r2 / (N)
        ll =  -0.5 * (logdetK + N*np.log(2.0*np.pi*sigma2) + r2/sigma2)
        F = 0
        
    return (ll, beta, F)


class KinshipCreator(Module):
    def __init__(self, kinship_data, is_normalized = False):
        """
        creates a kinship matrix from data
        
        params:
        kinship_data - the data to generate the kinship from
        is_normalized: True if the data supplied is normalized. False is the default must be normalized
                       Note that data is normalized according to axis = 0, so transopse the data before calling LMM if needed
        
        """
        self.data = kinship_data
        self.standardize = is_normalized

    
    def create_standard_kinship(self): #TODO change name
        """
        compute kinship matrix ( X * X.transpose() ) / (number of sites)
        returns matrix of dimensions nXn
        """
        logging.info("Creating a standard kinship...")
        X = self.data
        if not self.standardize:
            X = tools.standardize(X, axis = 0)

        # compute kinship matrix ( X * X.transpose() ) / (number of sites)
        return tools.symmetrize(blas.dsyrk(1.0, X, lower=1)) / X.shape[1]


class LMM(Module):
    def __init__(self, kinship_matrix):
        """
        initialize LMM with a kinship matrix
        Note - assumes kinship_data is of dimensions nXn where is number of samples
        """
        # preprocessing kinship
        self.s, self.U = tools.eigenDecompose(kinship_matrix)
    
    
    
    def run(self, data, pheno, covars, cpgnames, is_covars_normalized = False, logdelta = None, reml=True):
        """
        preprocess data and run lmm. 
        
        params:
        data - the methylation data to test (matrix of n sampels by m sites)
        pheno - the phenotype    (a 1D vector of size n (sampels) )
        covars - the covariates.(matrix of n sampels by x covariates or empty)
        is_covars_normalized - is the covariates matrix (supplied with param 'covar') is normalized (default is False).
                                if False - covars will be normalized.
                               Note that covariates will be normalized according to asix=0. transpose before calling this function if needed.
        
        data returned is sorted by pvalues (and is all of type ndarray)
        """
        # preprocess covariates and add a column of ones - before calc logdelta
        number_of_samples = pheno.shape[0]
        if covars is None:
            covars = np.empty((number_of_samples, 0))

        if not is_covars_normalized:
            covars = tools.standardize(covars, axis = 0)
        covars = np.concatenate((covars,np.ones((number_of_samples, 1))), axis=1)
        
        # if log delta is not supplied - calculate ir
        if logdelta is None:
            logdelta = findLogDelta(self.U, self.s, pheno, covars, reml=(reml>0))


        sorted_cpgnames, sorted_cpg_indices, p_vals, beta_est, sigma_e_est, sigma_g_est, statistics = \
                             self.lmm(data, pheno, covars, cpgnames, logdelta, reml)

        beta_est = np.array(beta_est)
        intercept_beta = beta_est[:,-1]     # interception coeff
        site_beta = beta_est[:,0]           # site coeff
        covariates_betas = beta_est[:,1:-1] # coeff for each covariate

        return  np.array(sorted_cpgnames),  \
                p_vals,                     \
                intercept_beta,             \
                covariates_betas,           \
                site_beta,                  \
                np.array(sigma_e_est),      \
                np.array(sigma_g_est),      \
                np.array(statistics)


    def lmm(self, data, phe, covars, cpgnames, logdelta, reml=True):
        """
        returns output sorted by pvalues:
        sorted_cpgnames, sorted_cpg_indices, p_vals, beta_est, sigma_e_est, sigma_g_est, statistics
        where beta_est is 2d array where beta_est[i] is the coefficients of site i 
            beta_est[i][0] is the coefficient of the interception
            beta_est[i][-1] is the coefficient of site i
            beta_est[i][1:-1] is the coefficient of the covariates
        """
        logging.info('Running LMM...')
        number_of_samples = phe.shape[0]
        t0 = time.time()
        
        #Prepare required matrices  
        Uy = np.dot(self.U.T, phe).flatten()
        UX = self.U.T.dot(covars)
        Sd = self.s + np.exp(logdelta)
        UyS = Uy / Sd
        yKy = UyS.T.dot(Uy) 
        logdetK = np.log(Sd).sum()

        num_of_non_zero_eigenvalues = len(Sd)
        num_of_zero_eigenvalues = number_of_samples - num_of_non_zero_eigenvalues
        logging.debug("found %d zero eigenvalue" % num_of_zero_eigenvalues)
        #Compute null LL
        if (covars.shape[1]>0):
            XX = covars.T.dot(covars)       
            [Sxx,Uxx]= la.eigh(XX)
            logdetXX  = np.log(Sxx).sum()
            null_ll, beta_0, null_F = lleval(Uy, UX, Sd, yKy, logdetK, logdetXX, reml=reml)
            logging.debug('null LL: %s' %null_ll)

        #Add an extra column to UX, that will hold UX for the tested site
        UX = np.concatenate((np.zeros((UX.shape[0], 1)), UX), axis=1)
        
        UX_all = self.U.T.dot(data)
        
        #Compute logdetXX - we assume it is the same for all sites because they are standardized 
        covars = np.concatenate((np.zeros((number_of_samples, 1)), covars), axis=1)
        covars[:,0] = data[:, 0] 
        XX = covars.T.dot(covars)
        [Sxx,Uxx]= la.eigh(XX)
        logdetXX  = np.log(Sxx).sum()

        #perform GWAS
        results = []
        for site_i, site_name in enumerate(cpgnames):  
            UX[:,0] = UX_all[:, site_i]
            
            ll, beta, F = lleval(Uy, UX, Sd, yKy, logdetK, logdetXX, reml=reml) # Note that the order of coefficient in beta is: site under test, covaraites, intercept
            # Calculate sigms_g, sigms_e
            sigma_g = np.sum([ ((Uy[i] - np.dot(UX[i,:],beta))**2) / Sd[i] for i in range(num_of_non_zero_eigenvalues)])
            sigma_g += np.sum([ ((Uy[i] - np.dot(UX[i,:],beta))**2) / np.exp(logdelta) for i in range(num_of_zero_eigenvalues)])
            if reml:
                sigma_g = (sigma_g/(number_of_samples-UX.shape[1]))**0.5
            else:
                sigma_g = (sigma_g/number_of_samples)**0.5

            sigma_e = (np.exp(logdelta) * (sigma_g**2))**0.5
            
            results.append((site_i, site_name, ll, F, beta, sigma_g, sigma_e))
                
        #sort and print results

        if reml:
            results.sort(key = lambda t: t[3], reverse=True)    
            fDist = stats.f(1, number_of_samples - 1)
            p_vals = fDist.sf([t[3] for t in results])
        else:
            results.sort(key = lambda t: t[2], reverse=True)
            chi2 = stats.chi2(1)
            p_vals = chi2.sf(2*(np.array([t[2] for t in results]) - null_ll))

        logging.info("LMM is done in %0.2f seconds" %(time.time()-t0))


        sorted_cpg_indices = [res[0] for res in results]    # sorted_cpg_indices[i] is the index of sorted_cpgnames[i] in cpgnames. i.e cpgnames[sorted_cpg_indices[i]] == sorted_cpgnames[i]
        sorted_cpgnames = [res[1] for res in results]

        beta_est = [res[4] for res in results]
        sigma_g_est = [res[5] for res in results]
        sigma_e_est = [res[6] for res in results]

        statistics=[]
        if reml:
            statistics = [res[3] for res in results]
        else:
            statistics = [res[2] for res in results]

        return sorted_cpgnames, sorted_cpg_indices, p_vals, beta_est, sigma_e_est, sigma_g_est, statistics




 