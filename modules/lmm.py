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

def compute_kinship(X):
    """
    param X: data of dimensions nXm where n is number of sampels and m is number of sites
             data must be normalized
    returns matrix of dimensions nXn
    """
    # compute kinship matrix ( X * X.transpose() ) / (number of sites)
    return tools.symmetrize(blas.dsyrk(1.0, X, lower=1)) / X.shape[1]

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


class LMM(Module):
    def __init__(self, kinship_data):
        """
        initialize LMM with the data for kinship matrix
        Note - assumes kinship_data is n sampels by m sites
        """
        X = tools.standardize(kinship_data, axis = 0)
        kinship = compute_kinship(X)
        # preprocessing kinship
        self.s, self.U = tools.eigenDecompose(kinship)
    
    
    
    def run(self, data, pheno, covars, cpgnames, logdelta = None, reml=True):
        """
        preprocess data and run lmm. 
        
        params:
        data - the methylation data to test (matrix of n sampels by m sites)
        pheno - the phenotype    (a 1D vector of size n (sampels) )
        covars - the covariates  (matrix of n sampels by x covariates or empty)
        """
        # preprocess covariates and add a column of ones - before calc logdelta
        number_of_samples = pheno.shape[0]
        if covars is None:
            covars = np.empty((number_of_samples, 0))

        covars = tools.standardize(covars, axis = 0)
        covars = np.concatenate((covars, np.ones((number_of_samples, 1))), axis=1)
        
        # if log delta is not supplied - calculate ir
        if logdelta is None:
            logdelta = findLogDelta(self.U, self.s, pheno, covars, reml=(reml>0))


        return self.lmm(data, pheno, covars, cpgnames, logdelta, reml)



    def lmm(self, data, phe, covars, cpgnames, logdelta, reml=True):

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
            
            ll, beta, F = lleval(Uy, UX, Sd, yKy, logdetK, logdetXX, reml=reml)
            results.append((site_i, site_name, ll, F))
                    
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
        return sorted_cpgnames, sorted_cpg_indices, p_vals 




 